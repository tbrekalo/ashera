#include "ashera/engine.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "biosoup/timer.hpp"
#include "detail/overlap.hpp"
#include "edlib.h"
#include "fmt/compile.h"
#include "fmt/core.h"
#include "racon/polisher.hpp"
#include "ram/minimizer_engine.hpp"

namespace ashera {

namespace detail {

auto constexpr kMinimizeBatchCap = 1ULL << 32ULL;  // ~4GB
auto constexpr kMapBatchCap = 1ULL << 30ULL;       // ~1GB

auto constexpr kAlignBatchCap = 1ULL << 29ULL;  // ~512MB

auto constexpr kFilterFrequency = 0.001;

[[nodiscard]] auto CreateMinimizerEngine(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MinimapParams params)
    -> ram::MinimizerEngine {
  return ram::MinimizerEngine(std::move(thread_pool), params.k, params.w,
                              params.bandwidth, params.chain, params.matches,
                              params.gap);
}

class BaseCoverage {
 public:
  auto AnoteMismatch(char const base, std::uint32_t const read_id) -> void {
    for (auto base_idx = 0U; base_idx < kNBases; ++base_idx) {
      if (base == kBaseForIndex[base_idx]) {
        ++mismatch_cnt_[base_idx];
        mismatch_ids_.emplace_back(base_idx, read_id);
      }
    }
  }

  [[nodiscard]] auto GetSnps() const -> std::vector<std::uint32_t> {
    auto dst = std::vector<std::uint32_t>();
    auto indices = std::array<std::uint8_t, kNBases>{0, 1, 2, 3};

    std::sort(indices.begin(), indices.end(),
              [&](std::uint8_t const lhs, std::uint8_t const rhs) -> bool {
                return mismatch_cnt_[lhs] > mismatch_cnt_[rhs];
              });

    auto const mismatches_cnt = std::accumulate(
        std::next(indices.begin()), indices.end(), 0U,
        [&](std::uint32_t const init, std::uint8_t idx) -> std::uint32_t {
          return init + mismatch_cnt_[idx];
        });

    auto const cand_cnt = mismatch_cnt_[indices.front()];
    if (cand_cnt >= 3) {
      for (auto const& it : mismatch_ids_) {
        if (it.first == indices.front()) {
          dst.push_back(it.second);
        }
      }
    }

    return dst;
  }

  auto ClearMismatchIds() -> void { mismatch_ids_.clear(); }

 private:
  /* clang-format off */
  static constexpr auto kNBases = 4;
  static constexpr std::array<char, kNBases> kBaseForIndex = {
    'A', 'T', 'C', 'G'
  };
  /* clang-format on */

  std::array<std::uint32_t, kNBases> mismatch_cnt_;
  std::vector<std::pair<std::uint8_t, std::uint32_t>> mismatch_ids_;
};

[[nodiscard]] auto MapSequences(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto dst = std::vector<std::vector<biosoup::Overlap>>(reads.size());
  auto timer = biosoup::Timer();

  auto minimizer_engine = CreateMinimizerEngine(thread_pool, MinimapParams());

  auto const store_overlap = [&](biosoup::Overlap const& ovlp) -> void {
    dst[ovlp.lhs_id].push_back(ovlp);
    dst[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
  };

  auto const find_batch_end_idx =
      [&reads](std::size_t const begin_idx, std::size_t const end_idx,
               std::uint64_t const kBackCapacity) -> std::size_t {
    auto sum = 0UL;
    auto curr_idx = begin_idx;
    while (sum < kBackCapacity && curr_idx < end_idx) {
      sum += reads[curr_idx++]->inflated_len;
    }

    return curr_idx;
  };

  auto map_batch_futures =
      std::vector<std::future<std::vector<biosoup::Overlap>>>();
  auto const map_batch = [&](std::size_t const batch_begin,
                             std::size_t const batch_end,
                             auto&& ovlp_sink) -> void {
    static_assert(
        detail::IsOverlapSinkV<decltype(ovlp_sink)>,
        "ovlp_filer must satisfy ovlp_sink(biosoup::overlap) -> void");

    for (auto i = batch_begin; i < batch_end; ++i) {
      map_batch_futures.emplace_back(thread_pool->Submit(
          [&](std::size_t const idx) -> std::vector<biosoup::Overlap> {
            return minimizer_engine.Map(reads[idx], true, true);
          },
          i));
    }

    for (auto& future : map_batch_futures) {
      for (auto&& ovlp : future.get()) {
        ovlp_sink(ovlp);
      }
    }

    map_batch_futures.clear();
  };

  for (auto minimize_batch_begin = 0U; minimize_batch_begin < reads.size();) {
    timer.Start();
    auto const minimize_batch_end = find_batch_end_idx(
        minimize_batch_begin, reads.size(), detail::kMinimizeBatchCap);

    minimizer_engine.Minimize(std::next(reads.cbegin(), minimize_batch_begin),
                              std::next(reads.cbegin(), minimize_batch_end),
                              true);
    minimizer_engine.Filter(detail::kFilterFrequency);

    fmt::print(stderr,
               FMT_COMPILE("[ashera::detail::MapSequences]({:12.3f}) : "
                           "minimized {} sequences\n"),
               timer.Stop(), (minimize_batch_end - minimize_batch_begin));

    for (auto map_batch_begin = 0; map_batch_begin < minimize_batch_end;) {
      auto const map_batch_end = find_batch_end_idx(
          map_batch_begin, minimize_batch_end, detail::kMapBatchCap);

      timer.Start();
      map_batch(minimize_batch_begin, minimize_batch_end, store_overlap);

      fmt::print(stderr,
                 FMT_COMPILE("[ashera::detail::MapSequences]({:12.3f}) : "
                             "mapped {} sequences\n"),
                 timer.Stop(), (map_batch_end - map_batch_begin));

      map_batch_begin = map_batch_end;
    }

    minimize_batch_begin = minimize_batch_end;
  }

  return dst;
}

// NOTE: overlap indices must match read indices
[[nodiscard]] auto OverlapsToSnpFreeAligments(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::vector<biosoup::Overlap>> overlaps)
    -> std::vector<std::vector<std::pair<biosoup::Overlap, EdlibAlignResult>>> {
  // TODO: remove after loging
  auto snp_annotations = std::vector<std::vector<std::uint32_t>>(reads.size());

  auto const find_batch_end_idx =
      [&reads](std::size_t const begin_idx, std::size_t const end_idx,
               std::uint64_t const kBackCapacity) -> std::size_t {
    auto sum = 0UL;
    auto curr_idx = begin_idx;
    while (sum < kBackCapacity && curr_idx < end_idx) {
      sum += reads[curr_idx++]->inflated_len;
    }

    return curr_idx;
  };

  auto const overlap_strings =
      [&reads](
          biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
    auto const query_str = reads[ovlp.lhs_id]->InflateData(
        ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

    auto target_str = reads[ovlp.rhs_id]->InflateData(
        ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

    if (!ovlp.strand) {
      auto rc = biosoup::NucleicAcid("", target_str);
      rc.ReverseAndComplement();

      target_str = rc.InflateData();
    }

    return std::make_pair(query_str, target_str);
  };

  auto const align_strings =
      [](std::string const& query_str,
         std::string const& target_str) -> EdlibAlignResult {
    return edlibAlign(
        query_str.c_str(), query_str.size(), target_str.c_str(),
        target_str.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  };

  auto const annote_mismatches = [](biosoup::Overlap const& ovlp,
                                    std::string const& target_str,
                                    EdlibAlignResult const& edlib_align_res)
      -> std::unordered_map<std::uint32_t, detail::BaseCoverage> {
    auto dst = std::unordered_map<std::uint32_t, detail::BaseCoverage>();

    auto const target_id = ovlp.rhs_id;

    auto query_pos = ovlp.lhs_begin;
    auto target_pos = 0;

    for (auto i = 0U; i < edlib_align_res.alignmentLength; ++i) {
      switch (edlib_align_res.alignment[i]) {
        case 0: {  // match
          ++target_pos;
          ++query_pos;

          break;
        }

        case 1: {  // insertion on the target
          ++query_pos;
          break;
        }

        case 2: {  // insertion on the query
          ++target_pos;
          break;
        }

        case 3: {  // mismatch
          ++target_pos;
          ++query_pos;

          dst[query_pos].AnoteMismatch(target_str[target_pos], target_id);

          break;
        }

        default: {
          break;
        }
      }
    }

    return dst;
  };

  auto const overlaps_to_alignment_coverage =
      [&overlap_strings, &align_strings,
       &annote_mismatches](std::vector<biosoup::Overlap> const overlaps)
      -> std::pair<std::vector<std::pair<biosoup::Overlap, EdlibAlignResult>>,
                   std::unordered_map<std::uint32_t, detail::BaseCoverage>> {
    auto overlaps_alignments =
        std::vector<std::pair<biosoup::Overlap, EdlibAlignResult>>();

    auto coverage = std::unordered_map<std::uint32_t, detail::BaseCoverage>();

    for (auto const& ovlp : overlaps) {
      auto const [query_str, target_str] = overlap_strings(ovlp);
      auto const alignment = align_strings(query_str, target_str);
      coverage.merge(annote_mismatches(ovlp, target_str, alignment));

      overlaps_alignments.emplace_back(ovlp, alignment);
    }

    return std::make_pair(overlaps_alignments, coverage);
  };

  auto const detect_snp_variants =
      [&snp_annotations](
          [[maybe_unused]] std::uint32_t query_id,
          std::unordered_map<std::uint32_t, detail::BaseCoverage> coverage_map)
      -> std::unordered_set<std::uint32_t> {
    auto dst = std::unordered_set<std::uint32_t>();

    for (auto& [pos, coverage] : coverage_map) {
      auto const snps = coverage.GetSnps();
      if (!snps.empty()) {
        // TODO: clean up after loggin
        snp_annotations[query_id].push_back(pos);
      }

      std::copy(snps.cbegin(), snps.cend(), std::inserter(dst, dst.end()));
      coverage.ClearMismatchIds();
    }

    return dst;
  };

  auto const transform_and_filter_overlaps =
      [&overlaps_to_alignment_coverage, &detect_snp_variants](
          [[maybe_unused]] std::uint32_t const query_id,
          std::vector<biosoup::Overlap> overlaps)
      -> std::vector<std::pair<biosoup::Overlap, EdlibAlignResult>> {
    if (overlaps.empty()) {
      return {};
    }

    auto [ovlp_align, coverage_map] =
        overlaps_to_alignment_coverage(std::move(overlaps));

    auto const snp_variants =
        detect_snp_variants(query_id, std::move(coverage_map));

    auto const is_snp_variant =
        [&snp_variants](std::uint32_t const target_id) -> bool {
      return !snp_variants.empty() &&
             snp_variants.find(target_id) != snp_variants.end();
    };

    auto const remove_begin = std::remove_if(
        ovlp_align.begin(), ovlp_align.end(),
        [&is_snp_variant](decltype(ovlp_align)::const_reference ref) -> bool {
          return is_snp_variant(ref.first.rhs_id);
        });

    ovlp_align.erase(remove_begin, ovlp_align.end());
    ovlp_align.shrink_to_fit();

    return ovlp_align;
  };

  auto timer = biosoup::Timer();

  auto dst =
      std::vector<std::vector<std::pair<biosoup::Overlap, EdlibAlignResult>>>();
  dst.reserve(overlaps.size());

  auto align_futures = std::vector<std::future<
      std::vector<std::pair<biosoup::Overlap, EdlibAlignResult>>>>();

  for (auto align_batch_begin = 0U; align_batch_begin < overlaps.size();) {
    timer.Start();
    auto const align_batch_end =
        find_batch_end_idx(align_batch_begin, overlaps.size(), kAlignBatchCap);

    for (auto id = align_batch_begin; id < align_batch_end; ++id) {
      align_futures.emplace_back(thread_pool->Submit(
          transform_and_filter_overlaps, id, std::move(overlaps[id]))

      );
    }

    for (auto& it : align_futures) {
      dst.emplace_back(it.get());
    }

    fmt::print(
        stderr,
        FMT_COMPILE("[ashera::detail::OverlapsToSnpFreeAligments]({:12.3f}) :"
                    " aligned and filtered {} overlaps\n"),
        timer.Stop(), align_batch_end - align_batch_begin);

    align_batch_begin = align_batch_end;
    align_futures.clear();
  }

  auto snp_ostrm = std::ofstream("./snp_annotations.txt");
  for (auto read_id = 0U; read_id < reads.size(); ++read_id) {
    if (!snp_annotations[read_id].empty()) {
      std::sort(snp_annotations[read_id].begin(),
                snp_annotations[read_id].end());

      snp_ostrm << read_id << ' ';
      for (auto const it : snp_annotations[read_id]) {
        snp_ostrm << it << ' ';
      }

      snp_ostrm << '\n';
    }
  }

  std::flush(snp_ostrm);
  snp_ostrm.close();

  return dst;
}

}  // namespace detail

MinimapParams::MinimapParams(std::uint32_t k, std::uint32_t w,
                             std::uint32_t bandwidth, std::uint32_t chain,
                             std::uint32_t matches, std::uint32_t gap)
    : k(k),
      w(w),
      bandwidth(bandwidth),
      chain(chain),
      matches(matches),
      gap(gap) {}

Engine::Engine(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
               std::uint32_t win_size)
    : thread_pool_(std::move(thread_pool)), win_size_(win_size) {}

auto Engine::Correct(std::vector<std::unique_ptr<biosoup::NucleicAcid>>&& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto timer = biosoup::Timer();

  timer.Start();
  auto overlaps = detail::MapSequences(thread_pool_, reads);

  auto const n_ovlps_before = std::accumulate(
      overlaps.cbegin(), overlaps.cend(), 0U,
      [](std::uint32_t const init, std::vector<biosoup::Overlap> const& vec)
          -> std::uint32_t { return init + vec.size(); });

  fmt::print(
      stderr,
      FMT_COMPILE("[ashera::Engine::Correct]({:12.3f}) : found {} overlaps\n"),
      timer.Stop(), n_ovlps_before);
  timer.Start();

  // for (auto& ovlp_vec : overlaps) {
  //   std::sort(
  //       ovlp_vec.begin(), ovlp_vec.end(),
  //       [&](biosoup::Overlap const& lhs, biosoup::Overlap const& rhs) -> bool
  //       {
  //         return detail::OverlapScore(lhs) > detail::OverlapScore(rhs);
  //       });

  //   if (ovlp_vec.size() > 16) {
  //     ovlp_vec.resize(16);
  //     ovlp_vec.shrink_to_fit();
  //   }
  // }

  // fmt::print(stderr,
  //            FMT_COMPILE("[ashera::Engine::Correct]({:12.3f}) : discarded low
  //            "
  //                        "quality overlaps\n"),
  //            timer.Stop());

  // timer.Start();

  auto overlaps_alignmetns = detail::OverlapsToSnpFreeAligments(
      thread_pool_, reads, std::move(overlaps));

  auto const n_ovlps_after = std::accumulate(
      overlaps_alignmetns.cbegin(), overlaps_alignmetns.cend(), 0U,
      [](std::uint32_t const init,
         std::vector<std::pair<biosoup::Overlap, EdlibAlignResult>> const& vec)
          -> std::uint32_t { return init + vec.size(); });

  fmt::print(stderr,
             FMT_COMPILE("[ashera::Engine::Correct]({:12.3f}) : "
                         "generated {} overlap alignment pairs\n"),
             timer.Stop(), n_ovlps_after);

  return {};
}

}  // namespace ashera
