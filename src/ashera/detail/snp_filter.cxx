#include "snp_filter.hpp"

#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "biosoup/timer.hpp"
#include "edlib.h"
#include "fmt/compile.h"
#include "fmt/core.h"

namespace ashera::detail {

auto constexpr kAlignSmallBatchCap = 1ULL << 29ULL;  // ~512MB
auto constexpr kAlignBigBatchCap = 1ULL < 31ULL;     // ~2GB

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
    auto const strongest_mismatch_idx = std::distance(
        mismatch_cnt_.cbegin(),
        std::max_element(mismatch_cnt_.cbegin(), mismatch_cnt_.cend()));

    auto const cand_cnt = mismatch_cnt_[strongest_mismatch_idx];
    if (cand_cnt >= 3) {
      for (auto const& it : mismatch_ids_) {
        if (it.first == strongest_mismatch_idx) {
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

// NOTE: overlap indices must match read indices
[[nodiscard]] auto OverlapsFilterSnps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::vector<biosoup::Overlap>> overlaps,
    std::size_t keep_best_n) -> std::vector<std::vector<biosoup::Overlap>> {
  auto snp_annotations = std::vector<std::vector<std::uint32_t>>(reads.size());

  auto const find_batch_end_idx =
      [&reads](std::size_t const begin_idx, std::size_t const end_idx,
               std::uint64_t const batch_capacity) -> std::size_t {
    auto sum = 0UL;
    auto curr_idx = begin_idx;
    while (sum < batch_capacity && curr_idx < end_idx) {
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

  auto const calculate_coverage =
      [&overlap_strings, &align_strings,
       &annote_mismatches](std::vector<biosoup::Overlap> const& overlaps)
      -> std::unordered_map<std::uint32_t, detail::BaseCoverage> {
    auto coverage = std::unordered_map<std::uint32_t, detail::BaseCoverage>();

    for (auto const& ovlp : overlaps) {
      auto const [query_str, target_str] = overlap_strings(ovlp);
      auto const alignment = align_strings(query_str, target_str);
      coverage.merge(annote_mismatches(ovlp, target_str, alignment));
    }

    return coverage;
  };

  auto const detect_snp_variants =
      [&snp_annotations](
          std::uint32_t query_id,
          std::unordered_map<std::uint32_t, detail::BaseCoverage> coverage_map)
      -> std::unordered_set<std::uint32_t> {
    auto dst = std::unordered_set<std::uint32_t>();

    for (auto& [pos, coverage] : coverage_map) {
      auto const snps = coverage.GetSnps();
      if (!snps.empty()) {
        snp_annotations[query_id].push_back(pos);
      }

      std::copy(snps.cbegin(), snps.cend(), std::inserter(dst, dst.end()));
      coverage.ClearMismatchIds();
    }

    return dst;
  };

  auto const filter_overlaps =
      [&calculate_coverage, &detect_snp_variants](
          std::uint32_t const query_id, std::vector<biosoup::Overlap> overlaps,
          std::size_t const keep_best_n) -> std::vector<biosoup::Overlap> {
    if (overlaps.empty()) {
      return {};
    }

    auto coverage_map = calculate_coverage(overlaps);

    auto const snp_variants =
        detect_snp_variants(query_id, std::move(coverage_map));

    if (!snp_variants.empty()) {
      auto const is_snp_variant =
          [&snp_variants](std::uint32_t const target_id) -> bool {
        return snp_variants.find(target_id) != snp_variants.end();
      };

      auto const remove_begin = std::remove_if(
          overlaps.begin(), overlaps.end(),
          [&is_snp_variant](biosoup::Overlap const& ovlp) -> bool {
            return is_snp_variant(ovlp.rhs_id);
          });

      overlaps.erase(remove_begin, overlaps.end());
    }

    if (overlaps.size() > keep_best_n) {
      overlaps.resize(keep_best_n);
    }

    overlaps.shrink_to_fit();

    return overlaps;
  };

  auto timer = biosoup::Timer();

  auto dst = std::vector<std::vector<biosoup::Overlap>>();
  dst.reserve(overlaps.size());

  auto align_futures =
      std::vector<std::future<std::vector<biosoup::Overlap>>>();

  for (auto align_batch_begin = 0U; align_batch_begin < overlaps.size();) {
    timer.Start();
    auto const align_batch_end = find_batch_end_idx(
        align_batch_begin, overlaps.size(), kAlignSmallBatchCap);

    auto const n_ovlps_in_batch = std::accumulate(
        std::next(overlaps.cbegin(), align_batch_begin),
        std::next(overlaps.cbegin(), align_batch_end), 0ULL,
        [](std::size_t const init,
           std::vector<biosoup::Overlap> const& ovlp_vec) -> std::size_t {
          return init + ovlp_vec.size();
        });

    for (auto id = align_batch_begin; id < align_batch_end; ++id) {
      align_futures.emplace_back(thread_pool->Submit(
          filter_overlaps, id, std::move(overlaps[id]), keep_best_n)

      );
    }

    for (auto& it : align_futures) {
      dst.emplace_back(it.get());
    }

    fmt::print(stderr,
               FMT_COMPILE("[ashera::detail::OverlapsFilterSnps]({:12.3f}) :"
                           " aligned and filtered {} overlaps\n"),
               timer.Stop(), n_ovlps_in_batch);

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
}  // namespace ashera::detail
