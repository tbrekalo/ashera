#include "ashera/engine.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <type_traits>
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

auto constexpr kCovgBatchCap = 1ULL << 29ULL;  // ~512MB

auto constexpr kFilterFrequency = 0.001;

auto CreateMinimizerEngine(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                           MinimapParams params) -> ram::MinimizerEngine {
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

  auto GetSnps() const -> std::vector<std::uint32_t> {
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

    auto const& cand = mismatch_cnt_[indices.front()];
    if (cand >= 3 + mismatches_cnt) {
      for (auto const& it : mismatch_ids_) {
        dst.push_back(it.second);
      }
    }

    return dst;
  }

  auto ClearMismatchIds() -> void { mismatch_ids_.clear(); }

 private:
  static constexpr auto kNBases = 4;
  static constexpr std::array<char, kNBases> kBaseForIndex = {'A', 'T', 'C',
                                                              'G'};

  std::array<std::uint16_t, kNBases> mismatch_cnt_;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> mismatch_ids_;
};

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
  auto overlaps = MapSequences(reads);
  fmt::print(
      stderr,
      FMT_COMPILE("[ashera::Engine::Correct]({:12.3f}) : found {} overlaps\n"),
      timer.Stop(), overlaps.size());

  auto overlap_alignment =
      std::vector<std::pair<biosoup::Overlap, EdlibAlignResult>>(reads.size());

  timer.Start();

  for (auto& ovlp_vec : overlaps) {
    std::sort(
        ovlp_vec.begin(), ovlp_vec.end(),
        [&](biosoup::Overlap const& lhs, biosoup::Overlap const& rhs) -> bool {
          return detail::OverlapScore(lhs) > detail::OverlapScore(rhs);
        });

    // ovlp_vec.resize((ovlp_vec.size() * 2) / 3);
    if (ovlp_vec.size() > 16) {
      ovlp_vec.resize(16);
      ovlp_vec.shrink_to_fit();
    }
  }

  fmt::print(
      stderr,
      FMT_COMPILE("[ashera::Engine::Correct]({:12.3f}) : updated overlaps\n"),
      timer.Stop());

  timer.Start();

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

  {
    auto const align_sequences = [&](biosoup::Overlap const& ovlp)
        -> std::pair<std::string, EdlibAlignResult> {
      auto const get_target = [&]() {
        auto const inflated_data = reads[ovlp.rhs_id]->InflateData(
            ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);
        if (ovlp.strand) {
          return inflated_data;
        } else {
          auto rc = biosoup::NucleicAcid("", inflated_data);
          rc.ReverseAndComplement();

          return rc.InflateData();
        }
      };

      auto const query = reads[ovlp.rhs_id]->InflateData(
          ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

      auto target = get_target();
      auto result = edlibAlign(
          query.c_str(), query.size(), target.c_str(), target.size(),
          edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

      return std::make_pair(std::move(target), result);
    };

    auto const calculate_coverage = [&](std::uint32_t const query_id)
        -> std::unordered_map<std::uint32_t, detail::BaseCoverage> {
      auto dst = std::unordered_map<std::uint32_t, detail::BaseCoverage>();

      for (auto const& ovlp : overlaps[query_id]) {
        auto const [target, edlibResult] = align_sequences(ovlp);
        auto const target_id = ovlp.rhs_id;

        auto query_pos = ovlp.lhs_begin;
        auto target_pos = 0U;

        for (auto i = 0U; i < edlibResult.alignmentLength; ++i) {
          switch (edlibResult.alignment[i]) {
            case 0: {  // match
              ++query_pos;
              ++target_pos;

              break;
            }
            case 3: {  // mismatch
              dst[query_pos].AnoteMismatch(target[target_pos], target_id);
              ++query_pos;
              ++target_pos;

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
          }
        }

        edlibFreeAlignResult(edlibResult);
      }

      return dst;
    };

    auto detect_snp_variants =
        [&](std::uint32_t const query_id) -> std::unordered_set<std::uint32_t> {
      auto coverage_map = calculate_coverage(query_id);
      auto snp_reads = std::unordered_set<std::uint32_t>();

      for (auto& [pos, coverage] : coverage_map) {
        auto const snps = coverage.GetSnps();

        std::copy(snps.cbegin(), snps.cend(),
                  std::inserter(snp_reads, snp_reads.end()));
        coverage.ClearMismatchIds();
      }

      return snp_reads;
    };

    timer.Start();

    auto const ovlp_cnt = std::accumulate(
        overlaps.cbegin(), overlaps.cend(), 0ULL,
        [](std::uint64_t const init,
           std::vector<biosoup::Overlap> const& ovlp_vec) -> std::uint64_t {
          return init + ovlp_vec.size();
        });

    auto ovlp_discard_cnt = 0ULL;
    auto covg_futures = std::vector<std::future<void>>();
    for (auto covg_batch_begin = 0UL; covg_batch_begin < reads.size();) {
      auto const covg_batch_end = find_batch_end_idx(
          covg_batch_begin, reads.size(), detail::kCovgBatchCap);

      timer.Start();
      for (auto id = covg_batch_begin; id < covg_batch_end; ++id) {
        covg_futures.emplace_back(thread_pool_->Submit(
            [&](std::uint32_t const& query_id) -> void {
              auto const snp_variants = detect_snp_variants(query_id);
              auto const is_snp_varaint =
                  [&snp_variants](std::uint32_t const target_id) -> bool {
                return snp_variants.find(target_id) != snp_variants.end();
              };

              auto const remove_begin = std::remove_if(
                  overlaps[query_id].begin(), overlaps[query_id].end(),
                  [&](biosoup::Overlap const& ovlp) -> bool {
                    return is_snp_varaint(ovlp.rhs_id);
                  });

              ovlp_discard_cnt +=
                  std::distance(remove_begin, overlaps[query_id].end());

              overlaps[query_id].erase(remove_begin, overlaps[query_id].end());
              overlaps[query_id].shrink_to_fit();
            },
            id));
      }

      for (auto& it : covg_futures) {
        it.get();
      }

      covg_futures.clear();

      fmt::print(stderr,
                 FMT_COMPILE("[ashera::Engine::Correct]({:12.3f}) : purged "
                             "snp variants for {} reads\n"),
                 timer.Stop(), (covg_batch_end - covg_batch_begin));

      covg_batch_begin = covg_batch_end;
    }

    fmt::print(
        stderr,
        FMT_COMPILE(
            "[ashera::Engine::Correct] discarded {} out of {} overlaps\n"),
        ovlp_discard_cnt, ovlp_cnt);
  }

  return {};
}

auto Engine::MapSequences(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto dst = std::vector<std::vector<biosoup::Overlap>>(reads.size());
  auto timer = biosoup::Timer();

  {
    auto minimizer_engine =
        detail::CreateMinimizerEngine(thread_pool_, MinimapParams());

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
        map_batch_futures.emplace_back(thread_pool_->Submit(
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

    for (auto minimize_batch_begin = 0UL;
         minimize_batch_begin < reads.size();) {
      auto const minimize_batch_end = find_batch_end_idx(
          minimize_batch_begin, reads.size(), detail::kMinimizeBatchCap);

      minimizer_engine.Minimize(std::next(reads.cbegin(), minimize_batch_begin),
                                std::next(reads.cbegin(), minimize_batch_end),
                                true);
      minimizer_engine.Filter(detail::kFilterFrequency);

      fmt::print(stderr,
                 FMT_COMPILE("[ashera::Engine::MapSequences]({:12.3f}) : "
                             "minimized {} sequences\n"),
                 timer.Stop(), (minimize_batch_end - minimize_batch_begin));

      for (auto map_batch_begin = 0; map_batch_begin < minimize_batch_end;) {
        auto const map_batch_end = find_batch_end_idx(
            map_batch_begin, minimize_batch_end, detail::kMapBatchCap);

        timer.Start();
        map_batch(minimize_batch_begin, minimize_batch_end, store_overlap);

        fmt::print(stderr,
                   FMT_COMPILE("[ashera::Engine::MapSequences]({:12.3f}) : "
                               "mapped {} sequences\n"),
                   timer.Stop(), (map_batch_end - map_batch_begin));

        map_batch_begin = map_batch_end;
      }

      minimize_batch_begin = minimize_batch_end;
    }
  }

  return dst;
}

// auto Engine::PruneLowQualityOvlps(
//     std::vector<std::vector<biosoup::Overlap>>&& overlaps)
//     -> std::vector<std::vector<biosoup::Overlap>> {
//
//   return overlaps;
// }

}  // namespace ashera
