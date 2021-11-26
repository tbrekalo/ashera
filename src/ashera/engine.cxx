#include "ashera/engine.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <type_traits>
#include <unordered_map>

#include "ashera/thread_pool.hpp"
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
auto constexpr kMapBatchCap = 1ULL << 30LLU;       // ~1GB

auto constexpr kAlignBatchCap = 1ULL < 28ULL;  // ~250MB

auto constexpr kFilterFrequency = 0.001;

auto CreateMinimizerEngine(MinimapParams params) -> ram::MinimizerEngine {
  return ram::MinimizerEngine(GetThreadPoolPtr(), params.k, params.w,
                              params.bandwidth, params.chain, params.matches,
                              params.gap);
}

struct Mismatch {
  std::uint32_t seq_id;
  std::uint32_t pos;

  std::uint32_t code;
};

class MismatchCnt {
 public:
  MismatchCnt() : cnts_{0U, 0U, 0U, 0U} {}

  auto Increment(std::uint8_t const code) noexcept -> void { ++cnts_[code]; }

  auto CountFor(std::uint8_t const code) const noexcept -> std::uint32_t {
    return cnts_[code];
  }

  auto MostFrequent() const noexcept -> std::pair<std::uint8_t, std::uint32_t> {
    auto const max_iter = std::max_element(cnts_.begin(), cnts_.end());
    auto const max_elem_code = std::distance(cnts_.begin(), max_iter);

    return std::make_pair(max_elem_code, *max_iter);
  };

 private:
  std::array<std::uint32_t, 4> cnts_;
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

Engine::Engine(float top_len_percentile, std::uint32_t win_size)
    : top_len_percentile_(top_len_percentile), win_size_(win_size) {}

auto Engine::Correct(std::vector<std::unique_ptr<biosoup::NucleicAcid>>&& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto timer = biosoup::Timer();

  auto overlaps = std::vector<std::vector<biosoup::Overlap>>(
      reads.size(), std::vector<biosoup::Overlap>());

  auto overlaps_begin_idx = std::vector<std::uint32_t>(reads.size(), 0);

  auto mismatches =
      std::vector<std::unordered_map<std::uint32_t, detail::MismatchCnt>>(
          reads.size());

  auto minimizer_engine = detail::CreateMinimizerEngine(MinimapParams());
  auto const thread_pool = GetThreadPoolPtr();

  timer.Start();

  {
    auto const align_sequences =
        [&](biosoup::Overlap const& ovlp) -> EdlibAlignResult {
      auto const query = reads[ovlp.lhs_id]->InflateData(
          ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

      auto const target = reads[ovlp.rhs_id]->InflateData(
          ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

      return edlibAlign(
          query.c_str(), query.size(), target.c_str(), target.size(),
          edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
    };

    // auto annotation_futures = std::vector<std::future<void>>();
    auto annotation_futures =
        std::vector<std::future<std::vector<detail::Mismatch>>>();

    auto const submit_annotate_mismatches =
        [&](biosoup::Overlap const& ovlp) -> std::vector<detail::Mismatch> {
      auto dst = std::vector<detail::Mismatch>();

      auto const edlib_result = align_sequences(ovlp);

      auto const query_id = ovlp.lhs_id;
      auto const target_id = ovlp.rhs_id;

      auto query_pos = ovlp.lhs_begin;
      auto target_pos = ovlp.rhs_begin;

      for (auto i = 0U; i < edlib_result.alignmentLength; ++i) {
        switch (edlib_result.alignment[i]) {
          case 0:
          case 3: {
            if (edlib_result.alignment[i] == 3) {
              dst.push_back(
                  detail::Mismatch{.seq_id = query_id,
                                   .pos = query_pos,
                                   .code = static_cast<std::uint8_t>(
                                       reads[target_id]->Code(target_pos))});

              dst.push_back(
                  detail::Mismatch{.seq_id = target_id,
                                   .pos = target_pos,
                                   .code = static_cast<std::uint8_t>(
                                       reads[query_id]->Code(query_pos))});
            }

            ++query_pos;
            ++target_pos;

            break;
          }
          case 1: {
            ++query_pos;
            break;
          }
          case 2: {
            ++target_pos;
            break;
          }
          default: {
            break;
          }
        }
      }

      return dst;
    };

    auto const store_overlap = [&](biosoup::Overlap const& ovlp) -> void {
      overlaps[ovlp.lhs_id].push_back(ovlp);
      overlaps[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
    };

    auto const find_batch_end_idx =
        [&reads](std::size_t const begin_idx,
                 std::size_t const end_idx) -> std::size_t {
      auto sum = 0UL;
      auto curr_idx = begin_idx;
      while (sum < detail::kMinimizeBatchCap && curr_idx < end_idx) {
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
          annotation_futures.emplace_back(
              thread_pool->Submit(submit_annotate_mismatches, ovlp));
        }
      }

      map_batch_futures.clear();

      // collect annotations
      for (auto& future : annotation_futures) {
        for (auto&& [seq_id, pos, code] : future.get()) {
          mismatches[seq_id][pos].Increment(code);
        }
      }

      annotation_futures.clear();
    };

    for (auto minimize_batch_begin = 0UL;
         minimize_batch_begin < reads.size();) {
      auto const minimize_batch_end =
          find_batch_end_idx(minimize_batch_begin, reads.size());

      minimizer_engine.Minimize(std::next(reads.cbegin(), minimize_batch_begin),
                                std::next(reads.cbegin(), minimize_batch_end),
                                true);
      minimizer_engine.Filter(detail::kFilterFrequency);

      fmt::print(
          stderr,
          FMT_COMPILE(
              "[ashera::Engine::Correct]({:12.3f}) : minimized {} sequences\n"),
          timer.Stop(), (minimize_batch_end - minimize_batch_begin));

      timer.Start();
      for (auto map_batch_begin = 0; map_batch_begin < minimize_batch_end;) {
        auto const map_batch_end =
            find_batch_end_idx(map_batch_begin, minimize_batch_end);

        map_batch(minimize_batch_begin, minimize_batch_end, store_overlap);

        fmt::print(
            stderr,
            FMT_COMPILE(
                "[ashera::Engine::Correct]({:12.3f}) : mapped {} sequences\n"),
            timer.Stop(), (map_batch_end - map_batch_begin));

        map_batch_begin = map_batch_end;
      }
      minimize_batch_begin = minimize_batch_end;
    }
  }

  timer.Start();

  return {};
}

}  // namespace ashera
