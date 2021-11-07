#include "ashera/engine.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <type_traits>

#include "ashera/thread_pool.hpp"
#include "biosoup/timer.hpp"
#include "detail/assembly.hpp"
#include "fmt/compile.h"
#include "fmt/core.h"
#include "ram/minimizer_engine.hpp"

namespace ashera {

namespace detail {

// ~2GB
auto constexpr kMinimapBatchCap = 1U << 31U;
auto constexpr kFilterFrequency = 0.001;

auto CreateMinimizerEngine(MinimapParams params) -> ram::MinimizerEngine {
  return ram::MinimizerEngine(GetThreadPoolPtr(), params.k, params.w,
                              params.bandwidth, params.chain, params.matches,
                              params.gap);
}
auto constexpr kEdgeCandiatePoolSz = 8;

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
  auto edge_candidates = std::vector<detail::EdgeCandidate>();

  auto timer = biosoup::Timer();
  timer.Start();

  // sort reads by lenght reads
  std::sort(reads.begin(), reads.end(),
            [](std::unique_ptr<biosoup::NucleicAcid> const& lhs,
               std::unique_ptr<biosoup::NucleicAcid> const& rhs) -> bool {
              return lhs->inflated_len > rhs->inflated_len;
            });

  // reindex ids
  for (auto i = 0U; i < reads.size(); ++i) {
    reads[i]->id = i;
  }

  fmt::print(
      stderr,
      FMT_COMPILE(
          "[ashera::Engine::Correct]({:12.3f}s) : reordered sequences\n"),
      timer.Stop());

  auto minimizer_engine = detail::CreateMinimizerEngine(MinimapParams());
  auto const thread_pool = GetThreadPoolPtr();

  auto const len_pivot =
      static_cast<decltype(reads.size())>(reads.size() * top_len_percentile_);

  auto const find_batch_end_idx =
      [&reads](std::size_t const begin_idx,
               std::size_t const end_idx) -> std::size_t {
    auto sum = 0UL;
    auto curr_idx = begin_idx;
    while (sum < detail::kMinimapBatchCap && curr_idx < end_idx) {
      sum += reads[curr_idx++]->inflated_len;
    }

    return curr_idx;
  };

  timer.Start();

  {
    auto const is_edge_ovlp = [&reads](biosoup::Overlap const& ovlp) -> bool {
      auto const ovlp_type =
          detail::DetermineOverlapType(ovlp, reads[ovlp.lhs_id]->inflated_len,
                                       reads[ovlp.rhs_id]->inflated_len);

      return ovlp_type == detail::OverlapType::kLhsToRhs ||
             ovlp_type == detail::OverlapType::kRhsToLhs;
    };

    std::vector<std::future<std::vector<biosoup::Overlap>>> batch_futures;
    auto const map_batch = [&](std::size_t const batch_begin,
                               std::size_t const batch_end,
                               auto&& ovlp_filer) -> void {
      static_assert(
          detail::IsOverlapFilerV<decltype(ovlp_filer)>,
          "ovlp_filer must satisfy ovlp_filter(biosoup::overlap) -> bool");

      for (auto i = batch_begin; i < batch_end; ++i) {
        batch_futures.emplace_back(thread_pool->Submit(
            [&](std::size_t const idx) -> std::vector<biosoup::Overlap> {
              return minimizer_engine.Map(reads[idx], true, true);
            },
            i));
      }

      for (auto& future : batch_futures) {
        for (auto&& ovlp : future.get()) {
          if (ovlp_filer(ovlp)) {
            edge_candidates.emplace_back(ovlp);
          }
        }
      }

      batch_futures.clear();
    };

    for (auto minimize_batch_begin = 0UL; minimize_batch_begin < len_pivot;) {
      auto const minimize_batch_end =
          find_batch_end_idx(minimize_batch_begin, len_pivot);

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
      for (auto map_batch_begin = minimize_batch_begin;
           map_batch_begin < reads.size();) {
        auto const map_batch_end =
            find_batch_end_idx(map_batch_begin, reads.size());
        map_batch(minimize_batch_begin, minimize_batch_end, is_edge_ovlp);

        fmt::print(
            stderr,
            FMT_COMPILE(
                "[ashera::Engine::Correct]({:12.3f}) : mapped {} sequences\n"),
            timer.Stop(), (map_batch_end - map_batch_begin));

        map_batch_begin = map_batch_end;
      }

      minimize_batch_begin = minimize_batch_end;
    }

    timer.Start();

    auto const n_half_edges = edge_candidates.size() / 2;
    std::partial_sort(edge_candidates.begin(),
                      std::next(edge_candidates.begin(), n_half_edges),
                      edge_candidates.end(),
                      [](detail::EdgeCandidate const& lhs,
                         detail::EdgeCandidate const& rhs) -> bool {
                        return lhs.confidence_ > rhs.confidence_;
                      });

    edge_candidates.erase(std::next(edge_candidates.begin(), n_half_edges),
                          edge_candidates.end());
    edge_candidates.shrink_to_fit();

    fmt::print(
        stderr,
        FMT_COMPILE(
            "[ashera::Engine::Correct]({:12.3f}) : {} edge candidates\n"),
        timer.Stop(), edge_candidates.size());
  }

  auto back_bones = detail::AssembleBackbones(reads, edge_candidates);

  return {};
}

}  // namespace ashera
