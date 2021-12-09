#include "ashera/engine.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <optional>
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
auto constexpr kMapBatchCap = 1ULL << 30ULL;       // ~1GB

auto constexpr kPileBatchCap = 1ULL << 29ULL;  // ~512MB

auto constexpr kFilterFrequency = 0.001;

auto CreateMinimizerEngine(MinimapParams params) -> ram::MinimizerEngine {
  return ram::MinimizerEngine(GetThreadPoolPtr(), params.k, params.w,
                              params.bandwidth, params.chain, params.matches,
                              params.gap);
}

struct Coverage {
  std::uint32_t a = 0U;
  std::uint32_t t = 0U;
  std::uint32_t c = 0U;
  std::uint32_t g = 0U;
  std::uint32_t i = 0U;
};

using Pile = std::vector<Coverage>;

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

Engine::Engine(std::uint32_t win_size) : win_size_(win_size) {}

auto Engine::Correct(std::vector<std::unique_ptr<biosoup::NucleicAcid>>&& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto timer = biosoup::Timer();

  auto overlaps = std::vector<std::vector<biosoup::Overlap>>(
      reads.size(), std::vector<biosoup::Overlap>());

  auto overlaps_begin_idx = std::vector<std::uint32_t>(reads.size(), 0);
  auto piles = std::vector<detail::Pile>();

  piles.reserve(reads.size());
  std::transform(
      reads.begin(), reads.end(), std::back_inserter(piles),
      [](std::unique_ptr<biosoup::NucleicAcid> const& read) -> detail::Pile {
        return detail::Pile(read->inflated_len);
      });

  auto minimizer_engine = detail::CreateMinimizerEngine(MinimapParams());
  auto const thread_pool = GetThreadPoolPtr();

  timer.Start();

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

      auto const target = get_target();

      return std::make_pair(
          std::move(query),
          edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(),
                     edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH,
                                         nullptr, 0)));
    };

    auto const update_pile = [&](std::uint32_t const pile_id) -> void {
      for (auto const& ovlp : overlaps[pile_id]) {
        auto const [target, edlibResult] = align_sequences(ovlp);

        auto query_pos = ovlp.lhs_begin;
        auto target_pos = 0U;

        for (auto i = 0U; i < edlibResult.alignmentLength; ++i) {
          switch (edlibResult.alignment[i]) {
            case 0:
            case 3: {
              switch (target[target_pos]) {
                /* clang-format off */
                case 'A': ++piles[pile_id][query_pos].a; break;
                case 'T': ++piles[pile_id][query_pos].t; break;
                case 'C': ++piles[pile_id][query_pos].c; break;
                case 'G': ++piles[pile_id][query_pos].t; break;
                  /* clang-format on */
              }

              ++query_pos;
              ++target_pos;

              break;
            }
            case 1: {
              ++piles[pile_id][query_pos].i;
              ++query_pos;
            }
            case 2: {
              ++target_pos;
            }
          }
        }

        edlibFreeAlignResult(edlibResult);
      }
    };

    auto pile_futures = std::vector<std::future<void>>();
    auto const update_piles = [&](std::size_t const batch_begin,
                                  std::size_t const batch_end) -> void {
      for (auto pile_id = batch_begin; pile_id < batch_end; ++pile_id) {
        pile_futures.emplace_back(thread_pool->Submit(update_pile, pile_id));
      }
    };

    // auto annotation_futures = std::vector<std::future<void>>();
    auto const store_overlap = [&](biosoup::Overlap const& ovlp) -> void {
      overlaps[ovlp.lhs_id].push_back(ovlp);
      overlaps[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
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

    for (auto minimize_batch_begin = 0UL;
         minimize_batch_begin < reads.size();) {
      auto const minimize_batch_end = find_batch_end_idx(
          minimize_batch_begin, reads.size(), detail::kMinimizeBatchCap);

      minimizer_engine.Minimize(std::next(reads.cbegin(), minimize_batch_begin),
                                std::next(reads.cbegin(), minimize_batch_end),
                                true);
      minimizer_engine.Filter(detail::kFilterFrequency);

      fmt::print(
          stderr,
          FMT_COMPILE(
              "[ashera::Engine::Correct]({:12.3f}) : minimized {} sequences\n"),
          timer.Stop(), (minimize_batch_end - minimize_batch_begin));

      for (auto map_batch_begin = 0; map_batch_begin < minimize_batch_end;) {
        auto const map_batch_end = find_batch_end_idx(
            map_batch_begin, minimize_batch_end, detail::kMapBatchCap);

        timer.Start();
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

    for (auto pile_batch_begin = 0UL; pile_batch_begin < reads.size();) {
      auto const pile_batch_end = find_batch_end_idx(
          pile_batch_begin, reads.size(), detail::kPileBatchCap);

      timer.Start();
      update_piles(pile_batch_begin, pile_batch_end);

      for (auto& it : pile_futures) {
        it.wait();
      }

      pile_futures.clear();

      fmt::print(stderr,
                 FMT_COMPILE(
                     "[ashera::Engine::Correct]({:12.3f}): updated {} piles\n"),
                 timer.Stop(), (pile_batch_end - pile_batch_begin));

      pile_batch_begin = pile_batch_end;
    }
  }

  auto mstrm = std::ofstream("mismatches");
  auto log_vec = std::vector<std::uint32_t>();

  std::flush(mstrm);

  timer.Start();

  return {};
}

}  // namespace ashera
