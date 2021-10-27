#include "ashera/engine.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>

#include "ashera/thread_pool.hpp"
#include "biosoup/timer.hpp"
#include "ram/minimizer_engine.hpp"

namespace ashera {

namespace detail {

auto constexpr kStampLen = 12UL;
auto constexpr kStampPrecision = 3UL;

// ~2GB
auto constexpr kMinimapBatchCap = 1U << 31U;
auto constexpr kFilterFrequency = 0.001;

auto CreateMinimizerEngine(MinimapParams params) -> ram::MinimizerEngine {
  return ram::MinimizerEngine(GetThreadPoolPtr(), params.k, params.w,
                              params.bandwidth, params.chain, params.matches,
                              params.gap);
};

enum class OverlapType : std::uint8_t {
  kInternal,
  kLhsContained,
  kRhsContained,
  kLhsToRhs,
  kRhsToLhs,

  kUnclassified
};

auto ClassifyOverlap(biosoup::Overlap const ovlp,
                     std::uint32_t const lhs_seq_size,
                     std::uint32_t const rhs_seq_size) -> OverlapType {
  auto const lhs_begin = ovlp.lhs_begin;
  auto const lhs_end = ovlp.lhs_end;

  auto const rhs_begin =
      ovlp.strand ? ovlp.rhs_begin : rhs_seq_size - ovlp.rhs_end;
  auto const rhs_end =
      ovlp.strand ? ovlp.rhs_end : rhs_seq_size - ovlp.rhs_begin;

  auto const lhs_ovlp_size = lhs_end - lhs_begin;
  auto const rhs_ovlp_size = rhs_end - rhs_begin;

  auto const ovlp_size = std::max(lhs_ovlp_size, rhs_ovlp_size);
  auto const overhang =
      std::min(lhs_begin, rhs_begin) +
      std::max(lhs_seq_size - lhs_end, rhs_seq_size - rhs_end);

  /*
    NOTE: order of overlap type evaluation is important

    if we do not check for internal overlaps
    before checking if one side is contained in another,
    we risk wrong classification
  */

  if (lhs_ovlp_size < (lhs_ovlp_size + overhang) * 0.875 ||
      rhs_ovlp_size < (rhs_ovlp_size + overhang) * 0.875) {
    return OverlapType::kInternal;

  } else if (lhs_begin < rhs_begin &&
             lhs_seq_size - lhs_end < rhs_seq_size - rhs_end) {
    return OverlapType::kLhsContained;

  } else if (rhs_begin < lhs_begin &&
             rhs_seq_size - rhs_end < lhs_seq_size - lhs_end) {
    return OverlapType::kRhsContained;

  } else if (lhs_begin > rhs_begin) {
    return OverlapType::kLhsToRhs;

  } else if (rhs_begin > lhs_begin) {
    return OverlapType::kRhsToLhs;
  }

  return OverlapType::kUnclassified;
}

auto ReverseOverlap(biosoup::Overlap const& ovlp) -> biosoup::Overlap {
  /* clang-format off */
  return biosoup::Overlap(
    ovlp.rhs_id, ovlp.rhs_begin, ovlp.rhs_end, 
    ovlp.lhs_id, ovlp.lhs_begin, ovlp.lhs_end,

    ovlp.score, ovlp.strand
  );
  /* clang-format on */
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

Engine::Engine(std::uint32_t win_size) : win_size_(win_size) {}

auto Engine::Correct(std::vector<std::unique_ptr<biosoup::NucleicAcid>>&& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto overlaps = std::vector<std::vector<biosoup::Overlap>>(reads.size());

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

  std::cerr << "[ashera::Engine::Correct]("
            /* clang-format off */
            << std::setw(detail::kStampLen)
            << std::setprecision(detail::kStampPrecision)
            << timer.Stop()
            /* clang-format on */
            << ") reordered sequences" << std::endl;

  auto minimizer_engine = detail::CreateMinimizerEngine(MinimapParams());
  auto const thread_pool = GetThreadPoolPtr();

  auto const mid_index = reads.size() / 2;

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
    std::vector<std::future<std::vector<biosoup::Overlap>>> batch_futures;
    auto const map_batch = [&](std::size_t const batch_begin,
                               std::size_t const batch_end) -> void {
      // TODO: think about rewrite with mutex
      for (auto i = batch_begin; i < batch_end; ++i) {
        batch_futures.emplace_back(thread_pool->Submit(
            [&](std::size_t const idx) -> std::vector<biosoup::Overlap> {
              return minimizer_engine.Map(reads[idx], true, true);
            },
            i));
      }

      for (auto& future : batch_futures) {
        for (auto&& ovlp : future.get()) {
          overlaps[ovlp.lhs_id].push_back(ovlp);
          overlaps[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
        }
      }

      batch_futures.clear();
    };

    for (auto batch_begin = 0UL; batch_begin < mid_index;) {
      auto const batch_end = find_batch_end_idx(batch_begin, mid_index);

      minimizer_engine.Minimize(std::next(reads.cbegin(), batch_begin),
                                std::next(reads.cbegin(), batch_end), true);
      minimizer_engine.Filter(detail::kFilterFrequency);

      std::cerr << "[ashera::Engine::Correct]("
                /* clang-format off */
                << std::setw(detail::kStampLen)
                << std::setprecision(detail::kStampPrecision)
                << timer.Stop()
                /* clang-format on */
                << ") minimized " << (batch_end - batch_begin) << " sequences"
                << std::endl;

      timer.Start();
      map_batch(batch_begin, batch_end);

      std::cerr << "[ashera::Engine::Correct]("
                /* clang-format off */
                << std::setw(detail::kStampLen)
                << std::setprecision(detail::kStampPrecision)
                << timer.Stop()
                /* clang-format on */
                << ") mapped " << (batch_end - batch_begin) << " sequences"
                << std::endl;

      batch_begin = batch_end;
    }

    for (auto batch_begin = 0UL; batch_begin < reads.size();) {
      auto const batch_end = find_batch_end_idx(batch_begin, reads.size());
      map_batch(batch_begin, batch_end);

      std::cerr << "[ashera::Engine::Correct]("
                /* clang-format off */
                << std::setw(detail::kStampLen)
                << std::setprecision(detail::kStampPrecision)
                << timer.Stop()
                /* clang-format on */
                << ") mapped " << (batch_end - batch_begin) << " sequences"
                << std::endl;

      batch_begin = batch_end;
    }
  }

  // TODO: change
  return {};
}

}  // namespace ashera
