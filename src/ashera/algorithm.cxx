#include "ashera/algorithm.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <numeric>

#include "biosoup/timer.hpp"
#include "detail/mapping_algo.hpp"
#include "detail/snp_filter.hpp"
#include "edlib.h"
#include "fmt/compile.h"
#include "fmt/core.h"

namespace ashera {

namespace detail {

auto constexpr kAlignSmallBatchCap = 1ULL << 29ULL;  // ~512MB
auto constexpr kAlignBigBatchCap = 1ULL < 31ULL;     // ~2GB

}  // namespace detail

auto FindSnpFreeOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    RamConfig const ram_config,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads) -> std::vector<std::vector<biosoup::Overlap>> {
  auto timer = biosoup::Timer();

  timer.Start();
  auto overlaps = detail::MapSequences(thread_pool, reads);

  auto const n_ovlps_before = std::accumulate(
      overlaps.cbegin(), overlaps.cend(), 0U,
      [](std::uint32_t const init, std::vector<biosoup::Overlap> const& vec)
          -> std::uint32_t { return init + vec.size(); });

  fmt::print(
      stderr,
      FMT_COMPILE(
          "[ashera::FindSnpFreeOverlaps]({:12.3f}) : found {} overlaps\n"),
      timer.Stop(), n_ovlps_before);
  timer.Start();

  auto overlaps_filtered =  // TODO: think about exposing last param
      detail::OverlapsFilterSnps(thread_pool, reads, std::move(overlaps), 16U);

  auto const n_ovlps_after = std::accumulate(
      overlaps_filtered.cbegin(), overlaps_filtered.cend(), 0U,
      [](std::uint32_t const init, std::vector<biosoup::Overlap> const& vec)
          -> std::uint32_t { return init + vec.size(); });

  fmt::print(stderr,
             FMT_COMPILE("[ashera::Engine::FindSnpFreeOverlaps]({:12.3f}) : "
                         "found {} solid overlaps\n"),
             timer.Stop(), n_ovlps_after);

  return overlaps_filtered;
}

[[nodiscard]] auto PolishReads(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, PolishConfig config,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>&& reads,
    std::vector<std::vector<biosoup::Overlap>>&& overlaps)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  /*
      form batches
        -> put a memory limit on how much you gonna work with

      you need a concept of a window and an interval
        -> slice up the target sequence into windows of specified length
          -> take an overlap, calculate the alignment
              -> from alignment, split up the reference according to windows
                  -> windows are based on the query (NOTE: different than racon
     impl)
  */

  auto const find_batch_end_idx =
      [&reads](std::size_t const begin_idx, std::size_t const end_idx,
               std::size_t const batch_cap) -> std::size_t {
    auto sum = 0ULL;
    auto curr_idx = 0ULL;
    while (curr_idx < end_idx && sum < batch_cap) {
      sum += reads[curr_idx]->inflated_len;
      ++curr_idx;
    }

    return curr_idx;
  };
}

}  // namespace ashera
