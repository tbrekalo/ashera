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
#include "detail/polish.hpp"
#include "detail/snp_filter.hpp"
#include "fmt/compile.h"
#include "fmt/core.h"

namespace ashera {

namespace detail {

auto constexpr kAlignBigBatchCap = 1ULL < 30ULL;  // ~1GB

}  // namespace detail

auto FindSnpFreeOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    RamConfig const ram_config,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
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
  auto timer = biosoup::Timer();

  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  dst.reserve(reads.size());
  /*
      form batches
        -> put a memory limit on how much you gonna work with

      you need a concept of a window and an interval
        -> slice up the target sequence into windows of specified length
          -> take an overlap, calculate the alignment
              -> from alignment, split up the reference according to windows
                  -> windows are based on the query
                    (NOTE: different than racon impl)
  */

  auto const find_batch_end_idx =
      [&reads, &overlaps](std::size_t const begin_idx,
                          std::size_t const end_idx,
                          std::size_t const batch_cap) -> std::size_t {
    auto sum = 0ULL;
    auto curr_idx = 0ULL;
    while (curr_idx < end_idx && sum < batch_cap) {
      sum += std::accumulate(
          overlaps[curr_idx].cbegin(), overlaps[curr_idx].cend(),
          reads[curr_idx]->inflated_len,
          [&reads](std::size_t const init,
                   biosoup::Overlap const& ovlp) -> std::size_t {
            return init + reads[ovlp.rhs_id]->inflated_len;
          });

      ++curr_idx;
    }

    return curr_idx;
  };

  auto const package_for_polish =
      [&reads](std::uint32_t const read_id,
               std::vector<biosoup::Overlap> ovlp_vec) -> detail::PolishPack {
    auto target = std::make_unique<biosoup::NucleicAcid>(*reads[read_id]);
    auto references = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
    references.reserve(ovlp_vec.size());

    std::transform(
        ovlp_vec.cbegin(), ovlp_vec.cend(), std::back_inserter(references),
        [&reads](biosoup::Overlap const ovlp) mutable
        -> std::unique_ptr<biosoup::NucleicAcid> {
          auto dst =
              std::make_unique<biosoup::NucleicAcid>(*reads[ovlp.rhs_id]);

          if (!ovlp.strand) {
            dst->ReverseAndComplement();
          }

          return dst;
        });

    // reindex things...
    auto index_lookup = std::vector<std::pair<std::uint32_t, std::uint32_t>>();
    index_lookup.reserve(ovlp_vec.size());  // TODO: swap out for inlined vector

    std::transform(
        references.cbegin(), references.cend(),
        std::back_inserter(index_lookup),
        [nxt_free_idx =
             0U](std::unique_ptr<biosoup::NucleicAcid> const& ref_ptr) mutable
        -> std::pair<std::uint32_t, std::uint32_t> {
          return std::make_pair(ref_ptr->id, nxt_free_idx++);
        });

    for (auto& it : ovlp_vec) {
      auto reindex_val =
          std::find_if(
              index_lookup.cbegin(), index_lookup.cend(),
              [ref_id = it.rhs_id](std::pair<std::uint32_t, std::uint32_t> puu)
                  -> bool { return puu.first == ref_id; })
              ->second;

      it.rhs_id = reindex_val;
    }

    return detail::PolishPack{
        .target = std::make_unique<biosoup::NucleicAcid>(*reads[read_id]),

        .overlaps = std::move(ovlp_vec),
        .references = std::move(references)};
  };

  timer.Start();

  auto polish_futures =
      std::vector<std::future<std::unique_ptr<biosoup::NucleicAcid>>>();
  for (auto polish_batch_begin = 0; polish_batch_begin < reads.size();) {
    auto const polish_batch_end = find_batch_end_idx(
        polish_batch_begin, reads.size(), detail::kAlignBigBatchCap);

    for (auto id = polish_batch_begin; id < polish_batch_end; ++id) {
      thread_pool->Submit(
          [&package_for_polish, &overlaps](std::uint32_t const target_id)
              -> std::unique_ptr<biosoup::NucleicAcid> {
            return detail::Polish(  // TODO: think of a better API design
                package_for_polish(target_id, std::move(overlaps[target_id])));
          },
          id);
    }

    for (auto& it : polish_futures) {
      dst.emplace_back(it.get());
    }

    fmt::print(
        stderr,
        FMT_COMPILE("[ashera::PolishReads]({:12.3f}) : polished {} reads\n"),
        timer.Stop(), polish_batch_end - polish_batch_begin);

    polish_batch_begin = polish_batch_end;
    polish_futures.clear();
  }

  return dst;
}

}  // namespace ashera
