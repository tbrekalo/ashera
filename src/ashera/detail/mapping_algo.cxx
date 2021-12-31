#include "mapping_algo.hpp"

#include "ashera/configs.hpp"
#include "biosoup/timer.hpp"
#include "fmt/compile.h"
#include "fmt/core.h"
#include "overlap.hpp"
#include "ram/minimizer_engine.hpp"

namespace ashera::detail {

auto constexpr kMinimizeBatchCap = 1ULL << 32ULL;  // ~4GB
auto constexpr kMapBatchCap = 1ULL << 30ULL;       // ~1GB

auto constexpr kFilterFrequency = 0.001;

[[nodiscard]] auto CreateMinimizerEngine(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, RamConfig params)
    -> ram::MinimizerEngine {
  return ram::MinimizerEngine(std::move(thread_pool), params.k, params.w,
                              params.bandwidth, params.chain, params.matches,
                              params.gap);
}

[[nodiscard]] auto MapSequences(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto dst = std::vector<std::vector<biosoup::Overlap>>(reads.size());
  auto timer = biosoup::Timer();

  auto minimizer_engine = CreateMinimizerEngine(thread_pool, RamConfig());

  auto const store_overlap = [&](biosoup::Overlap const& ovlp) -> void {
    dst[ovlp.lhs_id].push_back(ovlp);
    dst[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
  };

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
}  // namespace ashera::detail
