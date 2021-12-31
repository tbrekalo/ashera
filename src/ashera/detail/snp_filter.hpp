#ifndef ASHERA_DETAIL_SNP_FILTER_HPP_
#define ASHERA_DETAIL_SNP_FILTER_HPP_

#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ashera::detail {

// NOTE: overlap indices must match read indices
[[nodiscard]] auto OverlapsFilterSnps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::vector<biosoup::Overlap>> overlaps,
    std::size_t keep_best_n) -> std::vector<std::vector<biosoup::Overlap>>;

}  // namespace ashera::detail

#endif /* ASHERA_DETAIL_SNP_FILTER_HPP_ */
