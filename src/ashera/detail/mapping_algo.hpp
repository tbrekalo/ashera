#ifndef ASHERA_DETAIL_MAPPING_ALGO_HPP_
#define ASHERA_DETAIL_MAPPING_ALGO_HPP_

#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ashera::detail {

[[nodiscard]] auto MapSequences(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>>;

}

#endif /* ASHERA_DETAIL_MAPPING_ALGO_HPP_ */
