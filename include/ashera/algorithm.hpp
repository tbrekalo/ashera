#ifndef ASHERA_ALGORITHM_HPP_
#define ASHERA_ALGORITHM_HPP_

#include <cstdint>
#include <memory>
#include <vector>

#include "ashera/configs.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ashera {

[[nodiscard]] auto FindSnpFreeOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    RamConfig const ram_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>>;

}  // namespace ashera

#endif /* ASHERA_ALGORITHM_HPP */
