#ifndef ASHERA_DETAIL_POLISH_HPP_
#define ASHERA_DETAIL_POLISH_HPP_

#include <memory>
#include <vector>

#include "ashera/configs.hpp"

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"

namespace ashera::detail {

[[nodiscard]] auto Polish(
    PolishConfig config, 
    std::unique_ptr<biosoup::NucleicAcid> target,
    std::vector<biosoup::Overlap> overlaps,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> references)
    -> std::unique_ptr<biosoup::NucleicAcid>;

}  // namespace ashera::detail

#endif /* ASHERA_DETAIL_POLISH_HPP_ */
