#ifndef ASHERA_IO_HPP_
#define ASHERA_IO_HPP_

#include <memory>
#include <string_view>
#include <vector>

#include "biosoup/nucleic_acid.hpp"

namespace ashera {

auto LoadReads(std::string const& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

auto LoadReads(std::vector<std::string> const& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

}  // namespace ashera

#endif /* ASHERA_IO_HPP_ */
