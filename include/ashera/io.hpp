#ifndef ASHERA_IO_HPP_
#define ASHERA_IO_HPP_

#include <memory>
#include <string_view>
#include <vector>

#include "biosoup/nucleic_acid.hpp"

namespace ashera {

/**
* @brief Load reads from the fasta/fastq file
*/
auto LoadReads(std::string const& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

/**
* @brief Load reads from multiple fasta/fastq files
*/
auto LoadReads(std::vector<std::string> const& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

}  // namespace ashera

#endif /* ASHERA_IO_HPP_ */
