#ifndef ASHERA_DETAIL_ASSEMBLY_HPP_
#define ASHERA_DETAIL_ASSEMBLY_HPP_

#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "overlap.hpp"

namespace ashera::detail {

struct EdgeCandidate {
  EdgeCandidate(biosoup::Overlap ovlp) : ovlp(ovlp) {
    auto const len =
        std::max(ovlp.lhs_end - ovlp.lhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

    confidence_ = static_cast<float>(ovlp.score) / len;
  }

  biosoup::Overlap ovlp;
  float confidence_;
};

auto AssembleBackbones(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::vector<EdgeCandidate> const& edge_candidates)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

}  // namespace ashera::detail

#endif /* ASHERA_DETAIL_ASSEMBLY_HPP_ */
