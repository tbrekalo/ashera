#include "overlap.hpp"

namespace ashera::detail {

auto DetermineOverlapType(biosoup::Overlap const ovlp,
                          std::uint32_t const lhs_seq_size,
                          std::uint32_t const rhs_seq_size) -> OverlapType {
  auto const lhs_begin = ovlp.lhs_begin;
  auto const lhs_end = ovlp.lhs_end;

  auto const rhs_begin =
      ovlp.strand ? ovlp.rhs_begin : rhs_seq_size - ovlp.rhs_end;
  auto const rhs_end =
      ovlp.strand ? ovlp.rhs_end : rhs_seq_size - ovlp.rhs_begin;

  auto const lhs_ovlp_size = lhs_end - lhs_begin;
  auto const rhs_ovlp_size = rhs_end - rhs_begin;

  auto const ovlp_size = std::max(lhs_ovlp_size, rhs_ovlp_size);
  auto const overhang =
      std::min(lhs_begin, rhs_begin) +
      std::max(lhs_seq_size - lhs_end, rhs_seq_size - rhs_end);

  /*
    NOTE: order of overlap type evaluation is important

    if we do not check for internal overlaps
    before checking if one side is contained in another,
    we risk wrong classification
  */

  if (lhs_ovlp_size < (lhs_ovlp_size + overhang) * 0.875 ||
      rhs_ovlp_size < (rhs_ovlp_size + overhang) * 0.875) {
    return OverlapType::kInternal;

  } else if (lhs_begin < rhs_begin &&
             lhs_seq_size - lhs_end < rhs_seq_size - rhs_end) {
    return OverlapType::kLhsContained;

  } else if (rhs_begin < lhs_begin &&
             rhs_seq_size - rhs_end < lhs_seq_size - lhs_end) {
    return OverlapType::kRhsContained;

  } else if (lhs_begin > rhs_begin) {
    return OverlapType::kLhsToRhs;

  } else if (rhs_begin > lhs_begin) {
    return OverlapType::kRhsToLhs;
  }

  return OverlapType::kUnclassified;
}

auto ReverseOverlap(biosoup::Overlap const& ovlp) -> biosoup::Overlap {
  /* clang-format off */
  return biosoup::Overlap(
    ovlp.rhs_id, ovlp.rhs_begin, ovlp.rhs_end, 
    ovlp.lhs_id, ovlp.lhs_begin, ovlp.lhs_end,

    ovlp.score, ovlp.strand
  );
  /* clang-format on */
}

auto OverlapLength(biosoup::Overlap const& ovlp) -> std::uint32_t {
  return std::max(ovlp.rhs_end - ovlp.rhs_begin, ovlp.lhs_end - ovlp.lhs_begin);
}

auto OverlapScore(biosoup::Overlap const& ovlp) -> double {
  return ((2. * OverlapLength(ovlp)) * ovlp.score) /
         (OverlapLength(ovlp) + ovlp.score);
}

}  // namespace ashera::detail
