#ifndef ASHERA_DETAIL_OVERLAP_HPP_
#define ASEHRA_DETAIL_OVERLAP_HPP_

#include <cstdint>
#include <type_traits>

#include "biosoup/overlap.hpp"

namespace ashera::detail {

enum class OverlapType : std::uint8_t {
  kInternal,
  kLhsContained,
  kRhsContained,
  kLhsToRhs,
  kRhsToLhs,

  kUnclassified
};

template <class Impl, class = std::void_t<>>
struct OverlapFilterConcept : std::false_type {};

template <class Impl>
struct OverlapFilterConcept<Impl,
                            std::void_t<decltype(std::declval<Impl>()(
                                std::declval<biosoup::Overlap const&>()))>>
    : std::is_invocable_r<bool, Impl, biosoup::Overlap const&> {};

template <class T>
bool constexpr IsOverlapFiler = OverlapFilterConcept<T>::value;

auto DetermineOverlapType(biosoup::Overlap const ovlp,
                          std::uint32_t const lhs_seq_size,
                          std::uint32_t const rhs_seq_size) -> OverlapType;

auto ReverseOverlap(biosoup::Overlap const& ovlp) -> biosoup::Overlap;

}  // namespace ashera::detail

#endif /* ASHERA_DETAIL_OVERLAP_HPP_ */
