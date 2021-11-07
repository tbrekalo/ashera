#ifndef ASHERA_DETAIL_OVERLAP_HPP_
#define ASEHRA_DETAIL_OVERLAP_HPP_

#include <cstdint>
#include <type_traits>

#include "biosoup/overlap.hpp"

namespace ashera::detail {

/**
 * @brief specifies overlap types
 */
enum class OverlapType : std::uint8_t {
  kInternal,      ///< overlap has a deep location in one of the sequnces
  kLhsContained,  ///< rhs sequence spans over lhs sequence
  kRhsContained,  ///< lhs sequence spans over rhs sequence
  kLhsToRhs,      ///< lhs sequence connects to rhs sequence
  kRhsToLhs,      ///< rhs sequence connects to lhs sequence

  kUnclassified
};

/**
 * @brief concept for overlap filter callable:
 *         Filter(biosoup::Overlap const&) -> bool;
 */
template <class Impl, class = std::void_t<>>
struct OverlapFilterConcept : std::false_type {};

template <class Impl>
struct OverlapFilterConcept<Impl,
                            std::void_t<decltype(std::declval<Impl>()(
                                std::declval<biosoup::Overlap const&>()))>>
    : std::is_invocable_r<bool, Impl, biosoup::Overlap const&> {};

/**
 * @brief helper for @ref OverlapFilterConcept
 */
template <class T>
bool constexpr IsOverlapFilerV = OverlapFilterConcept<T>::value;

/**
 * @brief Determine overlap type provided origin sequence length information
 */
auto DetermineOverlapType(biosoup::Overlap const ovlp,
                          std::uint32_t const lhs_seq_size,
                          std::uint32_t const rhs_seq_size) -> OverlapType;

/**
 * @brief swap rhs and lhs
 */
auto ReverseOverlap(biosoup::Overlap const& ovlp) -> biosoup::Overlap;

}  // namespace ashera::detail

#endif /* ASHERA_DETAIL_OVERLAP_HPP_ */
