#include "ashera/algorithm.hpp"

namespace ashera {

RamConfig::RamConfig(std::uint32_t k, std::uint32_t w,
                             std::uint32_t bandwidth, std::uint32_t chain,
                             std::uint32_t matches, std::uint32_t gap)
    : k(k),
      w(w),
      bandwidth(bandwidth),
      chain(chain),
      matches(matches),
      gap(gap) {}

}  // namespace ashera
