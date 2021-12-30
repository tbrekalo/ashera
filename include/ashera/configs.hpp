#ifndef ASHERA_CONFIGS_HPP_
#define ASHERA_CONFIGS_HPP_

#include <cstdint>

namespace ashera {

struct RamConfig {
  RamConfig() = default;
  RamConfig(std::uint32_t k, std::uint32_t w, std::uint32_t bandwidth,
                std::uint32_t chain, std::uint32_t matches, std::uint32_t gap);

  std::uint32_t k = 15;
  std::uint32_t w = 5;
  std::uint32_t bandwidth = 500;
  std::uint32_t chain = 4;
  std::uint32_t matches = 100;
  std::uint32_t gap = 10000;
};

}  // namespace ashera

#endif /* ASHERA_CONFIGS_HPP_ */
