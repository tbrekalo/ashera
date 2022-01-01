#ifndef ASHERA_CONFIGS_HPP_
#define ASHERA_CONFIGS_HPP_

#include <cstdint>

namespace ashera {

struct RamConfig {
  std::uint32_t k = 15;
  std::uint32_t w = 5;
  std::uint32_t bandwidth = 500;
  std::uint32_t chain = 4;
  std::uint32_t matches = 100;
  std::uint32_t gap = 10000;
};

struct SpoaConfig {
  std::int8_t match = 3;
  std::int8_t mismatch = -5;
  std::int8_t gap = -4;
};

struct PolishConfig {
  double quality_threashold = 10.0;
  double error_threshold = 0.3;

  std::uint32_t window_len = 500U;
  SpoaConfig spoa_config;
};

}  // namespace ashera

#endif /* ASHERA_CONFIGS_HPP_ */
