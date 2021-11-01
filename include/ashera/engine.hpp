#ifndef ASHERA_ENGINE_HPP_
#define ASHERA_ENGINE_HPP_

#include <cstdint>
#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"

namespace ashera {

struct MinimapParams {
  MinimapParams() = default;
  MinimapParams(std::uint32_t k, std::uint32_t w, std::uint32_t bandwidth,
                std::uint32_t chain, std::uint32_t matches, std::uint32_t gap);

  std::uint32_t k = 15;
  std::uint32_t w = 5;
  std::uint32_t bandwidth = 500;
  std::uint32_t chain = 4;
  std::uint32_t matches = 100;
  std::uint32_t gap = 10000;
};

/**
 * @brief read correction engine
 */
class Engine {
 public:
  /**
   * @brief construct read correction engine
   */
  Engine(std::uint32_t win_size_);

  // disable copying
  Engine(Engine const&) = delete;
  Engine& operator=(Engine const&) = delete;

  // enable owndership transfer
  Engine(Engine&&) = default;
  Engine& operator=(Engine&&) = default;

  /**
   * @brief Generate corrected reads from the input ones
   */
  auto Correct(std::vector<std::unique_ptr<biosoup::NucleicAcid>>&& reads)
      -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

 private:
  std::uint32_t win_size_;
};

}  // namespace ashera

#endif /* ASHERA_ENGINE_HPP */
