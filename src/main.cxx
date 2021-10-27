#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "ashera/engine.hpp"
#include "ashera/io.hpp"
#include "ashera/thread_pool.hpp"
#include "biosoup/timer.hpp"
#include "cxxopts.hpp"

int main(int argc, char** argv) {
  try {
    auto options =
        cxxopts::Options("ashera", "Ashera is a read correction tool");

    /* clang-format off */
    options.add_options()
      ("t,threads", "number of threads given to Ashera", cxxopts::value<std::uint32_t>())
      ("w,window", "window size used for correction", cxxopts::value<std::uint32_t>())
      ("reads", "input reads", cxxopts::value<std::vector<std::string>>());
    /* clang-format on */

    options.parse_positional({"reads"});
    auto const result = options.parse(argc, argv);

    auto const n_threads = result["threads"].as<std::uint32_t>();
    ashera::InitThreadPool(n_threads);

    auto const win_size = result["window"].as<std::uint32_t>();

    auto engine = ashera::Engine(win_size);

    auto timer = biosoup::Timer();
    timer.Start();

    auto const reads = result["reads"].as<std::vector<std::string>>();
    auto sequences = ashera::LoadReads(reads);

    std::cerr << "[ashera](" << std::setw(12) << std::setprecision(3)
              << timer.Stop() << "s) : loaded " << sequences.size()
              << " sequences" << std::endl;

    timer.Start();
    auto ans = engine.Correct(std::move(sequences));

    std::cerr << "[ashera](" << std::setw(12) << std::setprecision(3)
              << timer.Stop() << "s) : generated corrected reads " << std::endl;

  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}
