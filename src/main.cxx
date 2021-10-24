#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ashera/io.hpp"
#include "biosoup/timer.hpp"
#include "cxxopts.hpp"

int main(int argc, char** argv) {
  try {
    auto options =
        cxxopts::Options("ashera", "Ashera is a read correction tool");
    options.add_options()("reads", "input reads",
                          cxxopts::value<std::vector<std::string>>());

    options.parse_positional({"reads"});
    auto result = options.parse(argc, argv);

    auto timer = biosoup::Timer();

    timer.Start();

    auto const reads = result["reads"].as<std::vector<std::string>>();
    auto const sequences = ashera::LoadReads(reads);

    std::cerr << "[ashera](" << std::setw(12) << std::setprecision(3)
              << timer.Stop() << "s) : loaded " << sequences.size()
              << " sequences" << std::endl;

  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}
