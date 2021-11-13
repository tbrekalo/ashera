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
#include "fmt/compile.h"
#include "fmt/core.h"

int main(int argc, char** argv) {
  try {
    auto options =
        cxxopts::Options("ashera", "Ashera is a read correction tool");

    /* clang-format off */
    options.add_options()
      ("t,threads", "number of threads given to Ashera", cxxopts::value<std::uint32_t>())
      ("w,window", "window size used for correction", cxxopts::value<std::uint32_t>())
      ("l,len", "percentile of longest reads for backbone sequences", 
        cxxopts::value<float>())
      ("reads", "input reads", cxxopts::value<std::vector<std::string>>());
    /* clang-format on */

    options.parse_positional({"reads"});
    auto const result = options.parse(argc, argv);

    auto const n_threads = result["threads"].as<std::uint32_t>();
    ashera::InitThreadPool(n_threads);

    auto const win_size = result["window"].as<std::uint32_t>();
    auto const len_percentile = result["len"].as<float>();

    auto engine = ashera::Engine(len_percentile, win_size);

    auto timer = biosoup::Timer();
    timer.Start();

    auto const reads = result["reads"].as<std::vector<std::string>>();
    auto sequences = ashera::LoadReads(reads);

    fmt::print(stderr,
               FMT_COMPILE("[ashera]({:12.3f}s) : loaded {} sequences\n"),
               timer.Stop(), sequences.size());

    timer.Start();
    auto ans = engine.Correct(std::move(sequences));

    fmt::print(stderr,
               FMT_COMPILE("[ashera]({:12.3f}s) : generated corrected reads\n"),
               timer.Stop());

    for (auto const& it : ans) {
      fmt::print(stdout, FMT_COMPILE(">{}\n{}\n"), it->name, it->InflateData());
    }

  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}
