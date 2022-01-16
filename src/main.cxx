#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "ashera/algorithm.hpp"
#include "ashera/configs.hpp"
#include "ashera/io.hpp"
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
      ("reads", "input reads", cxxopts::value<std::vector<std::string>>());
    /* clang-format on */

    options.parse_positional({"reads"});
    auto const result = options.parse(argc, argv);

    auto const n_threads = result["threads"].as<std::uint32_t>();

    auto const win_size = result["window"].as<std::uint32_t>();
    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(n_threads);

    auto timer = biosoup::Timer();
    timer.Start();

    auto const reads_path = result["reads"].as<std::vector<std::string>>();
    auto raw_reads = ashera::LoadReads(reads_path);

    fmt::print(stderr,
               FMT_COMPILE("[ashera]({:12.3f}s) : loaded {} sequences\n"),
               timer.Stop(), raw_reads.size());

    timer.Start();
    auto solid_overlaps = ashera::FindSnpFreeOverlaps(
        thread_pool, ashera::RamConfig(), raw_reads);

    fmt::print(stderr,
               FMT_COMPILE("[ashera]({:12.3f}s) : generated solid"
                           "overlaps\n"),
               timer.Stop());

    auto polish_cfg = ashera::PolishConfig();
    polish_cfg.window_len = win_size;

    timer.Start();
    auto const polished_reads =
        ashera::PolishReads(thread_pool, polish_cfg, std::move(raw_reads),
                            std::move(solid_overlaps));

    fmt::print(stderr,
               FMT_COMPILE("[ashera]({:12.3f}s) : polished {} sequences\n"),
               timer.Stop(), polished_reads.size());

    for (auto const& it : polished_reads) {
      fmt::print(stdout, FMT_COMPILE(">{}\n{}\n"), it->name, it->InflateData());
    }

  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}
