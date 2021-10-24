#include <string>
#include <vector>

#include "ashera/io.hpp"
#include "catch2/catch_test_macros.hpp"

TEST_CASE("Load fasta", "[fasta][gz][io]") {
  using namespace std::string_literals;

  /* clang-format off */
  auto const correct_ans = std::vector<std::string>{
    "AACCAATTCCAGAACCAATTCCAGAACCAATTCCAGAACCAATTCCAGAACCAATTCCAGAACCAATTCCAG",
    "GACCTTAACCAAGACCTTAACCAAGACCTTAACCAAGACCTTAACCAAGACCTTAACCAAGACCTTAACCAA"
  };
  /* clang-format on */

  auto const file_path = "data/small.fasta.gz"s;
  auto reads = ashera::LoadReads(file_path);

  REQUIRE(reads.size() == 2UL);
  for (auto i = 0UL; i < correct_ans.size(); ++i) {
    CHECK(reads[i]->InflateData() == correct_ans[i]);
  }
}

TEST_CASE("Load fastq", "[fasta][gz][io]") {
  using namespace std::string_literals;

  /* clang-format off */
  auto const correct_ans_data = std::vector<std::string>{
    "GGGAAGGTGGTAACCAACCACAA", "GACCTTTACCAAGAACTTAACCAA"
  };
  
  // Takes the average
  auto const correct_ans_quality = std::vector<std::string>{
    "AAAAAAAAAAAAAAAAAAAAAAA", "------------------------"
  };
  /* clang-format on */

  auto const file_path = "data/small.fastq.gz"s;
  auto reads = ashera::LoadReads(file_path);

  REQUIRE(reads.size() == 2UL);
  for (auto i = 0UL; i < correct_ans_data.size(); ++i) {
    CHECK(reads[i]->InflateData() == correct_ans_data[i]);
    CHECK(reads[i]->InflateQuality() == correct_ans_quality[i]);
  }
}
