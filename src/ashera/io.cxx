#include "ashera/io.hpp"

#include <iterator>
#include <limits>
#include <stdexcept>
#include <string_view> 
#include <utility>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

std::atomic_uint32_t biosoup::NucleicAcid::num_objects{0U};

namespace ashera {

namespace detail {

auto constexpr kFastaExts =
    std::array<char const*, 4>{".fasta", ".fa", "fasta.gz", ".fa.gz"};
auto constexpr kFastqExts =
    std::array<char const*, 4>{".fastq", ".fq", "fastq.gz", ".fq.gz"};

auto IsSuffix(std::string_view path, std::string_view target_suffix) -> bool {
  if (path.size() < target_suffix.size()) {
    return false;
  } else {
    auto const path_suffix = path.substr(path.size() - target_suffix.length());
    return path_suffix == target_suffix;
  }
}

template <std::size_t N>
auto IsFileType(std::string_view path_str,
                std::array<char const*, N> const& exts) -> bool {
  for (auto const& ext : exts) {
    if (IsSuffix(path_str, ext)) {
      return true;
    }
  }

  return false;
}

auto CreateParser(std::string const& file_path)
    -> std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> {
  auto parser =
      bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(
          file_path);

  if (IsFileType(file_path, kFastaExts)) {
    return bioparser::Parser<biosoup::NucleicAcid>::Create<
        bioparser::FastaParser>(file_path);
  } else if (IsFileType(file_path, kFastqExts)) {
    return bioparser::Parser<biosoup::NucleicAcid>::Create<
        bioparser::FastqParser>(file_path);
  } else {
    throw std::invalid_argument(
        "[ashera::detail::CreateParser] Invalid file type");
  }
}
}  // namespace detail

auto LoadReads(std::string const& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto const parser = detail::CreateParser(reads);
  return parser->Parse(std::numeric_limits<std::uint64_t>::max());
}

auto LoadReads(std::vector<std::string> const& reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  for (auto const& read : reads) {
    auto seqs = LoadReads(read);
    std::move(seqs.begin(), seqs.end(), std::back_inserter(dst));
  }

  return dst;
}

}  // namespace ashera
