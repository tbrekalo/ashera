#include "polish.hpp"

#include <cmath>
#include <string>

#include "edlib.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/spoa.hpp"

/* clang-format off */
#include "fmt/core.h"  // TODO: remove
#include "fmt/compile.h"
#include "fmt/ranges.h"
/* clang-format on */

namespace ashera::detail {

auto constexpr kAllowedFuzzPercent = 0.01;
auto constexpr kSmallWindowPercent = 0.02;

struct Interval {
  std::uint32_t begin_idx;
  std::uint32_t end_idx;
};

auto NormalizeToWindow(Interval const intv, std::uint32_t const window_idx,
                       std::uint32_t const window_len) -> Interval {
  return {.begin_idx = intv.begin_idx - window_idx * window_len,
          .end_idx = intv.end_idx - window_idx * window_len};
}

auto IntervalLength(Interval const intv) -> std::uint32_t {
  return intv.end_idx - intv.begin_idx;
}

[[nodiscard]] auto Polish(
    PolishConfig config, std::unique_ptr<biosoup::NucleicAcid> target,
    std::vector<biosoup::Overlap> overlaps,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> references)
    -> std::unique_ptr<biosoup::NucleicAcid> {
  /* clang-format off */
  auto alignment_engine = 
    spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 
      config.spoa_config.match, 
      config.spoa_config.mismatch,
      config.spoa_config.gap);
  /* clang-format on */

  auto target_str = target->InflateData();

  auto window_graphs = std::vector<spoa::Graph>();
  window_graphs.resize(static_cast<std::size_t>(std::ceil(
      static_cast<double>(target->inflated_len) / config.window_len)));

  auto window_ends = std::vector<std::uint32_t>();
  window_ends.reserve(window_graphs.size());

  for (auto window_id = 0; window_id < window_graphs.size(); ++window_id) {
    auto const window_begin = window_id * config.window_len;
    auto const window_end =
        std::min(window_begin + config.window_len, target->inflated_len);

    auto const backbone_str =
        target_str.substr(window_begin, window_end - window_begin);

    window_graphs[window_id].AddAlignment(spoa::Alignment(), backbone_str, 0U);
    window_ends.push_back(window_end);
  }

  std::sort(overlaps.begin(), overlaps.end(),
            [](biosoup::Overlap const& a, biosoup::Overlap const& b) -> bool {
              return a.lhs_begin < b.lhs_begin;
            });

  auto const align_to_window = [&alignment_engine, &references, &window_graphs,
                                window_len = config.window_len](
                                   std::uint32_t const window_idx,
                                   std::uint32_t const reference_idx,
                                   Interval const target_interval,
                                   Interval const reference_interval) -> void {
    if (IntervalLength(reference_interval) < window_len * kSmallWindowPercent) {
      return;
    }

    auto const target_local_intv =
        NormalizeToWindow(target_interval, window_idx, window_len);

    auto const kBeginWhitin = window_len * kAllowedFuzzPercent;
    auto const kEndWithin = window_len - kBeginWhitin;

    auto const ref_substring = references[reference_idx]->InflateData(
        reference_interval.begin_idx, IntervalLength(reference_interval));

    auto alignment = spoa::Alignment();
    if (target_local_intv.begin_idx < kBeginWhitin &&
        target_local_intv.end_idx > kEndWithin) {
      alignment =
          alignment_engine->Align(ref_substring, window_graphs[window_idx]);
    } else {
      auto mapping = std::vector<spoa::Graph::Node const*>();
      auto subgraph = window_graphs[window_idx].Subgraph(
          target_local_intv.begin_idx, target_local_intv.end_idx - 1, &mapping);
      alignment = alignment_engine->Align(ref_substring, subgraph);
      subgraph.UpdateAlignment(mapping, &alignment);
    }

    window_graphs[window_idx].AddAlignment(alignment, ref_substring);
  };

  auto const overlap_strings =
      [&target_str, &references](
          biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
    auto target_substr =
        target_str.substr(ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

    auto ref_substr = references[ovlp.rhs_id]->InflateData(
        ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

    return std::make_pair(std::move(target_substr), std::move(ref_substr));
  };

  for (auto const& ovlp : overlaps) {
    auto const [target_substr, ref_substr] = overlap_strings(ovlp);

    auto edlib_align_res = edlibAlign(
        target_substr.c_str(), target_substr.size(), ref_substr.c_str(),
        ref_substr.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

    auto target_pos = ovlp.lhs_begin - 1U;
    auto reference_pos = ovlp.rhs_begin - 1U;

    bool active_interval = false;
    auto active_target_interval = Interval();
    auto active_reference_interval = Interval();

    auto window_idx = ovlp.lhs_begin / config.window_len;

    for (auto i = 0U; i < edlib_align_res.alignmentLength; ++i) {
      switch (edlib_align_res.alignment[i]) {
        case 3:    // mismatch
        case 0: {  // match
          ++reference_pos;
          ++target_pos;

          if (!active_interval) {
            active_target_interval.begin_idx = target_pos;
            active_reference_interval.begin_idx = reference_pos;
            active_interval = true;
          }

          active_target_interval.end_idx = target_pos;
          active_reference_interval.end_idx = reference_pos;
          if (active_target_interval.end_idx + 1U == window_ends[window_idx]) {
            ++active_target_interval.end_idx;
            ++active_reference_interval.end_idx;

            align_to_window(window_idx, ovlp.rhs_id, active_target_interval,
                            active_reference_interval);

            active_interval = false;
            ++window_idx;
          }

          break;
        }

        case 1: {  // insertion on the polishing reference
          ++target_pos;

          active_target_interval.end_idx = target_pos;
          if (active_target_interval.end_idx + 1U == window_ends[window_idx]) {
            if (active_interval) {
              ++active_target_interval.end_idx;
              ++active_reference_interval.end_idx;

              align_to_window(window_idx, ovlp.rhs_id, active_target_interval,
                              active_reference_interval);

              active_interval = false;
            }

            ++window_idx;
          }

          break;
        }

        case 2: {  // insertion on the target
          ++reference_pos;
          break;
        }

        default: {
          break;
        }
      }
    }
  }

  auto polished_data = std::string();
  for (auto& it : window_graphs) {
    auto const consensus = it.GenerateConsensus();
    polished_data += it.GenerateConsensus();
  }

  return std::make_unique<biosoup::NucleicAcid>(target->name, polished_data);
}

}  // namespace ashera::detail
