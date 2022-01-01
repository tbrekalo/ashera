#include "polish.hpp"

#include <cmath>
#include <string>

#include "edlib.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/spoa.hpp"

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

  auto windows = std::vector<spoa::Graph>();
  windows.resize(static_cast<std::size_t>(std::ceil(
      static_cast<double>(target->inflated_len) / config.window_len)));

  auto window_ends = std::vector<std::uint32_t>(windows.size());
  for (auto window_id = 0; window_id < windows.size(); ++window_id) {
    auto const window_begin = window_id * config.window_len;
    auto const window_end =
        std::min(window_begin + config.window_len, target->inflated_len);

    auto const backbone_str =
        target_str.substr(window_begin, window_end - window_begin);

    windows[window_id].AddAlignment(spoa::Alignment(), backbone_str, 0U);
    window_ends.push_back(window_end);
  }

  std::sort(overlaps.cbegin(), overlaps.cend(),
            [](biosoup::Overlap const& a, biosoup::Overlap const& b) -> bool {
              return a.lhs_begin < b.lhs_begin;
            });

  auto const align_to_window = [&alignment_engine, &references, &windows,
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
      alignment_engine->Align(ref_substring, windows[window_idx]);
    } else {
      auto mapping = std::vector<spoa::Graph::Node const*>();
      auto subgraph = windows[window_idx].Subgraph(
          target_local_intv.begin_idx, target_local_intv.end_idx - 1, &mapping);
      alignment = alignment_engine->Align(ref_substring, subgraph);
      subgraph.UpdateAlignment(mapping, &alignment);
    }
  };

  for (auto const& ovlp : overlaps) {
    auto const& ref = references[ovlp.rhs_id];

    auto const ref_str = ref->InflateData();
    auto edlib_align_res = edlibAlign(
        target_str.c_str(), target_str.size(), ref_str.c_str(), ref_str.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

    auto target_pos = ovlp.lhs_id - 1U;
    auto reference_pos = ovlp.rhs_id - 1U;

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
          active_reference_interval.begin_idx = reference_pos;
          if (active_target_interval.end_idx == window_ends[window_idx]) {
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
          ++reference_pos;
          break;
        }

        case 2: {  // insertion on the target
          ++target_pos;

          // TODO: target logic

          break;
        }

        default: {
          break;
        }
      }
    }
  }

  // TODO: generate consensus

  return std::unique_ptr<biosoup::NucleicAcid>();
}

}  // namespace ashera::detail
