#include "assembly.hpp"

#include <algorithm>
#include <iterator>
#include <unordered_set>

#include "fmt/compile.h"
#include "fmt/core.h"

namespace ashera::detail {

struct Node;
struct Edge;

std::atomic_uint32_t node_cnt;
std::atomic_uint32_t edge_cnt;
std::atomic_uint32_t bb_seq_cnt;

auto ResetObjCounters() -> void {
  node_cnt = 0;
  edge_cnt = 0;
  bb_seq_cnt = 0;
}

struct Node {
  Node(biosoup::NucleicAcid seq) : id(node_cnt++), seq(std::move(seq)) {}

  std::uint32_t const id;
  biosoup::NucleicAcid seq;

  Edge* in_edge = nullptr;
  Edge* out_edge = nullptr;
};

struct NodePairing {
  std::unique_ptr<Node> original;
  std::unique_ptr<Node> rev_comp;
};

struct Edge {
  Edge(Node* tail_ptr, Node* head_ptr, std::uint32_t length)
      : tail_ptr(tail_ptr), head_ptr(head_ptr), length(length) {}

  auto Data() const -> std::string { return tail_ptr->seq.InflateData(length); }
  auto Quality() const -> std::string {
    return tail_ptr->seq.InflateQuality(length);
  }

  Node* tail_ptr = nullptr;
  Node* head_ptr = nullptr;
  std::uint32_t length;
};

auto CreateNodePairing(biosoup::NucleicAcid seq) -> NodePairing {
  auto dst = NodePairing();

  seq.ReverseAndComplement();
  dst.rev_comp = std::make_unique<Node>(seq);

  seq.ReverseAndComplement();
  dst.original = std::make_unique<Node>(std::move(seq));

  return dst;
}

auto NodeChainToSeq(Node const* begin, Node const* end)
    -> std::unique_ptr<biosoup::NucleicAcid> {
  auto data = std::string();
  auto quality = std::string();

  auto curr_node = begin;
  auto const traverse = [&]() -> void {
    data += curr_node->out_edge->Data();
    quality += curr_node->out_edge->Quality();

    curr_node = curr_node->out_edge->head_ptr;
  };

  traverse();
  while (true) {
    if (curr_node != end) {
      traverse();
    } else {
      if (curr_node != begin) {
        data += curr_node->seq.InflateData();
        quality += curr_node->seq.InflateQuality();
      }

      break;
    }
  }

  auto name = fmt::format(FMT_COMPILE("BbSeq_{:7d}"), bb_seq_cnt++);

  return std::make_unique<biosoup::NucleicAcid>(
      std::move(name), std::move(data), std::move(quality));
}

auto SeqsToNodes(std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs)
    -> std::vector<NodePairing> {
  auto dst = std::vector<NodePairing>();
  dst.reserve(seqs.size());

  std::transform(
      seqs.cbegin(), seqs.cend(), std::back_inserter(dst),
      [](std::unique_ptr<biosoup::NucleicAcid> const& seq) -> NodePairing {
        return CreateNodePairing(*seq);
      });

  return dst;
}

/**
 * @brief assemble backbone sequences for consensus tool
 */
auto AssembleBackbones(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::vector<EdgeCandidate> const& edge_candidates)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  ResetObjCounters();

  auto nodes = SeqsToNodes(seqs);
  auto edges = std::vector<std::unique_ptr<Edge>>();

  auto const can_connect = [](Node* tail_ptr, Node* head_ptr) -> bool {
    return tail_ptr->out_edge == nullptr && head_ptr->in_edge == nullptr;
  };

  auto const connect_via = [](Node* tail_ptr, Node* head_ptr,
                              std::uint32_t edge_len) -> std::unique_ptr<Edge> {
    auto dst = std::make_unique<Edge>(tail_ptr, head_ptr, edge_len);
    tail_ptr->out_edge = dst.get();
    head_ptr->in_edge = dst.get();

    return dst;
  };

  // TODO: thing about adding reverse complement in both cases
  for (auto const& ec : edge_candidates) {
    auto const& ovlp = ec.ovlp;
    auto const ovlp_type = DetermineOverlapType(
        ovlp, seqs[ovlp.lhs_id]->inflated_len, seqs[ovlp.rhs_id]->inflated_len);

    auto& lhs_node = nodes[ovlp.lhs_id].original;
    auto& rhs_node =
        ovlp.strand ? nodes[ovlp.rhs_id].original : nodes[ovlp.rhs_id].rev_comp;

    if (ovlp_type == OverlapType::kLhsToRhs) {
      if (can_connect(lhs_node.get(), rhs_node.get())) {
        auto const edge_len = ovlp.lhs_begin - ovlp.rhs_begin;
        edges.emplace_back(
            connect_via(lhs_node.get(), rhs_node.get(), edge_len));
      }
    } else {  // OverlapType::kRhsToLhs
      if (can_connect(rhs_node.get(), lhs_node.get())) {
        auto const edge_len = ovlp.rhs_begin - ovlp.lhs_begin;
        edges.emplace_back(
            connect_via(rhs_node.get(), lhs_node.get(), edge_len));
      }
    }
  }

  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto visited = std::vector<std::uint8_t>(node_cnt, 0);

  for (auto const& e : edges) {
    auto begin = e->tail_ptr;
    auto end = e->head_ptr;

    auto chain_len = 2UL;

    if (visited[begin->id] || visited[end->id]) {
      continue;
    }

    visited[begin->id] = 1;
    visited[end->id] = 1;

    // expand left
    while (begin->in_edge != nullptr) {
      visited[begin->id] = 1;
      auto const new_begin = begin->in_edge->tail_ptr;
      if (new_begin == nullptr || visited[new_begin->id]) {
        break;
      } else {
        begin = new_begin;
        ++chain_len;
      }
    }

    // expand right
    while (end->out_edge != nullptr) {
      visited[end->id] = 1;
      auto const new_end = end->out_edge->head_ptr;
      if (new_end == nullptr || visited[new_end->id]) {
        break;
      } else {
        end = new_end;
        ++chain_len;
      }
    }

    if (chain_len >= 5) {
      dst.emplace_back(NodeChainToSeq(begin, end));
      dst.back()->id = dst.size() - 1U;
    }
  }

  return dst;
}

}  // namespace ashera::detail
