import networkx as nx
import numpy as np
from scipy.optimize import nnls
from .de_bruijn import get_leaves, get_roots


def add_path(path, ksize, DG):
    haplo = ""
    for edge in path:
        contig = DG.edges[edge]["kmer"]
        if haplo == "":
            haplo = contig
        else:
            idx = haplo.index(contig[: ksize - 1])
            haplo = haplo[:idx] + contig

    return haplo


def find_projected_paths(subgraph, graph):
    start_nodes = get_roots(subgraph)
    paths = [
        [(parent_node, child_node, edge_idx)]
        for parent_node in start_nodes
        for child_node in subgraph[parent_node]
        for edge_idx in subgraph[parent_node][child_node]
    ]

    ii = 0

    while ii < len(paths):
        path = paths[ii]

        while True:

            parent_node = path[-1][1]

            edges = [
                (parent_node, child_node, edge_idx)
                for child_node in subgraph[parent_node]
                for edge_idx in subgraph[parent_node][child_node]
            ]

            if len(edges) == 0:
                break
            elif len(subgraph[parent_node]) != len(graph[parent_node]):
                paths.extend(path + [edge, ] for edge in edges)
                break
            else:
                path.append(edges[0])
                paths.extend(path[:-1] + [edge, ] for edge in edges[1:])

        ii += 1

    return paths


def get_path_matrix(subgraph, graph):

    all_paths = find_projected_paths(subgraph, graph)

    edge_list = list(subgraph.edges)
    edge_list_lookup = dict(zip(edge_list, range(len(edge_list))))

    weights = np.array(list(([e[1]["weight"] for e in subgraph.edges.items()])))
    A = np.zeros((len(edge_list), len(all_paths)))

    for i in range(len(all_paths)):
        path_idxs = [edge_list_lookup[edge] for edge in all_paths[i]]
        A[path_idxs, i] = 1

    return A, weights, all_paths


def solve_path_matrix(A, weights, all_paths, cutoff, g, ksize):

    fnnl, res = nnls(A, weights)
    predicted_paths = [all_paths[i] for i in np.where(fnnl > cutoff)[0]]
    predicted_haplotypes = [add_path(path, ksize, g) for path in predicted_paths]

    pred_freqs = fnnl[fnnl > cutoff]

    f_sum = pred_freqs.sum()
    # pred_freqs = pred_freqs / f_sum

    return predicted_paths, predicted_haplotypes, pred_freqs, f_sum


def partition_graph(graph, window_size):
    partitions = []

    node_dists = max_dists(graph, forward=False)
    max_dist = max(dist for dist in node_dists.values())
    node_idx_pairs = [(max_dist - dist, node) for node, dist in node_dists.items()]
    idx_node_table = dict((idx, []) for idx in set(idx_ for idx_, _ in node_idx_pairs))
    for idx, node in node_idx_pairs:
        idx_node_table[idx].append(node)

    idx = 0
    max_idx = max(idx for idx, _ in node_idx_pairs)
    while idx <= max_idx:
        partition = [node for i in range(idx, idx + window_size) for
                     node in idx_node_table.get(i, [])]

        parents = set(parent for child in partition for parent in graph.pred[child])
        partition.extend(parents)
        partitions.append(partition)
        idx += window_size

    return partitions


def resolve_paths(partition, graph, cutoff, ksize):
    subgraph = graph.subgraph(partition)
    A, weights, all_paths = get_path_matrix(subgraph, graph)

    # note g can also be subgraph
    predicted_paths, predicted_haplotypes, pred_freqs, f_sum = solve_path_matrix(
        A, weights, all_paths, cutoff, subgraph, ksize=ksize
    )

    edges = [
        (path[0][0], path[-1][1], {"weight": freq, "kmer": h})
        for path, freq, h in zip(predicted_paths, pred_freqs, predicted_haplotypes)
    ]

    return edges, f_sum


def merge_bubbles(graph, cutoff, ksize, window_size=50):

    partitions = partition_graph(graph, window_size)
    solutions = [resolve_paths(partition, graph, cutoff, ksize) for partition in partitions]

    edges = [edge for edges, _ in solutions for edge in edges]
    freq_sums = [f_sum for _, f_sum in solutions]
    new_graph = nx.MultiDiGraph()

    new_graph.add_edges_from(edges)
    new_cutoff = cutoff / np.mean(freq_sums)

    return new_graph, cutoff


def max_dists(graph, forward=True):
    """
    Find the max distance from a node to either a source node or a sink node.

    Args:
        graph (nx.MultiDiGraph): the graph to traverse
        forward (bool): whether to go forward or backward through the graph
    Returns:
        node_dists (dict): a dictionary of (node, max_dist) pairs
    """

    if forward:
        next_nodes_lookup = graph.succ
        node_order = list(nx.topological_sort(graph))
        start_nodes = get_roots(graph)
    else:
        next_nodes_lookup = graph.pred
        node_order = list(nx.topological_sort(graph))[-1::-1]
        start_nodes = get_leaves(graph)

    node_dists = dict((node, 0) for node in graph.nodes)

    for start_node in start_nodes:

        node_dists[start_node] = 1
        idx = node_order.index(start_node)

        for node in node_order[idx:]:
            if node_dists[node] > 0:
                for next_node in next_nodes_lookup[node]:
                    node_dists[next_node] = max(
                        node_dists[next_node], node_dists[node] + 1
                    )

    return node_dists


# TODO: remove this
def find_source_sink_paths(DAG, start_nodes=None, max_len=np.inf):

    if start_nodes is None:
        start_nodes = get_roots(DAG)

    paths = [
        [(parent_node, child_node, edge_idx)]
        for parent_node in start_nodes
        for child_node in DAG[parent_node]
        for edge_idx in DAG[parent_node][child_node]
    ]

    ii = 0

    while ii < len(paths):
        path = paths[ii]

        while len(path) < max_len:

            parent_node = path[-1][1]

            edges = [
                (parent_node, child_node, edge_idx)
                for child_node in DAG[parent_node]
                for edge_idx in DAG[parent_node][child_node]
            ]

            if len(edges) == 0:
                break
            else:
                path.append(edges[0])
                paths.extend(path[:-1] + [edge,] for edge in edges[1:])

        ii += 1

    return paths

