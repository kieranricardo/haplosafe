import networkx as nx
from collections import Counter


def build_de_bruijn(reads, ksize=61, cutoff=10):

    kmer_counter = count_kmers(reads, ksize=ksize)

    graph = nx.MultiDiGraph()
    edges = (
        (kmer[:-1], kmer[1:], {"weight": w, "kmer": kmer})
        for kmer, w in kmer_counter.items()
    )
    graph.add_edges_from(edge for edge in edges if (edge[2]["weight"] > cutoff))
    components = list(nx.weakly_connected_components(graph))
    giant_component = components[np.argmax([len(c) for c in components])]
    graph = graph.subgraph(giant_component).copy()

    return graph


def count_kmers(reads, ksize=61):

    kmer_generator = (read[i: (i+ksize)] for read in reads for i in range(len(read)+1-ksize))
    kmer_counter = Counter(kmer_generator)
    return kmer_counter


def get_leaves(DG):
    leaves = [node for node in DG.nodes if (len(list(DG.succ[node])) == 0)]
    leaves.sort()
    return leaves


def get_roots(DG):
    roots = [node for node in DG.nodes if (len(list(DG.pred[node])) == 0)]
    roots.sort()
    return roots


def fast_trim(graph, forward=True):

    if forward:
        top_sorted_nodes = list(nx.topological_sort(graph))
        start_nodes = get_roots(graph)
        next_nodes_lookup = graph.succ
        prev_nodes_lookup = graph.pred
    else:
        top_sorted_nodes = list(nx.topological_sort(graph))[-1::-1]
        start_nodes = get_leaves(graph)
        next_nodes_lookup = graph.pred
        prev_nodes_lookup = graph.succ

    max_node_path_len = []

    for start_node in start_nodes:
        node_dists = dict((node, 0) for node in graph.nodes)
        node_dists[start_node] = 1
        idx = top_sorted_nodes.index(start_node)

        for node in top_sorted_nodes[idx:]:
            if node_dists[node] > 0:
                for next_node in next_nodes_lookup[node]:
                    node_dists[next_node] = max(
                        node_dists[next_node], node_dists[node] + 1
                    )

        max_node_path_len.append(max(node_dists.values()) - 1)

    max_path_len = max(max_node_path_len)

    bad_root_chains = [
        [node]
        for node, mpl in zip(start_nodes, max_node_path_len)
        if mpl < (0.8 * max_path_len)
    ]

    for chain in bad_root_chains:
        while len(next_nodes_lookup[chain[-1]]) == 1:
            next_node = list(next_nodes_lookup[chain[-1]])[0]
            if len(prev_nodes_lookup[next_node]) > 1:
                break
            else:
                chain.append(next_node)
        graph.remove_nodes_from(chain)

    return True
