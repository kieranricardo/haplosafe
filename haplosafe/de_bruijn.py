import networkx as nx
from collections import Counter


def build_de_bruijn(reads, ksize=61, cutoff=10, collapse=True):

    kmer_counter = count_kmers(reads, ksize=ksize)

    DG = nx.DiGraph()
    edges = ((kmer[:-1], kmer[1:], {"weight": w, "kmer": kmer}) for kmer, w in kmer_counter.items())
    DG.add_edges_from(edge for edge in edges if (edge[2]['weight'] > cutoff))

    connected_components = list(nx.weakly_connected_components(DG))
    DG = DG.subgraph(max(connected_components)).copy()

    del kmer_counter

    # assert (len(connected_components) == 2) or (len(connected_components) == 1)

    assert (nx.is_directed_acyclic_graph(DG))

    if len(connected_components) == 2:
        DG.remove_nodes_from(connected_components[0])

    if collapse:
        raise NotImplementedError
        # de_bruijn_graph = collapse_chains(DG.copy())
        DG.clear()
    else:
        de_bruijn_graph = DG

    return de_bruijn_graph


def count_kmers(reads, ksize=61):

    reads = list(reads)

    kmer_counter = Counter()
    buff_size = 10000
    idxs = list(range(0, len(reads) + buff_size, buff_size))

    for idx_start, idx_end in zip(idxs[:-1], idxs[1:]):

        kmers = []
        for read in reads[idx_start:idx_end]:
            kmers.extend(kmer for kmer in (read[i:(i + ksize)] for
                                           i in range(len(read) + 1 - ksize)))
        kmer_counter.update(kmers)

    return kmer_counter


# def collapse_chains(g, ksize):
#     to_delete = []
#     is_chain = [node for node in g.nodes if ((len(g.succ[node]) == 1) and (len(g.pred[node]) == 1))]
#     chains = g.subgraph(is_chain)
#
#     components = list(nx.components.weakly_connected_components(chains))
#     components = [chains.subgraph(component) for component in components]
#     for component in components:
#
#         nodes = list(component.nodes)
#         ordered_nodes = list(filter(lambda n: is_start(n, component), nodes))[0]
#
#         while len(component.succ[ordered_nodes[-1]]) == 1:
#             ordered_nodes.append(list(component.succ[ordered_nodes[-1]][0]))
#
#         substring = add_path(ordered_nodes, ksize, component)
#         g.add_edge((ordered_nodes[0], ordered_nodes[-1], {"weight": 0.0, "kmer": substring}))
#
#         to_delete.extend(ordered_nodes[1:-1])
#
#     g.remove_nodes_from(to_delete)
#
#     return g


def get_leaves(DG):
    leaves = [node for node in DG.nodes if
              (len(list(DG.succ[node])) == 0)]
    leaves.sort()
    return leaves


def get_roots(DG):
    roots = [node for node in DG.nodes if
             (len(list(DG.pred[node])) == 0)]
    roots.sort()
    return roots