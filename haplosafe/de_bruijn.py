import networkx as nx
from collections import Counter


def build_de_bruijn(reads, ksize=61, cutoff=10):

    kmer_counter = count_kmers(reads, ksize=ksize)

    DG = nx.MultiDiGraph()
    edges = ((kmer[:-1], kmer[1:], {"weight": w, "kmer": kmer}) for kmer, w in kmer_counter.items())
    DG.add_edges_from(edge for edge in edges if (edge[2]['weight'] > cutoff))

    # connected_components = list(nx.weakly_connected_components(DG))
    # DG = DG.subgraph(max(connected_components)).copy()
    #
    # del kmer_counter
    #
    # # assert (len(connected_components) == 2) or (len(connected_components) == 1)
    #
    # assert (nx.is_directed_acyclic_graph(DG))
    #
    # if len(connected_components) == 2:
    #     DG.remove_nodes_from(connected_components[0])
    #
    # if collapse:
    #     raise NotImplementedError
    #     # de_bruijn_graph = collapse_chains(DG.copy())
    #     DG.clear()
    # else:
    #     de_bruijn_graph = nx.MultiDiGraph()
    #     de_bruijn_graph.add_edges_from(DG.edges.items())

    return DG


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