import numpy as np
from .utils import reverse_complement
import networkx as nx
from collections import Counter
from itertools import chain


def construct_debruijn(reads, trim_depth=500, ksize=61, cutoff=10):

    kmer_counter = Counter()

    for read in reads:
        revcomp = reverse_complement(read)
        kmer_counter.update(kmer for kmer in (read[i:(i + ksize)] for i in range(len(read) + 1 - ksize)))
        kmer_counter.update(kmer for kmer in (revcomp[i:(i + ksize)] for i in range(len(revcomp) + 1 - ksize)))

    DG = nx.DiGraph()
    DG.add_weighted_edges_from((kmer[:-1], kmer[1:], w) for kmer, w in kmer_counter.items() if (w > cutoff))
    connected_components = list(nx.weakly_connected_components(DG))

    #assert ((len(connected_components) % 2) == 0)
    assert (nx.is_directed_acyclic_graph(DG))

    matched_components = []
    for i, component_1 in enumerate(connected_components[:-1]):
        reverse_string = ''.join(sorted(reverse_complement(kmer) for kmer in component_1))
        for j, component_2 in enumerate(connected_components[i + 1:]):
            idx = i + j + 1
            forward_string = ''.join(sorted(component_2))
            if (forward_string == reverse_string):
                matched_components.append((i, idx))

    duplicate_nodes = chain(*[connected_components[matched[0]] for matched in matched_components])
    DG.remove_nodes_from(duplicate_nodes)
    #de_bruijn_graph.remove_nodes_from(connected_components[0])

    DG_trimmed = trim_graph(DG.copy(), trim_depth=trim_depth)
    DG.clear()
    de_bruijn_graph, attrs = collapse_chains(DG_trimmed.copy())
    DG_trimmed.clear()
    return de_bruijn_graph


def traverse_graph(DG, node, depth, forward=True):
    if depth == 0:
        return []
    else:
        nodes = [node]
        if forward:
            for next_node in DG.succ[node]:
                nodes += traverse_graph(DG, next_node, depth-1, forward=True)
        else:
            for next_node in DG.pred[node]:
                nodes += traverse_graph(DG, next_node, depth-1, forward=False)
        return nodes


def trim_graph(DG, trim_depth=1000):

    roots = [node for node in DG.nodes if ((len(DG.pred[node]) == 0))]
    leaves = [node for node in DG.nodes if ((len(DG.succ[node]) == 0))]

    excess_nodes = []
    for root in roots:
        excess_nodes.extend(traverse_graph(DG, root, trim_depth, forward=True))
    for leaf in leaves:
        excess_nodes.extend(traverse_graph(DG, leaf, trim_depth, forward=False))

    DG.remove_nodes_from(excess_nodes)
    return DG


def decompose_graph(DG, ordered_antichains):
    graph_decompositions = []

    for i in range(len(ordered_antichains) - 1):
        roots = ordered_antichains[i]
        leaves = ordered_antichains[i + 1]
        descendants = set()
        descendants = descendants.union(*[nx.algorithms.dag.descendants(DG, n) for n in roots])
        ancestors = set()
        ancestors = ancestors.union(*[nx.algorithms.dag.ancestors(DG, n) for n in leaves])
        subgraph_nodes = list(ancestors.intersection(descendants))
        subgraph_nodes.extend(roots + leaves)
        graph_decompositions.append((subgraph_nodes, roots, leaves))

    return graph_decompositions


def is_well_connected(DG, nodes):
    descendants = set()
    descendants = descendants.union(*[nx.algorithms.dag.descendants(DG, n) for n in nodes])
    ancestors = set()
    ancestors = ancestors.union(*[nx.algorithms.dag.ancestors(DG, n) for n in nodes])
    related_nodes = list(ancestors.union(descendants))

    return ((len(related_nodes) + len(nodes)) == len(DG.nodes))


def unresolvable_node(DG, node):
    return ((len(node) > 100) and (len(DG.pred[node]) > 0) and (len(DG.succ[node]) > 0) and
            ((len(DG.pred[node]) > 1) or (len(DG.succ[node]) > 1)))


def resolvable_node(DG, node):
    return ((len(node) <= 100) and (len(DG.pred[node]) > 0) and (len(DG.succ[node]) > 0) and
            ((len(DG.pred[node]) > 1) or (len(DG.succ[node]) > 1)))


def is_small_chain(DG, node):
    return (len(DG.succ[node]) == 1) and (len(DG.pred[node]) >= 2)


def is_start(node, chain):
    return (len(chain.pred[node]) == 0) and (len(chain.succ[node]) == 1)


def is_end(node, chain):
    return (len(chain.pred[node]) == 1 and (len(chain.succ[node]) == 0))


def collapse_chains(g):
    to_delete = []
    is_chain = [node for node in g.nodes if ((len(g.succ[node]) == 1) and (len(g.pred[node]) == 1))]
    chains = g.subgraph(is_chain)

    components = list(nx.components.weakly_connected_components(chains))
    components = [chains.subgraph(component) for component in components]
    attrs = dict()
    for component in components:
        nodes = list(component.nodes)

        start_nodes = list(filter(lambda n: is_start(n, component), nodes))
        end_nodes = list(filter(lambda n: is_end(n, component), nodes))
        weights = [g[u][v]['weight'] for (u, v) in component.edges]

        num_start_end = (len(start_nodes), len(end_nodes))
        if (num_start_end == (1, 1)):
            n1, n2 = start_nodes[0], end_nodes[0]
            super_node = ''
            node = n1
            for i in range(len(nodes) - 1):
                super_node = super_node + node[0]
                node = list(component.succ[node])[0]
            super_node = super_node + n2
            g.add_edge(super_node, list(g.succ[n2])[0], weight=np.mean(weights), std=np.std(weights))
            g.add_edge(list(g.pred[n1])[0], super_node, weight=np.mean(weights), std=np.std(weights))

            to_delete.extend(nodes)

    g.remove_nodes_from(to_delete)

    return g, attrs