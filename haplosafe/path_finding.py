import networkx as nx
import numpy as np
from scipy.optimize import nnls
from .de_bruijn import get_leaves, get_roots


def add_path(path, ksize, DG):
    haplo = ''
    for edge in zip(path[:-1], path[1:]):
        contig = DG.edges[edge]['kmer']
        if haplo == '':
            haplo = contig
        else:
            idx = haplo.index(contig[:ksize-1])
            haplo = haplo[:idx] + contig

    return haplo


def find_source_sink_paths(DAG, start_nodes=None, max_len=np.inf):

    if start_nodes is None:
        start_nodes = get_roots(DAG)

    paths = [[node] for node in start_nodes]

    ii = 0

    while ii < len(paths):
        path = paths[ii]
        next_nodes = list(DAG.succ[path[-1]])

        while (len(next_nodes)) > 0 and (len(path) < max_len):
            paths.extend(path + [node, ] for node in next_nodes[1:])
            path.append(next_nodes[0])

            next_nodes = list(DAG.succ[path[-1]])

        ii += 1

    return paths


def get_subgraph(nodes, g):

    sg = g.subgraph(nodes)
    node_counter = 0

    # TODO: check if only one iteration is needed
    # TODO: check how the false path is being created --> tip??!?!?!?!?!!?!?!?!?
    while len(sg.nodes) > node_counter:

        leaves = [node for node in sg.nodes if (len(list(sg.succ[node])) == 0)]
        roots = [node for node in sg.nodes if (len(list(sg.pred[node])) == 0)]

        nodes_to_add = []
        for node in sg.nodes:
            if not node in leaves:
                nodes_to_add.extend(g.succ[node])
            if not node in roots:
                nodes_to_add.extend(g.pred[node])
                #print("backtracked!")

        node_counter = len(sg.nodes)
        sg = g.subgraph(list(sg.nodes) + nodes_to_add)

    return sg


def get_path_matrix(g):

    all_paths = find_source_sink_paths(g)

    edge_list = list(g.edges)
    edge_list_lookup = dict(zip(edge_list, range(len(edge_list))))

    weights = np.array(list(([e[1]['weight'] for e in g.edges.items()])))
    A = np.zeros((len(edge_list), len(all_paths)))

    for i in range(len(all_paths)):
        path_idxs = [edge_list_lookup[(n1, n2)] for n1, n2 in
                     zip(all_paths[i][:-1], all_paths[i][1:])]
        A[path_idxs, i] = 1

    return A, weights, all_paths


# code to iterate this procedure


def solve_path_matrix(A, weights, all_paths, cutoff, g, ksize):

    fnnl, res = nnls(A, weights)

    predicted_paths = [all_paths[i] for i in np.where(fnnl > cutoff)[0]]
    predicted_haplotypes = [add_path(path, ksize, g) for path in predicted_paths]

    pred_freqs = fnnl[fnnl > cutoff]

    f_sum = pred_freqs.sum()
    pred_freqs = pred_freqs / f_sum

    return predicted_paths, predicted_haplotypes, pred_freqs, f_sum


def merge_bubbles(g, cutoff, ksize, window_size=50):

    start_nodes = get_roots(g)
    leaves = get_leaves(g)

    new_edges = []
    freq_sums = []

    while (start_nodes != leaves):

        paths = find_source_sink_paths(
            g, max_len=window_size, start_nodes=start_nodes
        )

        subgraph_nodes = sum(paths, [])
        subgraph = get_subgraph(subgraph_nodes, g)

        A, weights, all_paths = get_path_matrix(subgraph)

        # note g can also be subgraph
        predicted_paths, predicted_haplotypes, pred_freqs, f_sum = solve_path_matrix(
            A, weights, all_paths, cutoff, g, ksize=ksize
        )

        new_edges.extend(
            (path[0], path[-1], {"weight": freq, "kmer": h})
            for path, freq, h in
            zip(predicted_paths, pred_freqs, predicted_haplotypes)
        )

        freq_sums.append(f_sum)
        # new roots are old leaves
        start_nodes = get_leaves(subgraph)

    new_graph =  nx.DiGraph() # nx.MultiDiGraph()

    new_graph.add_edges_from(new_edges)
    new_cutoff = cutoff / np.mean(freq_sums)

    return new_graph, new_cutoff, new_edges
