import networkx as nx
import numpy as np
from itertools import product


def get_antichains(DG):
    antichains = list(nx.algorithms.dag.antichains(DG))
    max_ac_card = max(len(ac) for ac in antichains)

    decomposition_antichains = [ac for ac in antichains if (len(ac) == max_ac_card)]

    ordered_antichains = []
    antichains_copy = decomposition_antichains.copy()

    while (len(antichains_copy) > 0):
        ordered_antichains.append(antichains_copy[-1])
        antichains_copy = [ac for ac in antichains_copy
                           if not any(a in ordered_antichains[-1] for a in ac)]

    return ordered_antichains, max_ac_card


def estimate_frequencies(DG, ordered_antichains):
    max_ac_counts = []
    for antichain in ordered_antichains:
        try:
            max_ac_counts.append(sorted(np.round(DG.edges[(n, list(DG.succ[n])[0])]['weight']) for n in antichain))
        except IndexError:
            max_ac_counts.append(sorted(np.round(DG.edges[(list(DG.pred[n])[0], n)]['weight']) for n in antichain))

    max_ac_counts = np.vstack(max_ac_counts)

    haplo_freqs = max_ac_counts.mean(0)

    return haplo_freqs


def fit_haplo(sub_graph, subgraph_haplotypes, haplo_freqs, m_ac):
    edge_list = list(sub_graph.edges)
    weights = np.array(list(([e[1]['weight'] for e in sub_graph.edges.items()])))

    haplo_vecs = [[] for x in subgraph_haplotypes]
    paths = [[] for i in range(m_ac)]
    for i, haplotype in enumerate(subgraph_haplotypes):
        for path in list(product(*haplotype)):
            path_nodes = []
            for sub_path_nodes in path:
                path_nodes += sub_path_nodes[:-1]
            path_nodes.append(sub_path_nodes[-1])
            paths[i].append(path_nodes)

    for i, haplotype_paths in enumerate(paths):
        for path in haplotype_paths:
            haplo_vec = np.zeros(weights.shape)
            for j, node in enumerate(path[:-1]):
                try:
                    idx = edge_list.index((node, path[j + 1]))
                    haplo_vec[idx] = 1
                except ValueError:
                    pass
            haplo_vecs[i].append(haplo_vec)

    combos = list(product(*[range(len(haplo_vec)) for haplo_vec in haplo_vecs]))
    As = [np.vstack([haplo_vecs[i][combo[i]] for i in range(m_ac)]) for combo in combos]
    combos = [combo for combo, A in zip(combos, As) if (A.sum(0).all())]
    As = [A for A in As if (A.sum(0).all())]
    As = np.dstack(As).transpose()
    errors = (((np.dot(As, haplo_freqs) - weights) ** 2) / weights).sum(-1)
    combo = combos[np.argmin(errors)]
    return [paths[i][combo[i]] for i in range(m_ac)]


def label_nodes(subgraph, DG, roots, leaves, haplo_freqs, max_ac_card):
    labelled_nodes = [[] for i in range(max_ac_card)]

    antichains = [ac for ac in nx.algorithms.dag.antichains(subgraph) if (len(ac) == (max_ac_card - 1))]

    for antichain in antichains:

        weight_sums = np.array(list(np.round(sum(DG.edges[(n, n_)]['weight'] for n_ in DG.succ[n])) for n in antichain))
        idxs = np.argsort(weight_sums)

        missing_weights = np.array([])
        missing_idxs = []
        is_good_fit = []
        for i, weight in enumerate(haplo_freqs):
            delta_w = np.abs(weight_sums - weight)
            is_close = (delta_w < (1.5 * np.sqrt(weight))) * (delta_w < (1.5 * np.sqrt(weight_sums)))
            is_good_fit.append(is_close)

            if (is_close.sum() == 0):
                missing_weights = np.append(missing_weights, weight)
                missing_idxs.append(i)

        if missing_weights.shape[0] == 2:
            weight = missing_weights.sum()
            delta_w = np.abs(weight_sums - weight)
            is_close = (delta_w < (1.5 * np.sqrt(weight))) * (delta_w < (1.5 * np.sqrt(weight_sums)))
            is_good_fit[missing_idxs[0]] = is_close
            is_good_fit[missing_idxs[1]] = is_close
        if missing_weights.shape[0] > 2:
            continue

        is_good_fit = (np.vstack(is_good_fit))
        for i, row in enumerate(is_good_fit):
            if row.sum() == 1:
                if is_good_fit[:, row].sum() == 1:
                    labelled_nodes[i].append(antichain[np.where(row)[0][0]])
    return labelled_nodes


def add_path(path, ksize):
    h = ''
    for i, node in enumerate(path):
        h += node[:-(ksize-2)]
    h += node[-(ksize-2):]
    return h