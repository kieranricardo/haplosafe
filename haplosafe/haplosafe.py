from .de_bruijn import *
from .path_finding import *
from itertools import chain


def predict_haplotypes(fq_1_filepath=None, fq_2_filepath=None, trim_depth=500, ksize=61, cutoff=10):

    if not (fq_1_filepath is None):
        with open(fq_1_filepath) as fq:
            reads_1 = (line.strip() for i, line in enumerate(fq.readlines()) if ((i % 4) == 1))

        if not (fq_2_filepath is None):
            with open(fq_2_filepath) as fq:
                reads_2 = (line.strip() for i, line in enumerate(fq.readlines()) if ((i % 4) == 1))
            reads = chain(reads_1, reads_2)
        else:
            reads = reads_1

    else:
        raise ValueError("A fastq file must be passed.")

    de_bruijn_graph = construct_debruijn(reads, trim_depth=trim_depth, ksize=ksize, cutoff=cutoff)

    ordered_antichains, max_ac_card = get_antichains(de_bruijn_graph)

    graph_decompositions = decompose_graph(de_bruijn_graph, ordered_antichains)

    haplo_freqs = estimate_frequencies(de_bruijn_graph, ordered_antichains)

    potential_haplotypes = [[] for i in range(max_ac_card)]
    all_nodes_classified = [[] for i in range(max_ac_card)]

    for nodes, roots, leaves in graph_decompositions:
        subgraph = de_bruijn_graph.subgraph(nodes)

        roots = sorted(roots, key=lambda n: de_bruijn_graph.edges[(n, list(de_bruijn_graph.succ[n])[0])]['weight'])
        leaves = sorted(leaves, key=lambda n: de_bruijn_graph.edges[(list(de_bruijn_graph.pred[n])[0], n)]['weight'])

        labelled_nodes = label_nodes(subgraph, de_bruijn_graph, roots, leaves, haplo_freqs, max_ac_card)

        top_sort = list(nx.algorithms.dag.topological_sort(subgraph))
        subgraph_haplotypes = [[] for i in range(max_ac_card)]

        for k, node_list in enumerate(labelled_nodes):
            node_list.append(roots[k])
            node_list.append(leaves[k])
            all_nodes_classified[k].extend(node_list)
            node_list = sorted(set(node_list), key=lambda n: top_sort.index(n))

            for i in range(len(node_list) - 1):
                new_paths = list(nx.all_simple_paths(subgraph, node_list[i], node_list[i + 1]))
                subgraph_haplotypes[k].append(new_paths)

        subgraph_paths = fit_haplo(subgraph, subgraph_haplotypes, haplo_freqs, max_ac_card)
        for k in range(max_ac_card):
            potential_haplotypes[k].append(subgraph_paths[k])

    predicted_haplotypes = []

    for i, haplo in enumerate(potential_haplotypes):
        path = []
        for sub_path in haplo:
            path += sub_path[:-1]
        path.append(sub_path[-1])

        pred_haplo = add_path(path, ksize)
        predicted_haplotypes.append(pred_haplo)

    de_bruijn_graph.clear()

    return predicted_haplotypes






























