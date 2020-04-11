import numpy as np
import networkx as nx
from haplosafe.path_finding import find_source_sink_paths, merge_bubbles
from haplosafe.de_bruijn import build_de_bruijn


def random_graph(ksize):

    base_haplo = np.random.choice(["A", "C", "G", "T"], size=20 * ksize)
    haplotypes = [base_haplo.copy() for _ in range(5)]

    for _ in range(20):
        idx = np.random.randint(0, len(base_haplo) - 1)

        for hidx in np.random.randint(0, len(haplotypes) - 1, size=2):
            haplo = haplotypes[hidx]
            haplo[idx] = np.random.choice(["A", "C", "G", "T"])

    reads = ["".join(haplo) for haplo in haplotypes]

    g = build_de_bruijn(reads, ksize=ksize, cutoff=0)
    return g


def test_merge_bubbles_0():

    ksize = 10
    g = random_graph(ksize=ksize)
    cutoff = 1
    for _ in range(3):

        print("#### cutoff", cutoff)
        print("#### multi-edge g", any(i > 0 for _, _, i in g.edges))

        homo = nx.DiGraph(g)
        print("#### acyclic:", nx.algorithms.dag.is_directed_acyclic_graph(homo))

        new_g, cutoff, *_ = merge_bubbles(g, cutoff=cutoff, ksize=ksize, window_size=20)

        # as soon as there's a multi-edge in g this whole thing fails
        print("#### multi-edge g", any(i > 0 for _, _, i in new_g.edges))

        assert len(new_g.nodes) <= len(g.nodes)
        assert len(find_source_sink_paths(g)) >= len(find_source_sink_paths(new_g))

        g = new_g
        print()


def test_merge_bubbles_1():

    ksize = 10

    base_haplo = np.random.choice(["A", "C", "G", "T"], size=4 * ksize)

    edges = []

    for h in [base_haplo[: 2 * ksize], base_haplo[2 * ksize - ksize :]]:

        nodes = ["".join(h[:ksize]), "".join(h[-ksize:])]

        edited_kmer = h.copy()
        edited_kmer[ksize + 1] = "A" if edited_kmer[ksize + 1] != "A" else "T"

        edges.extend(
            [
                (nodes[0], nodes[1], {"weight": 1, "kmer": "".join(h)}),
                (nodes[0], nodes[1], {"weight": 1, "kmer": "".join(edited_kmer)}),
            ]
        )

    assert edges[0][1] == edges[-1][0]
    g = nx.MultiDiGraph()
    g.add_edges_from(edges)

    assert any(i > 0 for _, _, i in g.edges)
    #
    new_g, *_ = merge_bubbles(g, cutoff=0, ksize=ksize, window_size=2)
    # assert len(find_source_sink_paths(g)) != len(find_source_sink_paths(new_g))
