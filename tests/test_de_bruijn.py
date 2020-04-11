from haplosafe.de_bruijn import count_kmers, get_leaves, get_roots
import networkx as nx


def graph():
    test_g = nx.MultiDiGraph()
    test_g.add_edges_from([("a", "b"), ("a", "b"), ("b", "c"), ("b", "c"), ("c", "d")])
    test_g.add_edges_from([("A", "b")])
    return test_g


def test_count_kmers():
    def string_func(i):
        if i % 2 == 0:
            return "abcdefg"
        else:
            return "hijklmn"

    reads = (string_func(i) for i in range(20000))
    kmer_counts = count_kmers(reads, ksize=3)

    kmers = ["abc", "bcd", "cde", "def", "efg", "hij", "ijk", "jkl", "klm", "lmn"]
    assert sorted(kmer_counts.keys()) == kmers
    assert all(count == 10000 for count in kmer_counts.values())


def test_get_roots():
    g = graph()
    roots = get_roots(g)
    assert roots == ["A", "a"]


def test_get_leaves():
    g = graph()
    roots = get_leaves(g)
    assert roots == ["d"]
