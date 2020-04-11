from .de_bruijn import build_de_bruijn, get_roots, get_leaves, fast_trim
from .path_finding import merge_bubbles
from .utils import load_reads
import numpy as np


def predict_haplotypes(
    filepaths=None, reads=None, graph=None, ksize=61, cutoff=10, window_size=3,
):

    if graph is None:
        if reads is None:
            if filepaths is None:
                raise ValueError("A fastq file must be passed.")

            reads = load_reads(filepaths=filepaths)

        graph = build_de_bruijn(reads, ksize=ksize, cutoff=cutoff)
        del reads
        fast_trim(graph, forward=True)
        fast_trim(graph, forward=False)

    while True:
        graph, cutoff = merge_bubbles(graph, cutoff, ksize, window_size=window_size)

        terminal_nodes = set(get_roots(graph) + get_leaves(graph))
        if len(terminal_nodes) == len(graph.nodes):
            break

    haplotypes = [edge[1]["kmer"] for edge in graph.edges.items()]
    freqs = np.array([edge[1]["weight"] for edge in graph.edges.items()])

    return haplotypes, freqs
