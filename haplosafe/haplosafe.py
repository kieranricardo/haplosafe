from .de_bruijn import *
from .path_finding import *
from .utils import load_reads
from itertools import chain


def predict_haplotypes(
        fq_1_filepath=None,
        fq_2_filepath=None,
        reads=None, de_bruijn_graph=None,
        trim_depth=500,
        ksize=61, cutoff=10):

    if de_bruijn_graph is None:
        if reads is None:
            if (fq_1_filepath is None):
                raise ValueError("A fastq file must be passed.")

            reads = load_reads(fq_1_filepath=fq_1_filepath, fq_2_filepath=fq_2_filepath)

        de_bruijn_graph = build_de_bruijn(reads, trim_depth=trim_depth, ksize=ksize, cutoff=cutoff)
        del reads

    raise NotImplementedError

    return predicted_haplotypes, haplo_freqs






























