from haplosafe import predict_haplotypes
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from collections import Counter


if __name__ == "__main__":

    fq_1_filepath = None
    fq_2_filepath = None
    trim_depth = 500
    ksize = 61
    cutoff = 10

    for i, arg in enumerate(sys.argv):
        if arg == '-1':
            fq_1_filepath = sys.argv[i + 1]
        elif arg == '-2':
            fq_2_filepath = sys.argv[i + 1]
        elif arg == '-t':
            trim_depth = int(sys.argv[i + 1])
        elif arg == '-c':
            cutoff = int(sys.argv[i + 1])
        elif arg == '-k':
            ksize = int(sys.argv[i + 1])
        elif arg == '-o':
            outpath = sys.argv[i + 1]



    haplotypes, haplo_freqs = predict_haplotypes(
        fq_1_filepath=fq_1_filepath,
        fq_2_filepath=fq_2_filepath,
        trim_depth=trim_depth,
        ksize=ksize,
        cutoff=cutoff
    )

    seq_fp = os.path.split(fq_1_filepath)[0] + '/multi_seqgen.out'

    with open(seq_fp) as fasta:
        lines = fasta.readlines()
    true_haplotypes = [line.split()[1] for line in lines[1:]]

    true_haplo_freqs = Counter(true_haplotypes).values()

    sequences = [SeqRecord(Seq(haplo), str(i)) for i, haplo in enumerate(haplotypes)]
    SeqIO.write(sequences, outpath, 'fasta')

    freq_path = outpath + ".freqs.txt"
    with open(freq_path, 'w') as freq_file:
        freq_file.writelines(
            [','.join(haplo_freqs.astype(str)),
             '\n',
             ','.join(map(str, true_haplo_freqs))])
