from haplosafe import predict_haplotypes
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


    haplotypes = predict_haplotypes(fq_1_filepath=fq_1_filepath, fq_2_filepath=fq_2_filepath,
                            trim_depth=trim_depth, ksize=ksize, cutoff=cutoff)

    sequences = [SeqRecord(Seq(haplo), str(i)) for i, haplo in enumerate(haplotypes)]
    SeqIO.write(sequences, outpath, 'fasta')
