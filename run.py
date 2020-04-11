import os
import sys
from haplosafe import predict_haplotypes
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
        if arg == "-f":
            filepath = sys.argv[i + 1]
        elif arg == "-w":
            window_size = int(sys.argv[i + 1])
        elif arg == "-c":
            cutoff = int(sys.argv[i + 1])
        elif arg == "-k":
            ksize = int(sys.argv[i + 1])
        elif arg == "-o":
            outpath = sys.argv[i + 1]

    if os.path.isdir(filepath):
        filepaths = [os.path.join(filepath, fn) for fn in os.listdir(filepath)]
    else:
        filepaths = list(filepath.split(","))

    haplotypes, haplo_freqs = predict_haplotypes(
        filepaths=filepaths, window_size=window_size, ksize=ksize, cutoff=cutoff
    )

    sequences = [SeqRecord(Seq(haplo), str(i)) for i, haplo in enumerate(haplotypes)]
    SeqIO.write(sequences, outpath, "fasta")

    freq_path = outpath + ".freqs.txt"
    with open(freq_path, "w") as freq_file:
        freq_file.writelines([",".join(haplo_freqs.astype(str))])
