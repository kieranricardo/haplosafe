def reverse_complement(read):
    read = read[::-1]
    read = read.replace("A", "t").replace("C", "g")
    read = read.replace("T", "a").replace("G", "c")
    return read.upper()


def parse_fasta(filepath):

    with open(filepath) as fasta_file:
        reads = [
            line.strip()
            for i, line in enumerate(fasta_file.readlines())
            if ((i % 4) == 1)
        ]

    return reads


def parse_sam(filepath):
    def _record_check(line):
        return (len(line.split()) > 9) and (len(line.split()[9]) > 100)

    with open(filepath) as sam_file:
        reads = [
            line.split()[9].strip()
            for line in sam_file.readlines()
            if _record_check(line)
        ]

    return reads


def load_reads(filepaths):

    reads = []
    for filepath in filepaths:
        extension = filepath.split(".")[-1]
        if extension == "sam":
            reads.extend(parse_sam(filepath))
        elif extension == "fas" or extension == "fq":
            reads.extend(parse_fasta(filepath))
        else:
            raise ValueError(
                f"filepath: expected extension to one of '.sam', '.fas', found {extension}"
            )

    reads = [read for read in reads if 'N' not in read]
    rev_comp_reads = [reverse_complement(read) for read in reads]

    return reads + rev_comp_reads
