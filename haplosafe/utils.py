from itertools import chain


def complement(base):
    if (base == 'A'):
        return 'T'
    elif (base == 'C'):
        return 'G'
    elif (base == 'T'):
        return 'A'
    elif (base == 'G'):
        return 'C'
    else:
        raise ValueError("Base must be 'A', 'C', 'T', or 'G'. Recieved value {} instead.".format(base))


def reverse_complement(read):
    revcomp = ''.join(map(complement, read[::-1].upper()))
    return revcomp


def load_reads(fq_1_filepath=None, fq_2_filepath = None):

    with open(fq_1_filepath) as fq:
        reads_1 = (line.strip() for i, line in enumerate(fq.readlines()) if ((i % 4) == 1))

    rev_reads_1 = (map(reverse_complement, reads_1))

    reads = chain(reads_1, rev_reads_1)

    if not (fq_2_filepath is None):
        with open(fq_2_filepath) as fq:
            reads_2 = (line.strip() for i, line in enumerate(fq.readlines()) if ((i % 4) == 1))

        rev_reads_2 = (map(reverse_complement, reads_2))

        reads = chain(reads, reads_2, rev_reads_2)

    return reads