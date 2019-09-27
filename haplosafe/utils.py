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