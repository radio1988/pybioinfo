import pysam
from Bio.Seq import Seq

def contain_A(read, n=2):
    '''test if contain A in fastq read found in bam line
    :return: true/false
    '''
    seq = Seq(read.query_sequence)
    if read.flag & 16:
        # reverse
        seq = seq.reverse_complement()
    else:
        # forward alignment
        pass

    if seq.endswith("A"*n):
        return True
    else:
        return False


def contain_non_templated_A(read, n=1):
    '''test if contain non-genomic-templated ending A in fastq read found in bam line
    :return: true/false
    '''
    if read.flag & 16:
        # reverse
        seq = Seq(read.query_sequence)
        seq = seq.reverse_complement()
        three_prime_cigar = read.cigar[0]
    else:
        # forward
        seq = Seq(read.query_sequence)
        three_prime_cigar = read.cigar[-1]

    flag = False
    if three_prime_cigar[0] == 4:   # S: soft-clip
        soft_clipped_seq = seq[-three_prime_cigar[1]:]
        if soft_clipped_seq.endswith("A"*n):
            flag = True

    return flag


def high_mq(read, min=20):
    '''test if has large MQ for the current bam line
    :return: true/false
    '''
    if read.mapping_quality >= min:
        return True
    else:
        return False


def mapped_locations(read, max_locations = 1):
    '''
    :return: num mapping location (NH) for the current bam line
    '''
    print("TBA")
    return 0


def num_mismatch(read):
    '''
    :param read:
    :return: mis-match as measured by NM
    '''
    stats = read.get_cigar_stats()
    num = stats[0][10]  # NM, editing distance excluding clipping
    return num


def mapped_length(read):
    stats = read.get_cigar_stats()
    mapped_length = stats[0][0]  # M, matches
    return mapped_length
