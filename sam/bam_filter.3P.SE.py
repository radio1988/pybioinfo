import argparse

import pysam
from Bio.Seq import Seq

from sam import contain_A, contain_non_templated_A, high_mq, num_mismatch, mapped_length


# Parse args
parser = argparse.ArgumentParser(description="filter bam for 3P.SE data")
parser.add_argument('--bam',
                    help="input of bam file")
parser.add_argument('--genome',
                    default='NULL',
                    help="input of genome.fasta, optional")
parser.add_argument('--num_a',
                    default=2, type=int,
                    help="number of ending A's required in query read, default 2")
parser.add_argument('--num_untemp_a',
                    default=1, type=int,
                    help="number of un-genome-templated ending A's required in query read, default 1")
parser.add_argument('--min_mq',
                    default=20, type=int,
                    help="min MAPQ/MQ, default 20")
parser.add_argument('--out',
                    default="out.bam", type=str,
                    help="fname for output bam file, default: out.bam")
args = parser.parse_args()


# Read And Filter
samfile = pysam.AlignmentFile(args.bam, "rb")

with pysam.AlignmentFile(args.out, "wb", header=samfile.header) as outf:
    for i, read in enumerate(samfile):
        if not high_mq(read, args.min_mq):
            continue
        if num_mismatch(read) > 2:
            continue
        if mapped_length(read) < 15:
            continue
        if not contain_A(read, args.num_a):
            continue
        if not contain_non_templated_A(read, args.num_untemp_a):
            continue

        outf.write(read)
