# Get polyA sites from 3P.SE filtered bam
# pAs as the three prime mapped location
# output into BED format: location, direction
# todo: the UCSC 0 index issue?
import argparse
import gzip

import pysam

from sam import three_prime_mapped_location, mapping_direction

# Parse args
parser = argparse.ArgumentParser(
    description='''Get polyA sites from 3P.SE filtered bam
    pAs as the three prime mapped location
    output into BED format: location, direction
    ''')
parser.add_argument('--bam',
                    help="input bam file")
parser.add_argument('--out',
                    default="out.bed.gz", type=str,
                    help="fname for output file, default: out.bed")
args = parser.parse_args()


# Read bam, Get pAs, Collapsing pAs with same location and direction
samfile = pysam.AlignmentFile(args.bam, "rb")
loc2c = {}
for i, read in enumerate(samfile):
    chr = read.reference_name
    location = three_prime_mapped_location(read)
    direction = mapping_direction(read)
    loc = "\t".join((chr,
                          str(location),
                          direction
                          ))
    try:
        loc2c[loc] = loc2c[loc] + 1
    except KeyError as e:
        loc2c[loc] = 1

# bed sort
sorted_keys = sorted(loc2c.keys(),
                     key=lambda x:(x.split()[0],
                                   int(x.split()[1]),
                                   x.split()[2]))

# output pAs
with gzip.open(args.out, 'wt') as outf:
    for k in sorted_keys:
        v = loc2c[k]
        [chr, location, direction] = k.split("\t")
        out_line = "\t".join([
            chr,
            str(location),
            str(int(location)+1),
            'pAs',
            str(v),
            direction,
            "\n"
        ])
        outf.write(out_line)

# # pAs clusters
# pA_clusters = {}  # start-end: count
# for k_ in sorted_keys:
#     v_ = loc2c[k]


