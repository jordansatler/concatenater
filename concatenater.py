#!/usr/bin/env python
"""
Concatenate a set of loci. Assumes loci are in
fasta format and data are not interleaved.
Option to filter loci based on percentage of 
missing taxa in concatenate_loci function.

date: 24 Feb 2022
author: J. Satler
version: 1

usage:
    python concatenater.py /path/to/fasta/files
"""

import os
import sys
import glob
import random

def get_loci(files):
    """get list of nexus files"""
    return glob.glob(files + "/*.fa*")

def read_locus(locus):
    """read and parse locus"""
    d = {}
    with open(locus, "r") as l:
        ind = ""
        for line in l:
            if line.startswith(">"):
                ind = line.strip()[1:]
            else:
                d[ind] = line.strip()
    return d

def concatenate_loci(taxa, data):
    """concatenate loci and keep track of partitions"""
    res = {i:[] for i in sorted(taxa)}
    part = []
    count = 0
    for k, v in sorted(data.items()):
        seq_len = len(random.choice(list(v.values())))

        # uncomment if you want to filter loci by missing data
        # set missing data threshold - currently at 25%
        #if float(len(v)) / len(taxa) < 0.25:
        #    continue

        for ind in taxa:
            if ind in v:
                res[ind].append(v[ind])
            else:
                # add missing data for ind
                res[ind].append("-" * seq_len)
        part.append([k, (count + 1, count + seq_len)])
        count += seq_len
    return res, part, count

def write_concat_to_file(concat, seq_len):
    """write concatenated data to new file in phylip format"""
    with open("concatenated_data.phy", "w") as out:
        out.write("{0} {1}\n".format(len(concat), seq_len))
        for ind, seq in concat.items():
            out.write("{0}{1}{2}\n".format(ind,
                                           " " * (30 - len(ind)),
                                           ''.join(seq)))

def write_partition_to_file(partition):
    """write partiton information to file"""
    with open("concatenated_partition.txt", "w") as out:
        for i in partition:
            out.write("{0} = {1}-{2}\n".format(i[0],
                                               i[1][0],
                                               i[1][1]))

def main():
    if len(sys.argv) != 2:
        print("python concatenater.py /path/to/fasta/files")
        sys.exit()

    # store results
    loci = {}
    taxa = []

    data = get_loci(sys.argv[1])
    for i in data:
        r = read_locus(i)
        loci[os.path.basename(i).split(".")[0]] = r
        [taxa.append(t) for t in r.keys() if t not in taxa]
    concat, part, seq_len = concatenate_loci(taxa, loci)

    # write to data and partitioning scheme to file
    write_concat_to_file(concat, seq_len)
    write_partition_to_file(part)

if __name__ == '__main__':
    main()
