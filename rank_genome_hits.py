#!/usr/bin/env python2.7
# Fredrik Boulund 2015
# Rank hits to different genomes and plot barchart

from sys import argv, exit
from collections import Counter, defaultdict
from math import floor
from string import replace
import numpy as np
import matplotlib.pyplot as plt
import argparse
import fnmatch
import os
import re


def parse_commandline(argv):
    """ Parse commandline arguments.
    """

    desc = "Produce basic genome hit ranking plots. (c) Fredrik Boulund 2015."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("BLAST8", nargs="+",
            help="BLAT output file(s).")
    parser.add_argument("-d", "--refseq-dir", dest="refseq_dir", metavar="DIR",
            default="/shared/genomes/NCBI/bacterial/20150210/",
            help="Path to NCBI RefSeq directory structure [%(default)s].")
    parser.add_argument("-i", "--min-identity", dest="min_identity", metavar="i",
            type=float,
            default=90.0,
            help="Minimum identity [%(default)s].")
    parser.add_argument("-l", "--min-length", dest="min_length", metavar="l",
            type=int,
            default=6,
            help="Minimum aligned length [%(default)s].")
    parser.add_argument("-p", "--print-ranking", dest="print_ranking", action="store_true",
            default=False,
            help="Print genome rankings to stdout in addition to plot production [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit()

    return parser.parse_args(argv[1:])


def find_files(directory, pattern):
    """ Yield filenames by recursively searching dir with glob pattern.
    """
    for root, subfolders, files in os.walk(directory, followlinks=True):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


def build_genome_dictionary(d):
    """ Construct an accno:genome dictionary.
    """
    genomes = {}    
    for filename in find_files(d, "*.fna"):
        other, genome, accno = filename.rsplit("/", 2)
        genomes[accno.split(".")[0]] = genome
    return genomes


def parse_blast8(blast8_file, genomes, regex, min_identity, min_length):
    """ Parse blast8 file.
    """
    # BLAT blast8 output format:
    # query_id, subject_id, %_identity, alignment_length, mismatches, gap_openings, q.start, q.end, s.start, s.end, e-value, bitscore
    with open(blast8_file) as f:
        for line in f:
            splitline = line.split()
            query = splitline[0]
            target = splitline[1]
            identity = float(splitline[2])
            length = int(splitline[3])
            mismatches = int(splitline[4])
            if identity > min_identity and length > min_length:
                accno = parse_accno(target, regex)
                genome = genomes[accno]
                yield genome


def parse_accno(s, regex):
    """ Parse GI number from string.
    """
    hit = re.search(regex, s)
    if hit:
        return hit.group(1)
    else:
        return hit


def plot_rankings(genome_stats, blast8):
    """ Plot barchart of genomes in falling order.
    """

    N = 20
    indices = np.arange(N)
    width = 0.65

    genomes, ranks = zip(*genome_stats.most_common(20))
    genomes = [replace(genome, "_", " ") for genome in genomes]

    fig, ax = plt.subplots(figsize=(10,15))
    rects = ax.bar(indices, ranks, width)

    ax.set_title("Genome rankings for {}".format(blast8))
    ax.set_ylabel("Hits")
    ax.set_xticks(indices+width)
    ax.set_xticklabels(genomes, rotation="vertical")

    fig.tight_layout()
    fig.savefig("{}_rankings.pdf".format(blast8))
    fig.savefig("{}_rankings.png".format(blast8))


def print_rankings(genome_stats, blast8):
    """ Print rankings to stdout.
    """

    print blast8
    for genome, count in genome_stats.most_common(20):
        print "{:>10} {:<}".format(count, genome)


def write_rankings(genome_stats, blast8):
    """ Write rankings to file.
    """
    with open(blast8+"_rankings.txt", "w") as f:
        for genome, count in genome_stats.most_common():
            f.write("{:>10} {:<}\n".format(count, genome))


def main(directory, blast8_files, min_identity, min_length, print_ranking):
    genomes = build_genome_dictionary(directory)
    gi_accno_regex = re.compile(r'ref\|(\w{1,2}_[\d\w]+)\.\d{1,2}\|')

    for blast8_file in blast8_files:
        genome_stats = Counter(parse_blast8(blast8_file, 
                                            genomes, 
                                            gi_accno_regex,
                                            min_identity,
                                            min_length))
        
        plot_rankings(genome_stats, blast8_file)
        write_rankings(genome_stats, blast8_file)
        if print_ranking:
            print_rankings(genome_stats, blast8_file)


if __name__ == "__main__":
    options = parse_commandline(argv)
    main(options.refseq_dir, 
         options.BLAST8, 
         options.min_identity, 
         options.min_length,
         options.print_ranking)
