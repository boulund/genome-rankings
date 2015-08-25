#!/usr/bin/env python
# Fredrik Boulund 2015
# Rank hits to different genomes and plot barchart

from sys import argv, exit
from collections import Counter, defaultdict
from math import floor
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
            default="/c3se/users/boulund/Glenn/c3-c3se605-15-1/sequences/proteotyping/bacterial_20150210",
            help="Path to NCBI RefSeq directory structure [%(default)s].")
    parser.add_argument("-i", "--min-identity", dest="min_identity", metavar="i",
            type=float,
            default=90.0,
            help="Minimum identity [%(default)s].")
    parser.add_argument("-l", "--min-length", dest="min_length", metavar="l",
            type=int,
            default=6,
            help="Minimum aligned length [%(default)s].")

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
    #fragment_hits = defaultdict(list)
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
                #fragment_hits[query].append(genome)
                yield genome
            #yield query, target


def parse_accno(s, regex):
    """ Parse GI number from string.
    """
    hit = re.search(regex, s)
    if hit:
        return hit.group(1)
    else:
        return hit


def main(directory, blast8_files, min_identity, min_length):
    genomes = build_genome_dictionary(directory)
    gi_accno_regex = re.compile(r'ref\|(\w{1,2}_[\d\w]+)\.\d{1,2}\|')

    for blast8_file in blast8_files:
        genome_stats = Counter(parse_blast8(blast8_file, 
                                            genomes, 
                                            gi_accno_regex,
                                            min_identity,
                                            min_length))
        
        print blast8_file
        for genome, count in genome_stats.most_common(20):
            print count, genome
        


if __name__ == "__main__":
    options = parse_commandline(argv)
    main(options.refseq_directory, 
         options.BLAST8, 
         options.min_identity, 
         options.min_length)
