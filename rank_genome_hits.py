#!/usr/bin/env python
# Fredrik Boulund 2015
# Rank hits to different genomes and plot barchart

from sys import argv, exit
from collections import Counter
#import argparse
import fnmatch
import os
import re


if len(argv) < 2:
    print "usage: script.py BLAST8_FILE[s]..."
    exit()


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


def parse_blast8(blast8_file, genomes, regex):
    """ Parse blast8 file.
    """
    with open(blast8_file) as f:
        for line in f:
            splitline = line.split()
            #query = splitline[0]
            target = splitline[1]
            accno = parse_accno(target, regex)
            genome = genomes[accno]
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


def generate_genome_stats(blast8_file, genomes, regex):
    """ Generate counts for all seen genomes.
    """
    return Counter(parse_blast8(blast8_file, genomes, regex))


def main(directory, blast8_files):
    genomes = build_genome_dictionary(directory)
    gi_accno_regex = re.compile(r'ref\|(\w{1,2}_[\d\w]+)\.\d{1,2}\|')

    for blast8_file in blast8_files:
        genome_stats = generate_genome_stats(blast8_file, genomes, gi_accno_regex)
        
        print blast8_file
        for genome, count in genome_stats.most_common(10):
            print count, genome
        


if __name__ == "__main__":
    directory = "/c3se/users/boulund/Glenn/c3-c3se605-15-1/sequences/proteotyping/bacterial_20150210"
    blast8_files = argv[1:]
    main(directory, blast8_files)
