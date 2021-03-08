#!/usr/bin/env python3
'''
    Generate genomic coordinates of all promoters, 5'UTR, 3'UTR, CDS
'''

import os
import argparse

if not os.path.exists('annotation'):
    os.makedirs('annotation')

#initialize outputfiles
PSEUDOGFF = open('annotation/genes.gff', 'w')
PROMOTERSGFF = open('annotation/promoters.gff', 'w')
UPSTREAMGFF = open('annotation/upstream.gff', 'w')
DOWNSTREAMGFF = open('annotation/downstream.gff', 'w')

def parse_genelocations(chromz, results, flank):
    """ Parse genomic regions
    Args:
        chromz (dict) : Chromosomal sizes
        results (string) : Individual gene coordinates
        flank (int) : genomic distance from start/end site
    """
    lines = results.split("\t")
    lines[3] = int(lines[3])
    lines[4] = int(lines[4])
    if lines[6] == "+":
        end = lines[3] + flank
        start = lines[3] - flank
        upend = lines[3] - 1
        upstart = lines[3] - flank
        downstart = lines[4] + 1
        downend = lines[4] + flank
    elif lines[6] == "-":
        end = lines[4] + flank
        start = lines[4] - flank
        upend = lines[4] + flank
        upstart = lines[4] + 1
        downend = lines[3] - 1
        downstart = lines[3] - flank

    if downstart < 1:
        downstart = 1
    if upstart < 1:
        upstart = 1
    if start < 1:
        start = 1

    if upend > int(chromz[lines[0]]):
        upend = chromz[lines[0]]
    if downend > int(chromz[lines[0]]):
        downend = chromz[lines[0]]

    PROMOTERSGFF.write("{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]),
                                                     start, end, "\t".join(lines[5:])))
    UPSTREAMGFF.write("{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]),
                                                    upstart, upend, "\t".join(lines[5:])))
    DOWNSTREAMGFF.write("{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]),
                                                      downstart, downend, "\t".join(lines[5:])))

def main():
    usage = "usage: %prog -g [GTF/GFF file] -f [FEATURE TYPE (gene/transcript)] \
             -c [CHROMSIZES] -d [DISTANCE(bp) from start/end site]"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("-g", "--gtf", dest="gtf", required=True,
                        help="Enter .gtf/gff file to be processed.")
    parser.add_argument("-d", "--dist", dest="distance", default=2000, type=int,
                        help="Distance (bp) from site [TSS/TES].")
    parser.add_argument("-c", "--chrom", dest="chrom", required=True,
                        help="Enter ucsc chrom sizes file to be processed.")

    options = parser.parse_args()
    print(options)

    flank = options.distance #flank distance from TSS / TES

    if options.chrom and options.gtf:
        chromsizefile = open(options.chrom, 'r')
        chrom_sizes = {}
        for line in chromsizefile:
            line = line.split('\t')
            chrom_sizes[line[0]] = line[1].rstrip("\n")

        #feature = options.feature
        gtf_file = open(options.gtf, 'r')
        
        feature_dict = {}
        if not options.gtf.split('.')[-1] == 'gtf':
            sys.exit("ERROR :\tGTF required")
            
        #reading the feature type to be used
        for line in gtf_file:
            if not line.startswith('#'):
                lines = line.split("\t")
                feature_dict[lines[2]] = lines[2]

        if 'transcript' in feature_dict:
            feature = "transcript"
        elif 'gene' in feature_dict:
            feature = "gene"
        else:
            sys.exit("ERROR :\tGTF with either transcript/gene annotation is needed")
            
        print("NOTE :\tFeature type used is '%s'" %(feature))
        
        #reading the gff_file to parse regions
        gff_file = open(options.gtf, 'r')
    
        if options.gtf.split('.')[-1] == 'gff':
            for line in gff_file:
                if not line.startswith('#'):
                    lines = line.split("\t")
                if lines[2] == feature:
                    results = ("chr{0}\t{1}".format(lines[0], "\t".join(lines[1:])))
                    PSEUDOGFF.write(results+"\n")
                    parse_genelocations(chrom_sizes, results, flank)
        elif options.gtf.split('.')[-1] == 'gtf':
            for line in gff_file:
                if not line.startswith('#'):
                    lines = line.split("\t")
                    if lines[2] == feature:
                        newline = lines[8].split(' ')
                        results = ("chr{0}\t{1}\t{2}={3}".format(lines[0],
                                                                 "\t".join(lines[1:8]),
                                                                 newline[0], newline[1]))
                        PSEUDOGFF.write(results + "\n")
                        parse_genelocations(chrom_sizes, results, flank)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
