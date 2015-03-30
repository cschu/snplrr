#!/usr/bin/env python

import os
import sys
import time
import datetime

try:
    import argparse
except:
    sys.stderr.write('TESTEXIT.\n')
    sys.exit(0)


import vcf

GTYPE_HOMOZYGOUS_REF = 1
GTYPE_HOMOZYGOUS_ALT = 2
GTYPE_HETEROZYGOUS = 4
GTYPE_UNKNOWN = 8

GENOTYPES = {'0/0': GTYPE_HOMOZYGOUS_REF,
             '1/1': GTYPE_HOMOZYGOUS_ALT,
             '0/1': GTYPE_HETEROZYGOUS}

logfile = None

def anabl_getSeqsFromFASTA(fn):
    """
    Returns generator object to access sequences from a multi-FASTA file.
    Originates from 'anabl' - BLAST analysing tool, hence the prefix.
    """
    head, seq = None, ''
    for line in open(fn):
        if line[0] == '>':
            if head is not None:
                yield (head, seq) 
            head, seq = line.strip().strip('>'), ''
        else:
            seq += line.strip()
    yield (head, seq)

def getSequenceLengths(fn):
    """
    Returns sequence lengths from a FastX file
    """ 
    return dict([(head, len(seq)) 
                 for head, seq in anabl_getSeqsFromFASTA(fn)])

def getAverageContigCoverage(fn, contigLengths, fo):
    """
    Returns the average read coverage for a contig
    """
    out = open(fo, 'wb')
    coverage = dict([(cid, 0.0) for cid in contigLengths])
    for line in open(fn):
        line = line.strip().split()
        coverage[line[0]] += int(line[3])
    for cid in sorted(coverage):
        coverage[cid] = coverage[cid] / contigLengths[cid]
        out.write('%s\t%.05f\n' % (cid, coverage[cid]))
        out.flush()
    out.close()

    return coverage

def getTimeDelta(t):
    now = time.time()
    return now - t, now

def extractCoverage(pileup, refcontigs, out):

    global logfile

    logfile.write('Getting contig lengths...')
    logfile.flush()
    tstamp = time.time()
    seqLengths = getSequenceLengths(refcontigs)
    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))

    logfile.write('Calculating average contig coverage...')
    logfile.flush()
    tstamp = time.time()
    
    avgCov = getAverageContigCoverage(pileup, seqLengths, out)
    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))
    return avgCov
    # pass
    

def main(argv):

    #open(argv[-1], 'wb').write(str(sys.version) + '\n')
    #sys.exit(0)
    # print [][0]
    #sys.exit(1)

    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--refcontigs', help='Fasta file with reference sequences against which SNPs are called.')
    parser.add_argument('--pileup', help='Pileup file with coverage information for each position in each contig.')
    parser.add_argument('--out', help='The output file (tabular).')
    parser.add_argument('--logfile', help='A log file.', default='snplrr.log')

    args = parser.parse_args()

    try:
        input = [args.refcontigs, args.out, args.pileup, args.logfile]
    except:
        sys.stderr.write('Error: Invalid input parameters.\n')
        sys.exit(1)

    global logfile 
    logfile = open(args.logfile, 'wb')
    # logfile.write(str(input) + '\n')
    #coverage = extractCoverage(args.pileup, args.refcontigs)
    extractCoverage(args.pileup, args.refcontigs, args.out)
    #open(args.out, 'wb').write('\n'.join(['%s\t%.07f' % (cid, coverage[cid])
    #                                      for cid in sorted(coverage)]))
  
    logfile.close()

    pass



if __name__ == '__main__': main(sys.argv[1:])
