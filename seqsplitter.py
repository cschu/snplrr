#!/usr/bin/env python
import sys
import argparse

from snplrr import anabl_getSeqsFromFASTA as getSequences


def splitSequence(id_, seq, fragsize, overlap):
    p = 0
    while 1:
        if p >= len(seq) - 1:
            break
        yield '%s:%i-%i' % (id_, p + 1, min(len(seq), p + fragsize)), seq[p:p + fragsize]
        p += (fragsize - overlap)


def main():

    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--input', help='Input Fasta file.')
    parser.add_argument('--fragsize', type=int, default=10000, help='Size of the sequence fragments in bp.')
    parser.add_argument('--overlap', type=int, default=1000, help='Size of overlap between sequence fragments in bp.')
    parser.add_argument('--output', help='Output Fasta file.')
    args = parser.parse_args()

    with open(args.output, 'wb') as fo:
        for id_, seq in getSequences(args.input):
            for fragid, fragseq in splitSequence(id_, seq, args.fragsize, args.overlap):
                fo.write('>%s\n%s\n' % (fragid, fragseq))









    pass

if __name__ == '__main__': main()
