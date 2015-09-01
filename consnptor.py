#!/usr/bin/env python
import sys
import csv

try:
    import argparse
except:
    sys.stderr.write('TESTEXIT.\n')
    sys.exit(0)

from snplrr import readVCFQuick, anabl_getSeqsFromFASTA as getSequences

"""
ConSNPtor - context SNP extractor
Extracts sequence context around variant positions.
"""

def iupac(b1, b2):
    dic = {('A', 'C'): 'M', ('A', 'G'): 'R', ('A', 'T'): 'W', ('C', 'G'): 'S', ('C', 'T'): 'Y', ('G', 'T'): 'K'}
    return dic.get(tuple(sorted([b1, b2])), 'N')

def main():

    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--reference', help='Fasta file with reference sequences against which SNPs are called.')
    parser.add_argument('--vcf', help='Variant information in VCF format.')
    parser.add_argument('--flank', type=int, default=250, help='Size of flanking sequence context.')
    parser.add_argument('--out', help='Output in Fasta format.')
    args = parser.parse_args()



    flank = args.flank
    snps = dict(readVCFQuick(args.vcf))
    #    Returns dictionary over all variant positions: {(contig, pos): (count, alt, ref, genotype)}
    contigs = set([snp[0] for snp in snps])

    with open(args.out) as fo:
        for id_, seq in getSequences(args.reference):
            if id_ in contigs:
                for snp in snps:
                    if snp[0] != id_:
                        continue
                    pos = snp[1]
                    id_out = '%s:%i:%i-%i' % (id_, pos, max(1, pos - flank), min(len(seq), pos + flank))
                    pos -= 1
                    seq_out = '%s%c%s' % (seq[pos - flank:pos], iupac(snps[snp][1], snps[snp][2]), seq[pos+1:pos+flank+1])
                    fo.write('>%s\n%s\n' % (id_out, seq_out))

    pass
