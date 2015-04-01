#!/usr/bin/env python

import os
import sys
import argparse

#from Bio import Blast

from Bio.Blast import NCBIXML
# from Blast import NCBIXML


def doStuff(args):
    blast_records = NCBIXML.parse(open(args.blastXMLInput))
    
    with open(args.blastHitsTSV, 'wb') as out:
        for blast_record in blast_records:
            qlen = blast_record.query_length        
            qid = blast_record.query
            # hits[qid] = []
            for alignment in blast_record.alignments:            
                sid = alignment.title.split()[1]
                tid = sid.split('|')[0]
                for hsp in alignment.hsps:
                    if hsp.expect >= args.evalue:
                        continue
                    qcov = (hsp.query_end - hsp.query_start + 1.0) / qlen
                    if qcov < args.min_query_coverage:
                        continue
                    identity = hsp.identities / float(hsp.align_length)
                    if identity < args.min_identity:
                        continue
               
                    hit = (tid, sid, hsp.expect, qcov, qlen, hsp.align_length, identity)
                    out.write('\t'.join([qid] + map(str, hit)) + '\n')  
    pass


def main(argv):
    
    descr = ''
    parser = argparse.ArgumentParser(description=descr)        
    parser.add_argument('--evalue', type=float, default=1e-10)
    parser.add_argument('--min-identity', type=float, default=0.75)
    parser.add_argument('--min-query-coverage', type=float, default=0.75)
    parser.add_argument('blastXMLInput', type=str)
    parser.add_argument('blastHitsTSV', type=str)

    try:
        args = parser.parse_args()
    except:
        sys.exit(1)

    if not os.path.exists(args.blastXMLInput):
        sys.stderr.write('Input file (%s) is missing.\n' % args.blastXMLInput)
        sys.exit(1)

    doStuff(args)
      
    pass


if __name__ == '__main__': main(sys.argv[1:])
