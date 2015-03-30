#!/usr/bin/env python

import sys
#from Bio.Blast import NCBIXML
from Bio import Blast

E_VALUE_THRESHOLD = 1e-10
IDENTITY_THRESHOLD = 0.75
QCOV_THRESHOLD = 0.75

def main(argv):
    blast_records = NCBIXML.parse(open(argv[0]))
    hits = {}

    for blast_record in blast_records:
        qlen = blast_record.query_length        
        qid = blast_record.query
        hits[qid] = []
        for alignment in blast_record.alignments:            
            sid = alignment.title.split()[1]
            tid = sid.split('|')[0]
            for hsp in alignment.hsps:
                if hsp.expect >= E_VALUE_THRESHOLD:
                    continue
                qcov = (hsp.query_end - hsp.query_start + 1.0) / qlen
                if qcov < QCOV_THRESHOLD:
                    continue
                identity = hsp.identities / float(hsp.align_length)
                if identity < IDENTITY_THRESHOLD:
                    continue
                
                hit = (tid, sid, hsp.expect, qcov, qlen, hsp.align_length, identity)
                print '\t'.join([qid] + map(str, hit))
                
                hits[qid].append(hit)
                

                
    pass


if __name__ == '__main__': main(sys.argv[1:])
