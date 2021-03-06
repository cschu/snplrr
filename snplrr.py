#!/usr/bin/env python

import os
import csv
import sys
import time
import datetime

try:
    import argparse
except:
    sys.stderr.write('TESTEXIT.\n')
    sys.exit(0)


CONTIG_TABLEHEADER = ['contig', 'length',
                      'avg(coverage, susP)', 'avg(coverage, susBulk)',
                      '#SNPs(susP)', '#SNPs(susBulk)', '#SNPs(common)',
                      'SNP_freq(susP)', 'SNP_freq(susBulk)', 'SNP_freq(common)',
                      'Synteny_GM', 'NLR-motifs?']

SNP_TABLEHEADER = ['contig', 'pos',
                   # 'resP/Ref:base', 'resP/Ref:coverage',
                   'refBase',
                   'susP/Ref:base', 'susP/Ref:coverage',
                   'susBulk/Ref:base','susBulk/Ref:coverage',
                   'susP identical to susBulk?']

NO_FILTER = 0
GTYPE_HOMOZYGOUS_REF = 1
GTYPE_HOMOZYGOUS_ALT = 2
GTYPE_HETEROZYGOUS = 4
GTYPE_UNKNOWN = 8

GENOTYPES = {'0/0': GTYPE_HOMOZYGOUS_REF,
             '1/1': GTYPE_HOMOZYGOUS_ALT,
             '0/1': GTYPE_HETEROZYGOUS}

logfile = None

def getMASTInformation(fn):
    global logfile
    reader = csv.reader(open(fn), delimiter='\t', quotechar='"')
    mastInfo = {}
    for row in reader:
        logfile.write(str(row) + '\n')
        if row[0].startswith('#') or len(row) < 3:
            continue
        if row[0] not in mastInfo:
            mastInfo[row[0]] = set([])
        mastInfo[row[0]].add((row[1].replace('N/A', '???'), row[2]))
    return mastInfo


def getSyntenyInformation(fn):
    reader = csv.reader(open(fn), delimiter='\t', quotechar='"')
    syntenyInfo = {}
    for row in reader:
        if row[0].startswith('#'):
            continue
        if row[0] not in syntenyInfo:
            syntenyInfo[row[0]] = set()
        syntenyInfo[row[0]].add(row[1])
    return syntenyInfo


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

def getSNPPositions(vcf_handle, crit=None, sample='Sample1'):
    """
    Returns all variant positions according to the filter criteria.
    """
    positions = set()
    for record in vcf_handle:
        genotype = GENOTYPES.get(record.genotype(sample).data[0],
                                 GTYPE_UNKNOWN)
        if (genotype & crit) == crit:
            key = (record.CHROM, record.POS)
            positions.add(key)
    return positions

def getSNPPositionsFast(vcf_handle, crit=None, sample='Sample1'):
    return set(((record.CHROM, record.POS, record.genotype(sample).data[0]) for record in vcf_handle if (GENOTYPES.get(record.genotype(sample).data[0], GTYPE_UNKNOWN) & crit) == crit))

def readVCFQuick(fn, crit=0):
    """
    Returns dictionary over all variant positions: {(contig, pos): (count, alt, genotype)}
    """
    for line in open(fn):
        if line.startswith('#'): continue
        line = line.strip().split()
        if line[6].strip() != 'PASS': continue
        contig, pos = line[:2]
        ref, alt = line[3:5]
        genotype_fields = line[9].split(':')
        genotype = GENOTYPES.get(genotype_fields[0], GTYPE_UNKNOWN)
        yield ((contig, int(pos)),
               (int(genotype_fields[3]), alt, ref, genotype))
        #if (genotype & crit) == crit:
        #    yield ((contig, int(pos)), (int(genotype_fields[3]), baseCall))
    pass


def getSequenceLengths(fn):
    """
    Returns sequence lengths from a FastX file
    """
    return dict([(head, len(seq))
                 for head, seq in anabl_getSeqsFromFASTA(fn)])

def getAverageSNPCoverage(SNPpos, contigLengths):
    coverage = dict([(cid, []) for cid in contigLengths])
    global logfile
    logfile.write(str(coverage.keys()) + '\n')

    for snp in SNPpos:
        logfile.write(str(snp) + '\n')
        coverage[snp[0]].append(SNPpos[snp][0])
    for contig in coverage:
        if coverage[contig]:
            coverage[contig] = sum(coverage[contig])/float(len(coverage[contig]))
        else:
            coverage[contig] = 0
    return coverage

def getAverageContigCoverage(fn, contigLengths):
    """
    Returns the average read coverage for a contig
    """
    coverage = dict([(cid, 0.0) for cid in contigLengths])
    for line in open(fn):
        line = line.strip().split()
        coverage[line[0]] += int(line[3])
    for cid in coverage:
        coverage[cid] = coverage[cid] / contigLengths[cid]

    return coverage

def getAverageContigCoverageQuick(fn, contigLengths):
    coverage = {}
    currentContig, count = None, 0
    for line in open(fn):
        line = line.strip().split()
        if line[0] != currentContig:
            if currentContig is not None:
                coverage[currentContig] = float(count) / contigLengths[currentContig]
            currentContig, count = line[0], 0
        count += int(line[3])
    coverage[currentContig] = float(count) / contigLengths[currentContig]
    return coverage


def countSNPsPerContig(positions):
    snpCount = {}
    for pos in positions:
        snpCount[pos[0]] = snpCount.get(pos[0], 0) + 1
    return snpCount

def getTimeDelta(t):
    now = time.time()
    return now - t, now

def getContigLengths_(refContigs, logfile=sys.stdout):
    logfile.write('Getting contig lengths...')
    logfile.flush()
    tstamp = time.time()
    contigLengths = getSequenceLengths(refContigs)
    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))
    return contigLengths

def getVariantData_(vcf, setLabel='', crit=NO_FILTER, logfile=sys.stdout):
    logfile.write('Reading SNP data from %s...' % setLabel)
    logfile.flush()
    tstamp = time.time()
    snpsDict = dict(readVCFQuick(vcf, crit=NO_FILTER))
    filtered = set([snp for snp in snpsDict
                    if (snpsDict[snp][-1] & crit) == crit])
    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))
    return snpsDict, filtered

def run_snplrr(refContigs, contigSummary, snpTable,
               resP_vs_refMP, resP_vs_refVCF,
               susP_vs_refMP, susP_vs_refVCF,
               susB_vs_refMP, susB_vs_refVCF,
               syntenyTable, mastTable):

    global logfile
    contigLengths = getContigLengths_(refContigs, logfile=logfile)
    crit = GTYPE_HOMOZYGOUS_ALT|GTYPE_HOMOZYGOUS_REF
    resPSNPs_d, resPSNPs = getVariantData_(resP_vs_refVCF, crit=crit,
                                           setLabel='resP', logfile=logfile)
    logfile.write('resPSNPs_d/resPSNPs: %i/%i\n' % (len(resPSNPs), len(resPSNPs)))
    crit = GTYPE_HOMOZYGOUS_ALT
    susPSNPs_d, susPSNPs = getVariantData_(susP_vs_refVCF, crit=crit,
                                           setLabel='susP', logfile=logfile)
    logfile.write('susPSNPs_d/susPSNPs: %i/%i\n' % (len(susPSNPs), len(susPSNPs)))
    susBSNPs_d, susBSNPs = getVariantData_(susB_vs_refVCF, crit=crit,
                                           setLabel='susBulk', logfile=logfile)
    logfile.write('susBSNPs_d/susBSNPs: %i/%i\n' % (len(susBSNPs), len(susBSNPs)))

    syntenyInfo = getSyntenyInformation(syntenyTable)
    mastInfo = getMASTInformation(mastTable)

    # variant positions common to susceptible parents and bulk
    #susVarCommon = susPSNPs.intersection(susBSNPs)
    #logfile.write('susVarCommon: %i\n' % len(susVarCommon))
    # variant positions only occurring in either parents or bulk
    #susVarUnique = susPSNPs.union(susBSNPs).difference(susVarCommon)
    #logfile.write('susVarUnique: %i\n' % len(susVarUnique))
    # invariant positions common to susceptible parents and bulk
    # susInvCommon = set([(contig, pos + 1)
    #                     for contig in contigLengths
    #                     for pos in xrange(contigLengths[contig])
    #                     if (contig, pos + 1) not in susVarCommon.union(susVarUnique)])
    # logfile.write('susInvCommon: %i\n' % len(susInvCommon))
    # return None

    logfile.write('Calculating SNP coverage/contig...')
    logfile.flush()
    tstamp = time.time()
    coverage_susP = getAverageSNPCoverage(susPSNPs_d, contigLengths)
    coverage_susB = getAverageSNPCoverage(susBSNPs_d, contigLengths)
    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))

    logfile.write('Counting SNPs per contig...')
    logfile.flush()
    tstamp = time.time()
    snpCount_susP = countSNPsPerContig(susPSNPs)
    snpCount_susB = countSNPsPerContig(susBSNPs)
    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))

    # first step
    # any position in the control set that shows a variant
    # cannot be reliably used
    # remove those positions from the susceptible sets
    logfile.write('Getting rid off variant positions in reference...')
    logfile.flush()
    tstamp = time.time()
    susPSNPs.difference_update(resPSNPs)
    susBSNPs.difference_update(resPSNPs)
    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))
    #logfile.write('Getting rid off heterozygous positions in resP...')
    #logfile.flush()
    #tstamp = time.time()


    # second step
    # all homozygous SNPs that are common between susceptible parents
    # and bulk may indicate NB-LRR contigs
    logfile.write('Finding common SNPs between susP and susB...')
    logfile.flush()
    tstamp = time.time()
    commonSusSNPs = susPSNPs.intersection(susBSNPs)
    snpCount_common = countSNPsPerContig(commonSusSNPs)
    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))

    # third step
    # take all contigs that contain exclusively common bulk/susP SNPs
    # and no non-shared SNPs
    logfile.write('Gathering non-chimeric contigs...')
    logfile.flush()
    tstamp = time.time()
    susPOnlySNPs = susPSNPs.difference(susBSNPs)
    susPOnlyContigs = set([x[0] for x in susPOnlySNPs])

    commonSusContigs = set([x[0] for x in commonSusSNPs])
    susHomocontigs = commonSusContigs.difference(susPOnlyContigs)

    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))

    logfile.write('commonSusSNPs: %i\n' % len(commonSusSNPs))
    logfile.write('commonSusContigs: %i\n' % len(commonSusContigs))
    logfile.write('susPOnlySNPs: %i\n' % len(susPOnlySNPs))
    logfile.write('susPOnlyContigs: %i\n' % len(susPOnlyContigs))
    logfile.write('susHomocontigs: %i\n' % len(susHomocontigs))

    logfile.write('Writing contig information...')
    logfile.flush()
    tstamp = time.time()
    out_contigSummary = open(contigSummary, 'wb')

    out_contigSummary.write('\t'.join(CONTIG_TABLEHEADER) + '\n')





    for contig in sorted(susHomocontigs):
        row = [contig, contigLengths.get(contig, 0),
               '%.3f' % coverage_susP.get(contig, 0),
               '%.3f' % coverage_susB.get(contig, 0),
               int(snpCount_susP.get(contig, 0)),
               int(snpCount_susB.get(contig, 0)),
               int(snpCount_common.get(contig, 0))]
        if row[1] == 0:
            row.extend(['NA', 'NA', 'NA'])
        else:
            row.extend(['%.3f' % (row[-3]/float(row[1])),
                        '%.3f' % (row[-2]/float(row[1])),
                        '%.3f' % (row[-1]/float(row[1]))])
        row[4:7] = map(lambda x:'%.3f'%x, row[4:7])

        row.append(','.join(sorted(syntenyInfo.get(contig, ['NA']))))

        # don't like this so much
        def f(x):
            return ('%s/%s' % x) if x != 'NA' else x
        row.append(','.join(sorted(map(f, mastInfo.get(contig, ['NA'])))))

        out_contigSummary.write('\t'.join(map(str, row)) + '\n')
    out_contigSummary.close()

    out_snpTable = open(snpTable, 'wb')
    out_snpTable.write('\t'.join(SNP_TABLEHEADER) + '\n')

    #                (int(genotype_fields[3]), alt, ref, genotype))
    # ref, alt, cov, alt, cov

    NA = ('NA', 'NA')
    for contig, pos in sorted(commonSusSNPs, key=lambda x:(x[0],int(x[1]))):
        # (list(reversed(resPSNPs_d.get((contig, pos), NA)[:-1]))) + \
        row = [contig, pos,] + \
               (list(reversed(susPSNPs_d.get((contig, pos), NA)[:-1]))) + \
               (list(reversed(susBSNPs_d.get((contig, pos), NA)[:-2])))
        row.append('YES' if row[-2] == row[-4] else 'NO')
        out_snpTable.write('\t'.join(map(str, row)) + '\n')
    out_snpTable.close()

    logfile.write(' %is\n' % int((time.time() - tstamp) + 0.5))
    pass


def main(argv):

    #open(argv[-1], 'wb').write(str(sys.version) + '\n')
    #sys.exit(0)

    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--refcontigs', help='Fasta file with reference sequences against which SNPs are called.')
    parser.add_argument('--controlMP', help='Mpileup pileup file of resistant parent reads against reference (control data)')
    parser.add_argument('--controlVCF', help='VCF file of resistant parent reads against reference (control data)')
    parser.add_argument('--susP-vs-resMP', help='Mpileup pileup file of susceptible parent reads against reference')
    parser.add_argument('--susP-vs-resVCF', help='VCF file of susceptible parent reads against reference')
    parser.add_argument('--susBulk-vs-resMP', help='Mpileup pileup file of susceptible bulk reads against reference')
    parser.add_argument('--susBulk-vs-resVCF', help='VCF file of susceptible bulk reads against reference')
    parser.add_argument('--contig-summary', help='Contig summary (tabular).')
    parser.add_argument('--snp-table', help='SNP table.')
    parser.add_argument('--synteny-table')
    parser.add_argument('--mast-table')
    parser.add_argument('--logfile', help='A log file.', default='snplrr.log')


    args = parser.parse_args()

    try:
        input = [args.refcontigs, args.contig_summary, args.snp_table,
                 args.controlMP, args.controlVCF,
                 args.susP_vs_resMP, args.susP_vs_resVCF,
                 args.susBulk_vs_resMP, args.susBulk_vs_resVCF,
                 args.synteny_table, args.mast_table]
    except:
        sys.stderr.write('Error: Invalid input parameters.\n')
        sys.exit(1)

    global logfile
    logfile = open(args.logfile, 'wb')
    logfile.write(str(input) + '\n')
    #logfile.close()
    #sys.exit()
    run_snplrr(args.refcontigs, args.contig_summary, args.snp_table,
               args.controlMP, args.controlVCF, args.susP_vs_resMP,
               args.susP_vs_resVCF, args.susBulk_vs_resMP, args.susBulk_vs_resVCF,
               args.synteny_table, args.mast_table)
    logfile.close()

    pass



if __name__ == '__main__': main(sys.argv[1:])
