#!/usr/bin/env python3
'''
A very fast multi-thread ORF finder.

The script generates same results with NCBI ORFFinder, except this script:
    - Do name serially by coordinates
    - Do return ORFs in small scaffolds
    - Do not return truncated ORF at scaffold ends
    - Do not handle nested ORFs (since they are in different frames)
'''

import sys
import argparse
# from multiprocessing.dummy import Pool
from multiprocessing import Pool
from tqdm import tqdm
from Bio.Data import CodonTable
from Bio.SeqIO.FastaIO import SimpleFastaParser

def parse_args():
    parser = argparse.ArgumentParser(description='Multi-thread ORF finder. Masked regions and gaps (N) are skipped.')
    parser.add_argument('fasta',
                        help='fasta file')
    parser.add_argument('-l', type=int, default=100,
                        help='minimum nucleotide length (100)')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout,
                        help='output file (stdout)')
    parser.add_argument('-t', type=int, default=1,
                        help='threads (1)')
    parser.add_argument('--format', type=str, choices=['bed', 'gff'], default='bed',
                        help='output format (bed)')
    parser.add_argument('--table', type=int, default=1,
                        help='codon table (1)')
    parser.add_argument('--atg', action='store_true',
                        help='use ATG as start codon only')
    return parser.parse_args()

def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in seq[::-1])

def parseTable(table=1, ATGOnly=False):
    # utilize Bio.Data.CodonTable
    try:
        codons = CodonTable.unambiguous_dna_by_id[table]
        startCodon = set(codons.start_codons)
        stopCodon = set(codons.stop_codons)
    except KeyError:
        raise KeyError('Invalid table ID ({}).'.format(table))
    if ATGOnly:
        startCodon = set(['ATG'])
    startCodonRC = set([revcomp(x) for x in startCodon])
    stopCodonRC = set([revcomp(x) for x in stopCodon])
    return startCodon, stopCodon, startCodonRC, stopCodonRC

def SearchCodons(seqID, seq, frame=0, is_positive=True, minLen=100, startCodon=set(), stopCodon=set()):
    '''
    seqID: str
    seq: str
    startCodon: set
    stopCodon: set
    '''
    assert isinstance(seqID, str)
    assert isinstance(seq, str)
    # seq = seq.upper()

    orf_list = []

    if is_positive:
        startPos = None
        for n in range(frame, len(seq), 3):
            codon = seq[n:n+3]
            if codon in startCodon and not startPos:
                startPos = n
            elif codon in stopCodon and startPos:
                if n - startPos + 3 > minLen:
                    orf_list.append((seqID, startPos, n+3, frame, '+'))
                    startPos = None
                else:
                    startPos = None
                    continue
            elif codon.find('N') != -1:  # gap or masked
                startPos = None
    else:
        # assume codons are reverse complemented
        startPos = None
        stopPos = None
        for n in range(frame, len(seq), 3):
            anticodon = seq[n:n+3]
            if anticodon in stopCodon:
                if startPos and stopPos:
                    if startPos - stopPos + 3 > minLen:
                        orf_list.append((seqID, stopPos, startPos+3, frame, '-'))
                    startPos = None
                stopPos = n
            elif anticodon in startCodon:
                startPos = n
            elif anticodon.find('N') != -1:  # gap or masked
                startPos = None
                stopPos = None
        # last check
        if startPos and stopPos and startPos - stopPos > minLen:
            orf_list.append((seqID, stopPos, startPos+3, frame, '-'))

    return orf_list

def generate_output(results, handle, format, fasta=[], table=1):
    if format == 'bed':
        for n, rec in enumerate(results, 1):
            print('{}\t{}\t{}\t{}\t{}\t{}'.format(
                rec[0],
                rec[1],
                rec[2],
                'ORF{}'.format(n),
                '.',
                rec[4]), file=handle)
    elif format == 'gff':
        for n, rec in enumerate(results, 1):
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                rec[0],
                'orf_finder',
                'gene',
                rec[1] + 1,
                rec[2],
                '.',
                rec[4],
                '.',
                'ID=ORF{0};Name=ORF{0}'.format(n)
            ), file=handle)
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                rec[0],
                'orf_finder',
                'mRNA',
                rec[1] + 1,
                rec[2],
                '.',
                rec[4],
                '.',
                'ID=ORF{0}.t01;Parent=ORF{0}'.format(n)
            ), file=handle)
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                rec[0],
                'orf_finder',
                'exon',
                rec[1] + 1,
                rec[2],
                '.',
                rec[4],
                '.',
                'Parent=ORF{0}.t01;Name=ORF{0}'.format(n)
            ), file=handle)
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                rec[0],
                'orf_finder',
                'CDS',
                rec[1] + 1,
                rec[2],
                '.',
                rec[4],
                '0',
                'ID=ORF{0}.p01;Parent=ORF{0}.t01;Name=ORF{0}'.format(n)
            ), file=handle)
    elif format == 'fna':
        for n, rec in enumerate(results, 1):
            pass
    elif format == 'faa':
        for n, rec in enumerate(results, 1):
            pass

# hack tqdm to behave like progressbar
# tqdm >= 4.30.0, tested with 4.32.2
class tqdm(tqdm):
    def update_to_value(self, n):
        if n > self.last_print_n:
            self.n = n
            with self._lock:
                self.display()
            self.last_print_n = self.n

def main():
    args = parse_args()
    in_fa = args.fasta
    minLen = args.l
    threads = args.t
    OutHandle = args.o
    OutFormat = args.format
    tableID = args.table

    print('Codon table = {}.'.format(tableID), file=sys.stderr)

    startCodon, stopCodon, startCodonRC, stopCodonRC = parseTable(tableID, args.atg)
    with open(in_fa) as f:
        FA = [(t.split()[0], s.upper()) for (t, s) in SimpleFastaParser(f)]

    print('Found {} sequences.'.format(len(FA)), file=sys.stderr)
    print('Search in six frames...', file=sys.stderr)

    results = []
    pool = Pool(processes=threads)
    for seqID, seq in FA:
        pool.apply_async(SearchCodons, (seqID, seq, 0, True, minLen, startCodon, stopCodon), callback=results.append)
        pool.apply_async(SearchCodons, (seqID, seq, 1, True, minLen, startCodon, stopCodon), callback=results.append)
        pool.apply_async(SearchCodons, (seqID, seq, 2, True, minLen, startCodon, stopCodon), callback=results.append)
        pool.apply_async(SearchCodons, (seqID, seq, 0, False, minLen, startCodonRC, stopCodonRC), callback=results.append)
        pool.apply_async(SearchCodons, (seqID, seq, 1, False, minLen, startCodonRC, stopCodonRC), callback=results.append)
        pool.apply_async(SearchCodons, (seqID, seq, 2, False, minLen, startCodonRC, stopCodonRC), callback=results.append)
    pool.close()

    pb = tqdm(total=len(FA) * 6)
    while len(results) < len(FA) * 6:
        pb.update_to_value(len(results))
    pb.update_to_value(len(results))
    pb.close()

    pool.join()

    print('Generating output...', file=sys.stderr)
    # flatten list
    results = [item for sublist in results for item in sublist]
    results.sort(key=lambda x: (x[0], x[1], x[2]))
    generate_output(results, OutHandle, OutFormat, FA, tableID)
    if OutHandle != sys.stdout:
        OutHandle.close()

    print('Finish!', file=sys.stderr)

if __name__ == '__main__':
    main()
