#!/usr/bin/env python3

import sys
import re

# https://github.com/lh3/readfq/blob/master/readfq.py
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

FIN = open(sys.argv[1])
re_DNA = re.compile('[atcgATCG]+')
re_N = re.compile('[nN]+')
# case sensitive?

currentPOS = 0
for name, seq, qual in readfq(FIN):
    while currentPOS < len(seq):
        if seq[currentPOS] in set('atcgATCG'):
            r = re_DNA.match(seq, pos=currentPOS)
            print('{}\t{}\t{}\t{}\t{}'.format(
                name,
                currentPOS,
                r.end(0),
                r.end(0) - currentPOS,
                'DNA'
            ))
            currentPOS = r.end(0)
            # breakpoint()
        elif seq[currentPOS] in set('nN'):
            r = re_N.match(seq, pos=currentPOS)
            print('{}\t{}\t{}\t{}\t{}'.format(
                name,
                currentPOS,
                r.end(0),
                r.end(0) - currentPOS,
                'N'
            ))
            currentPOS = r.end(0)
        else:
            raise ValueError('Not DNA sequence?')
    currentPOS = 0

FIN.close()
