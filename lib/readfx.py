#!/usr/bin/env python3
'''
A collection of functions parsing fasta/q.
'''

'''
https://github.com/lh3/readfq/blob/master/readfq.py
https://github.com/lh3/readfq/pull/6
https://github.com/lh3/minimap2/tree/master/python#miscellaneous-functions
'''

def readfx(fp, binary=False): # this is a generator function
    if binary:
        fx_header1, fx_header2, fx_space = b'>@', b'@+>', b' '
        fq_3rd = ord('+') # subscription truns binary character into int
        seq_sep = b''
    else:
        fx_header1, fx_header2, fx_space = '>@', '@+>', ' '
        fq_3rd = '+'
        seq_sep = ''
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in fx_header1: # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, _, comment = last[1:].partition(fx_space)
        seqs, last = [], None
        for l in fp: # read the sequence
            if l[0] in fx_header2:
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != fq_3rd: # this is a fasta record
            yield name, seq_sep.join(seqs), None, comment # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = seq_sep.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, seq_sep.join(seqs), comment # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None, comment # yield a fasta record instead
                break

'''
https://github.com/ndaniel/readfq

Why python2 iterate file faster than python3?
https://stackoverflow.com/q/52545269/7859425
'''

'''
text version
'''

def readfq(fp):
    for title, seq, _, qual in zip(*[fp]*4):
        assert title[0] == '@'
        assert _[0] == '+'
        name, _, comment = title[1:-1].partition(' ')
        yield name, seq[:-1], qual[:-1], comment

def readfa(fp):
    for line in fp:
        if line[0] == '>':
            name, _, comment = line[1:-1].partition(' ')
            seqs = []
            break
    for line in fp:
        if line[0] == '>':
            yield name, ''.join(seqs), None, comment
            name, _, comment = line[1:-1].partition(' ')
            seqs = []
        else:
            seqs.append(line[:-1])
    yield name, ''.join(seqs), None, comment

'''
binary version
'''

def readfq_bin(fp):
    for title, seq, _, qual in zip(*[fp]*4):
        assert title[0] == 64  # b'@'
        assert _[0] == 43  # b'+'
        name, _, comment = title[1:-1].partition(b' ')
        yield name, seq[:-1], qual[:-1], comment

def readfa_bin(fp):
    for line in fp:
        if line[0] == 62:  # b'>'
            name, _, comment = line[1:-1].partition(b' ')
            seq = b''
            break
    for line in fp:
        if line[0] == 62:  # b'>'
            yield name, seq, None, comment
            name, _, comment = line[1:-1].partition(b' ')
            seq = b''
        else:
            seq += line[:-1]
    yield name, seq, None, comment

if __name__ == '__main__':
    import sys
    fn = sys.argv[1]

    fp = open(fn, 'rb')
    n, slen, qlen = 0, 0, 0
    for name, seq, qual, comment in readfx(fp, binary=True):
        n += 1
        slen += len(seq)
        qlen += qual and len(qual) or 0
    print(n, '\t', slen, '\t', qlen)
