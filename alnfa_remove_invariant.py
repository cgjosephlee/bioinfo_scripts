#!/usr/bin/env python3
'''
Also try https://github.com/sanger-pathogens/snp-sites,
which is faster and mem-efficient.
'''

import sys
import numpy as np

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
                    yield name, seq, ''.join(seqs) # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

print('Now reading...', file=sys.stderr)
titles = []
seqs = []
with open(sys.argv[1]) as f:
    for t, s, q in readfq(f):
        titles.append(t)
        seqs.append(s)

# https://stackoverflow.com/a/9493192/7859425
FA = np.array(seqs, dtype=np.string_)
FA = FA.view('S1').reshape((len(titles), -1), order='C')
print(FA.shape, file=sys.stderr)
# print(FA.flags, file=sys.stderr)
del seqs

print('Now filtering...', file=sys.stderr)
pass_bool = np.empty((FA.shape[0]-1, FA.shape[1]), dtype=np.bool)
r = FA[0]
for i in range(1, FA.shape[0]):
    q = FA[i]
    pass_bool[i-1] = r != q  # vectorized function, np.equal does not work
pass_indices = np.nonzero(np.sum(pass_bool, axis=0))[0]
# print(pass_indices.shape, file=sys.stderr)

out = '''\

No. samples:  {}
No. columns:  {}
Pass columns: {}
Pass rate:    {:.4f}'''.format(
    len(titles),
    FA.shape[1],
    len(pass_indices),
    len(pass_indices) / FA.shape[1]
)

print('Now printing...', file=sys.stderr)
FA_flt = np.ascontiguousarray(FA[:,pass_indices])
FA_flt = FA_flt.view(f'S{len(pass_indices)}').ravel()
for t, s in zip(titles, FA_flt):
    print('>{}\n{}'.format(t, s.decode()))
print('Finish.', file=sys.stderr)
print(out, file=sys.stderr)
