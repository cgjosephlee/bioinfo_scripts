#!/usr/bin/env python3

import sys
import numpy as np
from tqdm import tqdm

def FastaParser(handle):
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == '>':
            title = line[1:].rstrip()
            break
    else:   # no break encountered
        return  # Premature end of file, or just empty?

    lines = []
    for line in handle:
        if line[0] == '>':
            yield title, list(''.join(lines).replace(" ", "").replace("\r", ""))
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())
    yield title, list(''.join(lines).replace(" ", "").replace("\r", ""))

print('Now reading...', file=sys.stderr)
titles = []
seqs = []
with open(sys.argv[1]) as f:
    for t, s in FastaParser(f):
        titles.append(t)
        seqs.append(s)
FA = np.array(seqs, dtype='U1')
print(FA.shape, file=sys.stderr)
# peak memory usage here
del seqs

print('Now filtering...', file=sys.stderr)
pbar = tqdm(total=FA.shape[1])
pass_indices = []
for n, s in enumerate(FA.T):
    if not all(x == s[0] for x in s):
        pass_indices.append(n)
    pbar.update(1)
pbar.close()

out = '''\
All samples: {}
All columns: {}
Pass columns: {}
Pass rate: {:.4f}'''.format(
    len(titles),
    FA.shape[1],
    len(pass_indices),
    len(pass_indices) / FA.shape[1]
)

print('Now picking...', file=sys.stderr)
FA_flt = np.take(FA.T, pass_indices, axis=0).T
print('Now printing...', file=sys.stderr)
pbar = tqdm(total=FA_flt.shape[0])
for t, s in zip(titles, FA_flt):
    print('>{}\n{}'.format(t, ''.join(s)))
    pbar.update(1)
pbar.close()
print('Finish.', file=sys.stderr)
print(out, file=sys.stderr)
