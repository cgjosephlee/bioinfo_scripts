#!/usr/bin/env python3
'''
remove colonal sequence
'''

import sys

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
            yield title, ''.join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())
    yield title, ''.join(lines).replace(" ", "").replace("\r", "")

with open(sys.argv[1]) as f:
    FA = [x for x in FastaParser(f)]

dup_seq = {}
blacklist = []
for x in range(1, len(FA)):
    for k, v in FA[:x]:
        if k not in blacklist and FA[x][1] == v:
            blacklist.append(FA[x][0])
            dup_seq.setdefault(k, []).append(FA[x][0])
            break

if len(dup_seq) > 0:
    for k, v in sorted(dup_seq.items(), key=lambda x: len(x[1]), reverse=True):
        print('Found clones: {} <- {}'.format(k, ','.join(v)), file=sys.stderr)
    print('Total: {}, Retain: {}, Discard: {}.'.format(len(FA),
                                                       len(FA) - len(blacklist),
                                                       len(blacklist)), file=sys.stderr)

    for k, v in FA:
        if k not in blacklist:
            print('>{}\n{}'.format(k, v))
else:
    print('No clonality was found!', file=sys.stderr)
    print('Total: {}.'.format(len(FA)), file=sys.stderr)
