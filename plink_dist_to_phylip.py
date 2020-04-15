#!/usr/bin/env python3

import sys

in_dist = sys.argv[1]
in_id = in_dist + '.id'

# read id
ids = []
with open(in_id) as f:
    for line in f:
        line = line.strip().split()
        ids.append(line[0])

# assuming square output
distance = []
with open(in_dist) as f:
    for line in f:
        line = line.strip().split()
        distance.append(line)

# output
out = in_dist + '.phy'
out_handle = open(out, 'w')
if len(ids) != len(distance):
    raise ValueError('Something wrong.')
print(len(ids), file=out_handle)
maxLen = max(len(x) for x in ids)
if maxLen < 10:
    maxLen = 10
else:
    # trim?
    pass
for i, ds in zip(ids, distance):
    print('{: <{}} {}'.format(
        i, maxLen,
        ' '.join(['{:.6f}'.format(float(x)) for x in ds])
    ), file=out_handle)
