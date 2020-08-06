#!/usr/bin/env python3

import sys
import gzip

in1 = sys.argv[1]
in2 = sys.argv[2]
sep = sys.argv[3]

if in1.endswith('.gz'):
    handle1 = gzip.open(in1, 'rt')
else:
    handle1 = open(in1)
if in2.endswith('.gz'):
    handle2 = gzip.open(in2, 'rt')
else:
    handle2 = open(in2)

for F1, R1 in zip(handle1, handle2):
    nameF = F1.partition(' ')
    nameR = R1.partition(' ')

    if not nameF[0] == nameR[0]:
        print(nameF, nameR, file=sys.stderr)
        raise ValueError()

    sys.stdout.write(nameF[0]+sep+'1 '+nameF[2])
    sys.stdout.write(next(handle1))
    sys.stdout.write(next(handle1))
    sys.stdout.write(next(handle1))
    sys.stdout.write(nameR[0]+sep+'2 '+nameR[2])
    sys.stdout.write(next(handle2))
    sys.stdout.write(next(handle2))
    sys.stdout.write(next(handle2))

handle1.close()
handle2.close()
