#!/usr/bin/env python
'''
A qsub in-line command wrapper.
'''

from __future__ import print_function, division
import sys
import subprocess as sp
import argparse

parser = argparse.ArgumentParser(description='A qsub in-line command wrapper.')
parser.add_argument('cmd', type=str,
                    help='Commands in quote')
parser.add_argument('-q', type=str, default='serial',
                    help='Queue')
parser.add_argument('-P', type=str, default='MST108408',
                    help='Project ID')
parser.add_argument('-N', type=str,
                    help='Job name')
parser.add_argument('-l', type=str,
                    help='Resource (select=1:ncpus=40:ngpus=4:mpiprocs=4)')
parser.add_argument('-W', type=str,
                    help='Other attributes')
parser.add_argument('-j', type=str, choices=['oe', 'eo'],
                    help='Output direction')
parser.add_argument('--opt', type=str,
                    help='Other options pass to qsub')
parser.add_argument('--shell', type=str,
                    help='Shell')
parser.add_argument('--dryrun', action='store_true',
                    help='Do not run')
args = parser.parse_args()

# qsub options
qsub_CMD = 'qsub -q {} -P {} '.format(args.q, args.P)
if args.N:
    qsub_CMD += '-N {} '.format(args.N)
if args.l:
    qsub_CMD += '-l {} '.format(args.l)
if args.W:
    qsub_CMD += '-W {} '.format(args.W)
if args.j:
    qsub_CMD += '-j {} '.format(args.j)
if args.shell:
    qsub_CMD += '-S {} '.format(args.shell)
if args.opt:
    qsub_CMD += args.opt

# cmd preset
user_CMD = 'set -e; cd \\$PBS_O_WORKDIR; ' + args.cmd

# full_CMD = '{} -- {} -c "{}"'.format(qsub_CMD, args.shell, user_CMD)
full_CMD = '''echo "{}" | {}'''.format(user_CMD, qsub_CMD)

print(full_CMD, file=sys.stderr)
if not args.dryrun:
    proc = sp.run(full_CMD, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    sys.exit(proc.returncode)
