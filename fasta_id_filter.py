#!/usr/bin/env python

"""
Sequence ID filter for fasta files.
"""

import sys
import re
import argparse

def filter_or_not(id, filter_list, include_mode, filter_mode):
    if filter_mode == 'n':
        id = id.split()[0]
    if filter_mode in ('n', 'p') and id in filter_list:
        found = True
    elif filter_mode == 'r' and any(regex.search(id) for regex in filter_list):
        found = True
    else:
        found = False
    if found and include_mode:
        return True
    elif not found and not include_mode:
        return True
    else:
        return False

parser = argparse.ArgumentParser(description='Filter sequence ID from the list.')
parser.add_argument('fasta')
group1 = parser.add_mutually_exclusive_group()
group1.add_argument('-i', metavar='list', help='included list')
group1.add_argument('-e', metavar='list', help='excluded list')
group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-n', help='default mode (match ID before first space)', action='store_true')
group2.add_argument('-p', help='enable perfect match mode', action='store_true')
group2.add_argument('-r', help='enable regular expression mode', action='store_true')
parser.add_argument('--dry', help='don\'t run', action='store_true')

args = parser.parse_args()

fin = args.fasta

msg = 'Mode: '
if args.i:
    flist = args.i
    include_mode = True
    msg = msg + 'include, '
else:
    flist = args.e
    include_mode = False
    msg = msg + 'exclude, '

if args.p:
    filter_mode = 'p'
    msg = msg + 'perfect match.'
elif args.r:
    filter_mode = 'r'
    msg = msg + 'regular expression.'
else:
    filter_mode = 'n'
    msg = msg + 'default.'

print(msg, file=sys.stderr)

with open(flist) as f:
    filter_list = [x for x in f.read().split() if not None]

if filter_mode == 'r':
    filter_list = [re.compile(x) for x in filter_list]
else:
    filter_list = set(filter_list)

print('{} items in list.'.format(len(filter_list)), file=sys.stderr)

if args.dry:
    print('Matched items:', file=sys.stderr)
    with open(fin) as f:
        for line in f:
            if line.startswith('>'):
                p = filter_or_not(line.strip()[1:], filter_list, include_mode, filter_mode)
                if p:
                    print(line.strip()[1:])
else:
    with open(fin) as f:
        for line in f:
            if line.startswith('>'):
                p = filter_or_not(line.strip()[1:], filter_list, include_mode, filter_mode)
                if p:
                    print(line, end='')
            elif p:
                print(line, end='')
