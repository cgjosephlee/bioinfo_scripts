#!/usr/bin/env python3
'''
See also:
    https://www.ncbi.nlm.nih.gov/sites/batchentrez
    https://github.com/kblin/ncbi-acc-download
    https://github.com/kblin/ncbi-genome-download
ToDo:
    - a record a file?
'''

import sys
import argparse
from Bio import Entrez

def rm_dup(seq):
    # python >= 3.7
    # https://www.peterbe.com/plog/fastest-way-to-uniquify-a-list-in-python-3.6
    if not sys.version_info >= (3, 7):
        print('Order may not be preserved due to old python version (<3.7).', file=sys.stderr)
    return list(dict.fromkeys(seq))

parser = argparse.ArgumentParser(description='''\
Download fasta/genbank records from ncbi.''')
parser.add_argument('-a', type=str, nargs='+', default=[], metavar='acc',
                    help='accessions, separated by space')
parser.add_argument('-f', type=str, nargs='+', default=[], metavar='file',
                    help='files of accession list, separated by space')
parser.add_argument('-o', type=str, default=None, metavar='file',
                    help='output file')
parser.add_argument('--db', type=str, choices=('protein', 'nucleotide', 'genome'), default='nucleotide',
                    help='ncbi db (%(default)s)')
parser.add_argument('--format', type=str, choices=('fasta', 'genbank', 'gff3'), default='fasta',
                    help='download format (%(default)s)')
# if more than 200 accessions, one should enlarge --batch to use POST method
parser.add_argument('--batch', type=int, default=1, metavar='int',
                    help='batch download mode, use POST method when > 1')
# parser.add_argument('--api', type=str,
#                     help='a file containing api key')
# parser.add_argument('--silent', action='store_true',
#                     help='supress error and proceed')
args = parser.parse_args()

acc_list = args.a
acc_files = args.f
out_file = open(args.o, 'w') if args.o else sys.stdout
db = args.db
dl_format = args.format
batch_size = args.batch
# api_file = args.api
# silent_mode = args.silent

Entrez.email = ''
Entrez.api_key = ''

if len(acc_list) + len(acc_files) == 0:
    parser.print_usage()
    print('-a or -l is required.')
    sys.exit(2)

for file in acc_files:
    with open(file) as f:
        acc_list += f.read().strip().split()

# remove duplicates and keep order
acc_list = rm_dup(acc_list)

print('Total {} accessions.'.format(len(acc_list)), file=sys.stderr)
print('Fetch from NCBI {} database.'.format(db), file=sys.stderr)
print('Download in {} format.'.format(dl_format), file=sys.stderr)
print('Batch size = {}.'.format(batch_size), file=sys.stderr)

if dl_format == 'genbank':
    dl_format = 'gb'
if not batch_size > 1:
    # 3 requests/sec
    handle = Entrez.efetch(
        id=acc_list,
        db=db,
        rettype=dl_format,
        retmode='text'
    )
    print(handle.read(), file=out_file)
else:
    print('Use POST method for batch download.', file=sys.stderr)
    search_results = Entrez.read(Entrez.epost(db, id=','.join(acc_list)))
    webenv = search_results['WebEnv']
    query_key = search_results['QueryKey']
    for start in range(0, len(acc_list), batch_size):
        end = min(len(acc_list), start + batch_size)
        handle = Entrez.efetch(
            db=db,
            rettype=dl_format,
            retmode='text',
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key
        )
        print(handle.read(), file=out_file)

print('Finish!', file=sys.stderr)
out_file.close()
