#!/usr/bin/env python3
'''
ont-fast5-api v3.0.1
https://github.com/nanoporetech/ont_fast5_api
'''

import sys
import argparse
import json
import datetime
from ont_fast5_api.conversion_tools.conversion_utils import yield_fast5_files
from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.multi_fast5 import MultiFast5File
from ont_fast5_api.fast5_interface import is_multi_read

def return_bool(bool):
    if bool:
        return 1
    else:
        return 0

parser = argparse.ArgumentParser(description='Check nanopore fast5 file metadata. Output in json format.')
parser.add_argument('fast5',
                    help='input file or folder (read first file only)')
parser.add_argument('--fast5-prefix', action='store_true',
                    help='{flow_cell_id}_{run_id}')
parser.add_argument('--dir-prefix', action='store_true',
                    help='{yyyymmdd}_{hhmm}_{device_id}_{flow_cell_id}_{protocol_run_id:1:8} (timestamp might be skewed.)')
parser.add_argument('--is-multi', action='store_true',
                    help='is multi fast5')
parser.add_argument('--is-vbz', action='store_true',
                    help='is vbz compressed')
parser.add_argument('--read-id', action='store_true',
                    help='print read ids in fast5 file')
args = parser.parse_args()

in_f5 = next(yield_fast5_files(args.fast5, True, True))  # grab one file only, maybe the oldest one
print('INPUT: {}'.format(in_f5), file=sys.stderr)

MULTI = False
VBZ = False
READ_ID = []
RET = 0
n = 0

if is_multi_read(in_f5):
    MULTI = True
    with MultiFast5File(in_f5, 'r') as h:
        read = next(h.get_reads())
        # print(read.get_run_id().decode())
        # print(read.get_read_id())
        TRACK = read.get_tracking_id()
        CHANNEL = read.get_channel_info()
        CONTEXT = read.get_context_tags()
        FC = TRACK['flow_cell_id']
        RUN = TRACK['run_id']
        COMPRESS = read.raw_compression_filters
        if '32020' in COMPRESS.keys():
            VBZ = True
        WARNED_FC = False
        WARNED_RUN = False
        for r in h.get_reads():  # alway iter from first one even I have used next()?
            READ_ID.append(r.get_read_id())
            t = r.get_tracking_id()
            if t['flow_cell_id'] != FC and not WARNED_FC:
                print('WARN: Contain multiple flow-cells.', file=sys.stderr)
                WARNED_FC = True
                RET = 99
            if t['run_id'] != RUN and not WARNED_RUN:
                print('WARN: Contain multiple runs.', file=sys.stderr)
                WARNED_RUN = True
                RET = 99
            n += 1
else:
    with Fast5File(in_f5, 'r') as read:
        TRACK = read.get_tracking_id()
        CHANNEL = read.get_channel_info()
        CONTEXT = read.get_context_tags()
        COMPRESS = read.raw_compression_filters
        if '32020' in COMPRESS.keys():
            VBZ = True
        READ_ID.append(read.get_read_id())
        n = 1

METADATA = {
    'is_multi': MULTI,
    'read_per_fast5': n,
    'is_vbz': VBZ,
    'compression': COMPRESS,
    **TRACK,
    **CHANNEL,
    **CONTEXT
}

if args.fast5_prefix:
    print('{}_{}'.format(METADATA['flow_cell_id'], METADATA['run_id']))
elif args.dir_prefix:
    t = datetime.datetime.strptime(METADATA['exp_start_time'], '%Y-%m-%dT%H:%M:%SZ')
    print('{}_{}_{}_{}'.format(
        datetime.datetime.strftime(t, '%Y%m%d_%H%M'),
        METADATA['device_id'],
        METADATA['flow_cell_id'],
        METADATA['protocol_run_id'][0:8]
    ))
elif args.is_multi:
    print(return_bool(METADATA['is_multi']))
elif args.is_vbz:
    print(return_bool(METADATA['is_vbz']))
elif args.read_id:
    print('\n'.join(READ_ID))
else:
    print(json.dumps(METADATA, indent=2))

sys.exit(RET)
