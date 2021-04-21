#!/usr/bin/env python3

import sys
import os
import json
import pyotp
import pexpect
import struct, fcntl, termios, signal

# https://stackoverflow.com/questions/44369112/pexpect-and-terminal-resizing
def get_terminal_size():
    s = struct.pack("HHHH", 0, 0, 0, 0)
    a = struct.unpack('hhhh', fcntl.ioctl(sys.stdout.fileno(), termios.TIOCGWINSZ, s))
    return a[0], a[1]

def sigwinch_passthrough(sig, data):
    global p
    if not p.closed:
        p.setwinsize(*get_terminal_size())

with open(os.path.expanduser('~/twnia_key.json')) as f:
    config = json.load(f)

HOST = sys.argv[1]
PW = config['pw']
OTP = pyotp.TOTP(config['otpkey']).now()

p = pexpect.spawn(f'ssh {HOST}')
p.expect('assword:')
p.sendline(PW)
p.expect('OTP:')
p.sendline(OTP)
p.setwinsize(*get_terminal_size())
signal.signal(signal.SIGWINCH, sigwinch_passthrough)
p.interact()
