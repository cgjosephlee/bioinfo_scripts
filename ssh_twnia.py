#!/usr/bin/env python3

import sys
import os
import json
import pyotp
import pexpect

with open(os.path.expanduser('~/twnia_key.json')) as f:
    config = json.load(f)

HOST = sys.argv[1]
PW = config['pw']
OTP = pyotp.TOTP(config['otpkey']).now()

child = pexpect.spawn('ssh {}'.format(HOST))
child.expect('assword:')
child.sendline(PW)
child.expect('OTP:')
child.sendline(OTP)
child.interact()
