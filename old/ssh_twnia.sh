#!/bin/bash
set -e

source ~/twnia_key.sh
export SSHPASS="$pw$(python -c 'import pyotp,os; print(pyotp.TOTP(os.environ["otpkey"]).now())')"

sshpass -v -e ssh $1

unset pw
unset otpkey
unset SSHPASS
