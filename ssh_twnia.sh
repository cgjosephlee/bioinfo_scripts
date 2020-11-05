#!/bin/bash

source ~/twnia_key.sh
export SSHPASS="$pw$(python -c 'import pyotp,os; print(pyotp.TOTP(os.environ["otpkey"]).now())')"

if [[ $1 == 'twnia' ]]; then
    sshpass -v -e ssh twnia_2
elif [[ $1 == 'twcc' ]]; then
    sshpass -v -e ssh twcc_1
fi

unset pw
unset otpkey
unset SSHPASS
