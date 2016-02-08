#!/bin/bash

echo "removing log files from /private/var/log/asl/*asl"
read -n 1 -r -p "sure to remove log files? [y/n] " reply
echo

if [[ $reply =~ ^[Yy]$ ]]; then
    sudo rm -rf /private/var/log/asl/*asl
    echo "done!"
fi
