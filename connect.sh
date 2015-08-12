#!/bin/bash

# Hsin-Han Lee

case $1 in
"2.71" | "71")
    echo "connect to 140.112.2.71..."
    ssh cychen@140.112.2.71
;;
"206")
    echo "connect to 140.112.183.206..."
    ssh hsinhan@140.112.183.206
;;
*) # if no given parameter
    echo -e "Usage: $0 [alias of IP]"
    echo "Available IP: 140.112.2.71 (2.71 or 71)"
    echo "              140.112.183.206 (206)"
    exit 0 # 0 means success
;;
esac

# if style

# if [ "$1" == "2.71" ]; then
#     echo "connect to 140.112.2.71..."
#     ssh cychen@140.112.2.71

# elif [ "$1" == "206" ]; then
#     echo "connect to 140.112.183.206..."
#     ssh hsinhan@140.112.183.206

# fi
