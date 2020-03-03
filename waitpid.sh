#!/bin/bash
# https://stackoverflow.com/a/19396161/7859425
while [ -e /proc/$1 ]; do sleep ${2:-60}; done
