#!/bin/bash

set -e

ME=`dirname $0`

## Set some colors
RED='\e[31m'
GREEN='\e[32m'
NC='\e[0m' # No Color

if [ -z $TESTDIR ]; then
    echo "There is no TESTDIR! Exiting."
    exit 1
fi

success="true"
for i in `ls $TESTDIR/tools/*.py`; do
    echo -n "Testing tool $(basename $i .py) "
    if echo -ne '\n' | xvfb-run -a $i $TESTDIR/parameters.xml > /dev/null; then
        printf "${GREEN}[ DONE ]${NC}\n"
    else
        printf "${RED}[ FAIL ]${NC}\n"
        success="false"
    fi
done

if [ "$success" == "false" ]; then
    echo
    echo "Test failed. Some tools are not working properly"
    exit 1
else
    echo
    echo "Test passed."
    echo
    exit 0
fi
