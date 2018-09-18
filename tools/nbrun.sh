#!/bin/bash

cd $(dirname ${1})

TMPFILE=`mktemp .nbrunXXXXXXXX.py`
jupyter nbconvert $(basename ${1}) --stdout --to python | sed '/get_ipython/d' > $TMPFILE

echo  >> $TMPFILE
echo 'input("Please press RETURN to exit ")' >> $TMPFILE

python $TMPFILE

rm -f $TMPFILE
