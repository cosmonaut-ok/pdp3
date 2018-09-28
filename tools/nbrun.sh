#!/bin/bash

SCRIPT_PATH=${1}

shift

cd $(dirname ${SCRIPT_PATH})

TMPFILE=`mktemp .nbrunXXXXXXXX.ipynb`

cp -f $(basename ${SCRIPT_PATH}) $TMPFILE

for i in "$@"; do
    VAR=$(echo $i | cut -d'=' -f1)
    VAL=$(echo $i | cut -d'=' -f2)
    # echo "var:" $VAR and "val:" $VAL
    sed -i "s|\"${VAR}.*=.*\\\\n\"|\"${VAR}=${VAL}\\\\n\"|g" $TMPFILE
done

jupyter nbconvert ${TMPFILE} --stdout --to python | sed '/get_ipython/d' > ${TMPFILE}.py

echo  >> ${TMPFILE}.py
echo 'input("Please press RETURN to exit ")' >> ${TMPFILE}.py

python ${TMPFILE}.py

rm -f $TMPFILE ${TMPFILE}.py
