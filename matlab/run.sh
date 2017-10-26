#!/bin/bash

cd $(dirname $0)
CURRENT=${PWD}
DATA_DIR=$1
MOVIE_DIR=$2

matlab -nodisplay -r "try; cd('${CURRENT}'); rho_movie_create_light3('${DATA_DIR}', '${MOVIE_DIR}',100,254,2046,[0 1],[0 1],[-1e-7 0]); catch; end; quit"