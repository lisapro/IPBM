#!/bin/bash

mkdir -p build
rm -r build/*
CWD=`pwd`
cp data/* build &&
cd build && cmake $CWD/src -DFABM_BASE=$FABMDIR \
                           -DFABM_ERSEM_BASE=$CWD/../ERSEM \
-DFABM_NIVA_BASE=$CWD/../brom_niva_module -DCMAKE_BUILD_TYPE=Debug
                           #-DERSEM_USE_IRON=ON \
