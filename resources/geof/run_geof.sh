#!/bin/bash

/usr/local/bin/geof /home/ravi/git/tyler/resources/geof/createGLB.json --verbose \
--output_file=$3 \
--path_metadata=$4 \
--path_feature_input_file=$5 \
--min_x=${6} \
--min_y=${7} \
--min_z=${8} \
--max_x=${9} \
--max_y=${10} \
--max_z=${11} \
--cotypes=${12}