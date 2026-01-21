#!/bin/sh

ssh "node-$1" << EOF
cd ~/research/oc/cubes-bd2
nohup mpirun -np 112 ./parmain $1 &
EOF
