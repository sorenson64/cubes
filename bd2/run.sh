#!/bin/sh

ssh node-0 nohup  mpirun -np 112 ./pairmain 0 &
ssh node-1 nohup  mpirun -np 112 ./pairmain 1 &
ssh node-2 nohup  mpirun -np 112 ./pairmain 2 &
ssh node-3 nohup  mpirun -np 112 ./pairmain 3 &
ssh node-4 nohup  mpirun -np 112 ./pairmain 4 &
