#!/bin/bash

mpirun -np 128 -hostfile  /etc/hostfile -x OMP_NUM_THREADS=1   -x LD_LIBRARY_PATH  ./qr_test.x -p 8 -q 16 -r 40000 80000 2000 100 -nb 100
