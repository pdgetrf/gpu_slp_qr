#!/bin/bash

mpiexec -np 48 -genv OMP_NUM_THREADS 1 ./qr_test.x -p 6 -q 8 -r 2000 20000 2000 -nb 100
