#/bin/bash

echo "Parallel computing of pi number"
mpicc -O3  parallel_pi.c -lm -o  parallel_pi.out
mpirun -np 4 ./parallel_pi.out 10000