#/bin/bash

echo "Parallel computing of exponent"
mpicc -O3  parallel_exponent.c -o parallel_exponent.out
mpirun -np 4 ./parallel_exponent.out 10000