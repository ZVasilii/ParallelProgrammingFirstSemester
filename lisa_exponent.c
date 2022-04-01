#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

long double fact(long double x)
{
	if (x < 1)
		return 1;
	else
		return x * fact(x - 1);
}

long double inversed_fact(long double x)
{
	return 1.0 / fact(x);
}

int main(int argc, char** argv)
{
	double starttime_global = MPI_Wtime();
	long N = atoi(argv[1]);
	//printf("N = %ld\n", N);
	int size, rank = 0;
	int root = 0;
	long double result = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	long iterations = N / size + 1;
	long target = 0;
	long double local_sum = 0;

	double starttime_local = MPI_Wtime();
	for (int i = 0; i < iterations; ++i)
	{
		target = i * size + rank;
		if (target <= N)
		{
			//printf("I'm worker number %d / %d, sum is %ld\n", rank, size, target);
			local_sum += inversed_fact(target);
		}
	}
	MPI_Reduce(&local_sum, &result, 1, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	double endtime_local = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == root)
	{
		printf("I'm root, result is %0.30Lg\n", result);
	}

	printf("I'm worker number %d / %d, local time is %lg\n",
				 rank, size, endtime_local - starttime_local);
	double endtime_global = MPI_Wtime();
	if (rank == root)
	{
		printf("I'm root, time_whole_execution is %lg\n", endtime_global - starttime_global);
	}
	MPI_Finalize();
	return 0;
}	