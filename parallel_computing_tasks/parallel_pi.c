#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

long double inversed_n_squared(long double x)
{
	if (x > 0)
		return 1.0 / (x * x);
	else return 0;
}

int main(int argc, char** argv)
{
	double starttime_global = MPI_Wtime();
	long N = atoi(argv[1]);
	//printf("N = %ld\n", N);
	int size, rank = 0;
	int root = 0;  //Collector process
	long double result = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	long iterations = N / size + 1; //Number of terms for each process
	long target = 0; //Each term
	long double local_sum = 0;  //Sum of terms for each process

	double starttime_local = MPI_Wtime();
	for (int i = 0; i < iterations; ++i)
	{
		target = i * size + rank;
		if (target <= N)
		{
			//printf("I'm worker number %d / %d, sum is %ld\n", rank, size, target);
			local_sum += inversed_n_squared(target);
		}
	}
	//Sending all local sums to process root into result variable (Using MPI_SUM)
	MPI_Reduce(&local_sum, &result, 1, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	double endtime_local = MPI_Wtime();

	//Waiting for all processes to send local_sums
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == root)
	{
		printf("I'm root, result is %0.5Lg\n", sqrtl(result * 6)); // sum(1 / n^2) = sqrt(pi)/6
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