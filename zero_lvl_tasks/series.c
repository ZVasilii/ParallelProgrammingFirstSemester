#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char** argv)
{
	long N = atoi(argv[1]);
	//printf("N = %ld\n", N);
	MPI_Init(&argc, &argv);
	int size = 0;
	int rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	long each_sum = N / size;

	if (rank == 0)
	{	
		//Reciever
		long double sum = 0;
		long double tmp = 0;
		//Case when N < N_proc
		if (each_sum == 0)
		{
			for (int i = 1; i < N + 1; i++)
				sum += (long double) 1 / i;
			printf("WARNING! N < N_proc! Final sum is %0.5Lg\n", sum);
			return 0;
		}

		//Counting for 0 worker
		for (int i = 1; i <= each_sum; i++)
			sum += (long double)  1 / i;
		printf("Hello, World! I'm worker number %d out of %d, counting [%ld, %ld] my local_sum is %0.5Lg\n",
			rank, size, 1l, each_sum, sum);

		//Recieving messages
		for (int i = 0; i < size - 1; i++)
		{
			int errno = MPI_Recv(&tmp, 1, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			assert(errno == MPI_SUCCESS); 
			sum += tmp;
		}
		printf("Final sum is %0.5Lg\n", sum);

	}
	else if (rank == size - 1)
	{
		//Counting for last worker
		long double local_sum = 0;
		for (int i = rank * each_sum + 1; i <= N; i++)
			local_sum +=  (long double)  1 / i;
		printf("Hello, World! I'm worker number %d out of %d, counting [%ld, %ld] my local_sum is %0.5Lg\n",
			rank, size, rank * each_sum + 1, N, local_sum);
		int errno = MPI_Send(&local_sum, 1, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
		assert(errno == MPI_SUCCESS); 
	}
	else
	{
		//Counting for other workers
		long double local_sum = 0;
		for (int i = rank * each_sum + 1; i <= (rank + 1) * each_sum; i++)
			local_sum += (long double)  1 / i;
		printf("Hello, World! I'm worker number %d out of %d, counting [%ld, %ld] my local_sum is %0.5Lg\n",
			rank, size, rank * each_sum + 1, (rank + 1) * each_sum, local_sum);
		int errno = MPI_Send(&local_sum, 1 , MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
		assert(errno == MPI_SUCCESS); 
	}
	MPI_Finalize();
}