#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int size = 0;
	int rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		long int i = 1;
		printf("Hello, World! I'm worker number %d out of %d, initial value is %ld\n", rank, size, i);
		//Sending
		int errno = MPI_Send(&i, 1, MPI_LONG, rank + 1, 0, MPI_COMM_WORLD);
		assert(errno == MPI_SUCCESS);	
		//Recieving
		errno = MPI_Recv(&i, 1, MPI_LONG, size - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		assert(errno == MPI_SUCCESS);
		printf("Finish! I'm worker number %d out of %d, final value is %ld\n", rank, size, i);
	}
	else if (rank == size - 1)
	{
		long int i = 0;
		//Sending
		int errno = MPI_Recv(&i, 1, MPI_LONG, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		assert(errno == MPI_SUCCESS);
		i = i * 2;
		printf("Hello, World! I'm worker number %d out of %d, new value is %ld, finishing!\n",
			rank, size, i);
		//Sending
		errno = MPI_Send(&i, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		assert(errno == MPI_SUCCESS);	

	}
	else
	{
		long int i = 0;
		//Recieving
		int errno = MPI_Recv(&i, 1, MPI_LONG, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		assert(errno == MPI_SUCCESS);
		i = i * 2;
		printf("Hello, World! I'm worker number %d out of %d, new value is %ld, sending further!\n",
			rank, size, i);
		//Sending
		errno = MPI_Send(&i, 1, MPI_LONG, rank + 1, 0, MPI_COMM_WORLD);
		assert(errno == MPI_SUCCESS);	

	}

	MPI_Finalize();
}