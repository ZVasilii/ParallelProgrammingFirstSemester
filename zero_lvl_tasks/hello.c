#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int size = 0;
	int rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("Hello, World! I'm worker number %d out of %d\n", rank, size);
	MPI_Finalize();
}