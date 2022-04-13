#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

long change_value(long val)
{
	return val * 2;
}

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int size = 0;
	int rank = 0;
	long window_buf = 0;
	MPI_Win buf;
	int errno = MPI_Win_create(&window_buf, sizeof(long), sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &buf);
	assert(errno == MPI_SUCCESS);
	MPI_Win_fence(0, buf);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		long int initial = 1;
		printf("Hello, World! I'm worker number %d out of %d, initial value is %ld\n", rank, size, initial);
		MPI_Put(&initial, 1, MPI_LONG, 1, 0, 1, MPI_LONG, buf);
	}
	MPI_Win_fence(0, buf);

  for (int i = 1; i < size; ++i) {
      int next_rank = i + 1;
      long new_value = 0;
      if (i == size - 1)
      {
          next_rank = 0;
      }

      if (rank == i) {
          new_value = change_value(window_buf);
          MPI_Put(&new_value, 1, MPI_LONG, next_rank, 0, 1, MPI_LONG, buf);
          printf("Hello, World! I'm worker number %d out of %d, new value is %ld\n", rank, size, new_value);
      }
      
      MPI_Win_fence(0, buf);
  }
  if (rank == 0)
		printf("Hello, World! I'm worker number %d out of %d, final value is %ld\n", rank, size, window_buf);

	MPI_Win_free(&buf);
	MPI_Finalize();
}