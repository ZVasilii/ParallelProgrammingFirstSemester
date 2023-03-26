#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

void merge(int *a, int *b, int l, int m, int r);
void mergeSort(int *a, int *b, int l, int r);

int main(int argc, char** argv) 
{
	double starttime_global = MPI_Wtime();
	//Initialize MPI 
	int world_rank;
	int world_size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	//Create the array
	int N = atoi(argv[1]);

	//Create subarray for each process
	int size = N/world_size;
	int *sub_array = (int*) malloc(size * sizeof(int));
	int *original_array = NULL;

	if(world_rank == 0) 
	{
		original_array = (int*) malloc(N * sizeof(int));
		
		srand(time(NULL));
		#ifdef PRINT
		printf("I'm worker number %d / %d , ", world_rank, world_size);
		printf("this is the unsorted array: ");
		for(int c = 0; c < N; c++) 
		{
			original_array[c] = rand() % N;
			printf("%d ", original_array[c]);
		}

		printf("\n");
		printf("\n");
		#endif
		
	}
	//Send each subarray to each process
	MPI_Scatter(original_array, size, MPI_INT, sub_array, size, MPI_INT, 0, MPI_COMM_WORLD);

	
	// Perform the mergesort on each process 
	int *tmp_array = (int*) malloc(size * sizeof(int));
	mergeSort(sub_array, tmp_array, 0, (size - 1));
	
	//Gather the sorted subarrays into one
	int *sorted = NULL;
	if(world_rank == 0) 
		sorted = (int*) malloc(N * sizeof(int));
	
	MPI_Gather(sub_array, size, MPI_INT, sorted, size, MPI_INT, 0, MPI_COMM_WORLD);
	
	// Make the final mergeSort call 
	if(world_rank == 0) 
	{
		double endtime_global = MPI_Wtime();
		int *other_array = (int*) malloc(N * sizeof(int));
		mergeSort(sorted, other_array, 0, (N - 1));
		
		#ifdef PRINT
		// Display the sorted array 
		printf("I'm worker number %d / %d ,(root) ", world_rank, world_size);
		printf("this is the sorted array: ");
		for(int c = 0; c < N; c++) 
			printf("%d ", sorted[c]);
			
		printf("\n");
		printf("\n");
		#endif
		printf("I'm root, time_whole_execution is %lg\n", endtime_global - starttime_global);
			
		// Clean up root 
		free(sorted);
		free(other_array);
		free(original_array);	
	}

	// Clean up rest 
	free(sub_array);
	free(tmp_array);
	// Finalize MPI 
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

/********** Merge Function **********/
void merge(int *a, int *b, int l, int m, int r) 
{
	int h, i, j, k;
	h = l;
	i = l;
	j = m + 1;
	
	while((h <= m) && (j <= r))
	{	
		if(a[h] <= a[j]) 
		{
			b[i] = a[h];
			h++;
		}
		else 
		{
			b[i] = a[j];
			j++;
		}	
		i++;
	}
		
	if(m < h) 
	{
		for(k = j; k <= r; k++) 
		{
			b[i] = a[k];
			i++;
		}		
	}
	else 
	{
		for(k = h; k <= m; k++) 
		{
			b[i] = a[k];
			i++;
		}	
	}
	for(k = l; k <= r; k++) 
		a[k] = b[k];	
}

/********** Recursive Merge Function **********/
void mergeSort(int *a, int *b, int l, int r) 
{
	if(l < r) 
	{
		
		int m = (l + r)/2;
		mergeSort(a, b, l, m);
		mergeSort(a, b, (m + 1), r);
		merge(a, b, l, m, r);
		
	}		
}