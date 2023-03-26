#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double func(double x, double t);
double x_row_init(double x);
double t_row_init(double t);

double** create_matrix(size_t rows,  size_t cols);
void remove_matrix(double** m);
void print_matrix(double ** m, size_t rows,  size_t cols);
void initialize_rows(double** m,  size_t rows,  size_t cols, double tau, double h, double max_T, double max_X);

int main(int argc, char** argv)
{

	const double tau = strtod(argv[1], NULL);
	const double h = strtod(argv[2], NULL);
	const double max_T = strtod(argv[3], NULL);
	const double max_X = strtod(argv[4], NULL);
	size_t rows = round(max_T / tau);
	size_t cols = round(max_X / h);
	//u_matr[t][x], 0 < t < rows; o < x < cols


	MPI_Init(&argc, &argv);
	int size = 0;
	int rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double** u_matr = create_matrix(rows, cols);
	initialize_rows(u_matr, rows, cols, tau, h, max_T, max_X);

	size_t each_work = cols / size;

	#ifdef TIME
	double starttime_global = MPI_Wtime();
	#endif
	
	for(size_t row = 0; row < rows - 1; row++)
	{
		if (rank != size - 1)
		{
			if (rank != 0)
			{
				int errno = MPI_Recv(&u_matr[row + 1][rank * each_work], 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				assert(errno == MPI_SUCCESS);
			}
			for (size_t j = rank * each_work + 1; j < (rank + 1) * each_work + 1; j++)
			{
					u_matr[row + 1][j] = 2 * tau * h / (h + tau) * (func((row + 0.5) * tau / max_T, (j + 0.5) * h / max_X) - \
						(u_matr[row + 1][j - 1] - u_matr[row][j - 1] - u_matr[row][j])/ 2 / tau - (u_matr[row][j] - u_matr[row + 1][j - 1] - u_matr[row][j - 1])  / 2 / h);
			}
	
			int errno = MPI_Send(&u_matr[row + 1][(rank + 1) * each_work], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			assert(errno == MPI_SUCCESS);

			//MPI_Gather(&u_matr[row + 1][rank * each_work + 1], each_work, MPI_DOUBLE, &u_matr[row + 1][1], cols - 1, MPI_DOUBLE, size - 1, MPI_COMM_WORLD);
																																																																																													#ifdef PRINT
			errno = MPI_Send(&u_matr[row + 1][rank * each_work + 1], each_work, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD);
			assert(errno == MPI_SUCCESS);
																																																																																													#endif
		}

		if (rank == size - 1)
		{

			if (size != 1)
			{
				int errno = MPI_Recv(&u_matr[row + 1][rank * each_work], 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				assert(errno == MPI_SUCCESS);
			}
			for (size_t j = rank * each_work + 1; j < cols - 1; j++)
			{
					u_matr[row + 1][j] = 2 * tau * h / (h + tau) * (func((row + 0.5) * tau / max_T, (j + 0.5) * h / max_X) - \
						(u_matr[row + 1][j - 1] - u_matr[row][j - 1] - u_matr[row][j])/ 2 / tau - (u_matr[row][j] - u_matr[row + 1][j - 1] - u_matr[row][j - 1])  / 2 / h);

			}
			if (size != 1)
			{
				int errno = MPI_Send(&u_matr[row + 1][rank * each_work + 1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
				assert(errno == MPI_SUCCESS);
			} 
			//Upper left corner
			//u_matr[row + 1][cols - 1] = (func(row * tau / max_T, (cols - 1) * h / max_X) - (u_matr[row][cols - 1] - u_matr[row][cols - 2]) / h) * tau + u_matr[row][cols - 1];
																																																																																																																																																	#ifdef PRINT
			for (size_t i = 0; i < (size_t) size - 1; i++)
			{
				int errno = MPI_Recv(&u_matr[row + 1][i * each_work + 1], each_work, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				assert(errno == MPI_SUCCESS);
			}
																																																																																																																																																	#endif
			//MPI_Gather(&u_matr[row + 1][rank * each_work + 1], each_work, MPI_DOUBLE, &u_matr[row + 1][1], cols - 1, MPI_DOUBLE, size - 1, MPI_COMM_WORLD);	

		}
	}

	#ifdef TIME
	double endtime_global = MPI_Wtime();
	#endif

	#ifdef PRINT
	if (rank == size - 1)
	{
		print_matrix(u_matr, rows, cols);
		printf("\n");
	}
	#endif

	#ifdef TIME
	if (rank == size - 1)
	{
		printf("Time elapsed = %lgs\n", endtime_global - starttime_global);
	}
	#endif
	

	remove_matrix(u_matr);
	
	MPI_Finalize();
	return 0;

}


double func(double x, double t)
{
	return sin(x) * cos(t);
}

double t_row_init(double t)
{
	return sin(t);
}

double x_row_init(double x)
{
	return x * x;
}

double** create_matrix(size_t rows, size_t cols)
{
	double** m = (double**) calloc(rows, sizeof(int*));
	assert(m && "Create_matrix");
	m[0] = (double*) calloc(rows*cols, sizeof(double));
	for (size_t i = 1; i < rows; ++i)
		m[i] = m[i - 1] + cols;
	return m;
}

void remove_matrix(double** m)
{
	assert(m && "Remove matrix");
	free(m[0]);
	free(m);
}

void print_matrix(double ** m,  size_t rows, size_t cols)
{
	assert(m && "print matrix");
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
			printf("%.5lg ", m[rows - i - 1][j]);
		if (i != rows - 1)
			printf("\n");
	}
}

void initialize_rows(double** m,  size_t rows,  size_t cols, double tau, double h, double max_T, double max_X)
{
	assert(m && "initialize matrix");
	for (size_t i = 0; i < rows;  ++i)
		m[rows - 1 - i][0] = t_row_init(i * tau / max_T);
	for (size_t j = 0; j < cols; ++j)
		m[0][j] = x_row_init(j * h / max_X);
}