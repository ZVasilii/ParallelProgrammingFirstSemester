#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

struct args_t
{
	int worker;
	int number_of_thr;
	long number_to_count;
	long double* sum;
	pthread_mutex_t* mut;
};

void* series(void* data);

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		printf("Incorrect number of args, try: threads_number number_to_count\n");
		return 0;
	}
	int number_of_thr = atoi(argv[1]);
	long number_to_count = atol(argv[2]);
	pthread_t thr[number_of_thr];

	long double sum = 0.0;
	pthread_mutex_t mut;
	pthread_mutex_init(&mut, NULL); 
	struct args_t args[number_of_thr];
	for (int i = 0; i < number_of_thr; i++)
	{
		args[i].worker = i;
		args[i].number_of_thr = number_of_thr;
		args[i].number_to_count = number_to_count;
		args[i].sum = &sum;
		args[i].mut = &mut;
		pthread_create(thr + i, NULL, series, args + i);
	}
	for (int i = 0; i < number_of_thr; i++)
		pthread_join(thr[i], NULL);
	printf("Overall sum is %.9Lg\n", sum);
	return 0;
}

void* series(void* data)
{
	struct args_t* args = (struct args_t*) data;
	long double local_sum = 0.0;
	long from = args -> worker * (args -> number_to_count / args -> number_of_thr) + 1;
	long to = (args -> worker + 1) * (args -> number_to_count / args -> number_of_thr) + 1;
	if (to > (args -> number_to_count))
		to = args -> number_to_count;
	for (long i = from; i < to; i++)
		local_sum += (double)1 / i;
	printf("hello, world! I'm pthread with index #%d, my sum is %.9Lg\n", args -> worker, local_sum);
	pthread_mutex_lock(args -> mut);
	*(args -> sum) += local_sum;
	pthread_mutex_unlock(args -> mut);

	return NULL;

}