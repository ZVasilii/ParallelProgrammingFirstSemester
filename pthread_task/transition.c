#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

struct args_t
{
	int worker;
	int number_of_thr;
	pthread_mutex_t* mut;
};

int var = 2;

void* action(void* data);

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		printf("Incorrect number of args, try: threads_number \n");
		return 0;
	}
	int number_of_thr = atoi(argv[1]);
	pthread_t thr[number_of_thr];

	pthread_mutex_t mut;
	pthread_mutex_init(&mut, NULL); 
	struct args_t args[number_of_thr];
	for (int i = 0; i < number_of_thr; i++)
	{
		args[i].worker = i;
		args[i].number_of_thr = number_of_thr;
		args[i].mut = &mut;
		pthread_create(thr + i, NULL, action, args + i);
	}
	for (int i = 0; i < number_of_thr; i++)
		pthread_join(thr[i], NULL);
	printf("Final variable is %d\n", var);
	return 0;
}

void* action(void* data)
{
	struct args_t* args = (struct args_t*) data;
	pthread_mutex_lock(args -> mut);
	printf("hello, world! I'm pthread with index #%d/%d, variable is %d\n",\
		args -> worker, args -> number_of_thr, var);
	var *= 2;
	pthread_mutex_unlock(args -> mut);

	return NULL;

}