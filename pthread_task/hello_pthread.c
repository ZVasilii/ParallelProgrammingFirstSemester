#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void* hello(void* idx_);

int main(int argc, char** argv)
{
	int number_of_thr = atoi(argv[1]);
	pthread_t thr[number_of_thr];
	int index[number_of_thr];
	memset(index, 0, sizeof(index));
	for (int i = 0; i < number_of_thr; i++)
	{
		index[i] = i;
		pthread_create(thr + i, NULL, hello, index + i);
	}
	for (int i = 0; i < number_of_thr; i++)
		pthread_join(thr[i], NULL);
	return 0;
}

void* hello(void* idx_)
{
	int* idx = (int*) idx_;
	printf("hello, world! I'm pthread with index #%d\n", *idx);
	return NULL;
}