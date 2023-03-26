#include <pthread.h>
#include <stdio.h>

int main()
{
	pthread_attr_t tattr;
	int ret = -1;
	int scope = -1;
	size_t stack_size = -1;
	int policy = -1;
	int inheritsched = -1;
	int detachstate = -1;
	void* stackaddr = NULL;

	ret = pthread_attr_init(&tattr);
	if (ret != 0)
		perror("Init");

	pthread_attr_getscope(&tattr, &scope);
	printf("*****\n");
	printf("Scope = %d\n", scope);
	printf("Process == %d, system == %d\n", \
		PTHREAD_SCOPE_PROCESS, PTHREAD_SCOPE_SYSTEM);
	printf("*****\n");

	pthread_attr_getdetachstate(&tattr, &detachstate);
	printf("*****\n");
	printf("Detachstate = %d\n", detachstate);
	printf("DETACHED == %d, JOINABLE == %d\n", \
		 PTHREAD_CREATE_DETACHED,  PTHREAD_CREATE_JOINABLE);
	printf("*****\n");

	pthread_attr_getstack(&tattr, &stackaddr, &stack_size);
	printf("*****\n");
	printf("Stack_size = %lu bytes\n", stack_size);
	printf("Stack_addr = %p\n", stackaddr);
	printf("*****\n");

	pthread_attr_getinheritsched(&tattr, &inheritsched);
	printf("*****\n");
	printf("Inheritsched = %d\n", inheritsched);
	printf("INHERIT_SCHED == %d, EXPLICIT_SCHED == %d\n", \
		 PTHREAD_INHERIT_SCHED,  PTHREAD_EXPLICIT_SCHED);
	printf("*****\n");

	pthread_attr_getschedpolicy(&tattr, &policy);
	printf("*****\n");
	printf("Policy = %d\n", policy);
	printf("FIFO == %d, RR == %d, OTHER == %d\n", \
		 SCHED_FIFO,  SCHED_RR, SCHED_OTHER);
	printf("*****\n");


	ret = pthread_attr_destroy(&tattr); 
	if (ret != 0)
		perror("Destroy");
	return 0;
}
