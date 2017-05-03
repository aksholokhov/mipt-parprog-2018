#include "unistd.h"
#include "stdlib.h"
#include "pthread.h"
#include "stdio.h"
#include "math.h"
#include <sys/time.h>


const double PI = 3.141592653589793;
const double x_max = PI, y_max = 1, z_max = PI;


void* worker_routine(void* params) {
	int iter_num = ((int*)params)[0];
	int rank = ((int*)params)[1];
	int* counter = (int*)malloc(sizeof(int));
	double x, y, z, rnd_max = (double)RAND_MAX;
	unsigned int seed = time(NULL)*rank;
	for (int i = 0; i < iter_num; i++) {
		x = ((double)rand_r(&seed))/rnd_max*x_max;
		y = ((double)rand_r(&seed))/rnd_max*y_max;
		z = ((double)rand_r(&seed))/rnd_max*z_max;
		if ((x > 0 ) & (x < PI)) {
			if (y < sin(x)) {
				if (z < x*y) {
					(*counter)++;
				}
			}
		}
	}
	return (void*)counter;
}

int main(int argc, char const *argv[])
{
	int NUM_THREADS = atoi(argv[1]);
	int NUM_POINTS = atoi(argv[2]);
	div_t ppt = div(NUM_POINTS, NUM_THREADS);

	void* result;
	int* args[NUM_THREADS+1];
	
	pthread_t thread[NUM_THREADS];

	struct timeval stop, start;
	gettimeofday(&start, NULL);

	for (int i = 1; i <= NUM_THREADS; i++) {
		args[i] = (int*)malloc(sizeof(int)*2);
		args[i][0] = i == NUM_THREADS - 1 ? ppt.quot + ppt.rem : ppt.quot;
		args[i][1] = i+1;
		pthread_create(&thread[i], NULL, worker_routine, (void*)(args[i]));
	}

	double sum = 0;
	for (int i = 1; i <= NUM_THREADS; i++) {
		pthread_join(thread[i], &result);
		free(args[i]);
		sum += *(int*)result;
	}
	gettimeofday(&stop, NULL);
	double MP_perf = abs(stop.tv_usec - start.tv_usec);

	printf("MP result is: %.5f\n", sum/(double)NUM_POINTS*(x_max*y_max*z_max));

	gettimeofday(&start, NULL);
	args[0] = (int*)malloc(sizeof(int)*2);
	args[0][0] = NUM_POINTS;
	args[0][1] = 99;
	result = worker_routine(args[0]);
	gettimeofday(&stop, NULL);

	double SP_perf = abs(stop.tv_usec - start.tv_usec);

	printf("SP result is: %.5f\n", (double)(((int*)result)[0])/(double)NUM_POINTS*(x_max*y_max*z_max));	

	printf("Acceleration is: %.4f\n", SP_perf/MP_perf);
	return 0;
}