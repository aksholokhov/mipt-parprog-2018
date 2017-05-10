#include "unistd.h"
#include "stdlib.h"
#include "pthread.h"
#include "stdio.h"
#include "math.h"
#include <sys/time.h>


const double PI = 3.141592653589793;
const double x_max = PI, y_max = 1, z_max = PI;


void* worker_routine(void* params) {
	long iter_num = ((long*)params)[0];

	int rank = (int)((long*)params)[1];
	long* counter = (long*)malloc(sizeof(long));
	*counter = 0;
	double x, y, z, rnd_max = (double)RAND_MAX;
	unsigned int x_seed = 31*rank;
	for (int i = 0; i < iter_num; i++) {
		x = ((double)rand_r(&x_seed))/rnd_max*x_max;
		y = ((double)rand_r(&x_seed))/rnd_max*y_max;
		z = ((double)rand_r(&x_seed))/rnd_max*z_max;
		if ((x > 0 ) & (x < PI)) {
			if (y < sin(x)) {
				if (z < x*y) {
					(*counter)++;
				}
			}
		}
	}

	pthread_exit((void*)counter);
}

int main(int argc, char const *argv[])
{
	long NUM_THREADS = atol(argv[1]);
	long NUM_POINTS = atol(argv[2]);
	ldiv_t ppt = ldiv(NUM_POINTS, NUM_THREADS);

	void* result = malloc(sizeof(long));
	long* args[NUM_THREADS+1];
	
	pthread_t thread[NUM_THREADS+1];

	//clock_t begin = clock();
	struct timespec begin, end;

	clock_gettime(CLOCK_REALTIME, &begin);

	for (int i = 1; i <= NUM_THREADS; i++) {
		args[i] = (long*)malloc(sizeof(long)*2);
		args[i][0] = i == NUM_THREADS ? ppt.quot + ppt.rem : ppt.quot;
		args[i][1] = i;
		pthread_create(&thread[i], NULL, worker_routine, (void*)(args[i]));
	}

	long sum = 0;

	for (int i = 1; i <= NUM_THREADS; i++) {
		pthread_join(thread[i], &result);
		sum += *(long*)result;
	}

	clock_gettime(CLOCK_REALTIME, &end);

	double MP_perf = end.tv_sec - begin.tv_sec;
	MP_perf += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

	//printf("MP result is: %.5f\n", ((double)sum)/(double)NUM_POINTS*(x_max*y_max*z_max));

	clock_gettime(CLOCK_REALTIME, &begin);

	args[0] = (long*)malloc(sizeof(long)*2);
	args[0][0] = NUM_POINTS;
	args[0][1] = 999;
	pthread_create(&thread[0], NULL, worker_routine, (void*)(args[0]));
	pthread_join(thread[0], &result);

	clock_gettime(CLOCK_REALTIME, &end);

	double SP_perf = end.tv_sec - begin.tv_sec;
	SP_perf += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

	//printf("SP result is: %.5f\n",	(double)(((long*)result)[0])/(double)NUM_POINTS*(x_max*y_max*z_max));	

	//printf("Acceleration is: %.4f\n", (double)(SP_perf/MP_perf));
	
	printf("%ld,%ld,%.3f\n", NUM_THREADS, NUM_POINTS, SP_perf/MP_perf);
	return 0;
}