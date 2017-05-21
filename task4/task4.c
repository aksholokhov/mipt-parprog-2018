#include "unistd.h"
#include "stdlib.h"
#include "pthread.h"
#include "stdio.h"
#include "semaphore.h"
#include "math.h"
#include <fcntl.h>
#include <sys/time.h>

double c = 1, h = 0.1, dt = 0.05, T = 2, k;
int NUM_THREADS, NUM_POINTS;
double* points;
sem_t* sem[24];

void* worker_routine(void* params) {
	int rank = *(int*)params;
	div_t ppt = div(NUM_POINTS, NUM_THREADS);
	int N = ppt.quot;
	if (rank == NUM_THREADS - 1) {
		N += ppt.rem;
	}

	int start_point = rank*ppt.quot + 1;
	double t = 0;
	double buf[N];
	//printf("%d: %d, %d \n", rank, start_point, N);
	while (t < T) {
		if (rank != 0) {
			sem_wait(sem[rank]);
		}
		printf("B: %d-%.1f: %.3f \n", rank, t, points[start_point]);

		for (int i = start_point; i < start_point + N-1; i++) {
			buf[i-start_point] = k*(points[i-1] - points[i]) + points[i];
		}
		sem_post(sem[rank+12]);

		sem_wait(sem[(rank+1)+12]);
		for (int i = start_point; i < start_point + N-1; i++) {
			points[i] = buf[i - start_point];
		}
		printf("E: %d-%.1f: %.3f \n", rank, t, points[start_point + N]);
		sem_post(sem[(rank+1)]);

		t += dt;
	}

	return NULL;
}

double g(double x) {
	if ((x < 0) | (x > 2)) {
		return 0;
	}
	return x*(2-x);
}

int main(int argc, char** argv) {
	NUM_THREADS = atol(argv[1]);
	NUM_POINTS = (int)(10.0/h) + 1;
	k = c*dt/h;

	points = (double*)malloc(sizeof(double)*NUM_POINTS);
	points[0] = 0;
	for (int i = 1; i < NUM_POINTS; i++) {
		double x = (i-1)*10.0/NUM_POINTS;
		points[i] = g(x);
	}

	pthread_t thread[NUM_THREADS+1];
	int* args[NUM_THREADS+1];

	for (int i = 0; i < NUM_THREADS; i++) {
		sem[i] = sem_open("/my_sem", O_CREAT, 0777, 0);
		sem[i+12] = sem_open("/my_sem", O_CREAT, 0777, 1);
		args[i] = (int*)malloc(sizeof(int));
		args[i][0] = i;
		pthread_create(&thread[i], NULL, worker_routine, (void*)(args[i]));
	}

	struct timespec begin, end;
	int result;

	for (int i = 0; i < NUM_THREADS; i++) {
		pthread_join(thread[i], &result);
	}

	for (int i = 0; i < NUM_POINTS; i++) {
		printf("%.3f ", points[i]);
	}
	printf("\n");

	for (int i = 0; i < NUM_THREADS; i++) {
		sem_close(sem[i]);
		sem_close(sem[i+12]);
	}

	return 0;
}
