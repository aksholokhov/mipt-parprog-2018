#include "unistd.h"
#include "stdlib.h"
#include "pthread.h"
#include "stdio.h"
#include "semaphore.h"
#include "math.h"
#include <fcntl.h>
#include <sys/time.h>

double c = 1, h = 0.1, dt = 0.05, T = 4, k;
int NUM_THREADS, NUM_POINTS;
double* points;
double* buf;
sem_t* sem[24];

void* worker_routine(void* params) {
	int rank = *(int*)params;
	div_t p = div(NUM_POINTS, NUM_THREADS);
	int num_pts = rank == NUM_THREADS - 1 ? p.quot + p.rem : p.quot;
	int start_point = rank * p.quot;
	double t = 0;
	sem_post(sem[rank]);
	while (t < T) {
		if (rank != 0) 
			sem_wait(sem[rank + 12]);
		for (int i = start_point + 1; i < start_point + num_pts; i++) {
			buf[i] = (points[i] + points[i-1])/2;
		}
		if (rank != NUM_THREADS - 1) sem_wait(sem[rank+1]);
		buf[start_point + num_pts] = (points[start_point + num_pts] 
			+ points[start_point + num_pts-1])/2;
		if (rank != 0) sem_post(sem[rank]);
		for (int i = start_point + 1; i <= start_point + num_pts; i++) {
			points[i] = buf[i];
		}
		if (rank != NUM_THREADS - 1) sem_post(sem[rank+1+12]);
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

	points = (double*)malloc(sizeof(double)*NUM_POINTS + 1);
	buf = (double*)malloc(sizeof(double)*NUM_POINTS + 1);
	points[0] = 0;
	for (int i = 1; i < NUM_POINTS; i++) {
		double x = (i-1)*10.0/NUM_POINTS;
		points[i] = g(x);
		buf[i] = 0;
	}

	// single thread
	double t = 0;
	while (t < T) {
		for (int i = 1; i < NUM_POINTS; i++) {
			buf[i] = (points[i-1] + points[i])/2;
		}
		for (int i = 1; i < NUM_POINTS; i++) {
			points[i] = buf[i];
		}
		t += dt;
	}

	for (int i = 0; i < NUM_POINTS; i+=5) {
		printf("%.3f ", points[i]);
		if (i == 50) printf("\n");
	}
	printf("\n");

	//multi thread

	points[0] = 0;
	for (int i = 1; i < NUM_POINTS; i++) {
		double x = (i-1)*10.0/NUM_POINTS;
		points[i] = g(x);
		buf[i] = 0;
	}

	pthread_t thread[NUM_THREADS+1];
	int* args[NUM_THREADS+1];

	for (int i = 0; i < NUM_THREADS; i++) {
		char name[10] = "/my_sem";
		name[7] = i / 10 + '0';
		name[8] = i % 10 + '0';
		name[9] = 0;
		sem[i] = sem_open(name, O_CREAT, 0777, 0);
		name[7] = (i+12) / 10 + '0';
		name[8] = (i+12) % 10 + '0';
		name[9] = 0;
		sem[i+12] = sem_open(name, O_CREAT, 0777, 0);
		args[i] = (int*)malloc(sizeof(int));
		args[i][0] = i;
		pthread_create(&thread[i], NULL, worker_routine, (void*)(args[i]));
	}

	struct timespec begin, end;
	int result;

	for (int i = 0; i < NUM_THREADS; i++) {
		pthread_join(thread[i], &result);
	}

	for (int i = 0; i < NUM_POINTS; i+=5) {
		printf("%.3f ", points[i]);
		if (i == 50) printf("\n");
	}
	printf("\n");

	for (int i = 0; i < NUM_THREADS; i++) {
		sem_close(sem[i]);
		sem_close(sem[i+12]);
	}

	return 0;
}
