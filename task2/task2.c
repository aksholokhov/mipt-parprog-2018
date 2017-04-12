#include <unistd.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>

double recalculate(double* points, double k, double dt, double h) {
    return points[1] + (k*dt/(h*h))*(points[2] - 2*points[1] + points[0]);
}

double exact_solution(double x, double t, double k, double u0, double l) {
	double sum = 0;

	double PI = 3.141592653589793;
	for (int y = 0; y < 100500; y++) {
		double m = (double)y;
		sum += exp(-k*PI*PI*(2*m+1)*(2*m+1)*t/(l*l))*sin(PI*(2*m+1)*x/l)/(2*m+1);
	}

	return 4*u0/PI*sum;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int size, rank;


    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int points_num = atoi(argv[1]);
    double h = ((double)1)/(double)points_num;

    double T = 0.00000001, k = 1, dt = 2e-12, u0 = 0, u1 = 1;

    if (dt >= h*h/k) {
    	printf("Courant cond. fail:dt = %f, h*h/k = %f\n.", dt, h*h/k);
    	MPI_Finalize();   
    	return 0;
    }

    div_t p = div(points_num,size-1);
    int points_per_proc = p.quot;
    int steps = (int)(T/dt);

    MPI_Status status;
    MPI_Request request;

    double msg[2];

    if (rank != 0) {
        if (rank == size - 1) {
            points_per_proc += p.rem;
        } 

        double* points = (double*)malloc((points_per_proc+2)*sizeof(double));
        double* buf = (double*)malloc((points_per_proc+2)*sizeof(double));
        
        for (int i = 1; i < points_per_proc+1; i++) {
            points[i] = u1;
            buf[i] = u1;
        }

        
        for (int i = 0; i < steps; i++) {
        	if (rank != 1) {
        		MPI_Send(&points[1], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
            	MPI_Recv(&points[0], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
            } else {
            	points[0] = u0;
            }

            for (int j = 1; j < points_per_proc; j++) {
                buf[j] = recalculate(&points[j-1], k, dt, h);
            }
            
            if (rank != size - 1) {
	            MPI_Recv(&points[points_per_proc+1], 1, MPI_DOUBLE, (rank+1) % size, 0,
	                    MPI_COMM_WORLD, &status);
	            buf[points_per_proc] = recalculate(&points[points_per_proc - 1], k, dt, h);
	            MPI_Send(&points[points_per_proc], 1, MPI_DOUBLE, (rank+1) % size, 0,
	                    MPI_COMM_WORLD);
	        } else {
	        	buf[points_per_proc+1] = u0;	
	        	buf[points_per_proc] = recalculate(&points[points_per_proc - 1], k, dt, h);
	        }

            for (int j = 1; j < points_per_proc+1; j++) {
                points[j] = buf[j];
            }
        }

        MPI_Send(&points[1], points_per_proc, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        free(points);
        free(buf);
    }

    if (rank == 0) {

    	//Single-process solution:
    	double begin, end;

    	double* points = (double*)malloc((points_num+2)*sizeof(double));
        double* buf = (double*)malloc((points_num+2)*sizeof(double));
        
        points[0] = u0;
        points[points_num+1] = u0;

        for (int i = 1; i < points_num+1; i++) {
            points[i] = u1;
            buf[i] = u1;
        }

        // Multi-process solution
        begin = MPI_Wtime();

        for (int j = 1; j < size-1; j++) {
            MPI_Recv(points + 1 + (j-1)*points_per_proc, points_per_proc, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);    
        }
        MPI_Recv(points + 1 + (size-2)*points_per_proc, points_per_proc + p.rem, MPI_DOUBLE, size-1, 1, MPI_COMM_WORLD, &status);

        end = MPI_Wtime();

        points[points_num+1] = u0;
        printf("Multi-process time is: %f, with result below \n", end-begin);
        for (double j = 0; j<=1; j+=0.1) {
	        int point = (int)((points_num+1)*j);
	        if (j > 0.5) point++;
	            printf("%.3f ", points[point]);
	    }
	    printf("\n \n");

        begin = MPI_Wtime();

        for (int i = 1; i < points_num+1; i++) {
            points[i] = u1;
            buf[i] = u1;
        }
        
        for (int i = 0; i < steps; i++) {
            for (int j = 1; j < points_num+1; j++) {
                buf[j] = recalculate(&points[j-1], k, dt, h);
            }
            for (int j = 1; j < points_num+1; j++) {
                points[j] = buf[j];
            }
        }
        
        end = MPI_Wtime();
        printf("Single-process time is: %f, with result below \n", end-begin);

        for (double j = 0; j<=1; j+=0.1) {
            int point = (int)((points_num+1)*j);
            if (j > 0.5) point++;
            printf("%.3f ", points[point]);
        }
        printf("\n \n");
	
        printf("Exact solution is: \n");
        for (double j = 0; j<=1; j+=0.1) {
	        int point = (int)((points_num+1)*j);
	        if (j > 0.5) point++;
	            printf("%.3f ", exact_solution(j, T, k, u1, 1));
	    }
	    printf("\n");

    }

    MPI_Finalize();   
    return 0;
}