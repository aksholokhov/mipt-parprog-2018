#include <unistd.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

double recalculate(double* points, double k, double dt, double h) {
    return points[1] + (k*dt/(h*h))*(points[2] - 2*points[1] + points[0]);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    float T = 0.01, h = 0.02, k = 1, dt = 0.0002, step_log = 0.1, u0 = 0, u1 = 1;

    int points_num = (int)(1/h);
    int points_per_proc = points_num/(size-1);
    int steps = (int)(T/dt);

    MPI_Status status;

    if (rank != 0) {
        double* points = (double*)malloc((points_per_proc+2)*sizeof(double));
        double* buf = (double*)malloc((points_per_proc+2)*sizeof(double));
        
        for (int i = 0; i < points_per_proc; i++) {
            points[i] = u1;
            buf[i] = u1;
        }

        for (int i = 0; i < steps; i++) {
            MPI_Send(&points[0], 2, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
            MPI_Recv(&points[0], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);

            for (int j = 1; j < points_per_proc - 2; j++) {
                buf[j] = recalculate(&points[j-1], k, dt, h);
            }
            
            MPI_Recv(&points[points_per_proc-1], 2, MPI_DOUBLE, (rank+1) % size, 0,
                    MPI_COMM_WORLD, &status);
            buf[points_per_proc-2] = recalculate(&points[points_per_proc - 3], k, dt, h);
            buf[points_per_proc-1] = recalculate(&points[points_per_proc - 2], k, dt, h);

            for (int j = 1; j < points_per_proc; j++) {
                points[j] = buf[j];
            }

            if (rank != size - 1) MPI_Send(&points[points_per_proc-1], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                        
            MPI_Send(&points[points_per_proc-1], 1, MPI_DOUBLE, (rank+1) % size, 0,
                    MPI_COMM_WORLD);
        }
        free(points);
        free(buf);
    }

    if (rank == 0) {
        double buf[2], u[2];
        u[0] = u0;
        u[1] = u0;
        for (int i = 0; i < steps; i++) {
            MPI_Recv(buf, 2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
            MPI_Send(u, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            printf("%.2f ", u0);
            for (int j = 1; j < size - 1; j++) {
                MPI_Recv(buf, 1, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
                printf("%.2f ", buf[0]); 
            }
            printf("%.2f \n", u0);
            
            MPI_Send(u, 2, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD);
            MPI_Recv(buf, 1, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD, &status);
        }
    }

    MPI_Finalize();   
    return 0;
}
