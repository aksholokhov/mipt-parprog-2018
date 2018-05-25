/*
 * =====================================================================================
 *
 *       Filename:  calc_tau.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  25.05.2018 13:45:34
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double x = 20;
    MPI_Status status;
    MPI_Barrier(MPI_COMM_WORLD);
    double t = MPI_Wtime();
    if (rank == 0) {
        for (int i = 0; i < 10000; i++){
            x = x /1.01;
        }
        double arithmetic = MPI_Wtime() - t;
        x = 20.0;
        for (int i = 0; i < 10000; i++) {
            MPI_Send(&x, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double transmission = MPI_Wtime() - (t+arithmetic);
        printf("tau = %f \n", transmission/arithmetic);
    }
    else {
        double x = 0;
        for (int i = 0; i < 10000; i++) {
            MPI_Recv(&x, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
        }
        MPI_Barrier(MPI_COMM_WORLD);

    }
    MPI_Finalize();
    return 0;

}

