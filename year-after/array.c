/*
 * =====================================================================================
 *
 *       Filename:  array.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  24.05.2018 12:59:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>

float A(float i, float j) {
    return 1.0/(i+j+1);
}

int main(int argc, char **argv) {
    int N = atoi(argv[1]);
    double t = 0;

    MPI_Init(&argc, &argv);
    int size = 0;
    int rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;


    double* x = (double*)malloc(N*sizeof(double));
    int block_size = N/size;
    double* ptr = x + rank*block_size;
    for (int i = 0; i <= block_size; i++){
        *ptr = rank*block_size + i;
        ptr++;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    
    if (rank == 0) {
        ptr = x + block_size;
        for (int i = 1; i < size; i++) {
            MPI_Recv(ptr, block_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            ptr += block_size;
        }
        //printf("Got all the stuff! %.2f, %.2f", x[0], x[N-1]);
        for (int i = 1; i < size; i++) {
            MPI_Send(x, N, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        }
    }
    else {
        ptr = x + rank*block_size;
        MPI_Send(ptr, block_size, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
        MPI_Recv(x, N, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &status);
        //printf("%d: %.1f\n", rank, x[N-1]);
    }
    double res = 0;
    for (int i =  0; i < block_size; i++) {
        double row = 0;
        for (int j = 0; j < N; j++) {
            row += A(rank*block_size + i,j)*x[j];
        }
        res += row*row;
    }

    if (rank == 0) {
        double buf = 0;
        for (int i = 1; i < size; i++) {
            MPI_Recv(&buf, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            res += buf;
        }

        t = MPI_Wtime() - t;
        
        printf("%f \n", t);
        //printf("Answer: %.1f\n", sqrt(res));
    }
    else {
        MPI_Send(&res, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
    }


    free(x);
    MPI_Finalize();
    
    return 0;
}
