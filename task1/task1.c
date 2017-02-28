#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double f(double x) {
    return 1/(1+x*x);
}

float calculate_integral(double a, double b, double step) {
    if (a == b) return 0;
    double x1 = a, x2 = x1+step;
    double sum = 0;
    while(x2 <= b) {
        sum += (f(x1)+f(x2))/2*step;
        x1 = x2;
        x2 += step;
    }
    return sum;
}

int main(int argc, char const *argv[]) {
    int N = atoi(argv[1]);
    double a = 0, b = 1;
    MPI_Init(&argc, &argv);
    int size = 0;   
    int my_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
        double begin = MPI_Wtime();
        double result_single_thread = calculate_integral(a, b, (b-a)/N);
        double end = MPI_Wtime();
        printf("The single-process result is %.10f, with the time: %f secs", result_single_thread, end-begin);

        begin = MPI_Wtime();
        //TODO: FIX THIS 
        int N_per_proc = N/size > 0 ? N/size : 1;
        double step = (b-a)/N;
        double x1 = a+N_per_proc*step, x2 = x1+N_per_proc*step;
        
        double* buf = (double*)malloc(3*sizeof(double));
        for (int i = 1; i < size; i++) {
            buf[1] = i == size -1 ? b : x2;
            buf[0] = i < N? x1 : buf[1];
            buf[2] = step;
            MPI_Send(&buf[0], 3, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
            x1 = x2;
            x2 += step*N_per_proc;
        }

        double my_result = calculate_integral(a, a+N_per_proc*step, step);
        double result = my_result;
        MPI_Status Status;
        for (int i = 1; i < size; i++) {
            MPI_Recv(&buf[0], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Status);
            result += buf[0];
        }
        end = MPI_Wtime();
        printf("%d: %.10f\n", my_rank, result); 
        printf("%d: The multiple-process result is %.10f with the time: %f\n", my_rank, result, end-begin); 

    }
    if (my_rank != 0) {
        double* buf = (double*)malloc(3*sizeof(double));
        MPI_Status Status;
        MPI_Recv(&buf[0], 3, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &Status);
        double result = calculate_integral(buf[0], buf[1], buf[2]);
        printf("%d: %.10f\n", my_rank, result);
        buf[0] = result;
        MPI_Send(&buf[0], 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
