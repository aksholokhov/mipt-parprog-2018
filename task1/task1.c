#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double f(double x) {
    return 1/(1+x*x);
}

double calculate_integral2(double a, double b, double step) {
    if (a == b) return 0;
    double x1 = a, x2 = x1+step;
    double sum = 0;
    while(x2 <= b) {
        double result = (f(x1)+f(x2))/2*step;
      //  printf("%f\n", result);
        sum += result;
        x1 = x2;
        x2 += step;
    }
    //printf("sum: %.10f\n", sum);
    return sum;
}

double calculate_integral(double a, double b, double step)
{
    int i=0;
    double int_value=0;
    for(i=0; i<(b-a)/step; i++)
    {
        int_value += ( f( a + step * i ) + f( a + ( i + 1 ) * step ) ) * step /2.0;
    }
    return int_value;
}



int main(int argc, char const *argv[]) {
    int N = atoi(argv[1]);
    double a = 0, b = 1;
    double step = (b-a)/N;
    MPI_Init(&argc, &argv);
    int size = 0;   
    int my_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
        double begin = MPI_Wtime();
        double result_single_thread = calculate_integral(a, b, (b-a)/N);
        double end = MPI_Wtime();
        printf("The single-process result is %.10f, with the time: %f secs\n", result_single_thread, end-begin);

        begin = MPI_Wtime();
        
        int npp = N/size;


        double x1 = a, x2 = a+step*npp;
        for (int i = 1; i < size; i++) {
            double* buf = (double*)malloc(2*sizeof(double));
            buf[0] = x1;
            buf[1] = x2;
            MPI_Send(&buf[0], 2, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
            x1 = x2;
            x2 += step*npp;
        }

        double my_result = calculate_integral(x1, b, step);
        double result = my_result;

        double buf;
        MPI_Status status;
        for (int i = 1; i < size; i++) {
            MPI_Recv(&buf, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            result += buf;
        }

        end = MPI_Wtime();
        printf("%d: The multiple-process result is %.10f with the time: %f\n", my_rank, result, end-begin); 

    }
     if (my_rank != 0) {
        double* buf = (double*)malloc(2*sizeof(double));
        MPI_Status Status;
        MPI_Recv(&buf[0], 2, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &Status);
        double result = calculate_integral(buf[0], buf[1], step);
        MPI_Send(&result, 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
    } 

    MPI_Finalize();
    return 0;
}
