#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define pi 3.1415926535897932384626433832795028841971693993751058209749445923L

#define X 1.0
//#define T 4*pi
#define T 1.0

#define C 1 //constant C in equasion du/dt + (C)*du/dx = f

//#define ASYNC_MSG  // if defined uses async send and recv; else - uses sync counterparts

//#define TIMING  // if defined measures time from start of calculations to end. Disables error count, debug.

#define max(smth1, smth2) (smth1 > smth2 ? smth1 : smth2)

#define TRUNCATE // if defined uses M x 2 arrays insteaad of M x K

/*
K_PLACE:    placeholder for K
t_PLACE:    placeholder for t
t1_PLACE:   placeholder for t+1
*/
#ifndef TRUNCATE    
    #define K_PLACE K   
    #define t_PLACE t
    #define t1_PLACE (t+1)
#else
    #define t_PLACE t%2
    #define t1_PLACE (t+1)%2
    #define K_PLACE 2
#endif

#define DEBUG 0 // defines level of Debug

#ifdef TIMING   // disables debug
    #undef DEBUG
#endif

#ifdef DEBUG 
    #define d_printf(level, ...) if(DEBUG > level){printf(__VA_ARGS__);}
#else
    #define d_printf(...) 
#endif

long double u(long double x, long double t){
    return expl(x)*cosl(t);
}
// C = 1
long double f(long double x, long double t){
    return expl(x)*(cosl(t) - sinl(t));
}

// long double u(long double x, long double t){
//     if(2*t <= x)
//         return x*t-t*t/2 + cosl(pi*(2*t - x));
//     else
//         return x*t - t*t/2 + (2*t-x)*(2*t-x)/8 + expl(x/2-t);
// }
//// C = 2
// long double f(long double x, long double t){
//     return x+t;
// }

const char filename[] = "output.txt";

int main(int argc, char *argv[]){

    int commsize = -1, my_rank = -1, err = 0;
    int root = 0;

    long long int T_N = 10;
    long long int X_N = 10;

    if(argc >= 3){
        T_N = atoll(argv[1]);
        X_N = atoll(argv[2]);
    }


    const long long int K  = (T_N);
    const long long int M  = (X_N);
    const long double h = (1.0*X/(X_N));
    const long double tau = (1.0*T/(T_N));
    
    err = MPI_Init(&argc, &argv);
    if(err){
        fprintf(stderr, "Problem with MPI_Init\n");
    }
    
    long double res[M][K_PLACE];
    long double real_res[M][K_PLACE];

    for(long long int x = 0; x < M; x++){
        res[x][0] = u(x*h, 0);
    }

    for(long long int t = 0; t < K; t++){
        for(long long int x = 0; x < M; x++){
            real_res[x][t_PLACE] = u(x*h, t*tau);
        }
    }
    
    // res[x][t]
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    const long long int start_x  = my_rank == 0 ? 1 : M * my_rank / commsize;
    const long long int end_x    = M * (my_rank + 1) / commsize; //[start_x, end_x)

    FILE *file = NULL;
    if(my_rank == 0){
        file = fopen(filename, "a");
    }

    d_printf(1, "start_x = %d, end_x = %d\n", start_x, end_x);

    long double res_err = 0;

    clock_t start = clock(), diff;

    for(long long int t = 0; t < K-1; t++){

        res[0][t1_PLACE] = u(0, (t+1)*tau);

        MPI_Request request_send = NULL, request_recv = NULL;
        MPI_Status status = {};

        if(start_x != 1){
#ifdef ASYNC_MSG
            MPI_Irecv(&res[start_x - 1][t1_PLACE], 1, MPI_LONG_DOUBLE, my_rank - 1, t+1, MPI_COMM_WORLD, &request_recv);
#else
            MPI_Recv(&res[start_x - 1][t1_PLACE], 1, MPI_LONG_DOUBLE, my_rank - 1, t+1, MPI_COMM_WORLD, &status);
#endif
        }

        for(long long int x = start_x; x < end_x; x++){

#ifdef ASYNC_MSG
            if(x == start_x && x != 1){
                MPI_Wait(&request_recv, MPI_STATUS_IGNORE);
            }
#endif

            res[x][t1_PLACE] = (4/C*h*tau*f((x-0.5)*h, (t+0.5)*tau) + 2/C*h*(res[x-1][t_PLACE] + res[x][t_PLACE] - res[x-1][t1_PLACE]) 
                                + 2*tau*(res[x-1][t1_PLACE] + res[x-1][t_PLACE] - res[x][t_PLACE]))/(2/C*h+2*tau);
#ifndef TIMING
            res_err += fabsl(res[x][t1_PLACE]  - u(x*h, (t+1)*tau));
#endif       
            d_printf(10, "rank [%d] res[%d][%d] = %Lf\n", my_rank, x, t1_PLACE, res[x][t1_PLACE]);

            if(x == end_x - 1 && my_rank != commsize - 1){
#ifdef ASYNC_MSG
                MPI_Isend(&res[x][t1_PLACE],1, MPI_LONG_DOUBLE, my_rank + 1, t+1, MPI_COMM_WORLD, &request_send);
#else
                MPI_Send(&res[x][t1_PLACE],1, MPI_LONG_DOUBLE, my_rank + 1, t+1, MPI_COMM_WORLD);
#endif
            }

        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;

    if(my_rank == commsize - 1){

        for(long long int t = 0; t < K_PLACE; t++){

            d_printf(2, "t[%d]:\n", t);
            for(long long int x = max(start_x, end_x - 10); x < end_x; x++){
                d_printf(2, "%llf ", res[x][t_PLACE]);
            }
            d_printf(2, "\n");
            for(long long int x = max(start_x, end_x - 10); x < end_x; x++){
                d_printf(2, "%Lf ", real_res[x][t_PLACE]);
            }
            d_printf(2, "\n");

        }
    }

    d_printf(5, "[%d] Err = %.32Lf\n", my_rank, res_err);

    long double final_err = 0;
    MPI_Reduce(&res_err, &final_err, 1, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

    if(my_rank == 0){
        final_err /= (X_N-1)*(T_N-1);

        d_printf(0, "Err = %.32Lf\n", final_err);
        d_printf(0, "Time spent = %d ms\n", msec);

        fprintf(file, "%.32LF, %.32LF, %.32Lf, %d\n", h, tau, final_err, msec);
        fclose(file);
    }

    MPI_Finalize();
}