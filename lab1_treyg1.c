#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#define pi 3.1415926535897932384626433832795028841971693993751058209749445923L
#define X 1.0
#define T 4*pi
#define TRUNCATE 1

#ifndef TRUNCATE

    #define K_PLACE K
    #define K_LAST_PLACE (K-1)
    #define t_PLACE t
    #define t1_PLACE (t+1)
#else

    #define t_PLACE t%2
    #define t1_PLACE (t+1)%2
    #define K_LAST_PLACE (K-1)%2
    #define K_PLACE 2
#endif

#define DEBUG 2

#ifdef DEBUG
    #define d_printf(level, ...) if(DEBUG > level){printf(__VA_ARGS__);}
#else
    #define d_printf() 
#endif

long double u(long double x, long double t){
    return expl(x)*cosl(t);
}

long double f(long double x, long double t){
    return expl(x)*(0.1*cosl(t) - sinl(t));
}

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
    const long double h = (1.0*X/(X_N-1));
    const long double tau = (1.0*T/(T_N-1));
    
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

    const long long int start_x  = M * my_rank / commsize;
    const long long int end_x    = M * (my_rank + 1) / commsize; //[start_x, end_x)

    d_printf(1, "start_x = %d, end_x = %d\n", start_x, end_x);

    for(long long int t = 0; t < K-1; t++){

        res[0][t1_PLACE] = u(0, (t+1)*tau);

        MPI_Request request_send = NULL, request_recv = NULL;

        if(start_x != 0 && t != 0){

            MPI_Irecv(&res[start_x - 1][t_PLACE], 1, MPI_LONG_DOUBLE, my_rank - 1, t, MPI_COMM_WORLD, &request_recv);
        }

        for(long long int x = end_x - 1; x >= start_x; x--){

            if(x == start_x && x != 0 && t != 0){
                MPI_Wait(&request_recv, MPI_STATUS_IGNORE);
            }

            if(x == 0){
                continue;
            }
            else{
                res[x][t1_PLACE] = res[x][t_PLACE] + tau*f(x*h, t*tau) - 0.1*(res[x][t_PLACE] - res[x-1][t_PLACE])*tau/h;
            }

            if(x == end_x - 1 && my_rank != commsize - 1){
                MPI_Isend(&res[x][t1_PLACE],1, MPI_LONG_DOUBLE, my_rank + 1, t+1, MPI_COMM_WORLD, &request_send);
            }

        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    long double res_err = 0;

    if(my_rank == commsize - 1){

        for(long long int x = start_x; x < end_x; x++){
            res_err += fabsl(res[x][K_LAST_PLACE] - u(x*h, T));
        }
        printf("Err = %.32Lf\n", res_err);
        
        for(long long int t = 0; t < K_PLACE; t++){

            d_printf(2, "t[%d]:\n", t);
            for(long long int x = end_x - 10; x < end_x; x++){
                d_printf(2, "%llf ", res[x][t_PLACE]);
            }
            d_printf(2, "\n");
            for(long long int x = end_x - 10; x < end_x; x++){
                d_printf(2, "%Lf ", real_res[x][t_PLACE]);
            }
            d_printf(2, "\n");

        }
    }

    MPI_Finalize();
}