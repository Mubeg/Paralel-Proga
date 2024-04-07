#include <mpi.h>
#include <stdio.h>
#include <time.h>

int main(int argc, char *argv[]){
    int commsize = -1, my_rank = -1, err = 0;
    err = MPI_Init(&argc, &argv);
    if(err != MPI_SUCCESS){
        fprintf(stderr, "Problem with MPI_Init\n");
    }
    int res = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int root = 0;

#define BIG_N 100
    int N = 0;

    if(argc >= 2){
        N = atoi(argv[1]);
    }

    if(my_rank == root){
        MPI_Status status = {};
        clock_t start = clock(), diff;
        for(int i = 0; i < N; i++){
            MPI_Recv(&res, 1, MPI_INT, (my_rank - 1) % commsize, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&res, 1, MPI_INT, (my_rank + 1) % commsize, 0, MPI_COMM_WORLD);
        }
        diff = clock() - start;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("SEND: Time taken per iteration %d miliseconds\n", msec);
    }
    else{
        MPI_Status status = {};
        clock_t start = clock(), diff;

        for(int i = 0; i < N; i++){
            MPI_Send(&res, 1, MPI_INT, (my_rank + 1) % commsize, 0, MPI_COMM_WORLD);
            MPI_Recv(&res, 1, MPI_INT, (my_rank - 1) % commsize, 0, MPI_COMM_WORLD, &status);
        }
        diff = clock() - start;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("RECV: Time taken per iteration %d miliseconds\n", msec);

    }

    MPI_Finalize();
}