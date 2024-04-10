#include <mpi.h>
#include <stdio.h>

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

    if(my_rank == root){
        err = MPI_Send(&res, 1, MPI_INT, (my_rank + 1) % commsize, my_rank, MPI_COMM_WORLD);
        if(err != MPI_SUCCESS){
            fprintf(stderr, "Problem with MPI_Send root\n");
        }
    }

    MPI_Status status = {};
    MPI_Recv(&res, 1, MPI_INT, (my_rank - 1) % commsize, (my_rank-1)%commsize, MPI_COMM_WORLD, &status);
    if(err != MPI_SUCCESS){
        fprintf(stderr, "Problem with MPI_Send rank %d\n", my_rank);
    }
    res += 3;
    printf("My rank = %d, res = %d\n", my_rank, res);


    if(my_rank != root){
        MPI_Send(&res, 1, MPI_INT, (my_rank + 1) % commsize, my_rank, MPI_COMM_WORLD);
        if(err != MPI_SUCCESS){
            fprintf(stderr, "Problem with MPI_Send, rank %d\n", my_rank);
        }
    }


    MPI_Finalize();
}