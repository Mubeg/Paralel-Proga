#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]){
    int commsize = -1, my_rank = -1, err = 0;
    err = MPI_Init(&argc, &argv);
    if(err){
        fprintf(stderr, "Problem with MPI_Init\n");
    }
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    printf("Communicator size = %d My rank = %d\n", commsize, my_rank);
    MPI_Finalize();
}