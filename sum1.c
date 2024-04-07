#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]){
    int commsize = -1, my_rank = -1, err = 0;
    err = MPI_Init(&argc, &argv);
    if(err){
        fprintf(stderr, "Problem with MPI_Init\n");
    }
    int N = -1;
    if (argc >= 2){
        N =  atoi(argv[1]);
    }
    double res = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    N;

    int n = N * my_rank / commsize + 1;
    int m = N * (my_rank + 1) / commsize + 1; 

    printf("n = %d, m = %d\n", n, m);

    for(int i = n; i < m; i++){
        double k = i;
        res += 1/k;
    }

    double final = -1;
    int root = 0;

    printf("Communicator size = %d My rank = %d, res = %f\n", commsize, my_rank, res);

    MPI_Reduce(&res, &final, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    if(my_rank == root){
        printf("Sum is %f\n", final);
    }

    MPI_Finalize();
}