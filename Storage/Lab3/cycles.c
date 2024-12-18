#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <math.h>

#define NUM_THREADS_DEFAULT 4
#define ISIZE 5000
#define JSIZE 5000
#define INNER_BATCH_SIZE 500
#define PRECISION "3"


static const double _ = 1;

/*
arr1 - input
arr2 - output
i - cur_i
j - cur_j
init - init mode
*/
typedef struct function{
    void (*f)(void **arr1, void **arr2, int i, int j, int init);
    int D[2];
} function_t;

typedef struct shifts{
    int start_shifti;
    int start_shiftj;
    int end_shifti;
    int end_shiftj;
} shifts_t;

shifts_t get_shifts(function_t *f){
    shifts_t temp;
    temp.start_shifti = f->D[0] > 0 ? f->D[0] : 0;
    temp.start_shiftj = f->D[1] > 0 ? f->D[1] : 0;
    temp.end_shifti = f->D[0] < 0 ? f->D[0] : 0;
    temp.end_shiftj = f->D[1] < 0 ? f->D[1] : 0;
    return temp;
}


void etalon(void **arr1, void **arr2, int i, int j, int initialize){

    double **a = (double **) arr1;
    double **b = (double **) arr2;

    if(!initialize)
        b[i*ISIZE][j] = 10*i + j;
    else
        b[i*ISIZE][j] = sin(2*_*a[i*ISIZE][j]);
}

typedef struct _thread_data_t {
    int tid;
    function_t* f;
    void ** arr1;
    void ** arr2;

} thread_data_t;

typedef struct batch {
    int transpose;
    int barrier;
    int i_start;
    int i_end;
    int j_start;
    int j_end;
} batch_t;

batch_t *batch_array = NULL;
int n_batches = 0;
int next_batch = 0;
pthread_mutex_t lock_x;
pthread_barrier_t barrier; 

/* thread function */
void*thr_func(void* arg) {

    thread_data_t *data = (thread_data_t *)arg;
    while(1){
        pthread_mutex_lock(&lock_x);
        if(next_batch >= n_batches){
            pthread_mutex_unlock(&lock_x);
            break;
        }
        batch_t b = batch_array[next_batch];
        next_batch++;
        pthread_mutex_unlock(&lock_x);

        if(b.transpose){
            for(int j = b.j_start; j < b.j_end; j++){
                for(int i = b.i_start; i < b.i_end; i++){
                    data->f->f(data->arr1, data->arr2, i, j, 0);
                }
            }
        }
        else{
            for(int i = b.i_start; i < b.i_end; i++){
                for(int j = b.j_start; j < b.j_end; j++){
                    data->f->f(data->arr1, data->arr2, i, j, 0);
                }
            }
        }

        if(b.barrier){
            pthread_barrier_wait(&barrier);
        }
    }
    pthread_exit(NULL);
}

double run_threads(function_t *f, void **arr1, void**arr2, int n_runners){

    pthread_t* thr = calloc(n_runners, sizeof(pthread_t));
    thread_data_t* thr_data = calloc(n_runners, sizeof(thread_data_t));
    struct timespec start, finish;
    int rc;

    clock_gettime(CLOCK_MONOTONIC, &start);
    for(int n = 0; n < n_runners; n++) {
            
        thr_data[n].tid = n;
        thr_data[n].f = f;
        thr_data[n].arr1 = arr1;
        thr_data[n].arr2 = arr2;

        if((rc = pthread_create(&thr[n], NULL, thr_func, &thr_data[n]))) {
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
            return -1;
        }
    }

    for(int i = 0; i < n_runners; i++) {
        pthread_join(thr[i], NULL);
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);

    double elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    return elapsed;
}

double parallel_upper(function_t *f, void **arr1, void **arr2, int n_runners, shifts_t sh, int transpose){

    if(transpose){
        n_batches = JSIZE - sh.start_shiftj + sh.end_shiftj;
        batch_array = calloc(n_batches, sizeof(batch_t));
        for(int j = sh.start_shiftj; j < JSIZE + sh.end_shiftj; j++){
            batch_array[j].transpose = 1;
            batch_array[j].j_start = j;
            batch_array[j].j_end = j;
            batch_array[j].i_start = sh.start_shifti;
            batch_array[j].i_end = ISIZE + sh.end_shifti;
        }
    }
    else{
        n_batches = ISIZE - sh.start_shifti + sh.end_shifti;
        batch_array = calloc(n_batches, sizeof(batch_t));
        for(int i = sh.start_shifti; i < ISIZE + sh.end_shifti; i++){
            batch_array[i].i_start = i;
            batch_array[i].i_end = i;
            batch_array[i].j_start = sh.start_shiftj;
            batch_array[i].j_end = JSIZE + sh.end_shiftj;
        }
    }
    next_batch = 0;

    double time = run_threads(f, arr1, arr2, n_runners);

    free(batch_array);

    return time;
}

#define min(a, b) ((a) > (b) ? (b) : (a))

double parallel_inner_with_barrier(function_t *f, void **arr1, void **arr2, int n_runners, shifts_t sh){

    n_runners = min(n_runners, (int)((JSIZE - sh.start_shiftj + sh.end_shiftj)/INNER_BATCH_SIZE));

    n_batches = (int) ceil((ISIZE - sh.start_shifti + sh.end_shifti) * _ * (JSIZE - sh.start_shiftj + sh.end_shiftj));
    batch_array = calloc(n_batches, sizeof(batch_t));
    int curr_batch = 0;
    for(int i = sh.start_shifti; i < ISIZE + sh.end_shifti; i++){
         for(int j = sh.start_shiftj; j < JSIZE + sh.end_shifti; j += INNER_BATCH_SIZE){
                batch_array[curr_batch].i_start = i;
                batch_array[curr_batch].i_end = i;
                batch_array[curr_batch].j_start = j;
                batch_array[curr_batch].j_end = min(j + INNER_BATCH_SIZE, JSIZE + sh.end_shiftj);
                curr_batch++;
        }
        for(int n = 0; n < n_runners; n++){
            batch_array[curr_batch - n].barrier = 1;
        }
    }
    next_batch = 0;
    pthread_barrier_init(&barrier, NULL, n_runners);

    double time = run_threads(f, arr1, arr2, n_runners);

    free(batch_array);

    return time;
}


double run_double_loop(function_t *f, const char *out, int n_runners){

    if(n_runners <= 0 || out == NULL || f == NULL || f->f == NULL){
        return -1;
    }

    double arr[ISIZE][JSIZE];
    struct timespec start, finish;
    int solo = n_runners == 1;
    double elapsed;

    //init
    for(int i = 0; i < ISIZE; i++){
        for(int j = 0; j < JSIZE; j++){
            f->f(NULL, (void **)arr, i, j, 1);
        }
    }
    
    shifts_t sh = get_shifts(f);
    
    //solo
    if(solo){
        clock_gettime(CLOCK_MONOTONIC, &start);
        for(int i = sh.start_shifti; i < ISIZE + sh.end_shifti; i++){
            for(int j = sh.start_shiftj; j < JSIZE + sh.end_shifti; j++){
                f->f((void **)arr, (void **)arr, i, j, 0);
            }
        }
        clock_gettime(CLOCK_MONOTONIC, &finish);
    } else{

        /*paralleled
        D - solution
        D = [-1, 1] === arr[i][j] = f(arr[i+1][j-1])
        [0, 0] - parallel upper loop
        [0, 1] - parallel by upper loop
        [0, -1] - parallel by upper loop
        [1, 0] - switch pass order and parr upper
        [1, 1] - parr by inner with barriers
        [1, -1] - parr by inner with barriers
        [-1, 0] - switch pass order and parr by upper
        [-1, 1] - copy table and upper
        [-1, -1] - copy table and upper
        */
        if(f->D[0] == 0 || f->D[1] == 0){
            elapsed = parallel_upper(f, (void **)arr, (void **)arr, n_runners, sh, f->D[0] != 0);
        }
        else if(f->D[0] > 0){
            elapsed = parallel_inner_with_barrier(f, (void **)arr, (void **)arr, n_runners, sh);
        }
        else if(f->D[0] < 0){
            double arr2[ISIZE][JSIZE];
            for(int i = 0; i < ISIZE; i++){
                for(int j = 0; j < JSIZE; j++){
                    arr2[i][j] = arr[i][j];
                }
            }
            elapsed = parallel_upper(f, (void **)arr, (void **)arr2, n_runners, sh, 0);
        }
    }

    //finalize
    FILE * file = fopen(out, "w");
    for(int i = 0; i < ISIZE; i++){
        for(int j = 0; j < JSIZE; j++){
            fprintf(file, "%."PRECISION"lf", arr[i][j]);
        }
        fprintf(file,"\n");
    }
    fclose(file);

    if(solo){
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    }
    return elapsed;

}


int main(int argc, char **argv) {

    unsigned int n_thread = NUM_THREADS_DEFAULT;

    if(argc >= 2){
        n_thread = atoi(argv[1]);
    }

    function_t etalon_s;
    etalon_s.f = &etalon;
    etalon_s.D[0] = 0;
    etalon_s.D[1] = 0;

    double elapsed = run_double_loop(&etalon_s, "etalon_1.txt", 1);
    printf("Etalon ran in one thread for %.3lf", elapsed);
    elapsed = run_double_loop(&etalon_s, "etalon_mul.txt", n_thread);
    printf("Etalon ran in %d threads for %.3lf", n_thread, elapsed);


    // printf("Max_err_on_interval = %llf\nTemp:\n", thr_data[0].max_interval_err);
    // for(int i = 0; i < NUM_INTERVALS; i++){
    //     printf("[%d] = %llf, ", i, temp[i]);
    // }
    // printf("\n");

    // printf("Ans = %.16llf, Real ans = %.16llf\n"
    //        "Actual err = %.16llf, Max err = %.16lf\n"
    //        "%s\n"
    //        "Time spent = %lf msec; %lf msec\n", ans, real_ans, fabsl(ans - real_ans), max_error, max_error > fabsl(ans - real_ans) ? "GOOD": "NOT GOOD", elapsed*1000, elapsed1*1000);

    // FILE *file = NULL;
    // file = fopen(filename, "a");
    // fprintf(file, "%d, %lf\n", n_thread, elapsed*1000);
    // fclose(file);

    return EXIT_SUCCESS;

}