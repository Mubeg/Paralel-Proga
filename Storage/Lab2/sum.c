#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#define NUM_THREADS_DEFAULT 4


typedef struct _thread_data_t {
    int tid;
    unsigned long long int start;
    unsigned long long int end;
    long double sum;
} thread_data_t;

/* thread function */
void*thr_func(void* arg) {

    thread_data_t *data = (thread_data_t *)arg;

    for(unsigned long long int i = data->start; i < data->end; i++){
        data->sum += 1.0/i;
    }


    pthread_exit(NULL);
}

int main(int argc, char **argv) {

    unsigned int n_thread = NUM_THREADS_DEFAULT;
    unsigned long long int N = 0;
    long double ans = 0;

    if(argc >= 3){
        n_thread = atoi(argv[2]);
    }
    if(argc >= 2){
        N = atoll(argv[1]);
    }

    
    pthread_t* thr = calloc(n_thread, sizeof(pthread_t));
    int rc;
    thread_data_t* thr_data = calloc(n_thread, sizeof(thread_data_t));
    
    clock_t start = clock(), diff;
    for(int i = 0; i < n_thread; i++) {
        
        thr_data[i].tid = i;
        thr_data[i].start = N * i / n_thread + 1;
        thr_data[i].end = N * (i + 1) / n_thread + 1;
        thr_data[i].sum = 0;

        if((rc = pthread_create(&thr[i], NULL, thr_func, &thr_data[i]))) {
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
            return EXIT_FAILURE;
        }
    }
    clock_t start1 = clock(), diff1;

    /* block until all threads complete */
    for(int i = 0; i < n_thread; i++) {
        pthread_join(thr[i], NULL);
        ans += thr_data[i].sum;
    }

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;

    diff1 = clock() - start1;
    int msec1 = diff1 * 1000 / CLOCKS_PER_SEC;



    printf("Sum from 1 to %llu = %.16llf\nTime spent = %d msec; %d msec\n", N, ans, msec, msec1);

    return EXIT_SUCCESS;

}