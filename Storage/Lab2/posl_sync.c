#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#define NUM_THREADS_DEFAULT 4

int global_n = 0;
int last_access_id = -1;
int times_missed_lock = 0;
int times_locked = 0;
pthread_mutex_t lock_x;

typedef struct _thread_data_t {
    int tid;
    long double sum;
} thread_data_t;

/* thread function */
void*thr_func(void* arg) {

    thread_data_t *data = (thread_data_t *)arg;
    while(1){
        pthread_mutex_lock(&lock_x);
        times_locked++;
        if(last_access_id == data->tid - 1){
            global_n += (data->tid + 1)*3;
            last_access_id = data->tid;
            printf("[%d]: global_n = %d\n", data->tid, global_n);
            pthread_mutex_unlock(&lock_x);
            break;
        }
        times_missed_lock++;
        pthread_mutex_unlock(&lock_x);
    }
    pthread_exit(NULL);
}

int main(int argc, char **argv) {

    unsigned int n_thread = NUM_THREADS_DEFAULT;

    if(argc >= 2){
        n_thread = atoi(argv[1]);
    }

    
    pthread_t* thr = calloc(n_thread, sizeof(pthread_t));
    int rc;
    thread_data_t* thr_data = calloc(n_thread, sizeof(thread_data_t));
    
    clock_t start = clock(), diff;
    for(int i = 0; i < n_thread; i++) {
        
        thr_data[i].tid = i;
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
    }

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;

    diff1 = clock() - start1;
    int msec1 = diff1 * 1000 / CLOCKS_PER_SEC;



    printf("Time spent = %d msec; %d msec\nTimes took lock for nothing = %d out of %d total locks\n", msec, msec1, times_missed_lock, times_locked);

    return EXIT_SUCCESS;

}