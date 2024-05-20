#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#define NUM_THREADS_DEFAULT 4


typedef struct _thread_data_t {
    int tid;
    double stuff;
} thread_data_t;

/* thread function */
void*thr_func(void* arg) {

    thread_data_t *data = (thread_data_t *)arg;
    printf("[%d]: Hello world\n", data->tid);
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
    for(int i = 0; i < n_thread; i++) {
        thr_data[i].tid = i;
        if((rc = pthread_create(&thr[i], NULL, thr_func, &thr_data[i]))) {
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
            return EXIT_FAILURE;
        }
    }
    /* block until all threads complete */
    for(int i = 0; i < n_thread; i++) {
        pthread_join(thr[i], NULL);
    }
    return EXIT_SUCCESS;

}