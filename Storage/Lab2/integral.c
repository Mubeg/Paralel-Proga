#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <math.h>

#define NUM_THREADS_DEFAULT 4
#define NUM_POINTS_FOR_ANALYSYS_PER_INTERVAL 100
#define NUM_INTERVALS 1000
#define START_X -100
#define END_X 1000
#define INTERVAL_LENGTH (1.0 * ((END_X) - (START_X)) / NUM_INTERVALS)

// long double temp[NUM_INTERVALS] = {0};
long double ans = 0;
int next_calculated = 0;
pthread_mutex_t lock_x;


const char filename[] = "output.txt";


long double f(long double x){
    return x*x - 2*x + 1;
}

long double f_real(long double x){
    return x*x*x/3 - x*x + x;
}

typedef struct _thread_data_t {
    int tid;
    long double max_interval_err;
} thread_data_t;

/* thread function */
void*thr_func(void* arg) {

    thread_data_t *data = (thread_data_t *)arg;
    
    while(1){
        int my_interval_n = -1;
        long double my_ans = 0;
        long double my_max_f_d = -1;
        long double my_h = -1;
        long double my_left, my_right;

        pthread_mutex_lock(&lock_x);
        if(next_calculated == NUM_INTERVALS){
            pthread_mutex_unlock(&lock_x);
            break;
        }
        my_interval_n = next_calculated;
        next_calculated++;
        pthread_mutex_unlock(&lock_x);

        my_left = START_X + INTERVAL_LENGTH * my_interval_n;
        my_right = START_X + INTERVAL_LENGTH * (my_interval_n + 1);
        long double local_h = INTERVAL_LENGTH / NUM_POINTS_FOR_ANALYSYS_PER_INTERVAL;
        for(long double i = my_left; i < my_right; i += local_h){
            my_max_f_d = fmaxl(my_max_f_d, fabsl((f(i + local_h) - f(i)))/local_h);
        }
        
        long double n = (my_right - my_left)*(my_right - my_left)/2*my_max_f_d/data->max_interval_err;
        my_h = (my_right - my_left)/n;
        for(long double i = my_left; i < my_right; i += my_h){
            if(i + my_h > my_right){
                my_ans += f(i)*(my_right - i);
            }
            else{
                my_ans += f(i)*my_h;
            }
        }

        pthread_mutex_lock(&lock_x);
        ans += my_ans;
        pthread_mutex_unlock(&lock_x);

    }
    pthread_exit(NULL);
}

int main(int argc, char **argv) {

    struct timespec start, start1, finish;
    double elapsed, elapsed1;

    unsigned int n_thread = NUM_THREADS_DEFAULT;
    double max_error = 0;

    if(argc >= 3){
        n_thread = atoi(argv[2]);
    }
    if(argc >= 2){
        max_error = atof(argv[1]);
    }

    
    pthread_t* thr = calloc(n_thread, sizeof(pthread_t));
    int rc;
    thread_data_t* thr_data = calloc(n_thread, sizeof(thread_data_t));
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    for(int i = 0; i < n_thread; i++) {
        
        thr_data[i].tid = i;
        thr_data[i].max_interval_err = max_error / NUM_INTERVALS;

        if((rc = pthread_create(&thr[i], NULL, thr_func, &thr_data[i]))) {
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
            return EXIT_FAILURE;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &start1);

    /* block until all threads complete */
    for(int i = 0; i < n_thread; i++) {
        pthread_join(thr[i], NULL);
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    elapsed1 = (finish.tv_sec - start1.tv_sec);
    elapsed1 += (finish.tv_nsec - start1.tv_nsec) / 1000000000.0;


    // printf("Max_err_on_interval = %llf\nTemp:\n", thr_data[0].max_interval_err);
    // for(int i = 0; i < NUM_INTERVALS; i++){
    //     printf("[%d] = %llf, ", i, temp[i]);
    // }
    // printf("\n");

    long double real_ans = f_real(END_X) - f_real(START_X);
    // printf("Ans = %.16llf, Real ans = %.16llf\n"
    //        "Actual err = %.16llf, Max err = %.16lf\n"
    //        "%s\n"
    //        "Time spent = %lf msec; %lf msec\n", ans, real_ans, fabsl(ans - real_ans), max_error, max_error > fabsl(ans - real_ans) ? "GOOD": "NOT GOOD", elapsed*1000, elapsed1*1000);

    FILE *file = NULL;
    file = fopen(filename, "a");
    fprintf(file, "%d, %lf\n", n_thread, elapsed*1000);
    fclose(file);

    return EXIT_SUCCESS;

}