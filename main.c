#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>

#define INPUT_FILENAME "input.txt"
#define INF 99999

typedef struct {
    int n;
    int i;
    int k;
} args_t;

void initalize_graph(args_t* args);
void get_shortest_paths(args_t* args);
void* worker(void* args);
void* worker_st(void* args);
void display_result_mt(void);
void display_result_st(void);

int** graph;
int** dist;

// Number of nodes
int N;
// Number of edges
int M;
// Keeps track of the number of threads reading the dist array
int reader_count;

// Stores the run time for the single thread(st) and mulit-threaded(mt) FW
double st_time;
double mt_time;

// Semaphores for writing to the dist array and reading
sem_t sem_dist_lock;
sem_t sem_read_lock;


/*
 * Get the adjancy matrix from the input file and
 * initalizes the graph and dist arrays.
 */
void initalize_graph(args_t* args) {
    FILE* file = fopen(INPUT_FILENAME, "r");

    if (file == NULL) {
        perror("Input file cannot be opened.");
	exit(EXIT_FAILURE);
    }
    
    fscanf(file, "%d %d", &N, &M);

    graph = (int**)malloc(N * sizeof(int*));
    for (int i = 0; i < N; ++i) {
        graph[i] = (int*)malloc(N * sizeof(int));
	for (int j = 0; j < N; ++j) {
	    if (i != j) graph[i][j] = INF;
	    else graph[i][j] = 0;
	}
    }

    dist = (int**)malloc(N * sizeof(int*));
    for (int i = 0; i < N; ++i) {
        dist[i] = (int*)malloc(N * sizeof(int));
	for (int j = 0; j < N; ++j) {
	    if (i != j) dist[i][j] = INF;
	    else dist[i][j] = 0;
	}
    }

    int edge_count = 0;

    while(!feof(file)) {
	int u, v, w;
        fscanf(file, "%d %d %d", &u, &v, &w);

	if (w < 0) {
	    perror("Cannot have a neagtive weight!");
	    exit(EXIT_FAILURE);
	}
	
	if (u > N || u < 0 || v > N || v < 0) {
	    perror("A node is out of bounds. Check your input file.");
	    exit(EXIT_FAILURE);
	}

	if (edge_count > M) {
	    perror("The number of edges exceeds the maximum");
	    exit(EXIT_FAILURE);
	}

	graph[u - 1][v - 1] = w;
	graph[v - 1][u - 1] = w;
	
	dist[u - 1][v - 1] = w;
	dist[v - 1][u - 1] = w;

	edge_count++;
    }
}

/*
 * Worker for the multi-threaded floyd warshall
 */
void* worker(void* _args) {
    args_t* args = (args_t*) _args;
    
    for (int j = 0; j < args->n; ++j) {
        sem_wait(&sem_read_lock); 
	reader_count++;
	if (reader_count == 1) {
            sem_wait(&sem_dist_lock);
        }
	sem_post(&sem_read_lock);

	int update_dist = dist[args->i][args->k] + dist[args->k][j] < dist[args->i][j];

	sem_wait(&sem_read_lock);
        reader_count--;
	if (reader_count == 0) {
            sem_post(&sem_dist_lock);
        }
	sem_post(&sem_read_lock);

	if (update_dist) {	
            sem_wait(&sem_dist_lock);
	    dist[args->i][j] = dist[args->i][args->k] + dist[args->k][j];
	    sem_post(&sem_dist_lock);
	}
    }
    
    pthread_exit(NULL);
}

/*
 * Worker for the single thread floyd warshall
 */
void* worker_st(void* args) {
    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (graph[i][k] + graph[k][j] < graph[i][j]) graph[i][j] = graph[i][k] + graph[k][j];
	    }
        }
    }
    pthread_exit(NULL);
}

/*
 * Runs the multi-threaded and the single-threaded
 * implamentation of the floyd warshall algorthum and
 * get the run time of each.
 */
void get_shortest_paths(args_t* args) {
    sem_init(&sem_dist_lock, 0, 1);
    sem_init(&sem_read_lock, 0, 1);

    pthread_t* workers = (pthread_t*)malloc(N * sizeof(pthread_t));
    args = (args_t*)malloc(N * sizeof(args_t));

    // Mulit-Threaded Implementation
    mt_time = 0.0;    
    struct timeval start, end;
    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; ++i) {
	    args[i].n = N;
	    args[i].i = i;
	    args[i].k = k;
	    pthread_create(&workers[i], NULL, worker, (void*)(&args[i]));    
	}

        gettimeofday(&start, 0);	
        
	for (int i = 0; i < N; ++i) {
            pthread_join(workers[i], NULL);
	}

        gettimeofday(&end, 0);

        double time_iter = (double) (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.;
	mt_time += time_iter;
    }
    
    // Single-Threaded Implementation
    pthread_t worker;
    pthread_create(&worker, NULL, worker_st, (void*)NULL);
    gettimeofday(&start, 0);	
    pthread_join(worker, NULL);
    gettimeofday(&end, 0);
    st_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.;

    sem_destroy(&sem_dist_lock);
    sem_destroy(&sem_read_lock);
}

/*
 * Display the result from the Multi-Threaded Floyd Warshalls
 */
void display_result_mt(void) {
    if (N > 30) {
        printf("Matrix is too big to print.\n");
	return;
    }

    printf("\nOutput:\n");

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
	    if (dist[i][j] == INF) printf("INF ");
            else printf("%d ", dist[i][j]);
	}
	printf("\n");
    }
}

/*
 * Display the result from the Single-Thread Floyd Warshalls
 */
void display_result_st(void) {
    if (N > 30) {
        printf("Matrix is too big to print.\n");
	return;
    }

    printf("\nOutput:\n");

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
	    if (dist[i][j] == INF) printf("INF ");
            else printf("%d ", dist[i][j]);
	}
	printf("\n");
    }
}

int main(void) {
    args_t* args;
    initalize_graph(args);
    get_shortest_paths(args);

    printf("|-----------------------------------------------|\n");
    printf("|------- Mulit-Threaded Floyd Warshalls --------|\n");
    printf("|-----------------------------------------------|\n");
    printf("Elapsed Time: %f seconds\n", mt_time);
    display_result_mt();  

    printf("|-----------------------------------------------|\n");
    printf("|------- Single-Threaded Floyd Warshalls -------|\n");
    printf("|-----------------------------------------------|\n");
    printf("Elapsed Time: %f seconds\n", st_time);
    display_result_st();

    printf("\nSpeed Up: %f\n", st_time/mt_time);   

    return EXIT_SUCCESS;
}

