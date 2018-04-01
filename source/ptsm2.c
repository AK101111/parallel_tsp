#include "ptsm.h"
#include <omp.h>

int **finalPath;	//	best paths explored by different threads.
int **currentPath;	//	current path of each thread.
int **visited;		// visited array for each thread.
int *dirty; 		// array which stores dirty bit of each thread
					// denoting if that thread needs to update its 
					// final path from its local path.
int *threadInfo;
int *firstMin;
int *secondMin;

int NUM_THREADS;

volatile int bestLength = INT_MAX;	// Current Best, shared among all threads.
volatile int bestIndex =  -1;
/*
	Keeping bestLength volatile, however finalPath is not.
	Instead each thread has its own path array to minimize time
	spent in locked zone. (updating bestLength only).
*/
Graph* populate_data(char* fileName, unsigned int numCities){
	FILE *file;
	file=fopen(fileName, "r");
	if(!file)
		return NULL;
	
	Graph* ret = malloc(sizeof(Graph));
	if(!ret)
		return NULL;
	
	ret->numCities = numCities+1;
	ret->distance = malloc(sizeof(int*)*(numCities+1));
	for(int i=0; i<numCities+1; ++i){
		if(i==0){
			ret->distance[0] = calloc(numCities+1, (numCities+1)*sizeof(int));
			continue;
		}
		ret->distance[i] = malloc(sizeof(int)*(numCities+1));
		ret->distance[i][0] = 0;
		for(int j=1; j<numCities+1; ++j){
			if(!fscanf(file, "%d", &(ret->distance[i][j])))
				return NULL;
		}
	}
	fclose(file);
	return ret;
}
/*
int _get(int threadId){
	int k = -1;
	int j = -1;
	bool switch = -1;
	#pragma omp critical
	{
		for(int i=0;i<NUM_THREADS; ++i)
			if(threadInfo[i] == threadId)
				return i;
		int j=-1;
		for(int i=0; i<NUM_THREADS; ++i){
			if(threadInfo[i] == -1){
				threadInfo[i] = threadId;
				j = i;
				i = NUM_THREADS;
			}
		}
		return j;
	}
	return j;
}*/
__attribute__((always_inline)) static inline int _get(int threadId) { return threadId;}

/*
	utility min function, inlined for max perf. 
*/
__attribute__((always_inline)) static inline int _Min(int a, int b){
	return (a<b)?(a):(b);
}

void _precompute_bounds(Graph *G){
	#pragma omp parallel for
	for(int i=0; i<G->numCities; ++i){
		int min = INT_MAX;
		for(int j=0; j<G->numCities; ++j)
			if(j != i)
				min = _Min(min, G->distance[i][j]);
		firstMin[i] = min;
	}

	#pragma omp parallel for
	for(int i=0; i<G->numCities; ++i){
		int min = INT_MAX;
		int secondmin = INT_MAX;
		for(int j=0; j<G->numCities; ++j){
			if(j!=i){
				if(G->distance[i][j] <= min){
					secondmin = min;
					min = G->distance[i][j];
				}else if(G->distance[i][j] <= secondmin && 
					G->distance[i][j] != min){
					secondmin = G->distance[i][j];
				}	
			}
		}
		secondMin[i] = secondmin;
	}
	return;
}

/*
	Picks min weight outgoing edge from node
*/
__attribute__((always_inline)) static inline int _pick_min(Graph *G, int city){
	return firstMin[city];
}

/*
	Picks min/second min weight outgoing edge from node
*/
__attribute__((always_inline)) static inline int _pick_second(Graph* G, int city){
	return secondMin[city];
}

void _tsp_recursive_serial(Graph *G, float current_max, int current, 
	int length, int threadId);

void _tsp_recursive_parallel(Graph *G, float current_max, int current, 
	int length, int maxThreads){
	// length == 2
	if(length == G->numCities){
		
		if(current + G->distance[currentPath[maxThreads][length-1]][currentPath[maxThreads][0]] < bestLength){
			bestLength = current + G->distance[currentPath[maxThreads][length-1]][currentPath[maxThreads][0]];
			bestIndex = 0;
			dirty[maxThreads] = 1;
		}

		// hack to move heavy computation out of lock
		if(dirty[maxThreads]){
			memcpy(finalPath[maxThreads], currentPath[maxThreads], (G->numCities)*sizeof(int));
			dirty[maxThreads] = 0;
		}
		return;
	}
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
	for(int i=0; i<G->numCities; ++i){
		//#pragma omp task
		{
		int threadId = _get(omp_get_thread_num());
		if(threadId == -1)
		    perror("more threads than cities\n");
		if(!(visited[threadId][i])){
			float bound = current_max;
			int change = _pick_min(G,i);
			if(length == 1)
				change += _pick_min(G,currentPath[threadId][length-1]);
			else
				change += _pick_second(G,currentPath[threadId][length-1]);
			int wt = current + G->distance[currentPath[threadId][length-1]][i];
			bound -= ((float)change/2.0);
			//printf("done\n");
			if(bound + wt < bestLength){
				currentPath[threadId][length] = i;
				visited[threadId][i] = 1;
				_tsp_recursive_serial(G, bound, wt, 
					length+1, threadId);
				// reset visited.
				visited[threadId][i] =  0;
			}
		}
		}
	}
}
/*
	Recursive tsp with b&b.
*/
void _tsp_recursive_serial(Graph *G, float current_max, int current, 
	int length, int threadId){
	//printf("%lf %d %d\n", current_max, current, length);
	if(length == G->numCities){
		
		#pragma omp critical
		if(current + G->distance[currentPath[threadId][length-1]][currentPath[threadId][0]] < bestLength){
			bestLength = current + G->distance[currentPath[threadId][length-1]][currentPath[threadId][0]];
			bestIndex = threadId;
			dirty[threadId] = 1;
		}

		// hack to move heavy computation out of lock
		if(dirty[threadId]){
			memcpy(finalPath[threadId], currentPath[threadId], (G->numCities)*sizeof(int));
			dirty[threadId] = 0;
		}
		return;
	}

	
	for(int i=0; i<G->numCities; ++i){
		
		if(!(visited[threadId][i])){
			float bound = current_max;
			int change = _pick_min(G,i);
			if(length == 1)
				change += _pick_min(G,currentPath[threadId][length-1]);
			else
				change += _pick_second(G,currentPath[threadId][length-1]);
			int wt = current + G->distance[currentPath[threadId][length-1]][i];
			bound -= ((float)change/2.0);
			//printf("done\n");
			if(bound + wt < bestLength){
				currentPath[threadId][length] = i;
				visited[threadId][i] = 1;
				_tsp_recursive_serial(G, bound, wt, 
					length+1, threadId);
				// reset visited.
				visited[threadId][i] =  0;
			}
		}
	}
}

void solve_branch_bound(Graph *G, int numThreads){
	if(!G)
		return;
	// initialization of working memory.
	finalPath = malloc(sizeof(int*)*(numThreads+1));
	currentPath = malloc(sizeof(int*)*(numThreads+1));
	visited = malloc(sizeof(int*)*(numThreads+1));
	threadInfo = malloc(sizeof(int)*(numThreads+1));
	memset(threadInfo, -1, sizeof(int)*(numThreads+1));
	dirty = calloc(numThreads+1, sizeof(int));
	firstMin = calloc(G->numCities, sizeof(int));
	secondMin = calloc(G->numCities, sizeof(int));

	for(int i=0; i<numThreads+1; ++i){
		finalPath[i] = calloc(G->numCities, sizeof(int));
		currentPath[i] = calloc(G->numCities, sizeof(int));
		visited[i] = calloc(G->numCities, sizeof(int));
	}

	float current_max = 0.0;
	for(int i=0; i<G->numCities; ++i){
		current_max += _pick_min(G, i);
		current_max += _pick_second(G, i);
		//printf("%f\n", current_max);
	}
	current_max /= 2.0;

	NUM_THREADS = numThreads;
	// threadId = NUM_THREADS, master thread
	int threadId = NUM_THREADS;

	#pragma omp parallel for
	for(int i=0; i<NUM_THREADS; ++i){
		visited[i][0] = 1;
		currentPath[i][0] = 0;
		visited[i][1] = 1;
		currentPath[i][1] = 1;
	}
	_tsp_recursive_parallel(G, current_max, 0, 2, threadId);

	printf("Best path: ");
	for(int i=1; i<G->numCities; ++i){
		printf("%d ",finalPath[bestIndex][i]-1);
	}
	printf("\n Distance: %d\n",bestLength);
}
