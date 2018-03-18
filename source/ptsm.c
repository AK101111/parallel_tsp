#include "ptsm.h"
#include <omp.h>

int *finalPath;	//	best paths explored by different threads.
//int *currentPath;	//	current path of each thread.
//int *visited;		// visited array for each thread.
int *firstMin;
int *secondMin;

int OFFSET;
int NUM_THREADS;
int LEVELS = 3;

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

__attribute__((always_inline)) static inline int _get(int threadId, int p_level){ return (p_level*OFFSET) + threadId;}

/*
	utility min function, inlined for max perf. 
*/
int _Min(int a, int b){
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
		int second = INT_MAX;
		for(int j=0; j<G->numCities; ++j){
			if(j!=i){
				if(G->distance[i][j] <= min){
					second = min;
					min = G->distance[i][j];
				}else if(G->distance[i][j] <= second && 
					G->distance[i][j] != min){
					second = G->distance[i][j];
				}	
			}
		}
		secondMin[i] = second;
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
	int length, int *currentPath, int *visited);

void _tsp_recursive_parallel(Graph *G, float current_max, int current, 
	int length, int *currentPath, int *visited){
	if(length == G->numCities){
		#pragma omp critical		
		if(current + G->distance[currentPath[length-1]][currentPath[0]] < bestLength){
			bestLength = current + G->distance[currentPath[length-1]][currentPath[0]];
			memcpy(finalPath, currentPath, (G->numCities)*sizeof(int));
		}
		return;
	}
	
	//v#pragma omp parallel for
	//#pragma omp parallel for firstprivate(currentPath, visited)
	#pragma omp parallel for firstprivate(currentPath, visited) if(length<LEVELS)
	for(int i=2; i<G->numCities; ++i){
	//#pragma omp task untied if(length<LEVELS)
	{
		if(!(visited[i])){
			float bound = current_max;
			int change = _pick_min(G,i);
			if(length == 1)
				change += _pick_min(G,currentPath[length-1]);
			else
				change += _pick_second(G,currentPath[length-1]);
			int wt = current + G->distance[currentPath[length-1]][i];
			bound -= ((float)change/2.0);
			
			if(bound + wt < bestLength){
				//int *newPath = alloca(sizeof(int)*(G->numCities));
				//int *newVisited = alloca(sizeof(int)*(G->numCities));
			//	memcpy(newPath, currentPath, sizeof(int)*length);
			//	memcpy(newVisited, visited, sizeof(int)*(G->numCities));
				currentPath[length] = i;
				visited[i] = 1;
				
				_tsp_recursive_parallel(G, bound, wt, 
					length+1, currentPath, visited);
				// reset visited.
				visited[i] =  0;
			}
		}
	}
//	#pragma omp taskwait
	}
}
/*
	Recursive tsp with b&b.
*/

void solve_branch_bound(Graph *G, int numThreads){
	if(!G)
		return;
	// initialization of working memory.
	NUM_THREADS = numThreads;

	firstMin = calloc(G->numCities, sizeof(int));
	secondMin = calloc(G->numCities, sizeof(int));

	finalPath = calloc(G->numCities, sizeof(int));
	int *currentPath = calloc(G->numCities, sizeof(int));
	int *visited = calloc(G->numCities, sizeof(int));
	

	float current_max = 0.0;
	for(int i=0; i<G->numCities; ++i){
		current_max += _pick_min(G, i);
		current_max += _pick_second(G, i);
		//printf("%f\n", current_max);
	}
	current_max /= 2.0;

	omp_set_nested(1);
	omp_set_dynamic(1);
	omp_set_num_threads(NUM_THREADS);
	visited[0] = 1;
	currentPath[0] = 0;
	visited[1] = 1;
	currentPath[1] = 1;

	_tsp_recursive_parallel(G, current_max, 0, 2, currentPath, visited);
	
	printf("Best path: ");
	for(int i=1; i<G->numCities; ++i){
		printf("%d ",finalPath[i]-1);
	}
	printf("\n Distance: %d\n",bestLength);
}
