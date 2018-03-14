#include "ptsm.h"
//#include <omp.h>

int **finalPath;	//	best paths explored by different threads.
int **currentPath;	//	current path of each thread.
int **visited;		// visited array for each thread.
int *dirty; 		// array which stores dirty bit of each thread
					// denoting if that thread needs to update its 
					// final path from its local path.

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
	
	ret->numCities = numCities;
	ret->distance = malloc(sizeof(int*)*numCities);
	for(int i=0; i<numCities; ++i){
		ret->distance[i] = malloc(sizeof(int)*numCities);
		for(int j=0; j<numCities; ++j){
			if(!fscanf(file, "%d", &(ret->distance[i][j])))
				return NULL;
		}
	}
	fclose(file);
	return ret;
}

/*
	utility min function, inlined for max perf. 
*/
int _Min(int a, int b){
	return (a<b)?(a):(b);
}

/*
	Picks min weight outgoing edge from node
*/
int _pick_min(Graph *G, int city){
	int min = INT_MAX;
	for(int i=0; i<G->numCities; ++i)
		if(i != city)
			min = _Min(min, G->distance[city][i]);
	return min;
}

/*
	Picks min/second min weight outgoing edge from node
*/
int _pick_next(Graph* G, int city, int length){
	if(length == 1)
		return _pick_min(G,city);

	// fancy one pass second max of array
	int min = INT_MAX;
	int secondMin = INT_MAX;
	for(int i=0; i<G->numCities; ++i){
		if(i!=city){
			if(G->distance[city][i] <= min){
				secondMin = min;
				min = G->distance[city][i];
			}else if(G->distance[city][i] <= secondMin && 
				G->distance[city][i] != min){
				secondMin = min;
			}
		}
	}
	return secondMin;
}

/*
	Recursive tsp with b&b.
*/
void _tsp_recursive(Graph *G, float current_max, int current, 
	int length, int threadId){
	printf("%lf %d %d\n", current_max, current, length);
	if(length == G->numCities){
		// #pragma omp critical
		if(current < bestLength){
			// add lock safety
			bestLength = current;
			bestIndex = threadId;
			dirty[threadId] = 1;
		}

		// hack to move heavy computation out of lock
		if(dirty[threadId]){
			memcpy(finalPath[threadId], currentPath[threadId], (G->numCities)*sizeof(int));
			dirty[threadId] = 0;
		}
	}

	for(int i=0; i<G->numCities; ++i){
		if(!visited[threadId][i]){
			float bound = current_max;
			int change =  _pick_min(G,i) + 
						_pick_next(G,currentPath[threadId][length-1], length);
			int wt = current + G->distance[currentPath[threadId][length-1]][i];
			//printf("done\n");
			if(bound - change + wt < bestLength){
				currentPath[threadId][length] = i;
				visited[threadId][i] = 1;
				_tsp_recursive(G, bound - (float)change/2, wt, 
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
	finalPath = malloc(sizeof(int*)*numThreads);
	currentPath = malloc(sizeof(int*)*numThreads);
	visited = malloc(sizeof(int*)*numThreads);
	dirty = calloc(numThreads, sizeof(int)*numThreads);

	for(int i=0; i<numThreads; ++i){
		finalPath[i] = calloc(G->numCities, sizeof(int));
		currentPath[i] = calloc(G->numCities, sizeof(int));
		visited[i] = calloc(G->numCities, sizeof(int));
	}

	float current_max = 0.0;
	for(int i=0; i<G->numCities; ++i){
		current_max += (float)_pick_next(G, i, 1)/2;
		current_max += (float)_pick_next(G, i, 2)/2;
	}

	// threadId
	int threadId = 0;
	visited[threadId][0] = 1;
	currentPath[threadId][0] = 0;
	_tsp_recursive(G, current_max, 0, 1, 0);

	printf("Best path: ");
	for(int i=0; i<G->numCities; ++i){
		printf("%d ",finalPath[bestIndex][i]);
	}
	printf("\n Distance: %d\n",bestLength);
}