#include "ptsm.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int **finalPath;	//	best paths explored by different threads.
int **currentPath;	//	current path of each thread.
int **visited;
int *dirty;

volatile int bestLength = 0;	// Current Best, shared among all threads.

/*
	Keeping bestLength volatile, however finalPath is not.
	Instead each thread has its own path array to minimize time
	spent in locked zone. (updating bestLength only).
*/

Graph* populate_data(string fileName, unsigned int numCites){
	FILE *file;
	file=fopen(fileName, "r");
	if(!file)
		return NULL;
	
	Graph* ret = malloc(sizeof(Graph));
	if(!ret)
		return NULL;
	
	ret->numCites = numCites;
	for(int i=0; i<numCites; ++i){
		for(int j=0; j<numCites; ++j){
			if (!fscanf(file, "%d", &(ret->distance[i][j]))) 
				return NULL;
		}
	}
	return ret;
}

inline int Min(int a, int b){
	return (a<b)?(a):(b);
}

int _pick_next(Graph *G, int city){
	int min = INT_MAX;
	for(int i=0; i<G->numCites; ++i)
		if(i != city)
			min = Min(min, G->distance[city][i]);
	return min;
}

int _pick_next(Graph* G, int city, int length){
	if(length == 1)
		return _pick_next(G,city);

	// fancy one pass second max of array
	int min = INT_MAX;
	int secondMin = INT_MAX;
	for(int i=0; i<G->numCites; ++i){
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

void _tsp_recursive(graph *G, float current_max, int current, 
	int length, int threadId){

	if(length == G->numCites){
		// #pragma omp critical
		if(current < bestLength){
			// add lock safety
			bestLength = current;
			dirty[threadId] = 1;
		}

		// hack to move heavy computation out of lock
		if(dirty[threadId]){
			memcpy(finalPath[threadId], currentPath[threadId], 
				sizeof(currentPath[threadId]));
			dirty[threadId] = 0;
		}
	}

	for(int i=0; i<G->numCites; ++i){
		if(!visited[threadId][i]){
			float bound = current_max;
			int change =  _pick_next(G,i) + 
						_pick_next(G,currentPath[threadId[length-1]], length));
			int wt = current + G->distance[currentPath[threadId][length-1]][i];

			if(bound - change + wt < bestLength){
				currentPath[threadId][level] = i;
				visited[threadId][i] = 1;
				_tsp_recursive(G, bound - (float)change/2, wt, 
					length+1, threadId);
				// reset visited.
				visited[threadId][i] =  0;
			}
		}
	}
}

void tsp(graph *G, int numThreads){
	if(!G)
		return;
	// initialization of working memory.
	finalPath = malloc(sizeof(int*)*numThreads);
	currentPath = malloc(sizeof(int*)*numThreads);
	visited = malloc(sizeof(int*)*numThreads);
	dirty = calloc(numThreads, sizeof(int));
	for(int i=0; i<numThreads; ++i){
		finalPath[i] = calloc(G->numCites, sizeof(int));
		currentPath[i] = calloc(G->numCites, sizeof(int));
		visited[i] = calloc(G->numCites, sizeof(int));
	}

	float current_max = 0.0;
	for(int i=0; i<G->numCites; ++i){
		current_max += (float)_pick_next(G, i, 1)/2;
		current_max += (float)_pick_next(G, i, 2)/2;
	}

	visited[threadId][0] = 1;
	currentPath[threadId][0] = 0;
	_tsp_recursive(G, current_max, 0, 1, 0);
}