#include "ptsm.h"
#include <omp.h>

int *finalPath;	//	best paths explored by different threads.

volatile int bestLength = INT_MAX;	// Current Best, shared among all threads.
volatile int bestIndex =  -1;

unsigned long long int fact[25];

void init_facts(){
	unsigned long long int r = 1;
	for(unsigned long long int i=1; i<25; ++i){
		r *= i;
		facts[i] = r;
	}
}

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
	ret->distance = malloc(sizeof(int*)*(numCities));
	for(int i=0; i<numCities; ++i){
		ret->distance[i] = malloc(sizeof(int)*(numCities));
		for(int j=0; j<numCities; ++j){
			if(!fscanf(file, "%d", &(ret->distance[i][j])))
				return NULL;
		}
	}
	fclose(file);
	return ret;
}

void _tsp(Graph *G){
	for(unsigned long long int k= 0; k<facts[(G->numCities)-1]; ++k){
		int num[20];
		int l = -1;
		for (int i = 0; i < G->numCities; ++i)
        	num[i] = i;
        for (int i = 1; i < G->numCities; ++i) {
        	int facts = fact[(G->numCities)-i-2];
        	int incr = k / facts;
        	int t = num[i+incr];
        	for (int j = i+incr; j > i; j--)
            	num[j] = num[j-1];
        	num[i] = t;
        	k %= facts;
    	}

    	int curr = 0;

    	int bestSample = bestLength;

    	for(int i=0; i<G->numCities-1; ++i){
    		curr += G->distance[num[i]][num[i+1]];
    		if(curr > bestLength){
    			l = 1;
    			break;
    		}
    	}
    	if(l == -1){
    		#pragma omp critical
    		if(curr < bestLength){
    			bestLength = curr;
    			bestIndex = k;
    		}
    	}
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

	finalPath = calloc(G->numCities, sizeof(int));

	omp_set_dynamic(1);
	omp_set_num_threads(NUM_THREADS);
	
	_tsp(G);

	for (int i = 0; i < G->numCities; ++i)
        finalPath[i] = i;

    for (int i = 1; i < G->numCities; ++i) {
        int facts = fact[(G->numCities)-i-2];
        int incr = bestIndex / facts;
        int t = finalPath[i+incr];
        for (int j = i+incr; j > i; j--)
           	finalPath[j] = finalPath[j-1];
       	finalPath[i] = t;
       	bestIndex %= facts;
    }
	
	printf("Best path: ");
	for(int i=0; i<G->numCities; ++i){
		printf("%d ",finalPath[i]);
	}
	printf("\n Distance: %d\n",bestLength);
}
