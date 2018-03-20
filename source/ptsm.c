#include "ptsm.h"
#include <omp.h>

int *finalPath;	//	best paths explored by different threads.

volatile int bestLength = INT_MAX;	// Current Best, shared among all threads.
volatile long long int bestIndex =  -1;
int NUM_THREADS=1;
long long int fact[25];

void init_facts(){
	long long int r = 1;
	fact[0] = 1;
	for(long long int i=1; i<20; ++i){
		r *= i;
		fact[i] = r;
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

void print(int a[], int n){
	for(int i=0; i<n;++i)
		printf("%d ",a[i]);
	printf("\n");
}

void _tsp(Graph *G){
	#pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
	for(long long int k= 0; k<fact[(G->numCities)-1]; ++k){
		int num[20];
		int curr = 0;
		int bestSample = bestLength;
		int l = -1;
		for (int i = 0; i < G->numCities; ++i)
        		num[i] = i;
        	long long int kk = k;
		for (int i = 1; i < G->numCities; ++i) {
        		long long int facts = fact[(G->numCities)-i-1];
        		int incr = kk / facts;
	        	int t = num[i+incr];
        		for (int j = i+incr; j > i; j--)
            			num[j] = num[j-1];
	        	num[i] = t;
			curr += G->distance[num[i-1]][num[i]];
			if(curr > bestSample){
				l = 1;
				break;
			}
        		kk %= facts;
    		}
		if(l == 1)
			continue;
		//print(num,G->numCities);
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
	init_facts();
	omp_set_dynamic(1);
	omp_set_num_threads(numThreads);
	_tsp(G);
	
	for (int i = 0; i < G->numCities; ++i)
        	finalPath[i] = i;
	
	for (int i = 1; i < G->numCities; ++i) {
		int facts = fact[(G->numCities)-i-1];
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
