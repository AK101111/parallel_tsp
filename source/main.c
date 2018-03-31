#include "ptsm.h"

int main(int argc, char **argv){
	if(argc < 4){
		printf("Usage ./ptsm x t filename.txt \n Where /
			: \n x is the number of cities \n t is the /
			number of threads \n filename.txt is the /
			file that contains the distance matrix\n");
		return 0;
	}
	int numCities = atoi(argv[1]);
	int numThreads = atoi(argv[2]);
	Graph* G = populate_data(argv[3], numCities);
	if(G){
		solve_branch_bound(G, numThreads);
	}
	return 0;
}