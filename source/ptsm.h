#ifndef PTSM_H_
#define PTSM_H_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

typedef struct Graph{
  int numCities;
  int** distance;
}Graph;

/*
  Returns pointer to Graph on success, NULL otherwise.
  Populates Graph from file "filename".
  Assumes filename contains <numCities>*<numCities> numbers.
*/
Graph* populate_data(char* fileName, unsigned int numCities);

/*
  Outputs to stdout.
  Best path: <Bestpath>
  Distance: <Path length of Bestpath>
  Bestpath is the permutation of indexes of cities which 
  represents optimum path.
*/
void solve_branch_bound(Graph *G, int numThreads);
#endif // PTSM_H_
