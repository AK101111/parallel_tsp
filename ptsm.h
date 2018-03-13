#ifndef PTSM_H_
#define PTSM_H_

typedef struct Graph{
  int numCities;
  int** distance;
}Graph;

/*
  Returns pointer to Graph on success, NULL otherwise.
  Populates Graph from file "filename".
  Assumes filename contains <numCities>*<numCities> numbers.
*/
Graph* populate_data(string fileName, unsigned int numCites);

/*
  Outputs to stdout.
  Best path: <Bestpath>
  Distance: <Path length of Bestpath>
  Bestpath is the permutation of indexes of cities which 
  represents optimum path.
*/
void solve_branch_bound(graph *G);
#endif // PTSM_H_
