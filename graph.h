#ifndef GRAPH_H
#define GRAPH_H

#include "limit.h"
#include <vector>

using namespace std;

struct Edge 
{
	int u,v,c;	
	double w1,w2;

	// for uncertain information influence
	double truevalue1;
	double truevalue2;
	double observevalue1; // counter1/time1
	double delta1;   // for multiplicative chernoff bound    [p-exp(-delta1),p+exp(-delta1)]
	long int times1; // sample time
	long int counter1; // observe time
};

class Graph
{
private:
	static bool built;
	static int n;
	static int m;
	static vector<int> degree;
	static vector<int> index;
	static vector<Edge> edges;

	static void qsort_edges(int h, int t);

	// for uncertain information influence
	static vector<double> sigma;
	static vector<int> node;
	static int target[MAX_K];
	static int	targetSize;

public:
	static void	Build();
	static int	GetN();
	static int	GetM();
	static int	GetDegree(int node);
	static int	GetNeighbor(int node);
	static Edge	GetEdge(int node, int idx);
	static void Build2WC();
	static void Stats();

	// for uncertain information influence
	static void BuildUIC();

	// sample
	static void update_edge(Edge *);

	// set influence value
	static void set_edge_up(Edge *);
	static void set_edge_down(Edge *);
	static void set_all_edge_up();
	static void set_all_edge_down();
	static void set_edge_CLUCB(int size, int set[]);

	// update interval
	static void initial_interval();
	static void update_interval(Edge *);
	static void set_interval_up();
	static void set_interval_down();

	static void iteration_Wei(int t,double (*Run)(int num_iter, int size, int set[]));
	static void iteration_IC(int num_iter, double (*Run)(int num_iter, int size, int set[]));
	static int GetEdge_idx(int node, int idx);
	static void InformationCascadeSample(int num_iter, int size, int set[]);
	static void GeneralSample(int num_iter, int size, int set[]);	


	static double CalSampleNum();
	static void SampleEdge(Edge *e, int num);


	static void print_graph(char *);
};

#endif

