#include "graph.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <map>
#include <queue>
using namespace std;

bool Graph::built = false;
int  Graph::n = 0;
int  Graph::m = 0;
vector<int>  Graph::index;
vector<int> Graph::degree(MAX_NODE,0);
vector<Edge> Graph::edges;
vector<int> Graph::node;
vector<double> Graph::sigma;
int Graph::target[MAX_K]={0};
int	Graph::targetSize = 0;

void Graph::qsort_edges(int h, int t)
{
	if (h<t) 
	{
		int i = h, j = t;
		Edge mid = edges[(i+j)/2];
		edges[(i+j)/2] = edges[i];

		while (i<j) 
		{
			while ((i<j) && ((edges[j].u>mid.u)||((edges[j].u==mid.u)&&(edges[j].v>mid.v))))
				j--;
			if (i<j) {
				edges[i] = edges[j];
				i++;
			}
			while ((i<j) && ((edges[i].u<mid.u)||((edges[i].u==mid.u)&&(edges[i].v<mid.v))))
				i++;
			if (i<j) {
				edges[j] = edges[i];
				j--;
			}
		}

		edges[i] = mid;
		qsort_edges(h, i-1);
		qsort_edges(i+1, t);
	}
}

void Graph::Build()
{
	if (built) 
		return;
	built = true;

	scanf_s("%ld %ld", &n, &m);	
	degree.resize(n);
	edges.resize(2*m);

	for (int i=0; i<m; i++)
	{
		scanf_s("%ld %ld", &edges[i].u, &edges[i].v);
		edges[i+m].u = edges[i].v;
		edges[i+m].v = edges[i].u;
		edges[i].c	 = 1;
		edges[i+m].c = 1;
		degree[edges[i].u]++;
		degree[edges[i].v]++;
	}

	qsort_edges(0, 2*m-1);

	int m1 = 0;
	for (int i=1; i<2*m; i++)
	{
		if ((edges[i].u != edges[m1].u) || (edges[i].v != edges[m1].v))
		{
			m1++;
			edges[m1] = edges[i];
		}
		else 
		{
			edges[m1].c++;
		}
	}
	if (m!=0)
		m = m1+1;
	
	index.resize(n);
	for (int i=0; i<n; i++)
		index[i] = 0;
	for (int i=0; i<m; i++)
		index[edges[i].u] = i;
	for (int i=1; i<n; i++)
		if (index[i] < index[i-1]) // node i has no out-edge
			index[i] = index[i-1]; 
}

int Graph::GetN()
{
	if (!built)	Build();
	return n;
}

int Graph::GetM()
{
	if (!built)	Build();
	return m;
}

int Graph::GetDegree(int node)
{
	if (!built)	Build();
	return degree[node];
}

int Graph::GetNeighbor(int node)
{
	if (!built)	Build();
	if (node == 0)
		return index[node]+1;
	else 
		return index[node]-index[node-1];
}

Edge Graph::GetEdge(int node, int idx)
{
	if (!built)	Build();
	if (node == 0)
		return edges[idx];
	else
		return edges[index[node-1]+1+idx];
}

int Graph::GetEdge_idx(int node, int idx)
{
	if (!built)	Build();
	if (node == 0)
		return idx;
	else
		return index[node-1]+1+idx;
}

//----------------------------------start----------------------------------//
// for uncertain information influence

void Graph::BuildUIC()
{
		if (built) 
		return;
	built = true;

	scanf_s("%ld %ld", &n, &m);
	degree.resize(n);
	node.resize(n);
	sigma.resize(n);
	edges.resize(2*m);

	for (int i=0; i<2*m; i++)
	{
		scanf_s("%ld %ld %lg %lg", &edges[i].u, &edges[i].v,
			&edges[i].truevalue1, &edges[i].truevalue2);
		edges[i].u--;
		edges[i].v--;
		edges[i].truevalue1 = -log(edges[i].truevalue1);
		edges[i].truevalue2 = -log(edges[i].truevalue2);
		edges[i].observevalue1 = edges[i].truevalue1;
//		edges[i].delta1 = INIT_DELTA_PARA;
		edges[i].delta1 = -100;
		edges[i].times1 = 0;
		edges[i].counter1 = 0;
		edges[i].w1 = 0;
		edges[i].w2 = 0;
		edges[i].c = 1;
		degree[edges[i].u]++;
	}

	qsort_edges(0, 2*m-1);

	int m1 = 0;
	for (int i=1; i<2*m; i++)
	{
		if ((edges[i].u != edges[m1].u) || (edges[i].v != edges[m1].v))
		{
			m1++;
			edges[m1] = edges[i];
		}
		else 
		{
			edges[m1].c++;
		}
	}
	if (m!=0)
		m = m1+1;
	
	index.resize(n);
	for (int i=0; i<n; i++)
		index[i] = 0;
	for (int i=0; i<m; i++)
		index[edges[i].u] = i;
	for (int i=1; i<n; i++)
		if (index[i] < index[i-1])
			index[i] = index[i-1];


	//initial sample on edges
	for (int i=0;i<m;i++)
	{
		SampleEdge(&edges[i], INIT_SAMPLE_NUMBER);
		update_edge(&edges[i]);
	}

}

void Graph::update_edge(Edge *e)
{
	if(e->times1 == 0)
	{
		e->observevalue1 = 0;
	}
	else
	{
		// Chernoff bound
		if (e->counter1 == 0)
			e->observevalue1 = 9*log(10.0);
		else
			e->observevalue1 = -log((double)e->counter1/(double)e->times1);


//		double delta = 0.5*log((double)e->times1) - log(CHER_NUM);
		double delta = 0.5*(log((double)e->times1)-e->truevalue1-log(log((double)m))-log(4.0));
		if(delta > e->delta1)
			e->delta1 = delta;
	}
}

void Graph::set_edge_up(Edge *e)
{
		//update w1

		//multiplicative bound
		e->w1 = min(exp(-e->observevalue1)*(1+exp(-e->delta1)), 1.0);
		//additive bound
		e->w1 = min(exp(-e->observevalue1) + sqrt(log((double)m)/e->times1), e->w1);
		
		e->w1 = -log(e->w1);

		//update w2
		int d = degree[e->v];
		Edge e1;
		for(int j=0; j < d; j++)
			if (GetEdge(e->v, j).v == e->u)
				e1 = GetEdge(e->v, j);
		e->w2 = min(exp(-e1.observevalue1)*(1+exp(-e1.delta1)), 1.0);
		e->w2 = min(exp(-e1.observevalue1) + sqrt(log((double)m)/e1.times1), e->w2);
		
		e->w2 = -log(e->w2);

}

void Graph::set_all_edge_up()
{
	for (int i=0;i<m;i++)
		set_edge_up(&edges[i]);
}

void Graph::set_edge_down(Edge *e)
{
		//update w1

		//multiplicative bound
		e->w1 = max(exp(-e->observevalue1)*(1-exp(-e->delta1)), 1e-09);
		//additive bound
		e->w1 = max(exp(-e->observevalue1) - sqrt(log((double)m)/e->times1), e->w1);
		
		e->w1 = -log(e->w1);

		//update w2
		int d = degree[e->v];
		Edge e1;
		for(int j=0; j < d; j++)
			if (GetEdge(e->v, j).v == e->u)
				e1 = GetEdge(e->v, j);
		e->w2 = max(exp(-e1.observevalue1)*(1-exp(-e1.delta1)), 1e-09);
		e->w2 = max(exp(-e1.observevalue1) - sqrt(log((double)m)/e1.times1), e->w2);
		
		e->w2 = -log(e->w2);
}


void Graph::set_all_edge_down()
{
	for (int i=0;i<m;i++)
		set_edge_down(&edges[i]);
}


void Graph::set_edge_CLUCB(int size, int set[])
{
	set_all_edge_up();
	for(int i=0; i < size; i++)
	{
		int tempnode = set[i];
		int d = degree[tempnode];

		for(int j=0; j < d; j++)
			set_edge_down(&edges[GetEdge_idx(tempnode, j)]);
	}
}



void Graph::initial_interval()
{
	Edge *e;
	for(int i=0; i < m; i++)
	{
		e = &edges[i];
		e->delta1 = (e->truevalue1)+INIT_DELTA_PARA;
	}
}

void Graph::update_interval(Edge *e)
{
	e->delta1 = (e->delta1)+UPDATE_DELTA_PARA;
}

void Graph::set_interval_up()
{
	Edge* e;
	for(int i=0; i<m; i++)
	{
		e = &edges[i];
		e->w1 = min(exp(-e->truevalue1) + exp(-e->delta1), 1.0);
		e->w1 = -log(e->w1);
	}
}

void Graph::set_interval_down()
{
	Edge *e;
	for(int i=0; i<m; i++)
	{
		e = &edges[i];
		e->w1 = max(exp(-e->truevalue1) - exp(-e->delta1), 1e-09);
		e->w1 = -log(e->w1);
	}
}

void Graph::iteration_Wei(int threshold, double (*Run)(int num_iter, int size, int set[]))
{
/*
	FILE *fp;
	int set[1];
	int tempnode;
	printf("start calculating sigma for up-prob.\n");
	set_all_edge_up();
	for(int i=0; i < n; i++)
	{
		set[0] = i;
		sigma[i] = Run(NUM_ITER_1, 1, set);
		node[i] = i;
		if (i%(n/10)==0)
			printf("%d\n", i*10/n);
	}

	printf("start calculating sigma for low-prob.\n");
	set_all_edge_down();
	for(int i=0; i < n; i++)
	{
		set[0] = i;
		//difference
		sigma[i] = sigma[i] - Run(NUM_ITER_1, 1, set);
		//ratio
//		sigma[i] = (sigma_up[i] - sigma_down[i])/sigma_up[i];
		node[i] = i;
		if (i%(n/10)==0)
			printf("%d\n", i*10/n);
	}
	printf("finish calculating sigma.\n");


	int top[THRESHOLD];
	for(int i=0; i < THRESHOLD; i++)
		top[i] = i;

	for(int i=0; i < THRESHOLD; i++)
	{
		for(int j=i; j > 0; j--)
		{
			if(sigma[top[j]] > sigma[top[j-1]])
			{
				tempnode = top[j];
				top[j] = top[j-1];
				top[j-1] = tempnode;
			}
			else
				break;
		}
	}

	for(int i=THRESHOLD; i < n; i++)
	{
		if(sigma[i] > sigma[top[THRESHOLD-1]])
		{
			top[THRESHOLD-1] = i;
			for(int j=THRESHOLD-1; j > 0; j--)
			{
				if(sigma[top[j]] > sigma[top[j-1]])
				{
					tempnode = top[j];
					top[j] = top[j-1];
					top[j-1] = tempnode;
				}
				else
					break;
			}
		}
	}



	int d;
	Edge e;
	fopen_s(&fp, "update_node.txt", "a");

	for(int i=0; i < THRESHOLD; i++)
	{
		tempnode = top[i];
		d = degree[tempnode];

		printf("node: %d\t deltasigma: %lg\n", tempnode, sigma[tempnode]);
		fprintf(fp, "node: %d\t deltasigma: %lg\n", tempnode, sigma[tempnode]);

		for(int j=0; j < d; j++)
		{
			SampleEdge(&edges[GetEdge_idx(tempnode, j)], SAMPLE_PER_ITER);
			update_edge(&edges[GetEdge_idx(tempnode, j)]);
			e = GetEdge(tempnode, j);
			printf("edge:%d->%d\t %lg\n", e.u, e.v, e.delta1);//
		}
	}
	fprintf(fp, "\n");
	fclose(fp);



*/


// for uniform sampling


for (int i=0;i<m;i++)
{
	SampleEdge(&edges[i],SAMPLE_PER_ITER);
	update_edge(&edges[i]);
}




}


double Graph::CalSampleNum()   //calculate sample number, given all the edge[i].delta1
{
//	printf("Start Calculating Sample Number!");
	double out=0;
	for (int i=0;i<m;i++)
//		if ((exp(-edges[i].truevalue1)>0.0001) && (edges[i].delta1!=0))
//			out=out+4*log(double(m))/exp(-edges[i].truevalue1)/exp(-edges[i].delta1)/exp(-edges[i].delta1);
		out=out+edges[i].times1;
//	printf("Finish Calculating Sample Number!");
	return out;
}


void Graph::SampleEdge(Edge *e, int num)
{
	e->times1+=num;
	for (int i=0;i<num;i++)
		if (((double)rand()/(double)RAND_MAX) < exp(-e->truevalue1))
			e->counter1++;
}



void Graph::InformationCascadeSample(int num_iter, int size, int set[])
{
	targetSize = size;
	for (int i=0; i<size; i++)
		target[i] = set[i];

	int		h, t;
	int		*list  = new int[n];
	bool	*active= new bool[n];

	printf("start stimulate.\n");
	for (int it=0; it<num_iter; it++)
	{
//		printf("%d\n", it);
		memset(active, 0, sizeof(bool)*n);
		for (int i=0; i<targetSize; i++) 
		{
			list[i] = target[i];
			active[target[i]] = true;
		}

		h = 0; // current active node
		t = targetSize;

		while (h<t) 
		{
			int k = GetNeighbor(list[h]); // number of neighbors of h
			
			for (int i=0; i<k; i++)
			{
				int idx = GetEdge_idx(list[h], i); // Edge e = Graph::GetEdge(list[h], i);
				edges[idx].times1 ++;
				// we need to ative multi-times // if (active[edges[idx].v]) continue;
				if (((double)rand()/(double)RAND_MAX) < exp(-edges[idx].truevalue1))
				{
					edges[idx].counter1 ++;
					if(!active[edges[idx].v])
					{
						list[t] = edges[idx].v;
						active[edges[idx].v] = true;
						t++;
					}
				}

				// update edges[idx]
				update_edge(&edges[idx]);
			}
			h++;
		}
	}
	delete[] active;
	delete[] list;
}



void Graph::iteration_IC(int num_iter, double (*Run)(int num_iter, int size, int set[]))
{
	int set[1];
	int tempnode;
	printf("start running one node.\n");
	set_all_edge_up();
	for(int i=0; i < n; i++)
	{
		set[0] = i;
		sigma[i] = Run(NUM_ITER_1, 1, set);
		node[i] = i;
		if (i%(n/10)==0)
			printf("%d\n", i*10/n);
	}
	set_all_edge_down();
	for(int i=0; i < n; i++)
	{
		set[0] = i;
		//difference
		sigma[i] = sigma[i] - Run(NUM_ITER_1, 1, set);
		//ratio
//		sigma[i] = (sigma_up[i] - sigma_down[i])/sigma_up[i];
		node[i] = i;
		if (i%(n/10)==0)
			printf("%d\n", i*10/n);
	}
	printf("finish running one node.\n");


	int top[THRESHOLD];
	for(int i=0; i < THRESHOLD; i++)
		top[i] = i;

	for(int i=0; i < THRESHOLD; i++)
	{
		for(int j=i; j > 0; j--)
		{
			if(sigma[top[j]] > sigma[top[j-1]])
			{
				tempnode = top[j];
				top[j] = top[j-1];
				top[j-1] = tempnode;
			}
			else
				break;
		}
	}

	for(int i=THRESHOLD; i < n; i++)
	{
		if(sigma[i] > sigma[top[THRESHOLD-1]])
		{
			top[THRESHOLD-1] = i;
			for(int j=THRESHOLD-1; j > 0; j--)
			{
				if(sigma[top[j]] > sigma[top[j-1]])
				{
					tempnode = top[j];
					top[j] = top[j-1];
					top[j-1] = tempnode;
				}
				else
					break;
			}
		}
	}
	Graph::InformationCascadeSample(num_iter, THRESHOLD, top);
}




void Graph::GeneralSample(int num_iter, int size, int set[])
{
	int d;
	Edge e;
//	FILE *fp;

//	fopen_s(&fp, "update_node.txt", "a");

	printf("start general sample!\n");
	for(int i=0; i < size; i++)
	{
		int tempnode = set[i];
		d = degree[tempnode];

//		printf("node: %d\t deltasigma: %lg\n", tempnode, sigma[tempnode]);
//		fprintf(fp, "node: %d\t deltasigma: %lg\n", tempnode, sigma[tempnode]);

		for(int j=0; j < d; j++)
		{
			SampleEdge(&edges[GetEdge_idx(tempnode, j)],SAMPLE_PER_ITER);
			update_edge(&edges[GetEdge_idx(tempnode, j)]);
			e = GetEdge(tempnode, j);
//			printf("edge:%d->%d\t %lg\n", e.u, e.v, e.delta1);
		}
	}	


}



void Graph::print_graph(char *file)
{
	FILE *out;
	fopen_s(&out, file, "w");
	int u,v, counter, sample;
	Edge e;
	double down,up, truevalue, observe, w1, w2, delta;
	fprintf(out, "%d \t %d\n",n, m);
	for(int i=0; i < m; i++)
	{
		e = edges[i];
		u = e.u;
		v = e.v;
		down = max(exp(-e.observevalue1)*(1-exp(-e.delta1)), 1e-09);
		down = max(exp(-e.observevalue1) - sqrt(log((double)m)/e.times1), down);
		up = min(exp(-e.observevalue1)*(1+exp(-e.delta1)), 1.0);
		up = min(exp(-e.observevalue1) + sqrt(log((double)m)/e.times1), up);
		truevalue = exp(-e.truevalue1);
		counter = e.counter1;
		sample = e.times1;
		observe = exp(-e.observevalue1);
		w1 = exp(-e.w1);
		w2 = exp(-e.w2);
		delta = exp(-e.delta1);
		fprintf(out, "%d \t %d \t [%lg, %lg] \t\t %lg \t\t %lg \t\t %lg \t\t %lg \t %lg \t %d \t %d\n", u, v, down, up, truevalue, observe, w1, w2, delta, counter, sample);
	}
	fclose(out);
}


//----------------------------------end---------------------------------------//

void Graph::Build2WC()
{
	if (built) 
		return;
	built = true;

	scanf_s("%ld %ld", &n, &m);	
	degree.resize(n);
	edges.resize(2*m);

	for (int i=0; i<2*m; i++)
	{
		scanf_s("%ld %ld %lg %lg", &edges[i].u, &edges[i].v, &edges[i].w1, &edges[i].w2);
		edges[i].u--;
		edges[i].v--;
		edges[i].w1=-log(edges[i].w1);
		edges[i].w2=-log(edges[i].w2);
		//edges[i+m].u = edges[i].v;
		//edges[i+m].v = edges[i].u;
		edges[i].c	 = 1;
		//edges[i+m].c = 1;
		degree[edges[i].u]++;
		//degree[edges[i].v]++;
	}

	qsort_edges(0, 2*m-1);

	int m1 = 0;
	for (int i=1; i<2*m; i++)
	{
		if ((edges[i].u != edges[m1].u) || (edges[i].v != edges[m1].v))
		{
			m1++;
			edges[m1] = edges[i];
		}
		else 
		{
			edges[m1].c++;
		}
	}
	if (m!=0)
		m = m1+1;
	
	index.resize(n);
	for (int i=0; i<n; i++)
		index[i] = 0;
	for (int i=0; i<m; i++)
		index[edges[i].u] = i;
	for (int i=1; i<n; i++)
		if (index[i] < index[i-1])
			index[i] = index[i-1];
		//else
		//	printf("%d %d\n",edges[i].u, index[i]);
	//BuildWC();
}

void Graph::Stats()
{
	printf("number of vertices:\t%d\n",n);
	printf("number of edges:\t%d\n",m/2);
	printf("density:\t%lg\n",double(m)/n/(n-1));
	int maxdegree=0;
	double tdegree=0.0;
	int i,j,k;
	for (i=0;i<n;i++)
	{
		if (degree[i]>maxdegree) maxdegree=degree[i];
		tdegree+=degree[i];
		//if (degree[i]%2) printf("%d\n", i);
	}
	printf("average degree:\t%lg\n",tdegree/n);
	printf("maximal degree:\t%d\n",maxdegree);
	int maxcmp=0,ncmp=0;
	bool *used=new bool[n];
	memset(used,0,n);
	while (1)
	{
		queue<int> q;
		for (i=0;i<n;i++)
			if (!used[i]) break;
		if (i==n) break;
		ncmp++;
		int cmpsize=0;
		q.push(i);
		used[i]=true;
		while (!q.empty())
		{
			k=q.front();
			q.pop();
			cmpsize++;
			j=GetNeighbor(k);
			for (i=0;i<j;i++)
			{
				Edge e=GetEdge(k,i);
				if (used[e.v]) continue;
				q.push(e.v);
				used[e.v]=true;
			}
		}
		if (cmpsize>maxcmp) maxcmp=cmpsize;
	}
	printf("# of connected component:\t%d\n",ncmp);
	printf("largest component size:\t%d\n",maxcmp);
	printf("average component size:\t%lg\n",double(n)/ncmp);
}
