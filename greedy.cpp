#include "greedy.h"
#include "graph.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace std;

int Greedy::n = 0;
int Greedy::top = 0;
double Greedy::d[MAX_NODE] = {0};
int Greedy::list[MAX_NODE] = {0};
char Greedy::file[] = "greedy.txt";

void Greedy::Build(int num, double (*Run)(int num_iter, int size, int set[]))
{
	n = Graph::GetN();
	top = num;

	bool *used= new bool[n];
	memset(used, 0, sizeof(bool)*n);
	int set[SET_SIZE];

	double old = 0.0;

	double *improve=new double[n];
	int *lastupdate=new int[n];
	int *heap=new int[n];
	for (int i=0; i<n; i++)
	{
		heap[i] = i;
		lastupdate[i] = -1;
		improve[i] = (double)(n+1);//initialize largest
	}

	for (int i=0; i<top; i++)
	{
		int ccc = 0;
		//printf("%d\n",i);
		while (lastupdate[heap[0]] != i)
		{
			//printf("%d %d %d\n", i, heap[0], ccc);
			ccc++;
			lastupdate[heap[0]] = i;
			set[i] = heap[0];
			//printf("GreedyGC_SPM %d %d\n",heap[0],improve[heap[0]]);
			improve[heap[0]] = Run(NUM_ITER, i+1, set) - old;
			//printf("finish\n");

			char tmpfilename[200];
			sprintf_s(tmpfilename, "tmp/%02d%05d.txt", i, heap[0]);
			//FILE *tmpfile;
			//fopen_s(&tmpfile, tmpfilename,"w");
			//fprintf(tmpfile, "%lg\n", improve[heap[0]]); 
			//fclose(tmpfile);

			int x = 0;
			while (x*2+2<=n-i)
			{
				int newx=x*2+1;
				if ((newx+1<n-i) && (improve[heap[newx]]<improve[heap[newx+1]]))
					newx++;
				if (improve[heap[x]]<improve[heap[newx]])
				{
					int t=heap[x];
					heap[x] = heap[newx];
					heap[newx] = t;
					x = newx;
				}
				else
					break;
			}
		}

		used[heap[0]] = true;
		set[i] = heap[0];
		list[i] = heap[0];
		d[i] = improve[heap[0]];
		old+=d[i];

		heap[0] = heap[n-i-1];
		int x = 0;
		while (x*2+2<=n-i)//bug should-1
		{
			int newx=x*2+1;
			if ((newx+1<n-i) && (improve[heap[newx]]<improve[heap[newx+1]]))	//bug should-1
				newx++;
			if (improve[heap[x]]<improve[heap[newx]])
			{
				int t=heap[x];
				heap[x] = heap[newx];
				heap[newx] = t;
				x = newx;
			}
			else
				break;
		}
	}

	FILE *out;
	fopen_s(&out, file, "w");
	fprintf(out, "%d\n", top);
	for (int i=0; i<top; i++)
		fprintf(out, "%d\t%Lg\n", list[i], d[i]);	//the nodes we want!
	fclose(out);
	delete[] heap;
	delete[] lastupdate;
	delete[] improve;
	delete[] used;
}

void Greedy::BuildRanking(int num, double (*Run)(int num_iter, int size, int set[]))
{
	n = Graph::GetN();
	top = num;

	int* set = new int[top];//alert 100
	double* inf = new double[top];
	for(int i=0;i<top;i++)
	{
		set[i] =0;
		inf[i] =0;
	}

	double *improve=new double[n];
	int tmp[1];
	for (int i=0; i<n; i++)
	{
		tmp[0] = i;
		improve[i] = Run(NUM_ITER, 1, tmp);
		if(improve[i] > inf[top-1])
		{
			inf[top-1] = improve[i];
			set[top-1] = i;
			int j = top - 2;
			while(j>=0)
			{
				if(improve[i] > inf[j])
				{
					int int_tmp = set[j];
					double inf_tmp = inf[j];
					set[j] = i;
					inf[j] = improve[i];
					set[j+1] = int_tmp;
					inf[j+1] = inf_tmp;
				}
				else
					break;
				j--;
			}
		}
	}


	FILE *out;
	fopen_s(&out, file, "w");
	fprintf(out, "%d\n", top);
	for (int i=0; i<top; i++)
		fprintf(out, "%d\t%Lg\n", set[i], inf[i]);	//the nodes we want!
	fclose(out);
	//delete[] heap;
	//delete[] lastupdate;
	delete[] improve;
	delete[] set;
	delete[] inf;
	//delete[] used;
}


void Greedy::BuildFromFile(const char* name)
{
	n = Graph::GetN();
	FILE* in;
	fopen_s(&in, name, "r");
	fscanf_s(in, "%ld", &top);
	for (int i=0; i<top; i++)
		fscanf_s(in, "%ld %Lg", &list[i], &d[i]);
	fclose(in);
}

int  Greedy::GetNode(int i)
{
	if (i<0)
		return -1;
	if (i>=top) 
		return -1;
	return list[i];
}


