#include "limit.h"
#include "graph.h"
#include "greedy.h"
#include "pmia.h"
#include "independ_cascade.h"
#include "general_cascade.h"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <io.h>
#include <string>
#include <windows.h>
#include <set>
#include <crtdbg.h> 
#include <ctime>
using namespace std;

FILE* timetmpfile;
double timer;
LARGE_INTEGER start, Eend, freq;

void toSimulate(char *file, int (*GetNode)(int i), double (*Run)(int num_iter, int size, int set[]))
{
	FILE *out;
	fopen_s(&out, file, "w");
	int set[MAX_NODE];
	for (int t=0; t<SET_SIZE; t++)
	{
		set[t] = GetNode(t);
		fprintf(out, "%02d \t %10Lg\n", t+1, Run(NUM_ITER, t+1, set));
	}
	fclose(out);
}

double toSimulateOnce(int setsize, int (*GetNode)(int i), double (*Run)(int num_iter, int size, int set[]), int num_iter_test=NUM_ITER)
{
	int set[MAX_NODE];
	int t;
	for (t=0; t<setsize; t++)
	{
		set[t] = GetNode(t);
	}
	return Run(num_iter_test, t, set);
}

int main(int argc, char * argv[])
{
	srand((unsigned)time(NULL));
	QueryPerformanceFrequency(&freq);
	system("mkdir tmp");
	system("cd tmp");

	if (argc<=1) { 
		printf("-st : statistics of the weighted cascade graph\n");
		printf("-stg : statistics of the general independent cascade graph\n");
		printf("-b : baseline: random (-br), degree (-bd), degreediscount (-bdd), weighteddegree (-bw), pagerank (-bp) for general ic\n");
		printf("-g : greedy, SPM and SP1M for general ic\n");
		printf("-go : greedy with online bound for general ic\n");
		printf("-p bound1 bound2 : PMIA with 1/theta from bound1 to bound 2\n");
		printf("-m bound1 bound2 : MIA with 1/theta from bound1 to bound 2\n");
		printf("-gwc filename1 filename2... : c-greedy algorithm\n");
		printf("-gwcl seeds_list_file : c-greedy algorithm\n");
		printf("-mis seeds_list_file dist1 dist2 ... distd : MIS(marginal influence sort) algorithm \n");
		printf("-ts seeds_list_file dist1 dist2 ... distd : BTS(top-selection) algorithm \n");
		printf("-t seeds_file num_iter: use seeds and test on general ic model \n");
		printf("");
		printf("example: max_influence -p 20 2000 < hep.inf \n"); // read from file
		return 0;
	}
	system("del /Q _finished_.log"); 
	system("echo. 2> _running_.log");
	

	string s;

	s="-t";
	if(!s.compare(argv[1])){
		Graph::Build2WC();
		GeneralCascade::Build();
		char* default_file = "Test.txt";
		char* seeds_file = default_file;
		if (argc >= 3) {
			seeds_file = argv[2];
		}
		int num_iter = NUM_ITER;
		if (argc >= 4) {
			num_iter = atoi(argv[3]);
		}
		ifstream ft;
		ft.open(seeds_file);
		
		int temp[220]; //alert
		int nums;
		ft >> nums;//# of seeds
		double margin;//alert
		for(int i=0; i<nums; i++) {
			ft >> temp[i] >> margin;
		}

		FILE *out;
		char spreadfilename[]="spread.txt";
		fopen_s(&out, spreadfilename, "w");
		for (int i=1; i <= SET_SIZE; ++i) {
			double spread = 0;
			int set_size = i;
			if (set_size > nums) break;
			spread = GeneralCascade::Run(num_iter, set_size, temp);
			fprintf(out, "%d\t%f\n", set_size, spread);
		}
		fclose(out);
	}
	
		//----------------------------------start-------------------------------//
	// for uncertain information influence
	s="-w";
	if(!s.compare(argv[1])){

		int bound=320;
	
		double sigmaup;
		double sigmadown;
		double alpha;
		
		FILE *seed_out;

		Graph::BuildUIC();
//		Graph::initial_interval(); //	Graph::update_all_edge(INIT_TIMES);
		GeneralCascade::Build();

		FILE *alpha_out;
		fopen_s(&alpha_out, "uncertain_influence.txt", "w");

		Graph::set_all_edge_up();
		Graph::print_graph("uncertain_graph0.txt");
		
		printf("Running Pmia!\n");
		printf("%lg\n", SPT_new::Build(SET_SIZE,bound));
		printf("Pmia Complete!\n");
		fopen_s(&seed_out, "seed_set.txt", "a");
		fprintf(seed_out, "iteration 0:\n");
		for(int i=0; i < SET_SIZE; i++)
			fprintf(seed_out, "%d ", SPT_new::GetNode(i));
		fprintf(seed_out, "\n");
		
		sigmaup = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
		Graph::set_all_edge_down();
		sigmadown = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
		alpha = sigmadown/sigmaup;
		fprintf(alpha_out, "iteration 0:\n");
		fprintf(alpha_out, "%lg\n", alpha);
		fprintf(alpha_out, "%lg\n", Graph::CalSampleNum());

		printf("heihei\n");
		
		for(int i=0; i < ITER_TIMES; i++)
		{
			// intern update
			for(int j=0; j < ITER_TIMES_2; j++)
			{
				Graph::iteration_Wei(THRESHOLD, GeneralCascade::Run);
/*				Graph::set_all_edge_up();
				SPT_new::Build(SET_SIZE,bound);
				int greedy_set[SET_SIZE];
				for (int k=0;k<SET_SIZE;k++)
					greedy_set[k]=SPT_new::GetNode(k);
				Graph::GeneralSample(SAMPLE_PER_ITER, SET_SIZE, greedy_set);
*/				
				Graph::set_all_edge_up();
				sigmaup = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
				Graph::set_all_edge_down();
				sigmadown = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
				alpha = sigmadown/sigmaup;
				fprintf(alpha_out, "%lg\n", alpha);
				fprintf(alpha_out, "%lg\n", Graph::CalSampleNum());
				printf("Sub-iteration Complete! %d,%d\n",i,j);
			}
			fprintf(alpha_out, "\n");

			
			printf("%d\n", i);
			// greedy
			Graph::set_all_edge_up();
			printf("Running Pmia!\n");
			SPT_new::Build(SET_SIZE,bound);
			printf("Pmia Complete!\n");
			
			// print seed
			fprintf(seed_out, "iteration %d:\n", i+1);
			for(int k=0; k < SET_SIZE; k++)
				fprintf(seed_out, "%d ", SPT_new::GetNode(k));
			fprintf(seed_out, "\n");
		
			// calculate alpha
			sigmaup = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
			Graph::set_all_edge_down();
			sigmadown = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
			alpha = sigmadown/sigmaup;			
			fprintf(alpha_out, "iteration %d:\n", i+1);
			fprintf(alpha_out,"%lg\n", alpha);

				
		}
		
		fclose(seed_out);
		fclose(alpha_out);
		
	}

	// information cascade sample
	s="-c";
	if(!s.compare(argv[1])){

		int bound=320;
	
		double sigmaup;
		double sigmadown;
		double alpha;
		int seed_set[SET_SIZE];
		
		FILE *seed_out;

		Graph::BuildUIC();
//		Graph::initial_interval(); //	Graph::update_all_edge(INIT_TIMES);
		GeneralCascade::Build();

		FILE *alpha_out;
		fopen_s(&alpha_out, "uncertain_influence.txt", "w");

		Graph::set_all_edge_up();
		Graph::print_graph("uncertain_graph0.txt");
		
		SPT_new::Build(SET_SIZE,bound);
		
		printf("Pmia finished.\n");

		fopen_s(&seed_out, "seed_set.txt", "a");
		fprintf(seed_out, "iteration 0:\n");
		for(int i=0; i < SET_SIZE; i++)
			fprintf(seed_out, "%d ", SPT_new::GetNode(i));
		fprintf(seed_out, "\n");
		
		sigmaup = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
		Graph::set_all_edge_down();
		sigmadown = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
		alpha = sigmadown/sigmaup;
		fprintf(alpha_out, "iteration 0: ");
		fprintf(alpha_out, "%lg\n", alpha);
		fprintf(alpha_out, "%lg\n", Graph::CalSampleNum());

		for(int i=0; i < ITER_TIMES; i++)
		{
			
			for(int j=0; j < SET_SIZE; j++)
				seed_set[j] = SPT_new::GetNode(j);
		
			// intern update
			printf("Information Cascade Sample\n");
//			Graph::InformationCascadeSample(NUM_ITER, SET_SIZE, seed_set);
			Graph::iteration_IC(100*NUM_ITER, GeneralCascade::Run);

			// greedy
			Graph::set_all_edge_up();
			printf("Pmia start.\n");
			SPT_new::Build(SET_SIZE,bound);
			printf("Pmia finished.\n");
			
			// print seed
			fprintf(seed_out, "iteration %d:\n", i+1);
			for(int k=0; k < SET_SIZE; k++)
				fprintf(seed_out, "%d ", SPT_new::GetNode(k));
			fprintf(seed_out, "\n");
		
			// calculate alpha
			sigmaup = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
			Graph::set_all_edge_down();
			sigmadown = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
			alpha = sigmadown/sigmaup;			
			fprintf(alpha_out, "iteration %d: ", i+1);
			fprintf(alpha_out,"%lg\n", alpha);
			fprintf(alpha_out, "%lg\n", Graph::CalSampleNum());

			printf("%d complete!\n", i);
			
		}
		
		fclose(seed_out);
		fclose(alpha_out);
		
	}

	// for CLUCB algorithm
	s="-b";
	if(!s.compare(argv[1])){

		int bound=320;
	
		double sigmaup;
		double sigmadown;
		double alpha;

		int s1[SET_SIZE],s2[SET_SIZE]; //s1--bar{theta}, s2--adjusted theta
		int diff[2*SET_SIZE];
		int diff_num=0;

		
		FILE *seed_out;

		Graph::BuildUIC();
		GeneralCascade::Build();

		FILE *alpha_out;
		fopen_s(&alpha_out, "uncertain_influence.txt", "w");

		Graph::set_all_edge_up();
		Graph::print_graph("uncertain_graph0.txt");

		fopen_s(&seed_out, "seed_set.txt", "a");
		
		printf("Running Pmia!\n");
		SPT_new::Build(SET_SIZE,bound);
		printf("Pmia Complete!\n");
		
		sigmaup = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
		Graph::set_all_edge_down();
		sigmadown = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
		alpha = sigmadown/sigmaup;
		fprintf(alpha_out, "iteration 0:\n");
		fprintf(alpha_out, "%lg\n", alpha);
		fprintf(alpha_out, "%lg\n", Graph::CalSampleNum());

		printf("heihei\n");
		
		for(int i=0; i < ITER_TIMES; i++)
		{
			fprintf(seed_out, "iteration 0:\n");

			Graph::set_all_edge_up();
			printf("Running Pmia!\n");
			SPT_new::Build(SET_SIZE,bound);
			printf("Pmia Complete!\n");
			for(int k=0; k < SET_SIZE; k++)
			{
				s1[k]=SPT_new::GetNode(k);
				fprintf(seed_out, "%d ", s1[k]);
			}
			fprintf(seed_out, "\n");

			Graph::set_edge_CLUCB(SET_SIZE, s1);
			printf("Running Pmia!\n");
			SPT_new::Build(SET_SIZE,bound);
			printf("Pmia Complete!\n");
			for(int k=0; k < SET_SIZE; k++)
			{
				s2[k]=SPT_new::GetNode(k);
				fprintf(seed_out, "%d ", s2[k]);
			}
			fprintf(seed_out, "\n");



			//cal difference of s1 and s2
			diff_num=0;
			for (int k=0;k<SET_SIZE;k++)
			{
				bool flag=true;
				for (int j=0;j<SET_SIZE;j++)
					if (s1[k]==s2[j])
						flag=false;
				if (flag)
				{
					diff[diff_num]=s1[k];
					diff_num++;
				}
			}
			for (int k=0;k<SET_SIZE;k++)
			{
				bool flag=true;
				for (int j=0;j<SET_SIZE;j++)
					if (s2[k]==s1[j])
						flag=false;
				if (flag)
				{
					diff[diff_num]=s2[k];
					diff_num++;
				}
			}
			for(int k=0; k < diff_num; k++)
			{
				fprintf(seed_out, "%d ", diff[k]);
				printf("%d ", diff[k]);
			}
			printf("\n");

			Graph::GeneralSample(SAMPLE_PER_ITER, diff_num, diff);


			
			// greedy
			Graph::set_all_edge_up();
			printf("Running Pmia!\n");
			SPT_new::Build(SET_SIZE,bound);
			printf("Pmia Complete!\n");
			
			// print seed
			fprintf(seed_out, "iteration %d:\n", i+1);
			for(int k=0; k < SET_SIZE; k++)
				fprintf(seed_out, "%d ", SPT_new::GetNode(k));
			fprintf(seed_out, "\n");
		
			// calculate alpha
			sigmaup = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
			Graph::set_all_edge_down();
			sigmadown = toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
			alpha = sigmadown/sigmaup;			
			fprintf(alpha_out, "iteration %d:\n", i+1);
			fprintf(alpha_out,"%lg\n", alpha);
			fprintf(alpha_out, "%lg\n", Graph::CalSampleNum());

			printf("%d complete!\n", i);
				
		}
		
		fclose(seed_out);
		fclose(alpha_out);
		
	}


	//-----------------------------------end--------------------------------//

	// GreedyGC (improved by CELF)
	s="-g";
	if (!s.compare(argv[1])) {
		Graph::Build2WC();
		GeneralCascade::Build();	
		QueryPerformanceCounter(&start);
		Greedy::Build(SET_SIZE,GeneralCascade::Run);
		QueryPerformanceCounter(&Eend);
		timer = (double)(Eend.QuadPart - start.QuadPart) / freq.QuadPart;
		fopen_s(&timetmpfile, "time_greedy_gc.txt", "w");
		fprintf(timetmpfile,"%lg\n", timer);
		fclose(timetmpfile);
		system("copy greedy.txt greedy_gc.txt");
		toSimulate("GC_Greedy.txt", Greedy::GetNode, GeneralCascade::Run);
		system("del /Q tmp\\*");
	}

	// delete _running_.log to indicate finish
	system("del /Q _running_.log"); 
	system("echo. 2> _finished_.log"); 
	// cout << "_CrtDumpMemoryLeaks()" << _CrtDumpMemoryLeaks() << endl;
}

