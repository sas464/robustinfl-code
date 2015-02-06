#ifndef LIMIT_H
#define LIMIT_H

//#define MAX_NODE	660000
//#define	MAX_EDGE	4500000
#define MAX_NODE	66000
// #define	MAX_EDGE	450000 // not used
#define MAX_K		420

#define STR_LEN		200

#define SET_SIZE	50
#define EPS 1e-10

// for uncertain information influence
#define CHER_NUM 2.0
#define INIT_TIMES 100
#define SAMPLE_TIMES 100
#define ITER_TIMES 5
#define ITER_TIMES_2 3
#define THRESHOLD 50

// interval
#define INIT_DELTA_PARA 1;
#define UPDATE_DELTA_PARA 0.3;

//real sample
#define INIT_SAMPLE_NUMBER 200   //for every edge
#define SAMPLE_PER_ITER 500
#define NUM_ITER_1	50		 // and in calculating sigma function
#define NUM_ITER	10000    //number of cascades running each sub-iteration
	

#endif
