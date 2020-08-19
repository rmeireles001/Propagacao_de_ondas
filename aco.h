#ifndef aco_h
#define aco_h

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include "variable.h"
#include "values.h"
#include "parameters.h"
#include "ant.h"
#include "propagacao.h"


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876


using namespace std;
class lsaco;

struct trash{
	double *help_b;
	values tempval;
	trash();
	~trash();
	void reset(double *vec, int size);
};

class aco{
public:
	variable var;
	values itbestval, itworstval, bestval, *runbestval;
	ants *ant, bestant, itbestant, itworstant, *runbestant; 
	long int nevaluations, *nevalbest;
	trash t;
	long int seed;
	void get_data(double llimit, double ulimit, double increment);
	void ls_data(double llimit, double ulimit, double increment);
	void end();
	void ls_end();
	void initialize_ants_variables(int runno);
	void initialize_trail();
	void iteration_init(int iterationno, int nants);
	void find_values(int nants);
	int find_var_valueno();
	int find_no_of_max(int *noofmax);
	void find_prob();
	void decode_varvalues(double *varvalues, int antno);
	void analysis(int iterationno, int nants, int section, propagacao *p);
	//Objective function
	double objective_function(double *varvalue, int section, propagacao *p);
	void some_stats(int iterationno, int runno);
	void evaporate_trail();
	void update_trail_weighted(ants *antptr, int weight);
	void ras_update(double *help_b);
	void trail();
	void run(int section, propagacao *p);
	void ls_run(int section, propagacao *p, values *gvalptr, variable *varv);
	double get_var();
	#if localsearch
	lsaco *ls;
	#endif
};

class lsaco : public aco{
public:
	void init();
	void get_data(values *gvalptr, variable *varv);
	void run(values *gvalptr, variable *varv, int section, propagacao *p, long int *s);
	void some_stats(int iterationno);
	void evaporate_trail();
	void update_trail(ants *antptr);
	void trail();
	void end();
};

double drand(double low, double high);
double ran01( long *idum );

#endif
