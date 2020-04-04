#ifndef variable_h
#define variable_h

#include <iostream>
#include <string>
#include <math.h>
#include <fstream>

using namespace std;

class variable{
public:
	double ulimit;
	double llimit;
	double increment;
	int nvalues;
	double *trail, *prob;
	void init(double ulimit, double llimit, double increment);
	void free();
	void find_nvalue();
	double *allocate(int size);
	void reset(double *vetor, int size);
	void ls_nvalue();
};

#endif