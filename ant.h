#ifndef ants_h
#define ants_h

#include <iostream>
#include "parameters.h"

using namespace std;

class ants{
public:
	int variable;
	double ofn;
	double bestofn;
	int iter;
	ants();
	void init(int itno, double iofn);
	void copy_from(ants *from);
};

#endif