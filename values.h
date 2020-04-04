#ifndef valuesaco_h
#define valuesaco_h

#include <iostream>
#include "parameters.h"

using namespace std;

class values{
public:
	double var, ofn;
	int iter;
	values();
	void copy_from(values *from);
};

#endif