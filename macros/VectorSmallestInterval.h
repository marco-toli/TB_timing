
#ifndef FindSmallestInterval_h
#define FindSmallestInterval_h

#include <vector>
#include <iostream>


using namespace std;

double FindSmallestInterval(double mean, double meanErr, double min, double max,        
			    std::vector<double>* vals,                        
			    const double fraction, const bool verbosity);


#endif
