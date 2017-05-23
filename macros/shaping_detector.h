

#ifndef shaping_detector_h
#define shaping_detector_h

#include <vector>
#include <iostream>

using namespace std;


double shapeSiPM(double* x, double* par);
double shapeAPD(double* x, double* par);


double edepActivityDistr(double* x, double* par);
double generateScintiPulse(double* x, double* par);


#endif
