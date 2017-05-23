#include <vector>
#include <iostream>
#include "VectorSmallestInterval.h"
#include <algorithm>

using namespace std;

double FindSmallestInterval(double mean, double meanErr, double min, double max,                           std::vector<double>* vals,                          const double fraction, const bool verbosity)
{
  if( verbosity )
  std::cout << ">>>>>> FindSmallestInterval" << std::endl;


  std::sort(vals->begin(),vals->end());

  unsigned int nPoints = vals->size();
  unsigned int maxPoints = (unsigned int)(fraction * nPoints);

  unsigned int minPoint = 0;
  unsigned int maxPoint = 0;
  double delta = 999999.;
  for(unsigned int point = 0; point < nPoints-maxPoints; point++)
  {
    
//     std::cout << "maxPoints = " << maxPoints << " :: nPoints = " << nPoints <<  " :: nPoints -maxPoints = " << (nPoints-maxPoints) << std::endl;
    double tmpMin = vals->at(point);
    double tmpMax;
    if (point+maxPoints > 0) tmpMax = vals->at(point+maxPoints-1);
    else tmpMax = 1;
    
    if( tmpMax-tmpMin < delta )
    {
      delta = tmpMax - tmpMin;
      min = tmpMin;
      max = tmpMax;
      minPoint = point;
      maxPoint = point + maxPoints - 1;
    }
  }


/*
  TH1F* h_temp2 = new TH1F("h_temp2","",100,min-0.1*(max-min),max+0.1*(max-min));
  h_temp2 -> Sumw2();
  for(unsigned int point = minPoint; point <= maxPoint; ++point)
    h_temp2 -> Fill( vals.at(point) );

  mean = h_temp2 -> GetMean();
  meanErr = h_temp2 -> GetMeanError();

  delete h_temp2;*/
  
  return delta;
  
}