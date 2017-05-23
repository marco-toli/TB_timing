#include "shaping_detector.h"


#include <iostream>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TApplication.h"
#include "TFormula.h"
#include "TRandom.h"

#include <vector>

double shapeSiPM(double* x, double* par)
{
  
  double t = x[0];
  
  const double PI = 3.1415926;
  
  double A 		= par[0];	//cell signal amplitude
  double sigma_A	= par[1];	//cell amplitude jitter
  
  double input_t	= par[2];
  double tau_rise 	= par[3];
  double tau_decay 	= par[4];
  
  double b 		= tau_decay/tau_rise;
  
//   cout << " parameters are: [" << A << ", " << sigma_A << ", " << input_t << ", " << tau_rise << ", " << tau_decay << "]" << endl;

  //amplitude fluctuations
//   double A_bar = 1/(2*PI*sigma_A) * exp(-A_bar/(2*sigma_A));
  double A_bar = 1;
//   cout << " A_bar = " << A_bar << endl;
  
  //cell shaping
  double shaped_time = A_bar / (pow(b,(1/(1-b))) - pow(b,(1/(1/b-1)))) * (exp(-(t - input_t)/tau_decay) - exp(-(t-input_t)/tau_rise)) ;
  
  if (t < input_t) return 0;
  else 		   return shaped_time;
  
  
}

double shapeAPD(double* x, double* par)
{
  
  double t = x[0]/1000;
  
  const double PI = 3.1415926;
  
  double A 		= par[0];	//cell signal amplitude
  double sigma_A	= par[1];	//cell amplitude jitter
  
  double input_t	= par[2]/1000;
  double tau_rise 	= par[3];
  double tau_decay 	= par[4];
  
  /*
  TF1 *f = new TF1("f", "pow(exp(1.)*[1]*x/[0],[0])*exp(-[1]*x)",0.,20.);
  f->SetParameters(3.0, 0.8);
*/
  //amplitude fluctuations
//   double A_bar = 1/(2*PI*sigma_A) * exp(-A_bar/(2*sigma_A));
  double A_bar = 1;
//   cout << " A_bar = " << A_bar << endl;
  
  //apd shaping
  double shaped_time = A_bar * pow(exp(1)*tau_rise*(t-input_t)/tau_decay, tau_decay)* exp(-tau_rise*(t-input_t)) ;
  
//   cout << "exp(-tau_rise*(t-input_t)) = " << exp(-tau_rise*(t-input_t)) << " :: pow(exp(1)*tau_rise*(t-input_t)/tau_decay, tau_decay) = " << pow(exp(1)*tau_rise*(t-input_t)/tau_decay, tau_decay) << " :: shaped_time = " << shaped_time << endl;
  
  
  if (t < input_t) return 0;
  else 		   return shaped_time;
  
  
}


