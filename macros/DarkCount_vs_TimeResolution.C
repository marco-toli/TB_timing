// g++ -Wall -o DarkCount_vs_TimeResolution.exe  dict.cc setTDRStyle.cc ConfigParser.cc config_parser.cc ConfigFileLine.cc DarkCount_vs_TimeResolution.C  VectorSmallestInterval.cc shaping_detector.cc `root-config --cflags --glibs`

#include "VectorSmallestInterval.h"
#include "shaping_detector.h"
#include "dict.h"
#include "TLatex.h"
#include "TPaveText.h"

#include "setTDRStyle.h"
#include "ConfigParser.h"

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
#include <stdlib.h>
// using namespace std;

// double FindSmallestInterval(double mean, double meanErr, double min, double max, std::vector<double>* vals, const double fraction, const bool verbosity);


// void DarkCount_vs_TimeResolution()
// {
  
 
int main(int argc, char** argv)
{
  if( argc < 4 )
  {
    std::cerr << ">>>>> DarkCounts.C::usage:   " << argv[0] << " sipm_side_size [mm]    (maxEvents)  (PDE)  (input_file_name)  (output_file_name)" << std::endl;
    return -1;
  }
  
  float sipm_side_size;
  int length_cryst = 3;  

  if (argc > 1) 
  {
    sipm_side_size = atof(argv[1]);
    std::cout << " argv 1 = " << atof(argv[1]) << std::endl;
//     std::cout << "length crystal = " << length_cryst << std::endl;
  }
  
  float inputPDE = 0;
  
  int maxEvents = 1000;
  if (argc > 2) maxEvents = atoi(argv[2]);    
  if (argc > 3) inputPDE = atof(argv[3]);
  
  std::string inputFileName = "temp";
  if (argc > 4) inputFileName = argv[4];
  
  std::string outputFileName = "temp";
  if (argc > 5) outputFileName = argv[5];
  

//   
//   setTDRStyle();
//    TApplication* theApp = new TApplication("App", &argc, argv);
  
  TFile * RunFile;
//   RunFile = new TFile(Form("%s_%dmm_MERGED.root", inputFileName.c_str(), length_cryst),"READ");   //test_mu_3mm_noAbs
  RunFile = new TFile(Form("/afs/cern.ch/work/m/mlucchin/TB_timing/macros/ntuples/%s.root", inputFileName.c_str()),"READ");   //test the CMS tile conf
//   RunFile = new TFile("../data/CMS_tile_studies/out_cms_tile_3mm_12x12_wide_sipm5x5_grease_MERGED.root","READ");   //test the CMS tile conf
  
//   RunFile = new TFile(Form("../data/test_length/grease/merged/out_mu_grease_%dmm_MERGED_2.root", length_cryst),"READ");   //test_mu_3mm_noAbs

  
  const int NPH = 15;
  int nPhotons[NPH];
  nPhotons[0] = 1;
  nPhotons[1] = 2;
  nPhotons[2] = 3;
  nPhotons[3] = 4;
  nPhotons[4] = 5;
  nPhotons[5] = 7;
  nPhotons[6] = 10;
  nPhotons[7] = 15;
  nPhotons[8] = 20;
  nPhotons[9] = 30;
  nPhotons[10] = 40;
  nPhotons[11] = 50;
  nPhotons[12] = 75;
  nPhotons[13] = 100;
  nPhotons[14] = 200;
    
//   float thresh_amp = 0.003;
  float thresh_amp = 0.5;
  
  const int TIME_BINS = 550;
//   const int TIME_BINS = 1300;
//   const int TIME_BINS = 4202;
  int minTime = -200;
  int maxTime = 2000;
  

  
  const int nNoise = 8;
  float DCounts[nNoise];		//dark counts in 1 ns per mm^2 [ns-1 mm-2]
  DCounts[0] = 0;
  DCounts[1] = 1.;
  DCounts[2] = 2.5;     // GHz
  DCounts[3] = 5;
  DCounts[4] = 10;
  DCounts[5] = 25;
  DCounts[6] = 50;
  DCounts[7] = 100;
  

  
  
  double cher_PDE;
  double scint_PDE;
  
  if (inputPDE !=0) 
  {
      scint_PDE = inputPDE;
      cher_PDE  = inputPDE*0.8;
  }
  else              
  {
      scint_PDE = 0.25;
      cher_PDE = 0.21;
  }
//  cher_PDE = 0;
  
//   double read_out_eff = 0.36*0.7/0.3;	//6x6 mm² on 10x10 mm² + wrapping
  double read_out_eff = 1.;	//6x6 mm² on 10x10 mm² + wrapping
  scint_PDE = scint_PDE *read_out_eff;
  cher_PDE  = cher_PDE  *read_out_eff;
  
//   double read_out_eff = 1.;
//   scint_PDE = scint_PDE *0.7/0.3;
//   cher_PDE  = cher_PDE  *0.7/0.3;
  
  int SPTR = 66;
  float cell_size = 0.020;	//mm
  
  float cell_rise  = 200;  //ps
  float cell_decay = 10000; //ps
  float cell_amp   = 1;
  float cell_fluct = 0.1;
    
  int minTimeDC = -cell_decay*4;    //generating DC earlier as 2 times the recovery time of a cell, to obtain a proper baseline shift
  int maxTimeDC = maxTime;
  
  float beam_spot = 15;  //beam spot selection around tile center in mm
  
  
  TF1 * funcShapeSiPM = new TF1 ("funcShapeSiPM", shapeSiPM, minTimeDC, maxTime, 5);
  funcShapeSiPM->SetParameter(0, cell_amp);
  funcShapeSiPM->SetParameter(1, cell_fluct);
  funcShapeSiPM->SetParameter(3, cell_rise);	//in ps
  funcShapeSiPM->SetParameter(4, cell_decay);
  
  float sipm_section = sipm_side_size*sipm_side_size; //mm^2
  float tot_cells = 1./cell_size/cell_size; // per mm^2 
  
  //printing out configuration parameters
  std::cout << "************************** " << std::endl;
  std::cout << "*      CONFIGURATION     * " << std::endl;
  std::cout << "************************** " << std::endl;
  std::cout << "length     = " << length_cryst   << std::endl;  
  std::cout << "cher_PDE   = " << cher_PDE   << std::endl;
  std::cout << "scint_PDE  = " << scint_PDE  << std::endl;
  std::cout << "RO eff     = " << read_out_eff  << std::endl;
  std::cout << "SPTR       = " << SPTR       << std::endl;
  std::cout << "cell_size  = " << cell_size  << std::endl;
  std::cout << "cell_rise  = " << cell_rise  << std::endl;
  std::cout << "cell_decay = " << cell_decay << std::endl;
  std::cout << "cell_amp   = " << cell_amp   << std::endl;
  std::cout << "cell_fluct = " << cell_fluct << std::endl;
  std::cout << "tot cells  = " << tot_cells  << std::endl;
  std::cout << "beam spot  = " << beam_spot  << std::endl;
  std::cout << "************************** " << std::endl;
  
  
  float busy_frac[nNoise];
  for (int iNoise = 0; iNoise < nNoise; iNoise++) 
  {
    busy_frac[iNoise] =  DCounts[iNoise]*(5*cell_decay/1000) / tot_cells / sipm_section; //fraction of busy cells
    std::cout << "busy_frac [" << iNoise << "] = " << busy_frac[iNoise] << std::endl;
  }
  std::cout << "************************** " << std::endl;
  
  
  /// defining my histos and vectors
  TGraphErrors * gBaselineShift = new TGraphErrors ();
  TGraphErrors * gBaselineNoise = new TGraphErrors ();
  TGraphErrors * gPDE_drop = new TGraphErrors ();
  
  TGraphErrors * gE_dep_f = new TGraphErrors ();
  
  TGraphErrors * g_CTR_NPH [nNoise];
  TGraphErrors * g_CTR_NPH_corr [nNoise];    
  for (int iNoise = 0; iNoise < nNoise; iNoise++)
  {
	g_CTR_NPH[iNoise] = new TGraphErrors ();
	g_CTR_NPH_corr[iNoise] = new TGraphErrors ();	
  }
  
  std::vector<double> * ext_time_f[nNoise][NPH];
  std::vector<double> * ext_CTR[nNoise][NPH];
  std::vector<double> * ext_CTR_corr[nNoise][NPH];
  
  TH1F * hCTR[nNoise][NPH];
  TH1F * hCTR_corr[nNoise][NPH];
  
  //prompt photons counter
  TH1F * hCher_f[nNoise];
  TH1F * hScint_f[nNoise];
  //detected photons time distribution
  TH1F * hTime_cher_f[nNoise];
  TH1F * hTime_scint_f[nNoise];
  //photons survived after PDE deletion
  TH1F * hTime_cher_f_PDE[nNoise];
  TH1F * hTime_scint_f_PDE[nNoise];
  //total photons including cherenkov
  TH1F * hTime_scint_f_PDE_tot[nNoise];
  //adding thermal floor (DC+Crosstalk)
  TH1F * hTime_scint_f_PDE_tot_SPTR[nNoise];
  TH1F * hTime_scint_f_PDE_tot_SPTR_DRAW[nNoise];
  //smearing with SPTR
  TH1F * hTime_scint_f_PDE_tot_DC_SPTR[nNoise];
  //detector shaping
  TH1F * hTime_scint_f_PDE_tot_DC_SPTR_shaped[nNoise];
  
  TH1F * hTime_spike[nNoise];
  TH1F * hTime_spike_shaped[nNoise];
  
  TH1F * hDarkCount[nNoise];
  TH1F * hDarkCountRef[nNoise];
  TH1F * hBaseLineShift[nNoise];
  TH1F * hBaseLineShiftShaped[nNoise];
  
  TH2F * hCherCorr[nNoise];
  TH2F * hScatterTime[nNoise][NPH];
  TProfile * pTimeWalk[nNoise][NPH];
  TH2F * hTimeWalk[nNoise][NPH];
     
  TF1 *fitLin[nNoise][NPH];
     

  
  float maxAmp = 0.005 * length_cryst;
  
  for (int iNoise = 0; iNoise < nNoise; iNoise++) 
  { 
    for (int iNPH = 0; iNPH  < NPH; iNPH++)
    {
      hCTR[iNoise][iNPH] = new TH1F (Form("hCTR_%.1f_%dph", DCounts[iNoise], nPhotons[iNPH]), Form("hCTR_%.1f_%dph", DCounts[iNoise], nPhotons[iNPH]), 2000, 0, maxTime);  
      hCTR_corr[iNoise][iNPH] = new TH1F (Form("hCTR_corr_%.1f_%dph", DCounts[iNoise], nPhotons[iNPH]), Form("hCTR_corr_%.1f_%dph", DCounts[iNoise], nPhotons[iNPH]), 2000, 0, maxTime);  
      hScatterTime[iNoise][iNPH] = new TH2F (Form("hScatterTime_%.1f_%dph", DCounts[iNoise], nPhotons[iNPH]), Form("hCorr_%.1f_%dph", DCounts[iNoise], nPhotons[iNPH]), 2000, 0, maxTime, 2000, 0, maxTime);
      pTimeWalk[iNoise][iNPH] = new TProfile (Form("pTimeWalk_%.1f_%dph",  DCounts[iNoise], nPhotons[iNPH]), Form("pTimeWalk_%.1f_%dph",  DCounts[iNoise], nPhotons[iNPH]), 500, 0, maxAmp);
      hTimeWalk[iNoise][iNPH] = new TH2F (Form("hTimeWalk_%.1f_%dph",  DCounts[iNoise], nPhotons[iNPH]), Form("hTimeWalk_%.1f_%dph",  DCounts[iNoise], nPhotons[iNPH]), 500, 0, maxAmp, 2000, -2000, maxTime);
    }
    
    hCher_f[iNoise] = new TH1F (Form("hCher_f_%.1f", DCounts[iNoise]), Form("hCher_f_%.1f_", DCounts[iNoise]), 500, 0, 1000);
    hScint_f[iNoise] = new TH1F (Form("hScint_f_%.1f", DCounts[iNoise]), Form("hScint_f_%.1f_", DCounts[iNoise]), 500, 0, 10000);
        
    hTime_cher_f[iNoise] = new TH1F (Form("hTime_cher_f_%.1f", DCounts[iNoise]), Form("hTime_cher_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    hTime_scint_f[iNoise] = new TH1F (Form("hTime_scint_f_%.1f", DCounts[iNoise]), Form("hTime_scint_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    
    hTime_cher_f_PDE[iNoise] = new TH1F (Form("hTime_cher_f_PDE_%.1f", DCounts[iNoise]), Form("hTime_cher_PDE_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    hTime_scint_f_PDE[iNoise] = new TH1F (Form("hTime_scint_f_PDE_%.1f", DCounts[iNoise]), Form("hTime_scint_PDE_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    
    hTime_scint_f_PDE_tot[iNoise] = new TH1F (Form("hTime_scint_f_PDE_tot_%.1f", DCounts[iNoise]), Form("hTime_scint_PDE_tot_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    hTime_scint_f_PDE_tot_SPTR[iNoise] = new TH1F (Form("hTime_scint_f_PDE_tot_SPTR_%.1f", DCounts[iNoise]), Form("hTime_scint_PDE_tot_SPTR_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);    
    hTime_scint_f_PDE_tot_SPTR_DRAW[iNoise] = new TH1F (Form("hTime_scint_f_PDE_tot_SPTR_DRAW_%.1f", DCounts[iNoise]), Form("hTime_scint_PDE_tot_SPTR_DRAW_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    
    hTime_scint_f_PDE_tot_DC_SPTR[iNoise] = new TH1F (Form("hTime_scint_f_PDE_tot_DC_SPTR_%.1f", DCounts[iNoise]), Form("hTime_scint_PDE_tot_DC_SPTR_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise] = new TH1F (Form("hTime_scint_f_PDE_tot_DC_SPTR_shaped_%.1f", DCounts[iNoise]), Form("hTime_scint_PDE_tot_DC_SPTR_shaped_f_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    
    hTime_spike[iNoise] = new TH1F (Form("hTime_spike_%.1f", DCounts[iNoise]), Form("hTime_spike_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    hTime_spike_shaped[iNoise] = new TH1F (Form("hTime_spike_shaped_%.1f", DCounts[iNoise]), Form("hTime_spike_shaped_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    
    hDarkCount[iNoise] = new TH1F (Form("hDarkCount_%.1f", DCounts[iNoise]), Form("hDarkCount_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    hDarkCountRef[iNoise] = new TH1F (Form("hDarkCountRef_%.1f", DCounts[iNoise]), Form("hDarkCountRef_%.1f_", DCounts[iNoise]), TIME_BINS, minTime, maxTime);
    
    hBaseLineShift[iNoise] = new TH1F (Form("hBaseLineShift_%.1f", DCounts[iNoise]), Form("hBaseLineShift_%.1f_", DCounts[iNoise]), 1000, -20, 2000);
    hBaseLineShiftShaped[iNoise] = new TH1F (Form("hBaseLineShiftShaped_%.1f", DCounts[iNoise]), Form("hBaseLineShiftShaped_%.1f_", DCounts[iNoise]), 5000, -2, 10000);
    
    hCherCorr[iNoise] = new TH2F (Form("hCherCorr_%.1f", DCounts[iNoise]), Form("hCorr_%.1f_", DCounts[iNoise]), 500, 0, 500, 500, 0, 500);
  }
  
//   cout << " debug 0 " << endl;

  // get correction curves for TimeWalk   
//   TFile * storedCurves = new TFile(Form("./graphs_histos/mu_10mm_toyLSO_noise_lenght_%d.root", length_cryst), "READ");
//    TFile * storedCurves = new TFile(Form("./graphs_histos/12um_25PDE_10ns_decay_100RO_70LC/mu_10mm_toyLSO_noise_lenght_%d.root", length_cryst), "READ");
//   TFile * storedCurves = new TFile(Form("./graphs_histos/12um_25_PDE_10ns_decay_36RO_70LC/mu_10mm_toyLSO_noise_lenght_%d.root", length_cryst), "READ");
//   TFile * storedCurves = new TFile(Form("./graphs_histos/GLUE_12um_25PDE_10ns_decay_36RO_70LC/mu_10mm_toyLSO_noise_lenght_%d.root", length_cryst), "READ");
//   TFile * storedCurves = new TFile(Form("%s_%d.root", outputFileName.c_str(), length_cryst), "READ");
  TFile * storedCurves = new TFile(Form("/afs/cern.ch/work/m/mlucchin/TB_timing/macros/graphs_histos/PDE_scan/%s_PDE_%.0f.root", inputFileName.c_str(), scint_PDE*100), "READ");
  
  
  TF1 *fitTimeWalk[nNoise][NPH];
  for (int iNoise = 0; iNoise < nNoise; iNoise++)
  {
    for (int iNPH = 0; iNPH < NPH; iNPH++)    fitTimeWalk[iNoise][iNPH] = (TF1*) storedCurves->Get(Form("fitLin_%d_NPH_%d", iNoise, iNPH));
  }
  
  /// defining tree variablesù
  std::vector<float>* inputInitialPosition =  new vector<float>(3,0.); 
  Float_t depositedEnergyTotal;
  Float_t depositedEnergyCore_f;
  std::vector<float> * time_prod_scint = new std::vector<float>;
  std::vector<float> * time_prod_cher = new std::vector<float>;
  
  std::vector<float> * time_ext_scint = new std::vector<float>;
  std::vector<float> * time_ext_cher = new std::vector<float>;
  
  std::vector<float> * time_det_scint = new std::vector<float>;
  std::vector<float> * time_det_cher = new std::vector<float>;
  
  TTree* TreeRun = (TTree*) RunFile->Get("tree");
  TreeRun->SetBranchAddress("inputInitialPosition",&inputInitialPosition);
  TreeRun->SetBranchAddress("depositedEnergyCore_f",&depositedEnergyCore_f);
  TreeRun->SetBranchAddress("time_ext_scint",&time_ext_scint);
  TreeRun->SetBranchAddress("time_ext_cher",&time_ext_cher);
//  TreeRun->SetBranchAddress("time_det_scint",&time_det_scint);
//  TreeRun->SetBranchAddress("time_det_cher",&time_det_cher);
  
  TGraphErrors * gCTR_vs_Noise = new TGraphErrors();
  TGraphErrors * gCTR_vs_Noise_corr = new TGraphErrors();
    
  double min_ctr_length_noise[nNoise];
  double corr_min_ctr_length_noise[nNoise];
    
  for (int iNoise = 0; iNoise < nNoise; iNoise++)
  {
    
    for (int iNPH = 0; iNPH  < NPH; iNPH++)
    {
      ext_time_f[iNoise][iNPH] = new std::vector<double>;
      ext_CTR[iNoise][iNPH] = new std::vector<double>;
      ext_CTR_corr[iNoise][iNPH] = new std::vector<double>;
    }
 
    ///******************************************///
    ///		   Run over events	         ///
    ///******************************************///
    
    int NEVENTS = TreeRun->GetEntries();
    if (maxEvents < NEVENTS) NEVENTS = maxEvents;
    cout << "DCounts[" << iNoise << "] = " << DCounts[iNoise] << " :: nEvents = " << NEVENTS << endl;
    
    int count_weird = 0;
    
    for (Int_t iEvt= 0; iEvt < NEVENTS; iEvt++) 
    {

	TreeRun->GetEntry(iEvt);
//  	cout << " iEvt = " << iEvt << " :: energy_Dep = " << depositedEnergyCore_f << endl;	
 	if( iEvt%1 == 0 ) std::cout << ">>> reading entry " << iEvt << " / " << NEVENTS << "\r" << std::flush;
        
        float x = inputInitialPosition->at(0);
        float y = inputInitialPosition->at(1);
//         if ( x>6 || x<-6 || y > 6 || y<-6) continue;    // if particle is not contained
        if ( sqrt(pow(x,2) + pow(y,2)) > beam_spot) continue;    // if particle is not contained

	//photon counting
	int nChPh_f = 0;
	int nScintPh_f = 0;
	
// 	if (time_ext_scint->size() > 300000) 
	if (time_ext_scint->size() > 30000/5*length_cryst) 
	{
	  std::cout << std::endl;
	  count_weird++;
	  std::cout << "too many photons: pion conversion ?! --> scint ph[" << iEvt << "] = " << time_ext_scint->size() << " :: weird events = " << count_weird << " :: frac = " << (float) count_weird/iEvt <<  std::endl;
	  
	  continue;
	}
	
	//counting prompt photons
	int t_prompt = maxTime;
//  	std::cout << " scint ph[" << iEvt << "] = " << time_ext_scint->size() << " extracted = " << (int) time_ext_scint->size()*read_out_eff << std::endl;
	for (int iPhExt = 0; iPhExt < (int) time_ext_scint->size()*read_out_eff; iPhExt++)
	{
	  hTime_scint_f[iNoise]->Fill(time_ext_scint->at(iPhExt));
	  if (time_ext_scint->at(iPhExt) < t_prompt) nScintPh_f++;
	}	
	for (int iPhExt = 0; iPhExt < (int) time_ext_cher->size()*read_out_eff; iPhExt++)
	{
	  hTime_cher_f[iNoise]->Fill(time_ext_cher->at(iPhExt));
	  if (time_ext_cher->at(iPhExt) < t_prompt) nChPh_f++;
	}
	hCher_f[iNoise]->Fill(nChPh_f);
	hScint_f[iNoise]->Fill(nScintPh_f);
	
	
	 
	//0) cleaning up temp histos (filled and cleaned event-by-event)
	for (int iBin = 0; iBin < TIME_BINS; iBin++) 
	{	 
	  hDarkCount[iNoise]->SetBinContent(iBin+1, 0);	  
// 	  hTime_scint_f[iNoise]->SetBinContent(iBin+1, 0);	  
// 	  hTime_scint_f_PDE[iNoise]->SetBinContent(iBin+1, 0);	  
// 	  hTime_scint_f_PDE_tot[iNoise]->SetBinContent(iBin+1, 0);	  
	  hTime_scint_f_PDE_tot_SPTR[iNoise]->SetBinContent(iBin+1, 0);	  
	  hTime_scint_f_PDE_tot_DC_SPTR[iNoise]->SetBinContent(iBin+1, 0);	  
	  hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->SetBinContent(iBin+1, 0);	  
	  hTime_spike_shaped[iNoise]->SetBinContent(iBin+1, 0);
	}
	
	
	
	//1) PDE - random deletion according to PDE

	int IN_SIZE = time_ext_scint->size();
        float eff_PDE = scint_PDE * (1-busy_frac[iNoise]);
        gPDE_drop->SetPoint(iNoise, DCounts[iNoise], eff_PDE);
        
	int del_SIZE = (int) IN_SIZE*(1-eff_PDE);
	
	for (int iPhExt = 0; iPhExt < del_SIZE; iPhExt++)               time_ext_scint->erase(time_ext_scint->begin() +  (rand() % time_ext_scint->size()));
	
	for (int iPhExt = 0; iPhExt < time_ext_scint->size(); iPhExt++)	hTime_scint_f_PDE[iNoise]->Fill(time_ext_scint->at(iPhExt));
	
	
	
	
	//2) CHERENKOV - adding Cherenkov in good proportion wrt to PDE		
// 	cout << "cher phot = " << nChPh_f << endl;
	for (int iPhCh = 0; iPhCh <  time_ext_cher->size()*cher_PDE * (1-busy_frac[iNoise]); iPhCh++)
	{
	  time_ext_scint->push_back(time_ext_cher->at(iPhCh));	  
	}
	for (int iPhExt = 0; iPhExt < time_ext_scint->size(); iPhExt++)		hTime_scint_f_PDE_tot[iNoise]->Fill(time_ext_scint->at(iPhExt));

	
	
	
	// 3) SPTR - sostituisco ogni valore nel vettore con il suo valore smeareato per SPTR - 
	for (int iPhExt = 0; iPhExt < time_ext_scint->size(); iPhExt++)	  
	{
	  time_ext_scint->at(iPhExt) = gRandom->Gaus(time_ext_scint->at(iPhExt), SPTR);
	  hTime_scint_f_PDE_tot_SPTR[iNoise]->Fill(time_ext_scint->at(iPhExt));		
	  hTime_scint_f_PDE_tot_SPTR_DRAW[iNoise]->Fill(time_ext_scint->at(iPhExt));		
	}

// 	std::cout << "SPTR added ..." << std::endl;
//  	for (int iPhExt = 0; iPhExt < time_ext_scint->size(); iPhExt++)		hTime_scint_f_PDE_tot_DC_SPTR[iNoise]->Fill(time_ext_scint->at(iPhExt));
	
	
	// 5) detector shaping
	for (int iBin = 0; iBin< TIME_BINS; iBin++) 
	{		
  // 	cout 	 << " iBin = " << iBin << " :: time_ph = " << opPhoton_time_det->at(iPh) << " :: time_bin = " 
  // 		 << hTotDetAPD_shaped[iNoise]->GetBinCenter(iBin+1) << " :: shaping = " << funcShapeSiPM->Eval(hTotDetAPD_shaped[iNoise]->GetBinCenter(iBin+1)) << endl;
	  funcShapeSiPM->SetParameter(2, hTime_scint_f_PDE_tot_SPTR[iNoise]->GetBinCenter(iBin+1));
// 	  double time_raw = hTime_scint_f_PDE_tot_SPTR[iNoise]->GetBinContent(iBin+1);	
	  for (int jBin = 0; jBin< TIME_BINS; jBin++) 
	  {
	    double time_bin = hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->GetBinCenter(jBin+1);	    
	    hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->Fill(time_bin, funcShapeSiPM->Eval(time_bin)  * hTime_scint_f_PDE_tot_SPTR[iNoise]->GetBinContent(iBin+1)); //* 1/TIME_BINS);
	  }
	}
	
	
	//adding shaped dark counts ->  to save computing time
	//generate DC on longer gate but calculate their shaping only in interest region close to leading edge
	
	for (int iDC = 0; iDC < DCounts[iNoise]*(maxTimeDC-minTimeDC)/1000; iDC++)
	{
	  float t_DC = gRandom->Uniform(minTimeDC, maxTimeDC);
	  hDarkCount[iNoise]->Fill(t_DC);
// 	  time_ext_scint->push_back(t_DC);	 
  	  funcShapeSiPM->SetParameter(2, t_DC);
	  
	  for (int jBin = 0; jBin< TIME_BINS; jBin++) 
	  {
	    double time_bin = hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->GetBinCenter(jBin+1);	    
// 	    std::cout << " time_bin
	    hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->Fill(time_bin, funcShapeSiPM->Eval(time_bin) );//* 1/TIME_BINS);
	  }
	  
	}
	
	//baseline subtraction
	
	double ped_f = 0;
	int ped_bin  = (int) TIME_BINS*0.08 ;
// 	std::cout << "calculating event-by-event baseline subtraction based on first " << ped_bin << " bins :: i.e. a time interval of average = " << ped_bin*(maxTime-minTime)/TIME_BINS << std::endl;
	for (int iBin = 0; iBin< ped_bin; iBin++)
	{
	  ped_f += hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->GetBinContent (iBin+1);
	}
	ped_f /= ped_bin;
	hBaseLineShiftShaped[iNoise]->Fill(ped_f);
	
	//WARNING check if value is correct given the SiPM shaping time
	//	  use an average baseline subtraction instead of event by event
// 	float mean_baseline = 0.0354461*DCounts[iNoise];
	float mean_baseline = 0.266458 + 10.5757*DCounts[iNoise];
//  	ped_f = mean_baseline;
	
 	for (int iBin = 0; iBin< TIME_BINS; iBin++) hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->SetBinContent (iBin+1, hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->GetBinContent(iBin+1) - ped_f);
// 	std::cout << std::endl
// 	std::cout << " ped_f[" << iNoise << "] = " << ped_f << std::endl;

	
	//spike shaping for reference
	funcShapeSiPM->SetParameter(2, 0);
//	  double time_raw = hTime_scint_f_PDE_tot_DC_SPTR[iNoise]->GetBinContent(iBin+1);	
	for (int jBin = 0; jBin< TIME_BINS; jBin++) 
	{
	    double time_bin = hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->GetBinCenter(jBin+1);	    
	    hTime_spike_shaped[iNoise]->Fill(time_bin, funcShapeSiPM->Eval(time_bin) );//* (maxTime-minTime)/TIME_BINS);
	}
	
	
	// sort photons by time of arrival
//  	cout << "before sorting time vectors ... " << endl;
// 	sort(time_ext_scint->begin(), time_ext_scint->end());
// 	std::cout << "Photons sorted by time of arrival ..." << std::endl;
	
	
	for (int iNPH = 0; iNPH < NPH; iNPH++)
	{
 	  
	  double scint_mean_time_f = 0;
	  
	  for (int iBin = 0; iBin< TIME_BINS; iBin++) 
	  {
	    if (hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->GetBinContent(iBin+1) > thresh_amp*nPhotons[iNPH]) 
	    {
	      scint_mean_time_f = hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->GetBinCenter(iBin+1);
	      break;
	    }
	  }


	  /*
	  //use average time of arrival as estimator
	  for (int iPrSc = 0; iPrSc < nPhotons[iNPH]; iPrSc++)
	  {
	    if (iPrSc<time_ext_scint->size()) scint_mean_time_f += time_ext_scint->at(iPrSc);
	    else	
	    {
	      cout << "not enough front photons: max @ " << iPrSc << " ph " << " :: evt = " << iEvt << endl;
	      scint_mean_time_f = 0;
	      break;
	    }	    
//  	    cout << " here iNPH = " << iNPH << ":: iPrsc = " << iPrSc << endl;
	  }
	  
	  ext_time_f[iNoise][iNPH]->push_back(scint_mean_time_f/nPhotons[iNPH]);
	  */
	  
	  double CT = 0, CT_corr = 0;
	  if (scint_mean_time_f!=0)
	  {
	    CT = scint_mean_time_f;
	    if (fitTimeWalk[iNoise][iNPH]!= NULL) CT_corr = CT - fitTimeWalk[iNoise][iNPH]->Eval(depositedEnergyCore_f);
	  }
	  
//     	  double CT_corr = CT;
// 	  
//  	  if (iNoise == 7)
// 	  cout << "iNPH: " << iNPH << " :: time_f = " << scint_mean_time_f/nPhotons[iNPH] << " :: CT = " << CT << endl;
	    
	  if ( CT != 0   && depositedEnergyCore_f>0.0007 ) 
	  {
	    ext_CTR[iNoise][iNPH]->push_back(CT);
	    hCTR[iNoise][iNPH]->Fill(CT);
	    
	    ext_CTR_corr[iNoise][iNPH]->push_back(CT_corr);
	    hCTR_corr[iNoise][iNPH]->Fill(CT_corr);
	    
	    pTimeWalk[iNoise][iNPH]->Fill(depositedEnergyCore_f, CT);
	    hTimeWalk[iNoise][iNPH]->Fill(depositedEnergyCore_f, CT);
	  }
	}	//end of NPHE loop	
      }	//end of events loop         
      std::cout << std::endl;
      
//        std::cout << "end of events loop ..." << std::endl;
      
      for (int iBin = 0; iBin < TIME_BINS; iBin++) 
      {
	hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->SetBinError(iBin+1, 0);
	hTime_spike_shaped[iNoise]->SetBinError(iBin+1, 0);
      }
      
//       std::cout << "before CTR calculation ..." << std::endl;
      
      float min_ctr 	   = 9999;
      float min_ctr_corr = 9999;
        
      for (int iNPH = 0; iNPH < NPH; iNPH++)
      {
	double CL_ext_time_f 	= FindSmallestInterval((double) 200, 10, 0, 500, ext_time_f[iNoise][iNPH], 0.68, 0)/2.;
	double CL_ext_CTR    	= FindSmallestInterval((double) 200, 10, 0, 500, ext_CTR[iNoise][iNPH], 0.68, 0)/2.;	     
	double CL_ext_CTR_corr  = FindSmallestInterval((double) 200, 10, 0, 500, ext_CTR_corr[iNoise][iNPH], 0.68, 0)/2.;
	
	std::cout << "CL_ext_CTR = " << CL_ext_CTR << std::endl;
	
	g_CTR_NPH[iNoise] ->SetPoint(iNPH, nPhotons[iNPH], CL_ext_CTR);	
	g_CTR_NPH_corr[iNoise] ->SetPoint(iNPH, nPhotons[iNPH], CL_ext_CTR_corr);
        
        if (CL_ext_CTR < min_ctr)           min_ctr = CL_ext_CTR;
	if (CL_ext_CTR_corr < min_ctr_corr) min_ctr_corr = CL_ext_CTR_corr;	
	
	fitLin[iNoise][iNPH] = new TF1 (Form("fitLin_%d_NPH_%d",iNoise, iNPH), "[0] + [1]*x + [2]*x*x", 0., 2);
	fitLin[iNoise][iNPH]->SetLineColor(kRed);
	for (int i = 0; i<5; i++) pTimeWalk[iNoise][iNPH]->Fit(fitLin[iNoise][iNPH], "WQR");
      }
      
      min_ctr_length_noise[iNoise] = min_ctr;
      corr_min_ctr_length_noise[iNoise] = min_ctr_corr;
      
      
      // 	if (iNoise 
      gCTR_vs_Noise->SetPoint(iNoise, DCounts[iNoise]*1e9, min_ctr_length_noise[iNoise]);	
      gCTR_vs_Noise_corr->SetPoint(iNoise, DCounts[iNoise]*1e9, corr_min_ctr_length_noise[iNoise]);
      
    }  //end of iNoise loop
    std::cout << std::endl;
    
//     std::cout << "end of noise loop ..." << std::endl;

   //prompt photons
   cout << "prompt scint = " 		<< hScint_f[0]->GetMean() << endl;
   cout << "prompt cher = " 		<< hCher_f[0]->GetMean() << endl;

    
    
   ///******************** DRAWING **********************///
   TLegend *leg;
   int selNPH = 4;

   TPaveText *pt3 = new TPaveText(0.15, 0.84, 0.4, 0.88, "brNDC");
   pt3->SetShadowColor(0);
   pt3->SetFillColor(0);
   pt3->SetLineColor(0);   
//    pt3->AddText("Cherenkov + Scintillation");
   pt3->AddText(Form("Crystal length: %d mm", length_cryst));
   pt3->SetTextFont(42);
   pt3->SetTextAlign(12);
   pt3->SetTextSize(0.03);
   
   TPaveText *pt2 = new TPaveText(0.15, 0.78, 0.4, 0.83, "brNDC");
   pt2->SetShadowColor(0);
   pt2->SetFillColor(0);
   pt2->SetLineColor(0);
   pt2->AddText(Form("PDE = %.2f ", scint_PDE));
   pt2->SetTextFont(42);
   pt2->SetTextAlign(12);
   pt2->SetTextSize(0.03);
      
   TPaveText *pt = new TPaveText(0.15, 0.72, 0.5, 0.77, "brNDC");
   pt->SetShadowColor(0);
   pt->SetFillColor(0);
   pt->SetLineColor(0);
   pt->AddText(Form("#sigma_{SPTR} = %d ps", SPTR));
   pt->SetTextFont(42);
   pt->SetTextAlign(12);
   pt->SetTextSize(0.03);
   
   TLatex t_prel(0.59,.91,"SIMULATION: 10 mm LSO:Ce crystal"); 
   t_prel.SetTextFont(43);
   t_prel.SetTextSize(19);
   t_prel.SetNDC(kTRUE);
//    pt->Draw();


//    pt3->Draw();
 
    //scale histos
   
    for (int iNoise = 0; iNoise < nNoise; iNoise++)
    {
       hTime_scint_f[iNoise]->Scale(1./maxEvents);
       hTime_scint_f_PDE[iNoise]->Scale(1./maxEvents);
       hTime_scint_f_PDE_tot[iNoise]->Scale(1./maxEvents);
       hTime_scint_f_PDE_tot_SPTR[iNoise]->Scale(1./maxEvents);
       hTime_scint_f_PDE_tot_SPTR_DRAW[iNoise]->Scale(1./maxEvents);
       hTime_scint_f_PDE_tot_DC_SPTR[iNoise]->Scale(1./maxEvents);
       hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->Scale(1./maxEvents);
       hDarkCount[iNoise]->Scale(1./maxEvents);       
    }

    
   ///DRAW Simulation Workflow
   

   
   //step 1
   TCanvas * cScintillationTime = new TCanvas ("cScintillationTime", "cScintillationTime", 450, 450);
   hTime_scint_f[0]->SetLineColor(kBlack);
   hTime_scint_f[0]->SetStats(0);
   hTime_scint_f[0]->SetTitle("Scintillation photons at SiPM");
   hTime_scint_f[0]->GetXaxis()->SetTitle("Time [ps]");
   if (maxEvents > 1)   hTime_scint_f[0]->GetYaxis()->SetTitle("Average signal [photons / 4 ps]");
   else 	  	hTime_scint_f[0]->GetYaxis()->SetTitle("Signal [photons / 4 ps]");
   hTime_scint_f[0]->GetYaxis()->SetTitleOffset(1.3);
   hTime_scint_f[0]->GetXaxis()->SetRangeUser(-200, 800);
   hTime_scint_f[0]->Draw();
   cout << "mean arrival time for scint f = " << hTime_scint_f[0]->GetMean() << endl;
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
//      hTime_scint_f[iNoise]->SetLineColor(iNoise+1);
//      hTime_scint_f[iNoise]->Draw("same");
//      cout << "mean ch f = " << hTime_scint_f[iNoise]->GetMean() << endl;
   }
   gPad->SetGrid();
   
   leg = new TLegend(0.6,0.7,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   for (int iNoise = 0; iNoise < nNoise; iNoise++)      leg->AddEntry(hTime_scint_f[iNoise], Form("DCR = %.1e", DCounts[iNoise]*1e9), "lp");
   
//    leg->Draw();
   
   //step 2
   TCanvas * cScintillationTimePDE = new TCanvas ("cScintillationTimePDE", "cScintillationTimePDE", 450, 450);
   hTime_scint_f_PDE[0]->SetLineColor(kBlack);
   hTime_scint_f_PDE[0]->SetStats(0);
   hTime_scint_f_PDE[0]->SetTitle("Scintillation photons detected (after PDE)");
   hTime_scint_f_PDE[0]->GetXaxis()->SetTitle("Time [ps]");
   if (maxEvents > 1)   hTime_scint_f_PDE[0]->GetYaxis()->SetTitle("Average signal [photons / 4 ps]");
   else 	  	hTime_scint_f_PDE[0]->GetYaxis()->SetTitle("Signal [photons / 4 ps]");
   hTime_scint_f_PDE[0]->GetYaxis()->SetTitleOffset(1.3);
   hTime_scint_f_PDE[0]->GetXaxis()->SetRangeUser(-200, 800);
   hTime_scint_f_PDE[0]->Draw();
   cout << "mean arrival time for scint f PDE = " << hTime_scint_f_PDE[0]->GetMean() << endl;
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     hTime_scint_f_PDE[iNoise]->SetLineColor(iNoise+1);
//      hTime_scint_f_PDE[iNoise]->Draw("same");
   }
   gPad->SetGrid();
   //step 3
   TCanvas * cScintillationTimePDE_tot = new TCanvas ("cScintillationTimePDE_tot", "cScintillationTimePDE_tot", 450, 450);
   hTime_scint_f_PDE_tot[0]->SetLineColor(kBlack);
   hTime_scint_f_PDE_tot[0]->SetStats(0);
   hTime_scint_f_PDE_tot[0]->SetTitle("Scintillation + Cherenkov photons detected (after PDE)");
   hTime_scint_f_PDE_tot[0]->GetXaxis()->SetTitle("Time [ps]");
   if (maxEvents > 1)   hTime_scint_f_PDE_tot[0]->GetYaxis()->SetTitle("Average signal [photons / 4 ps]");
   else 	  	hTime_scint_f_PDE_tot[0]->GetYaxis()->SetTitle("Signal [photons / 4 ps]");
   hTime_scint_f_PDE_tot[0]->GetYaxis()->SetTitleOffset(1.3);
   hTime_scint_f_PDE_tot[0]->GetXaxis()->SetRangeUser(-200, 800);
   hTime_scint_f_PDE_tot[0]->Draw();
   
   cout << "mean arrival time for scint f PDE tot = " << hTime_scint_f_PDE_tot[0]->GetMean() << endl;
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     hTime_scint_f_PDE_tot[iNoise]->SetLineColor(iNoise+1);
//      hTime_scint_f_PDE_tot[iNoise]->Draw("same");
   }
   gPad->SetGrid();
   //step 4
   TCanvas * cScintillationTimePDE_tot_SPTR = new TCanvas ("cScintillationTimePDE_tot_SPTR", "cScintillationTimePDE_tot_SPTR", 450, 450);
   hTime_scint_f_PDE_tot_SPTR[0]->SetLineColor(kBlack);
   hTime_scint_f_PDE_tot_SPTR[0]->SetStats(0);
   hTime_scint_f_PDE_tot_SPTR[0]->SetTitle("Arrival time including Photodetector Time Jitter (SPTR)");
   hTime_scint_f_PDE_tot_SPTR[0]->GetXaxis()->SetTitle("Time [ps]");
   if (maxEvents > 1)   hTime_scint_f_PDE_tot_SPTR[0]->GetYaxis()->SetTitle("Average signal [photons / 4 ps]");
   else 	  	hTime_scint_f_PDE_tot_SPTR[0]->GetYaxis()->SetTitle("Signal [photons / 4 ps]");
   hTime_scint_f_PDE_tot_SPTR[0]->GetYaxis()->SetTitleOffset(1.3);
   hTime_scint_f_PDE_tot_SPTR[0]->GetXaxis()->SetRangeUser(-200, 800);
   hTime_scint_f_PDE_tot_SPTR[0]->Draw();
  
   cout << "mean arrival time for scint f PDE tot = " << hTime_scint_f_PDE_tot_SPTR[0]->GetMean() << endl;
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     hTime_scint_f_PDE_tot_SPTR[iNoise]->SetLineColor(iNoise+1);
//      hTime_scint_f_PDE_tot_SPTR[iNoise]->Draw("same");
   }
   gPad->SetGrid();
   /*
   //step 5
   TCanvas * cScintillationTimePDE_tot_DC_SPTR = new TCanvas ("cScintillationTimePDE_tot_DC_SPTR", "cScintillationTimePDE_tot_DC_SPTR", 450, 450);
   hTime_scint_f_PDE_tot_DC_SPTR[0]->SetLineColor(kBlack);
   hTime_scint_f_PDE_tot_DC_SPTR[0]->SetStats(0);
   hTime_scint_f_PDE_tot_DC_SPTR[0]->SetTitle("Time of pixels firing (including dark counts)");
   hTime_scint_f_PDE_tot_DC_SPTR[0]->GetXaxis()->SetTitle("Time [ps]");
   if (maxEvents > 1)   hTime_scint_f_PDE_tot_DC_SPTR[0]->GetYaxis()->SetTitle("Average signal [photons / 4 ps]");
   else 	  	hTime_scint_f_PDE_tot_DC_SPTR[0]->GetYaxis()->SetTitle("Signal [photons / 4 ps]");
   hTime_scint_f_PDE_tot_DC_SPTR[0]->GetYaxis()->SetTitleOffset(1.3);
   hTime_scint_f_PDE_tot_DC_SPTR[0]->Draw();
   cout << "mean arrival time for scint f PDE tot = " << hTime_scint_f_PDE_tot_DC_SPTR[0]->GetMean() << endl;
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     hTime_scint_f_PDE_tot_DC_SPTR[iNoise]->SetLineColor(iNoise+1);
     hTime_scint_f_PDE_tot_DC_SPTR[iNoise]->Draw("same");
   }
   gPad->SetGrid();
   */
   //step 6
   TCanvas * cScintillationTimePDE_tot_DC_SPTR_shaped = new TCanvas ("cScintillationTimePDE_tot_DC_SPTR_shaped", "cScintillationTimePDE_tot_DC_SPTR_shaped", 450, 450);
   hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->SetLineColor(kBlack);
   hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->SetStats(0);
   hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->SetTitle("SiPM signal (convolution with cell recovery time)");
   hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->GetXaxis()->SetTitle("Time [ps]");
   hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->GetYaxis()->SetTitle("Signal [a.u.]");
   hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->GetYaxis()->SetTitleOffset(1.3);
   hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->GetXaxis()->SetRangeUser(-200, 800);
   hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->Draw();
   
   cout << "mean arrival time for scint f PDE tot_shaped = " << hTime_scint_f_PDE_tot_DC_SPTR_shaped[0]->GetMean() << endl;
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     
     hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->SetLineColor(iNoise+1);
     hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise]->Draw("same");
   }
   TLine * l_thresh[NPH];
   for (int iNPH = 0; iNPH < NPH; iNPH++) 
   {
     l_thresh[iNPH] = new TLine(-200, thresh_amp*nPhotons[iNPH], 800, thresh_amp*nPhotons[iNPH]);
     l_thresh[iNPH]->SetLineColor(kBlack);
     l_thresh[iNPH]->SetLineStyle(7);
     l_thresh[iNPH]->SetLineWidth(2);
     l_thresh[iNPH]->Draw("same");
   }
   gPad->SetGrid();
      
   leg = new TLegend(0.6,0.7,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(hTime_scint_f_PDE_tot_DC_SPTR_shaped[iNoise], Form("DCR = %.1e", DCounts[iNoise]*1e9), "lp");
   }
   leg->Draw();
   
   //spike
   TCanvas * cScintillationSpike = new TCanvas ("cScintillationSpike", "cScintillationSpike", 450, 450);
   hTime_spike_shaped[0]->SetLineColor(kBlack);
   hTime_spike_shaped[0]->Draw();
//    cout << "mean arrival time for scint f PDE = " << hTime_spike_shaped[0]->GetMean() << endl;
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     hTime_spike_shaped[iNoise]->SetLineColor(iNoise+1);
     hTime_spike_shaped[iNoise]->Draw("same");
   }
   gPad->SetGrid();
   //baseline fluctuations
   /* 
   TCanvas * cBaseLineFluctuations = new TCanvas ("cBaseLineFluctuations", "cBaseLineFluctuations", 450, 450);
   hBaseLineShift[0]->SetLineColor(kBlack);
   hBaseLineShift[0]->Draw();   
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     cout << "rms baseline shift [" << iNoise << "] = " << hBaseLineShift[iNoise]->GetRMS() << endl;
     hBaseLineShift[iNoise]->SetLineColor(iNoise+1);
     hBaseLineShift[iNoise]->Draw("same");
   }
   gPad->SetGrid();
   */
   // shaped
   
   
   TCanvas * cBaseLineFluctuationsShaped = new TCanvas ("cBaseLineFluctuationsShaped", "cBaseLineFluctuationsShaped", 450, 450);
   hBaseLineShiftShaped[0]->SetLineColor(kBlack);
   hBaseLineShiftShaped[0]->Draw();
   hBaseLineShiftShaped[0]->GetXaxis()->SetTitle("Mean baseline [a.u.]");
   hBaseLineShiftShaped[0]->GetYaxis()->SetTitle("Counts");
   
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     cout << "*** CHECK *** rms shaped baseline shift[" << iNoise << "] = " << hBaseLineShiftShaped[iNoise]->GetMean() <<  " :: RMS = " << hBaseLineShiftShaped[iNoise]->GetMean() << endl;
     hBaseLineShiftShaped[iNoise]->SetLineColor(iNoise+1);
     hBaseLineShiftShaped[iNoise]->Draw("same");
     gBaselineShift->SetPoint(iNoise, DCounts[iNoise], hBaseLineShiftShaped[iNoise]->GetMean());
     gBaselineShift->SetPointError(iNoise, 0, hBaseLineShiftShaped[iNoise]->GetMeanError());
     
     gBaselineNoise->SetPoint(iNoise, DCounts[iNoise], hBaseLineShiftShaped[iNoise]->GetRMS());
     gBaselineNoise->SetPointError(iNoise, 0, hBaseLineShiftShaped[iNoise]->GetRMSError());
   }
   gPad->SetGrid();
   
   leg = new TLegend(0.6,0.7,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(hBaseLineShiftShaped[iNoise], Form("DCR = %.1e", DCounts[iNoise]*1e9), "lp");
   }
   leg->Draw();
   
   TCanvas * cGraphBaseline = new TCanvas ("cGraphBaseline", "cGraphBaseline", 450, 450);
   gBaselineShift->Draw("ALPE");
   gBaselineShift->GetXaxis()->SetTitle("DCR [counts per ns]");
   gBaselineShift->GetYaxis()->SetTitle("Mean baseline");
   TF1 * fitBaseLine = new TF1 ("fitBaseLine", "[0] + [1]*x", 0, 200);
   gBaselineShift->Fit(fitBaseLine, "QR");
   std::cout << "fitBaseLine->par (0) = " << fitBaseLine->GetParameter(0) << " :: fitBaseLine->par(1) = " << fitBaseLine->GetParameter(1) << std::endl;
   
   TCanvas * cGraphBaselineNoise = new TCanvas ("cGraphBaselineNoise", "cGraphBaselineNoise", 450, 450);
   gBaselineNoise->Draw("ALPE");
   gBaselineNoise->GetXaxis()->SetTitle("DCR [counts per ns]");
   gBaselineNoise->GetYaxis()->SetTitle("Pedestal noise");
   TF1 * fitBaseLineNoise = new TF1 ("fitBaseLineNoise", "[0]*sqrt(x)", 0, 200);
   gBaselineNoise->Fit(fitBaseLineNoise, "QR");
   std::cout << "fitBaseLineNoise->par (0) = " << fitBaseLineNoise->GetParameter(0) << std::endl;
   
   TCanvas * cPDE_drop = new TCanvas ("cPDE_drop", "cPDE_drop", 450, 450);
   gPDE_drop->Draw("ALPE");
   gPDE_drop->GetXaxis()->SetTitle("DCR [counts per ns]");
   gPDE_drop->GetYaxis()->SetTitle("PDE scint [%]");
//    TF1 * fitBaseLine = new TF1 ("fitBaseLine", "[0] + [1]*x", 0, 200);
//    gPDE_drop->Fit(fitBaseLine, "QR");
   
   // dark counts only
   TCanvas * cDarkCounts = new TCanvas ("cDarkCounts", "cDarkCounts", 450, 450);
   hDarkCount[0]->SetLineColor(kBlack);
   hDarkCount[0]->SetStats(0);
   hDarkCount[0]->SetTitle("Dark counts time distribution");
   hDarkCount[0]->GetXaxis()->SetTitle("Time [ps]");
   hDarkCount[0]->GetYaxis()->SetTitle("Number of dark counts");
   hDarkCount[0]->GetYaxis()->SetTitleOffset(1.3);   
   hDarkCount[0]->Draw();
   
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     hDarkCount[iNoise]->SetLineColor(iNoise+1);
     hDarkCount[iNoise]->Draw("same");
   }
   gPad->SetGrid();
   
   leg = new TLegend(0.6,0.7,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(hDarkCount[iNoise], Form("DCR = %.1e", DCounts[iNoise]*1e9), "lp");
   }
   leg->Draw();

   
   int sel_length = 0;
   
   TCanvas * cgCTR = new TCanvas ("cgCTR", "cgCTR", 500, 500);
   g_CTR_NPH[0]->Draw("ALPE");
   g_CTR_NPH[0]->GetXaxis()->SetTitle("Threshold [number of photons]");
   g_CTR_NPH[0]->GetYaxis()->SetTitle("#sigma_{t} [ps]");
   g_CTR_NPH[0]->GetYaxis()->SetTitleOffset(1.4);
   g_CTR_NPH[0]->GetXaxis()->SetLimits(0, nPhotons[NPH-1]);
   g_CTR_NPH[0]->GetYaxis()->SetRangeUser(0, 100);
   g_CTR_NPH[0]->SetLineWidth(2);
   
   for (int iNoise = 1; iNoise < nNoise; iNoise++)
   {
     g_CTR_NPH[iNoise]->SetLineWidth(2);
     g_CTR_NPH[iNoise]->SetLineColor(iNoise+1);
     g_CTR_NPH[iNoise]->Draw("LPE same");
   }
   
   gPad->SetGrid();
   
   leg = new TLegend(0.58,0.68,0.88,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(g_CTR_NPH[iNoise], Form("DCR = %.1e Hz", DCounts[iNoise]*1e9), "lp");
   }
   leg->Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();
   
   TCanvas * cgCTR_corr = new TCanvas ("cgCTR_corr", "cgCTR_corr", 500, 500);
   g_CTR_NPH_corr[0]->Draw("ALPE");
   g_CTR_NPH_corr[0]->SetTitle("SIMULATION: LYSO crystal 12x12 mm^{2} + SiPM 5x5 mm^{2}");
   g_CTR_NPH_corr[0]->GetXaxis()->SetTitle("Threshold [number of photons]");
   g_CTR_NPH_corr[0]->GetYaxis()->SetTitle("#sigma_{t} [ps]");
   g_CTR_NPH_corr[0]->GetYaxis()->SetTitleOffset(1.4);
   g_CTR_NPH_corr[0]->GetXaxis()->SetLimits(0, nPhotons[NPH-1]);
   g_CTR_NPH_corr[0]->GetYaxis()->SetRangeUser(0, 100);
   g_CTR_NPH_corr[0]->SetLineWidth(2);
   
   for (int iNoise = 1; iNoise < nNoise; iNoise++)
   {
     g_CTR_NPH_corr[iNoise]->SetLineWidth(2);
     g_CTR_NPH_corr[iNoise]->SetLineColor(iNoise+1);
     g_CTR_NPH_corr[iNoise]->Draw("LPE same");
   }
   
   gPad->SetGrid();
   
   leg = new TLegend(0.58,0.68,0.88,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(g_CTR_NPH_corr[iNoise], Form("DCR = %.1e Hz", DCounts[iNoise]*1e9), "lp");
   }
   leg->Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();
   
   TCanvas * cCTR = new TCanvas ("cCTR", "cCTR", 500, 500);
   hCTR[0][selNPH]->SetLineColor(kBlack);
   hCTR[0][selNPH]->Draw();
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     hCTR[iNoise][selNPH]->SetLineColor(iNoise+1);
     hCTR[iNoise][selNPH]->Draw("same");
   }
   gPad->SetGrid();
   
   leg = new TLegend(0.6,0.7,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(hCTR[iNoise][selNPH], Form("DCR = %.1e Hz", DCounts[iNoise]*1e9), "lp");
   }
   leg->Draw();
   
   int selNoise = 0;
   TCanvas * cCTR_NPH = new TCanvas ("cCTR_NPH", "cCTR_NPH", 500, 500);
   hCTR[selNoise][selNPH]->SetLineColor(kBlack);
   hCTR[selNoise][selNPH]->Draw();
   for (int iNPH = 0; iNPH < NPH; iNPH++)
   {
     hCTR[selNoise][iNPH]->SetLineColor(iNPH+1);
     hCTR[selNoise][iNPH]->Draw("same");
     hCTR[selNoise][iNPH]->GetXaxis()->SetTitle("CTR [ps]");
     hCTR[selNoise][iNPH]->GetYaxis()->SetTitle("Counts");
   }
   gPad->SetGrid();
   
   
   leg = new TLegend(0.6,0.7,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   leg->AddEntry(hCTR[1][1], "raw", "lp");
   leg->AddEntry(hCTR_corr[1][1], "timewalk corrected", "lp");
  
   
   TCanvas * cScintPrompt = new TCanvas ("cScintPrompt", "cScintPrompt", 500, 500);
   hScint_f[0]->SetLineColor(kBlack);
   hScint_f[0]->GetXaxis()->SetTitle(Form("Number of photons in first %d ps", maxTime));
   hScint_f[0]->Draw();
   gPad->SetGrid();
   
   TCanvas * cCherPrompt = new TCanvas ("cCherPrompt", "cCherPrompt", 500, 500);
   hCher_f[0]->SetLineColor(kBlack);
   hCher_f[0]->GetXaxis()->SetTitle(Form("Number of photons in first %d ps", maxTime));
   hCher_f[0]->Draw();
   gPad->SetGrid();
   
   TCanvas * cCTR_vs_Noise = new TCanvas ("cCTR_vs_Noise", "cCTR_vs_Noise", 500, 500);
   gCTR_vs_Noise->GetXaxis()->SetTitle("DCR [Hz]");
   gCTR_vs_Noise->GetYaxis()->SetTitle("#sigma_{t} [ps]");
   gCTR_vs_Noise->Draw("ALPE");
   gCTR_vs_Noise_corr->SetLineColor(kGreen+1);
   gCTR_vs_Noise_corr->SetMarkerColor(kGreen+1);
   gCTR_vs_Noise_corr->Draw("same LPE");
   
   /*
   TCanvas * cSingleCTR = new TCanvas ("cSingleCTR", "cSingleCTR", 500, 500);
   hCTR[1][1]->Draw("same");
   hCTR[1][1]->SetLineColor(kBlack);
   hCTR_corr[1][1]->Draw("same");
   hCTR_corr[1][1]->SetLineColor(kGreen+2);
   gPad->SetGrid();
   leg->Draw();
  

   ///time plots
   TCanvas * cCherenkovTime = new TCanvas ("cCherenkovTime", "cCherenkovTime", 600, 600);
   hTime_cher_f[0]->SetLineColor(kBlack);
   hTime_cher_f[0]->Draw();
   cout << "mean arrival time for ch f = " << hTime_cher_f[0]->GetMean() << endl;
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     hTime_cher_f[iNoise]->SetLineColor(iNoise+1);
     hTime_cher_f[iNoise]->Draw("same");
//      cout << "mean ch f = " << hTime_cher_f[iNoise]->GetMean() << endl;
   }
   gPad->SetGrid();
   
   leg = new TLegend(0.6,0.7,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
      leg->AddEntry(hTime_cher_f[iNoise], Form("%d", DCounts[iNoise]), "lp");
   }
   leg->Draw();
   


   TCanvas * cScatterTime = new TCanvas ("cScatterTime", "cScatterTime", 800, 800);
   cScatterTime->Divide(3,3);
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     cScatterTime->cd(iNoise+1);
     hScatterTime[iNoise][selNPH]->Draw("COLZ");
     hScatterTime[iNoise][selNPH]->GetXaxis()->SetTitle("time xtal_1 [ps]");
     hScatterTime[iNoise][selNPH]->GetYaxis()->SetTitle("time xtal_2 [ps]");
   }
   gPad->SetGrid();
*/
   
//    TCanvas * cBuffer = new TCanvas();

   
   TCanvas * cTimeWalk[nNoise];
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
//      fitLin[iNoise] = new TF1 (Form("fitLin_%d",iNoise), "[0] + [1]*x + [2]*x*x", 0.6+iNoise*0.05, 1.4-iNoise*0.025);

     cTimeWalk[iNoise] = new TCanvas (Form("cTimeWalk_%.1f", DCounts[iNoise]), Form("cTimeWalk_%.1f", DCounts[iNoise]), 500, 500);     
     pTimeWalk[iNoise][selNPH]->Draw();
//      pTimeWalk[iNoise][selNPH]->SetStats(0);
//      pTimeWalk[iNoise][selNPH]->SetFitParameters(1);
     pTimeWalk[iNoise][selNPH]->GetXaxis()->SetTitle("amp1/amp2");
     pTimeWalk[iNoise][selNPH]->GetYaxis()->SetTitle("time1-time2 [ps]");
     pTimeWalk[iNoise][selNPH]->GetYaxis()->SetTitleOffset(1.4);
     pTimeWalk[iNoise][selNPH]->GetXaxis()->SetRangeUser(0.,2);
     gPad->SetGrid();
   }
   /*
   TCanvas * cScatterTimeWalk[nNoise];
   for (int iNoise = 0; iNoise < nNoise; iNoise++)
   {
     cScatterTimeWalk[iNoise] = new TCanvas (Form("cScatterTimeWalk_%.1f", DCounts[iNoise]), Form("cScatterTimeWalk_%.1f", DCounts[iNoise]), 500, 500);
     hTimeWalk[iNoise][selNPH]->Fit(fitLin[iNoise], "QR");
     hTimeWalk[iNoise][selNPH]->Draw("COLZ");
     hTimeWalk[iNoise][selNPH]->GetXaxis()->SetTitle("amp1/amp2");
     hTimeWalk[iNoise][selNPH]->GetYaxis()->SetTitle("time1-time2 [ps]");
     hTimeWalk[iNoise][selNPH]->GetYaxis()->SetTitleOffset(1.4);
     hTimeWalk[iNoise][selNPH]->GetXaxis()->SetRangeUser(0.5,1.5);
     gPad->SetGrid();
   }*/
   
//    TCanvas * cEmission = new TCanvas ("cEmission", "cEmission", 500, 500);
   
    ///writing to output file
//      TFile * outputFile = new TFile(Form("./graphs_histos/12um_25PDE_10ns_decay_100RO_70LC/mu_10mm_toyLSO_noise_lenght_%d.root", length_cryst), "RECREATE");
    //     TFile * outputFile = new TFile(Form("./graphs_histos/12um_25_PDE_10ns_decay_36RO_70LC/mu_10mm_toyLSO_noise_lenght_%d.root", length_cryst), "RECREATE");
//     TFile * outputFile = new TFile(Form("./graphs_histos/GLUE_12um_25PDE_10ns_decay_36RO_70LC/mu_10mm_toyLSO_noise_lenght_%d.root", length_cryst), "RECREATE");
//     TFile * outputFile = new TFile(Form("%s_%d.root", outputFileName.c_str(), length_cryst), "RECREATE");



    TFile * outputFile = new TFile(Form("/afs/cern.ch/work/m/mlucchin/TB_timing/macros/graphs_histos/PDE_scan/%s_PDE_%.0f.root", inputFileName.c_str(), scint_PDE*100), "RECREATE");
    outputFile->cd();
   
    gBaselineShift->SetName("gBaselineShift");
    gBaselineNoise->SetName("gBaselineNoise");
    gPDE_drop->SetName("gPDE_drop");

    gBaselineShift->Write();
    gBaselineNoise->Write();
    gPDE_drop->Write();

    
    for (int iNoise = 0; iNoise <nNoise; iNoise++)
    {
      for (int iNPH = 0; iNPH < NPH; iNPH++)
      {	
	hTimeWalk[iNoise][iNPH]->Write();
	hScatterTime[iNoise][iNPH]->Write();	
	fitLin[iNoise][iNPH]->SetName(Form("fitLin_%d_NPH_%d", iNoise, iNPH));
        fitLin[iNoise][iNPH]->Write();

        pTimeWalk[iNoise][iNPH]->Write();
      }           
      g_CTR_NPH[iNoise]->SetName(Form("g_CTR_NPH_noise_%d", iNoise));
      g_CTR_NPH[iNoise]->Write();
      
      g_CTR_NPH_corr[iNoise]->SetName(Form("g_CTR_NPH_corr_noise_%d", iNoise));
      g_CTR_NPH_corr[iNoise]->Write();
      
      hBaseLineShiftShaped[iNoise]->Write();
      hTime_scint_f[iNoise]->Write();
      hTime_scint_f_PDE[iNoise]->Write();
      hTime_scint_f_PDE_tot[iNoise]->Write();
      hTime_scint_f_PDE_tot_SPTR_DRAW[iNoise]->Write();
      
    }
    outputFile->Close();    
    std::cout << "output file written" << std::endl;
//     std::cout << "\n" << std::endl;
    
//     theApp -> Run();
    return 0;
//     std::cin.ignore();
//     std::cin.get();
    

}

