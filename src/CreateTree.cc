#include "CreateTree.hh"
#include <algorithm>

using namespace std ;

CreateTree* CreateTree::fInstance = NULL ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::CreateTree (TString name)
{
  if( fInstance )
  {
    return ;
  }
  
  this -> fInstance = this ;
  this -> fname     = name ;
  this -> ftree     = new TTree (name,name) ;
  
  this -> GetTree ()->Branch ("Event", &this->Event, "Event/I") ;
  
  inputInitialPosition = new vector<float>(3,0.); 
  inputMomentum =        new vector<float>(4,0.); 
  this -> GetTree() -> Branch("inputInitialPosition", "vector<float>", &inputInitialPosition);
  this -> GetTree() -> Branch("inputMomentum",        "vector<float>", &inputMomentum);
  
  this -> GetTree() -> Branch("depositedEnergyTotal",     &this->depositedEnergyTotal,         "depositedEnergyTotal/F");
  this -> GetTree() -> Branch("depositedEnergyCore_f",      &this->depositedEnergyCore_f,           "depositedEnergyCore_f/F");
  this -> GetTree() -> Branch("depositedEnergyCore_r",      &this->depositedEnergyCore_r,           "depositedEnergyCore_r/F");
//   this -> GetTree() -> Branch("depositedEnergyCapillary", &this->depositedEnergyCapillary, "depositedEnergyCapillary/F");
//   this -> GetTree() -> Branch("depositedEnergyCladding",  &this->depositedEnergyCladding,   "depositedEnergyCladding/F");
  this -> GetTree() -> Branch("depositedEnergyWorld",     &this->depositedEnergyWorld,         "depositedEnergyWorld/F");
  
  this -> GetTree() -> Branch("time_prod_scint", &time_prod_scint);
  this -> GetTree() -> Branch("time_prod_cher",  &time_prod_cher);
  this -> GetTree() -> Branch("time_ext_scint",  &time_ext_scint);
  this -> GetTree() -> Branch("time_ext_cher", 	 &time_ext_cher);
  
  this -> GetTree() -> Branch("time_prod_scint_ref", &time_prod_scint_ref);
  this -> GetTree() -> Branch("time_prod_cher_ref",  &time_prod_cher_ref);
  this -> GetTree() -> Branch("time_ext_scint_ref",  &time_ext_scint_ref);
  this -> GetTree() -> Branch("time_ext_cher_ref",   &time_ext_cher_ref);
    
  this -> GetTree() -> Branch("tot_phot_sci",        &this->tot_phot_sci,               "tot_phot_sci/I");
  this -> GetTree() -> Branch("tot_phot_cer",        &this->tot_phot_cer,               "tot_phot_cer/I");
//   this -> GetTree() -> Branch("tot_latGap_phot_sci", &this->tot_latGap_phot_sci, "tot_latGap_phot_sci/I");
//   this -> GetTree() -> Branch("tot_latGap_phot_cer", &this->tot_latGap_phot_cer, "tot_latGap_phot_cer/I");
  this -> GetTree() -> Branch("tot_gap_phot_sci",    &this->tot_gap_phot_sci,       "tot_gap_phot_sci/I");
  this -> GetTree() -> Branch("tot_gap_phot_cer",    &this->tot_gap_phot_cer,       "tot_gap_phot_cer/I");
  this -> GetTree() -> Branch("tot_det_phot_sci",    &this->tot_det_phot_sci,       "tot_det_phot_sci/I");
  this -> GetTree() -> Branch("tot_det_phot_cer",    &this->tot_det_phot_cer,       "tot_det_phot_cer/I");
  
  h_phot_sci_lambda = new TH1F("h_phot_sci_lambda","",1000,250.,1250.);
//   h_phot_sci_E = new TH1F("h_phot_sci_E","",1000,0.,5.);
  h_phot_sci_time = new TH1F("h_phot_sci_time","",100000,0.,100000.);
//   h_phot_sci_angleAtProduction = new TH1F("h_phot_sci_angleAtProduction","",2000,-1.,1.);
  h_phot_cer_lambda = new TH1F("h_phot_cer_lambda","",1000,250.,1250.);
//   h_phot_cer_E = new TH1F("h_phot_cer_E","",1000,0.,5.);
  h_phot_cer_time = new TH1F("h_phot_cer_time","",100000,0.,100000.);
//   h_phot_cer_angleAtProduction = new TH1F("h_phot_cer_angleAtProduction","",2000,-1.,1.);
  
//   h_phot_sci_latGap_lambda = new TH1F("h_phot_sci_latGap_lambda","",1000,250.,1250.);
//   h_phot_sci_latGap_E = new TH1F("h_phot_sci_latGap_E","",1000,0.,5.);
//   h_phot_sci_latGap_time = new TH1F("h_phot_sci_latGap_time","",10000,0.,10000.);
//   h_phot_sci_latGap_angleAtProduction = new TH1F("h_phot_sci_latGap_angleAtProduction","",2000,-1.,1.);
//   h_phot_sci_latGap_angleWithSurfNormal = new TH1F("h_phot_sci_latGap_angleWithSurfNormal","",2000,-1.,1.);
//   h_phot_cer_latGap_lambda = new TH1F("h_phot_cer_latGap_lambda","",1000,250.,1250.);
//   h_phot_cer_latGap_E = new TH1F("h_phot_cer_latGap_E","",1000,0.,5.);
//   h_phot_cer_latGap_time = new TH1F("h_phot_cer_latGap_time","",10000,0.,10000.);
//   h_phot_cer_latGap_angleAtProduction = new TH1F("h_phot_cer_latGap_angleAtProduction","",2000,-1.,1.);
//   h_phot_cer_latGap_angleWithSurfNormal = new TH1F("h_phot_cer_latGap_angleWithSurfNormal","",2000,-1.,1.);
  
  h_phot_sci_gap_lambda = new TH1F("h_phot_sci_gap_lambda","",1000,250.,1250.);
//   h_phot_sci_gap_E = new TH1F("h_phot_sci_gap_E","",1000,0.,5.);
  h_phot_sci_gap_time = new TH1F("h_phot_sci_gap_time","",100000,0.,100000.);
//   h_phot_sci_gap_angleAtProduction = new TH1F("h_phot_sci_gap_angleAtProduction","",2000,-1.,1.);
//   h_phot_sci_gap_angleWithSurfNormal = new TH1F("h_phot_sci_gap_angleWithSurfNormal","",2000,-1.,1.);
  h_phot_cer_gap_lambda = new TH1F("h_phot_cer_gap_lambda","",1000,250.,1250.);
//   h_phot_cer_gap_E = new TH1F("h_phot_cer_gap_E","",1000,0.,5.);
  h_phot_cer_gap_time = new TH1F("h_phot_cer_gap_time","",100000,0.,100000.);
//   h_phot_cer_gap_angleAtProduction = new TH1F("h_phot_cer_gap_angleAtProduction","",2000,-1.,1.);
//   h_phot_cer_gap_angleWithSurfNormal = new TH1F("h_phot_cer_gap_angleWithSurfNormal","",2000,-1.,1.);
  
//   h_phot_sci_befgap_angleAtProduction = new TH1F("h_phot_sci_befgap_angleAtProduction","",2000,-1.,1.);
//   h_phot_cer_befgap_angleAtProduction = new TH1F("h_phot_cer_befgap_angleAtProduction","",2000,-1.,1.);
  
//   h_phot_sci_befdet_angleWithSurfNormal = new TH1F("h_phot_sci_befdet_angleWithSurfNormal","",2000,-1.,1.);
//   h_phot_cer_befdet_angleWithSurfNormal = new TH1F("h_phot_cer_befdet_angleWithSurfNormal","",2000,-1.,1.);
//   h_phot_sci_det_angleWithSurfNormal = new TH1F("h_phot_sci_det_angleWithSurfNormal","",2000,-1.,1.);
//   h_phot_cer_det_angleWithSurfNormal = new TH1F("h_phot_cer_det_angleWithSurfNormal","",2000,-1.,1.);
  
  this -> Clear() ;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



CreateTree::~CreateTree()
{}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



int CreateTree::Fill() 
{ 
  return this -> GetTree() -> Fill(); 
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



bool CreateTree::Write(TFile * outfile)
{
  outfile -> cd();
  ftree -> Write();
  
  h_phot_sci_lambda->Write();
//   h_phot_sci_E->Write();
  h_phot_sci_time->Write();
//   h_phot_sci_angleAtProduction -> Write();
  h_phot_cer_lambda->Write();
//   h_phot_cer_E->Write();
  h_phot_cer_time->Write();
//   h_phot_cer_angleAtProduction -> Write();
  /*
  h_phot_sci_latGap_lambda->Write();
  h_phot_sci_latGap_E->Write();
  h_phot_sci_latGap_time->Write();
  h_phot_sci_latGap_angleAtProduction->Write();
  h_phot_sci_latGap_angleWithSurfNormal->Write();
  h_phot_cer_latGap_lambda->Write();
  h_phot_cer_latGap_E->Write();
  h_phot_cer_latGap_time->Write();
  h_phot_cer_latGap_angleAtProduction->Write();
  h_phot_cer_latGap_angleWithSurfNormal->Write();*/
  
  h_phot_sci_gap_lambda->Write();
//   h_phot_sci_gap_E->Write();
  h_phot_sci_gap_time->Write();
//   h_phot_sci_gap_angleAtProduction->Write();
//   h_phot_sci_gap_angleWithSurfNormal->Write();
  h_phot_cer_gap_lambda->Write();
//   h_phot_cer_gap_E->Write();
  h_phot_cer_gap_time->Write();
//   h_phot_cer_gap_angleAtProduction->Write();
//   h_phot_cer_gap_angleWithSurfNormal->Write();
  /*
  h_phot_sci_befgap_angleAtProduction->Write();
  h_phot_cer_befgap_angleAtProduction->Write();
  
  h_phot_sci_befdet_angleWithSurfNormal->Write();
  h_phot_cer_befdet_angleWithSurfNormal->Write();
  
  h_phot_sci_det_angleWithSurfNormal->Write();
  h_phot_cer_det_angleWithSurfNormal->Write();*/
  
  return true ;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



void CreateTree::Clear()
{
  Event	= 0;
  
  depositedEnergyTotal = 0.;
  depositedEnergyCore_f = 0.;
  depositedEnergyCore_r = 0.;
//   depositedEnergyCapillary = 0.;
//   depositedEnergyCladding = 0.;
  depositedEnergyWorld = 0.;
  
  tot_phot_sci = 0;
  tot_phot_cer = 0;
//   tot_latGap_phot_sci = 0;
//   tot_latGap_phot_cer = 0;
  tot_gap_phot_sci = 0;
  tot_gap_phot_cer = 0;
  tot_det_phot_sci = 0;
  tot_det_phot_cer = 0;
  
  for (int i = 0 ; i < 3 ; ++i) 
  {
    inputInitialPosition -> at(i) = 0.;
  }
  for (int i = 0 ; i < 4 ; ++i) 
  {
    inputMomentum ->at(i) = 0.;
  }
  
  time_ext_cher.clear();
  time_ext_scint.clear();
  time_prod_cher.clear();
  time_prod_scint.clear();
  
  time_ext_cher_ref.clear();
  time_ext_scint_ref.clear();
  time_prod_cher_ref.clear();
  time_prod_scint_ref.clear();
}
