#ifndef CreateTree_H
#define CreateTree_H 1

#include <iostream>
#include <vector>
#include "TString.h"

#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

class CreateTree
{
private:
  
  TTree*  ftree ;
  TString fname ;
  
public:
  
  CreateTree (TString name);
  ~CreateTree () ;
  
  TTree*             GetTree() const { return ftree; };
  TString            GetName() const { return fname; };
  void               AddEnergyDeposit(int index, float deposit);
  void               AddScintillationPhoton(int index);
  void               AddCerenkovPhoton(int index);
  int                Fill();
  bool               Write(TFile *);
  void               Clear() ;
  
  static CreateTree* Instance() { return fInstance; } ;
  static CreateTree* fInstance;
  
  int Event;
  
  std::vector<float>* inputMomentum ; // Px Py Pz E
  std::vector<float>* inputInitialPosition ; // x, y, z
   
  std::vector<float> time_ext_scint;
  std::vector<float> time_ext_cher;
  std::vector<float> time_prod_scint;
  std::vector<float> time_prod_cher;
  
  std::vector<float> time_ext_scint_ref;
  std::vector<float> time_ext_cher_ref;
  std::vector<float> time_prod_scint_ref;
  std::vector<float> time_prod_cher_ref;
  
  float depositedEnergyTotal;
  float depositedEnergyCore_f;
  float depositedEnergyCore_r;
  float depositedEnergyCapillary;
  float depositedEnergyCladding;
  float depositedEnergyWorld;
  
  int tot_phot_sci;
  int tot_phot_cer;
  int tot_latGap_phot_sci;
  int tot_latGap_phot_cer;
  int tot_gap_phot_sci;
  int tot_gap_phot_cer;
  int tot_det_phot_sci;
  int tot_det_phot_cer;
  
  TH1F* h_phot_sci_lambda;
  TH1F* h_phot_sci_E;
  TH1F* h_phot_sci_time;
  TH1F* h_phot_sci_angleAtProduction;
  TH1F* h_phot_cer_lambda;
  TH1F* h_phot_cer_E;
  TH1F* h_phot_cer_time;
  TH1F* h_phot_cer_angleAtProduction;
  
  TH1F* h_phot_sci_latGap_lambda;
  TH1F* h_phot_sci_latGap_E;
  TH1F* h_phot_sci_latGap_time;
  TH1F* h_phot_sci_latGap_angleAtProduction;
  TH1F* h_phot_sci_latGap_angleWithSurfNormal;
  TH1F* h_phot_cer_latGap_lambda;
  TH1F* h_phot_cer_latGap_E;
  TH1F* h_phot_cer_latGap_time;
  TH1F* h_phot_cer_latGap_angleAtProduction;
  TH1F* h_phot_cer_latGap_angleWithSurfNormal;
  
  TH1F* h_phot_sci_gap_lambda;
  TH1F* h_phot_sci_gap_E;
  TH1F* h_phot_sci_gap_time;
  TH1F* h_phot_sci_gap_angleAtProduction;
  TH1F* h_phot_sci_gap_angleWithSurfNormal;
  TH1F* h_phot_cer_gap_lambda;
  TH1F* h_phot_cer_gap_E;
  TH1F* h_phot_cer_gap_time;
  TH1F* h_phot_cer_gap_angleAtProduction;
  TH1F* h_phot_cer_gap_angleWithSurfNormal;
  
  TH1F* h_phot_sci_befgap_angleAtProduction;
  TH1F* h_phot_cer_befgap_angleAtProduction;
  
  TH1F* h_phot_sci_befdet_angleWithSurfNormal;
  TH1F* h_phot_cer_befdet_angleWithSurfNormal;
  TH1F* h_phot_sci_det_angleWithSurfNormal;
  TH1F* h_phot_cer_det_angleWithSurfNormal;
};

#endif
