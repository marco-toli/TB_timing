#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "TString.h"
#include "TRandom3.h"
#include "TCint.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"

long int CreateSeed();



using namespace std;
using namespace CLHEP;



int to_int (string name)
{
  int Result ;             // int which will contain the result
  stringstream convert (name) ;
  string dummy ;           
  convert >> dummy ;       
  convert >> Result ;
  return Result ;
}


//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


SteppingAction::SteppingAction (DetectorConstruction* detectorConstruction,
                                const G4int& scint, const G4int& cher):
  fDetectorConstruction(detectorConstruction),
  propagateScintillation(scint),
  propagateCerenkov(cher)
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


SteppingAction::~SteppingAction ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void SteppingAction::UserSteppingAction (const G4Step * theStep)
{
 
  G4Track* theTrack = theStep->GetTrack () ;
  
  const G4ThreeVector& theTrackDirection = theTrack->GetMomentumDirection();
  const G4ThreeVector& theTrackVertexDirection = theTrack->GetVertexMomentumDirection();
  
  G4int trackID = theTrack->GetTrackID();
  TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformation());
  G4ParticleDefinition* particleType = theTrack->GetDefinition () ;
  
  G4StepPoint * thePrePoint  = theStep->GetPreStepPoint () ;
  G4StepPoint * thePostPoint = theStep->GetPostStepPoint () ;
  const G4ThreeVector & thePrePosition  = thePrePoint->GetPosition () ;
  G4VPhysicalVolume * thePrePV  = thePrePoint->GetPhysicalVolume () ;
  G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume () ;
  G4String thePrePVName  = "" ; if ( thePrePV )  thePrePVName  = thePrePV  -> GetName () ;
  G4String thePostPVName = "" ; if ( thePostPV ) thePostPVName = thePostPV -> GetName () ;
  G4VSolid* thePreS = thePrePV->GetLogicalVolume()->GetSolid();
    
  G4int nStep = theTrack -> GetCurrentStepNumber();
  
//        cout << " step length = " << theStep->GetStepLength() << endl;

  //-------------
  // get position
  G4double global_x = thePrePosition.x()/mm;
  G4double global_y = thePrePosition.y()/mm;
  G4double global_z = thePrePosition.z()/mm;
  
  
  // optical photon
  if( particleType == G4OpticalPhoton::OpticalPhotonDefinition() )
  {
    G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
    
    
    //----------------------------
    // count photons at production
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ) &&
        (nStep == 1) && (processName == "Scintillation") )
    {
      CreateTree::Instance()->tot_phot_sci += 1;
//       cout << " in log volume " << endl;
      //save only prompt photons
      if (thePrePoint->GetGlobalTime()/picosecond>1000) theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      else if (thePrePVName == "corePV") CreateTree::Instance()->time_prod_scint.push_back(thePrePoint->GetGlobalTime()/picosecond );
      else if (thePrePVName == "corePV_ref") CreateTree::Instance()->time_prod_scint_ref.push_back(thePrePoint->GetGlobalTime()/picosecond );
      
      if( !propagateScintillation ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      
      CreateTree::Instance()->h_phot_sci_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
//       CreateTree::Instance()->h_phot_sci_E      -> Fill( theTrack->GetTotalEnergy()/eV );
      CreateTree::Instance()->h_phot_sci_time   -> Fill( thePrePoint->GetGlobalTime()/ns );
//       CreateTree::Instance()->h_phot_sci_angleAtProduction -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection)) );
    }
    
        
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core")) &&
        (nStep == 1) && (processName == "Cerenkov") )
    {
      CreateTree::Instance()->tot_phot_cer += 1;
      if (thePrePVName == "corePV") CreateTree::Instance()->time_prod_cher.push_back(thePrePoint->GetGlobalTime()/picosecond );
      else if (thePrePVName == "corePV_ref") CreateTree::Instance()->time_prod_cher_ref.push_back(thePrePoint->GetGlobalTime()/picosecond );

      if( !propagateCerenkov ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
      
      CreateTree::Instance()->h_phot_cer_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
//       CreateTree::Instance()->h_phot_cer_E      -> Fill( theTrack->GetTotalEnergy()/eV );
      CreateTree::Instance()->h_phot_cer_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
//       CreateTree::Instance()->h_phot_cer_angleAtProduction -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection)) );
    }
    
    
    //-------------------------------------------
    // count photons exiting from lateral surface
    /*
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ) &&
        (processName == "Scintillation") &&
        (thePrePVName == "latGapLayerPV") && (thePostPVName == "latGapPV") )
    {
      CreateTree::Instance()->tot_latGap_phot_sci += 1;
      // if you do not want to kill a photon once it exits the fiber, comment here below
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      
      CreateTree::Instance()->h_phot_sci_latGap_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
//       CreateTree::Instance()->h_phot_sci_latGap_E      -> Fill( theTrack->GetTotalEnergy()/eV );
      CreateTree::Instance()->h_phot_sci_latGap_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
//       CreateTree::Instance()->h_phot_sci_latGap_angleAtProduction -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection)) );
//       CreateTree::Instance()->h_phot_sci_latGap_angleWithSurfNormal -> Fill( cos(thePreS->SurfaceNormal(thePrePosition).angle(theTrackDirection)) );
    }
    
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("capillary") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("cladding") ) &&
        (processName == "Cerenkov") &&
        (thePrePVName == "latGapLayerPV") && (thePostPVName == "latGapPV") )
    {
      CreateTree::Instance()->tot_latGap_phot_cer += 1;
      // if you do not want to kill a photon once it exits the fiber, comment here below
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      
      CreateTree::Instance()->h_phot_cer_latGap_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
//       CreateTree::Instance()->h_phot_cer_latGap_E      -> Fill( theTrack->GetTotalEnergy()/eV );
      CreateTree::Instance()->h_phot_cer_latGap_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
//       CreateTree::Instance()->h_phot_cer_latGap_angleAtProduction -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection)) );
//       CreateTree::Instance()->h_phot_cer_latGap_angleWithSurfNormal -> Fill( cos(thePreS->SurfaceNormal(thePrePosition).angle(theTrackDirection)) );
    }*/
    
    
    
    //----------------------------
    // count photons at fiber exit
    
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ) &&
        (processName == "Scintillation") &&
        ((thePrePVName == "gapLayerPV") && (thePostPVName == "gapPV") || (thePrePVName == "gapLayerPV_ref") && (thePostPVName == "gapPV_ref")) )
    {
      CreateTree::Instance()->tot_gap_phot_sci += 1;
      if (thePrePVName == "gapLayerPV") CreateTree::Instance()->time_ext_scint.push_back(thePrePoint->GetGlobalTime()/picosecond );
      else if (thePrePVName == "gapLayerPV_ref") CreateTree::Instance()->time_ext_scint_ref.push_back(thePrePoint->GetGlobalTime()/picosecond );
      // if you do not want to kill a photon once it exits the fiber, comment here below
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
      CreateTree::Instance()->h_phot_sci_gap_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
//       CreateTree::Instance()->h_phot_sci_gap_E      -> Fill( theTrack->GetTotalEnergy()/eV );
      CreateTree::Instance()->h_phot_sci_gap_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
//       CreateTree::Instance()->h_phot_sci_gap_angleAtProduction -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection)) );
//       CreateTree::Instance()->h_phot_sci_gap_angleWithSurfNormal -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackDirection)) );
    }
    
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("capillary") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("cladding") ) &&
        (processName == "Cerenkov") &&
        ((thePrePVName == "gapLayerPV" && thePostPVName == "gapPV") || (thePrePVName == "gapLayerPV_ref" && thePostPVName == "gapPV_ref")))
    {
      CreateTree::Instance()->tot_gap_phot_cer += 1;
      // if you do not want to kill a photon once it exits the fiber, comment here below
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      
      if (thePrePVName == "gapLayerPV") CreateTree::Instance()->time_ext_cher.push_back(thePrePoint->GetGlobalTime()/picosecond );
      else if (thePrePVName == "gapLayerPV_ref") CreateTree::Instance()->time_ext_cher_ref.push_back(thePrePoint->GetGlobalTime()/picosecond );
      
      CreateTree::Instance()->h_phot_cer_gap_lambda -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV) );
//       CreateTree::Instance()->h_phot_cer_gap_E      -> Fill( theTrack->GetTotalEnergy()/eV );
      CreateTree::Instance()->h_phot_cer_gap_time   -> Fill( thePrePoint->GetGlobalTime()/picosecond );
//       CreateTree::Instance()->h_phot_cer_gap_angleAtProduction -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackVertexDirection)) );
//       CreateTree::Instance()->h_phot_cer_gap_angleWithSurfNormal -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackDirection)) );
    }
    
   
    
    
    //------------------------------
    // count photons at the detector
    
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ) &&
        (processName == "Scintillation") &&
        (thePrePVName == "detLayerPV") && (thePostPVName == "detPV") )
    {
      CreateTree::Instance()->tot_det_phot_sci += 1;
//       CreateTree::Instance()->h_phot_sci_det_angleWithSurfNormal -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackDirection)) );
      // if you do not want to kill a photon once it enters the detector, comment here below
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }
    
    if( ( theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("capillary") ||
          theTrack->GetLogicalVolumeAtVertex()->GetName().contains("cladding") ) &&
        (processName == "Cerenkov") &&
        (thePrePVName == "detLayerPV") && (thePostPVName == "detPV") )
    {
      CreateTree::Instance()->tot_det_phot_cer += 1;
//       CreateTree::Instance()->h_phot_cer_det_angleWithSurfNormal -> Fill( cos(G4ThreeVector(0.,0.,1.).angle(theTrackDirection)) );
      // if you do not want to kill a photon once it enters the detector, comment here below
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }
    
    
    
    
    /*
    if( (theTrack->GetLogicalVolumeAtVertex()->GetName().contains("core")) && (nStep == 1) )
    {    
      //----------------------------------------------------------
      // storing time, energy and position at gap with fast timing
      Photon ph;
      ph.position.SetX(global_x);
      ph.position.SetY(global_y);
      ph.position.SetZ(global_z);
      ph.direction.SetX(theTrack->GetVertexMomentumDirection().x());
      ph.direction.SetY(theTrack->GetVertexMomentumDirection().y());
      ph.direction.SetZ(theTrack->GetVertexMomentumDirection().z());
      ph.dist = (global_z/(0.5*fiber_length));
      ph.energy = theTrack->GetTotalEnergy()/eV;
      
      Fiber* fib = fDetectorConstruction -> GetFiber();
      std::map<int,Travel> trc = GetTimeAndProbability(ph,fib,theTrackInfo->GetParticleProdTime());
      
      for(unsigned int it = 0; it < CreateTree::Instance()->attLengths->size(); ++it)
      {
        int attLength = int( CreateTree::Instance()->attLengths->at(it) );
        
        if( trc[attLength].prob[0] < 1.E-09 ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);      
        
        for(int it2 = 0; it2 < 3; ++it2)
        {
          CreateTree::Instance()->tot_gap_photFast_cer->at(it) += trc[attLength].prob[it2];
          
          //CreateTree::Instance()->h_photFast_cer_gap_lambda[attLength] -> Fill( MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV), trc[attLength].prob[it2] );
          //CreateTree::Instance()->h_photFast_cer_gap_E[attLength]      -> Fill( theTrack->GetTotalEnergy()/eV, trc[attLength].prob[it2] );
          //CreateTree::Instance()->h_photFast_cer_gap_time[attLength]   -> Fill( trc[attLength].time[it2], trc[attLength].prob[it2] );
        }
      }
    }
    */
  } // optical photon
  
  
  // non optical photon
  else
  {
    //G4cout << ">>> begin non optical photon" << G4endl;
    
    G4double energy = theStep->GetTotalEnergyDeposit() - theStep->GetNonIonizingEnergyDeposit();
    if ( energy == 0. ) return;
    
    CreateTree::Instance() -> depositedEnergyTotal += energy/GeV;
    
    if( thePrePVName == "corePV" )
    {
      CreateTree::Instance()->depositedEnergyCore_f += energy/GeV;
    }
    
    if( thePrePVName == "corePV_ref" )
    {
      CreateTree::Instance()->depositedEnergyCore_r += energy/GeV;
    }
    
    if( thePrePVName.contains("capillary") )
    {
//       CreateTree::Instance()->depositedEnergyCapillary += energy/GeV;
    }
    if( thePrePVName.contains("cladding") )
    {
//       CreateTree::Instance()->depositedEnergyCladding += energy/GeV;
    }
    if( thePrePVName.contains("world") )
    {
      CreateTree::Instance() -> depositedEnergyWorld += energy/GeV;
    }
    
    //G4cout << ">>> end non optical photon" << G4endl;
  } // non optical photon
  
  
  return ;
}

long int CreateSeed()
{
  TRandom3 rangen;
  
  long int sec = time(0);
  G4cout << "Time : " << sec << G4endl;
  
  sec += getpid();
  G4cout << "PID  : " << getpid() << G4endl;
  
  FILE* fp = fopen ("/proc/uptime", "r");
  int upsecs = 0;
  if( fp != NULL )
  {
    char buf[BUFSIZ];
    char *b = fgets(buf,BUFSIZ,fp);
    if( b == buf )
    {
      /* The following sscanf must use the C locale.  */
      setlocale(LC_NUMERIC, "C");
      setlocale(LC_NUMERIC, "");
    }
    fclose(fp);
  }
  G4cout << "Upsecs: " << upsecs << G4endl;
  sec += upsecs;
  
  G4cout << "Seed for srand: " << sec << G4endl;
  srand(sec);
  rangen.SetSeed(rand());
  long int seed = round(1000000000000*rangen.Uniform());
  return seed;
}