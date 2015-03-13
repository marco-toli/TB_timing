//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

#include "DetectorConstruction.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"


#include "DetectorConstruction.hh"
#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include "CreateTree.hh"
#include <algorithm>
#include <string>
#include <sstream>

using namespace CLHEP;



DetectorConstruction::DetectorConstruction (const string& configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------
  
  ConfigFile config (configFileName) ;
  
  config.readInto(checkOverlaps,"checkOverlaps");
  
  config.readInto(world_material,	"world_material");
  config.readInto(fibre_length,  	"fibre_length");
  config.readInto(fibre_isSquare,	"fibre_isSquare");
  config.readInto(detector,      	"detector");
  
  config.readInto(core_radius,   "core_radius");
  config.readInto(core_material, "core_material");
  config.readInto(core_rIndex,   "core_rIndex");
  config.readInto(core_absLength,"core_absLength");
  
  config.readInto(capillary_thickness,"capillary_thickness");
  config.readInto(capillary_material, "capillary_material");
  config.readInto(capillary_rIndex,   "capillary_rIndex");
  config.readInto(capillary_absLength,"capillary_absLength");
  
  config.readInto(cladding_thickness,"cladding_thickness");
  config.readInto(cladding_material, "cladding_material");
  config.readInto(cladding_rIndex,   "cladding_rIndex");
  config.readInto(cladding_absLength,"cladding_absLength");
  
  config.readInto(gap_l,       "gap_l");
  config.readInto(gap_material,"gap_material");
  
  config.readInto(det_l,       "det_l");
  config.readInto(det_material,"det_material");
  
  config.readInto(depth,	"depth");
  config.readInto(cryst_dist,	"cryst_dist");
  config.readInto(abs_thick,	"abs_thick");
  
  B_field_intensity = config.read<double>("B_field_intensity") * tesla ;
  
  expHall_x = 100.*cm;
  expHall_y = 100.*cm;
  expHall_z = 100.*cm;
  
  B_field_IsInitialized = false ;
  
  initializeMaterials();
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



DetectorConstruction::~DetectorConstruction ()
{
  delete stepLimit;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



G4VPhysicalVolume* DetectorConstruction::Construct ()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl ;
  
  
  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------
  
  
  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid * worldS = new G4Box ("worldS", 0.5 * expHall_x, 0.5 * expHall_y, 0.5 * expHall_z) ;
  G4LogicalVolume * worldLV = new G4LogicalVolume (worldS, WoMaterial, "worldLV", 0, 0, 0) ;
  G4VPhysicalVolume * worldPV = new G4PVPlacement (0, G4ThreeVector (), worldLV, "worldPV", 0, false, 0, checkOverlaps) ;
  
  // the pre-shower
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid* absorberS;
  absorberS = new G4Box  ("absorberS", core_radius*4, core_radius*4, 0.5*abs_thick);    
  G4LogicalVolume* absorberLV = new G4LogicalVolume (absorberS, MyMaterials::Lead(), "absorberLV") ;
  new G4PVPlacement(0, G4ThreeVector(0.,0., -0.5 * (fibre_length+abs_thick) - 3*mm), absorberLV, "absorberPV", worldLV, false, 0, checkOverlaps) ;
  
  
  
  // the first crystal
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid* coreS;
  if( !fibre_isSquare ) coreS = new G4Tubs ("coreS", 0., core_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;    
  else                  coreS = new G4Box  ("coreS", core_radius, core_radius, 0.5*fibre_length) ;    
  G4LogicalVolume* coreLV = new G4LogicalVolume (coreS, CoMaterial, "coreLV");
  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), coreLV, "corePV", worldLV, false, 0, checkOverlaps) ;
  new G4PVPlacement(0, G4ThreeVector(0.,0.,cryst_dist), coreLV, "corePV_ref", worldLV, false, 0, checkOverlaps) ;
  
  G4VSolid* capillaryS;
  G4LogicalVolume* capillaryLV;
  if( capillary_thickness > 0 )
  {
    if( !fibre_isSquare ) capillaryS = new G4Tubs ("capillaryS", core_radius, core_radius+capillary_thickness, 0.5*fibre_length, 0.*deg, 360.*deg) ;
    else
    {
      G4VSolid* dummyS = new G4Box ("dummyS", core_radius+capillary_thickness, core_radius+capillary_thickness, 0.5*fibre_length) ;
      G4VSolid* subS = new G4Box ("subS", core_radius, core_radius, 0.51*fibre_length);
      capillaryS = new G4SubtractionSolid ("capillaryS", dummyS, subS);    
    }
    capillaryLV = new G4LogicalVolume (capillaryS, CaMaterial, "capillaryLV") ;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), capillaryLV, "capillaryPV", worldLV, false, 0, checkOverlaps) ;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,cryst_dist), capillaryLV, "capillaryPV_ref", worldLV, false, 0, checkOverlaps) ;
  }
  
  G4VSolid* claddingS;
  G4LogicalVolume* claddingLV;
  if( cladding_thickness > 0 )
  {
    if( !fibre_isSquare ) claddingS = new G4Tubs ("claddingS", core_radius+capillary_thickness, core_radius+capillary_thickness+cladding_thickness, 0.5*fibre_length, 0.*deg, 360.*deg) ;    
    else
    {
      G4VSolid* dummyS = new G4Box ("dummyS", core_radius+capillary_thickness+cladding_thickness, core_radius+capillary_thickness+cladding_thickness, 0.5*fibre_length) ;
      G4VSolid* subS = new G4Box ("subS", core_radius+capillary_thickness, core_radius+capillary_thickness, 0.51*fibre_length);
      claddingS = new G4SubtractionSolid ("claddingS", dummyS, subS);
    }
    claddingLV = new G4LogicalVolume (claddingS, ClMaterial, "claddingLV") ;  
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), claddingLV, "claddingPV", worldLV, false, 0, checkOverlaps) ;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,cryst_dist), claddingLV, "claddingPV_ref", worldLV, false, 0, checkOverlaps) ;
  }
  
  
  G4double total_radius = core_radius + capillary_thickness + cladding_thickness;
  G4VSolid* totalS;
  G4VSolid* totalDepthS;
  if( fibre_isSquare ) totalS = new G4Box("totalS",total_radius,total_radius,0.51*fibre_length);
  if( fibre_isSquare ) totalDepthS = new G4Box("totalS",total_radius+depth,total_radius+depth,0.51*fibre_length);
  
  
  // fibre gap for photon counting
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  // lateral gaps
  G4VSolid* latGapLayerS;
  double lat_gap = 0.1*mm;
  if( !fibre_isSquare ) latGapLayerS = new G4Tubs ("latGapLayerS", total_radius, total_radius+depth, 0.5*fibre_length, 0.*deg, 360.*deg) ;
  else
  {
    G4VSolid* dummyS = new G4Box ("dummyS", total_radius+depth, total_radius+depth, 0.5*fibre_length) ;
    latGapLayerS = new G4SubtractionSolid("latGapLayerS", dummyS, totalS);
  }

  
  G4VSolid* latGapS;
  if( !fibre_isSquare ) latGapS = new G4Tubs ("latGapS", total_radius+depth, total_radius+lat_gap, 0.5*fibre_length, 0.*deg, 360.*deg) ;
  else
  {
    G4VSolid* dummyS = new G4Box ("dummyS", total_radius+lat_gap, total_radius+lat_gap, 0.5*fibre_length) ;
    latGapS = new G4SubtractionSolid ("latGapS", dummyS, totalDepthS);
  }
  G4LogicalVolume* latGapLayerLV = new G4LogicalVolume (latGapLayerS, GaMaterial, "latGapLayerLV") ;
  G4LogicalVolume* latGapLV      = new G4LogicalVolume (latGapS,      GaMaterial,      "latGapLV") ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.), latGapLayerLV, "latGapLayerPV", worldLV, false, 0, checkOverlaps) ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.),      latGapLV,      "latGapPV", worldLV, false, 0, checkOverlaps) ;
  
  new G4PVPlacement (0, G4ThreeVector (0., 0., cryst_dist), latGapLayerLV, "latGapLayerPV_ref", worldLV, false, 0, checkOverlaps) ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., cryst_dist),      latGapLV,      "latGapPV_ref", worldLV, false, 0, checkOverlaps) ;
  
  //end gaps
  G4VSolid* gapLayerS;
  if( !fibre_isSquare ) gapLayerS = new G4Tubs ("gapLayerS", 0., total_radius, 0.5*depth, 0.*deg, 360.*deg) ;
  else                  gapLayerS = new G4Box  ("gapLayerS", total_radius, total_radius, 0.5*depth) ;
  
  G4VSolid* gapS;
  if( !fibre_isSquare ) gapS      = new G4Tubs (     "gapS", 0., total_radius, 0.5*(gap_l-depth), 0.*deg, 360.*deg) ;
  else                  gapS      = new G4Box  (     "gapS", total_radius, total_radius, 0.5*(gap_l-depth)) ;
  
  G4LogicalVolume* gapLayerLV = new G4LogicalVolume (gapLayerS, GaMaterial, "gapLayerLV") ;
  G4LogicalVolume* gapLV      = new G4LogicalVolume (gapS,      GaMaterial,      "gapLV") ;
  
  
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*fibre_length+0.5*depth),          gapLayerLV, "gapLayerPV", worldLV, false, 0, checkOverlaps) ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*fibre_length+depth+0.5*(gap_l-depth)), gapLV,      "gapPV", worldLV, false, 0, checkOverlaps) ;
  
  new G4PVPlacement (0, G4ThreeVector (0., 0., cryst_dist + .5*fibre_length+0.5*depth),          gapLayerLV, "gapLayerPV_ref", worldLV, false, 0, checkOverlaps) ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., cryst_dist + .5*fibre_length+depth+0.5*(gap_l-depth)), gapLV,      "gapPV_ref", worldLV, false, 0, checkOverlaps) ;
  

  
  // Si detector for photon counting
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid * detLayerS = NULL;
  G4VSolid * detS      = NULL;
  G4LogicalVolume * detLayerLV = NULL;
  G4LogicalVolume * detLV      = NULL;
  if( detector )
  {
    detLayerS = new G4Box ("detLayerS", 0.5*3.*mm, 0.5*3.*mm, 0.5*depth) ;
    detS      = new G4Box (     "detS", 0.5*3.*mm, 0.5*3.*mm, 0.5*(det_l-depth)) ;
    detLayerLV = new G4LogicalVolume (detLayerS, DeMaterial, "detLayerLV") ;
    detLV      = new G4LogicalVolume (detS,      DeMaterial,      "detLV") ;
    //for standard conf
    new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*fibre_length+gap_l+0.5*depth),          detLayerLV, "detLayerPV", worldLV, false, 0, checkOverlaps) ;
    new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*fibre_length+gap_l+depth+0.5*(det_l-depth)), detLV,      "detPV", worldLV, false, 0, checkOverlaps) ;
    
    new G4PVPlacement (0, G4ThreeVector (0., 0., cryst_dist + 0.5*fibre_length+gap_l+0.5*depth),          detLayerLV, "detLayerPV_ref", worldLV, false, 0, checkOverlaps) ;
    new G4PVPlacement (0, G4ThreeVector (0., 0., cryst_dist + 0.5*fibre_length+gap_l+depth+0.5*(det_l-depth)), detLV,      "detPV_ref", worldLV, false, 0, checkOverlaps) ;
    

  }
  
  
  
  //-----------------------------------------------------
  //------------- Visualization attributes --------------
  //-----------------------------------------------------
  
  G4Colour white  (1.00, 1.00, 1.00);  // white
  G4Colour gray   (0.50, 0.50, 0.50);  // gray
  G4Colour black  (0.00, 0.00, 0.00);  // black
  G4Colour red    (1.00, 0.00, 0.00);  // red
  G4Colour green  (0.00, 1.00, 0.00);  // green
  G4Colour blue   (0.00, 0.00, 1.00);  // blue
  G4Colour cyan   (0.00, 1.00, 1.00);  // cyan
  G4Colour air    (0.90, 0.94, 1.00);  // cyan
  G4Colour magenta(1.00, 0.00, 1.00);  // magenta 
  G4Colour yellow (1.00, 1.00, 0.00);  // yellow
  G4Colour brass  (0.80, 0.60, 0.40);  // brass
  G4Colour brown  (0.70, 0.40, 0.10);  // brown
  
  G4VisAttributes* VisAttWorld = new G4VisAttributes(black);
  VisAttWorld -> SetVisibility(true) ;
  VisAttWorld -> SetForceWireframe(true) ;
  worldLV -> SetVisAttributes(VisAttWorld) ;
  
  G4VisAttributes* VisAttCore = new G4VisAttributes(green);
  VisAttCore -> SetVisibility(true);
  VisAttCore -> SetForceWireframe(true);
  coreLV -> SetVisAttributes(VisAttCore);
  
  if( capillary_thickness > 0 )
  {
    G4VisAttributes* VisAttCapillary = new G4VisAttributes(blue);
    VisAttCapillary -> SetVisibility(true);
    VisAttCapillary -> SetForceWireframe(true);
    capillaryLV -> SetVisAttributes(VisAttCapillary);
  }
  
  if( cladding_thickness > 0 )
  {
    G4VisAttributes* VisAttCladding = new G4VisAttributes(yellow);
    VisAttCladding -> SetVisibility(true);
    VisAttCladding -> SetForceWireframe(true);
    claddingLV -> SetVisAttributes(VisAttCladding);
  }
  
  G4VisAttributes* VisAttGapLayer = new G4VisAttributes(red);
  VisAttGapLayer -> SetVisibility(true);
  VisAttGapLayer -> SetForceWireframe(true);
  gapLayerLV -> SetVisAttributes(VisAttGapLayer);
  latGapLayerLV -> SetVisAttributes(VisAttGapLayer);
  
  G4VisAttributes* VisAttGap = new G4VisAttributes(gray);
  VisAttGap -> SetVisibility(true);
  VisAttGap -> SetForceWireframe(true);
  gapLV -> SetVisAttributes(VisAttGap);
  latGapLV -> SetVisAttributes(VisAttGap);
  
  
  if( detector )
  {
    G4VisAttributes* VisAttDetLayer = new G4VisAttributes(red);
    VisAttDetLayer -> SetVisibility(true);
    VisAttDetLayer -> SetForceWireframe(false);
    detLayerLV -> SetVisAttributes(VisAttDetLayer);
    
    G4VisAttributes* VisAttDet = new G4VisAttributes(gray);
    VisAttDet -> SetVisibility(true);
    VisAttDet -> SetForceWireframe(false);
    detLV -> SetVisAttributes(VisAttDet);
  }
  
  if (B_field_intensity > 0.1 * tesla) ConstructField () ; 
  
  
  
  // //-----------------------------------------------
  // //------------- Fast photon timing --------------
  // //-----------------------------------------------
  
  // std::vector<std::pair<double,double> > rIndVecCore;
  // std::vector<std::pair<double,double> > rIndVecClad;
  // std::vector<std::pair<double,double> > rIndVecAir;
  // std::vector<std::pair<double,double> > rIndVecGap;
  
  // G4MaterialPropertyVector* mpVec;
  
  // mpVec = ClMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   rIndVecCore.push_back(dummy);
  // }
  
  // mpVec = WoMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   std::pair<double,double> dummy2(mpVec->GetLowEdgeEnergy(it)/eV,fibre_cladRIndex);
  //   rIndVecAir.push_back(dummy);
  //   rIndVecClad.push_back(dummy2);
  // }
  
  // mpVec = GaMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   rIndVecGap.push_back(dummy);
  // }
  
  
  // fib = FiberInit(fibre_length,fibre_radius,CreateTree::Instance()->attLengths,rIndVecCore,rIndVecClad,rIndVecAir,rIndVecGap) ;
  
  
  
  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl ;
  return worldPV ;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



void DetectorConstruction::initializeMaterials ()
{
  //-----------------
  // define materials
  
  WoMaterial = NULL ;
  if      ( world_material == 1 ) WoMaterial = MyMaterials::Air () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Wo. material: "<< WoMaterial << G4endl ;
  
  
  CoMaterial = NULL ;
  if      ( core_material == 1 ) CoMaterial = MyMaterials::Quartz();
  else if ( core_material == 2 ) CoMaterial = MyMaterials::SiO2();
  else if ( core_material == 3 ) CoMaterial = MyMaterials::SiO2_Ce();
  else if ( core_material == 4 ) CoMaterial = MyMaterials::DSB_Ce();
  else if ( core_material == 5 ) CoMaterial = MyMaterials::LuAG_Ce();
  else if ( core_material == 6 ) CoMaterial = MyMaterials::YAG_Ce();
  else if ( core_material == 7 ) CoMaterial = MyMaterials::LSO();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << core_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Co. material: "<< CoMaterial << G4endl ;
  
  
  CaMaterial = NULL ;
  if      ( capillary_material == 1 ) CaMaterial = MyMaterials::Quartz();
  else if ( capillary_material == 2 ) CaMaterial = MyMaterials::SiO2();
  else if ( capillary_material == 3 ) CaMaterial = MyMaterials::SiO2_Ce();
  else if ( capillary_material == 4 ) CaMaterial = MyMaterials::DSB_Ce(); 
  else if ( capillary_material == 5 ) CaMaterial = MyMaterials::LuAG_Ce();
  else if ( capillary_material == 6 ) CaMaterial = MyMaterials::YAG_Ce();
  else if ( capillary_material == 7 ) CaMaterial = MyMaterials::LSO();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << capillary_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Ca. material: "<< CaMaterial << G4endl ;
  
  
  ClMaterial = NULL ;
  if      ( cladding_material == 1 ) ClMaterial = MyMaterials::Quartz();
  else if ( cladding_material == 2 ) ClMaterial = MyMaterials::SiO2();
  else if ( cladding_material == 3 ) ClMaterial = MyMaterials::SiO2_Ce();
  else if ( cladding_material == 4 ) ClMaterial = MyMaterials::DSB_Ce();
  else if ( cladding_material == 5 ) ClMaterial = MyMaterials::LuAG_Ce();
  else if ( cladding_material == 6 ) ClMaterial = MyMaterials::YAG_Ce();
  else if ( cladding_material == 7 ) ClMaterial = MyMaterials::LSO();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << cladding_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Cl. material: "<< ClMaterial << G4endl ;
  
  
  GaMaterial = NULL;
  if     ( gap_material == 1 ) GaMaterial = MyMaterials::Air();
  else if( gap_material == 2 ) GaMaterial = MyMaterials::OpticalGrease();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
    exit(-1);
  }
  G4cout << "Gap material: " << gap_material << G4endl;
  
  
  DeMaterial = NULL;
  if	 ( det_material == 1 ) DeMaterial = MyMaterials::Silicon();
  else if( det_material == 2 ) DeMaterial = MyMaterials::Quartz();
  else if( det_material == 3 ) DeMaterial = MyMaterials::Air();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
    exit(-1);
  }
  G4cout << "Detector material: " << det_material << G4endl;
  
  
  
  //------------------
  // change properties
  
  if( core_absLength > 0 )
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = { 1.*eV, 10.*eV };
    G4double Absorption[nEntries_ABS] = { core_absLength*mm, core_absLength*mm };
    
    CoMaterial -> GetMaterialPropertiesTable() -> RemoveProperty("ABSLENGTH");
    CoMaterial -> GetMaterialPropertiesTable() -> AddProperty("ABSLENGTH",PhotonEnergy_ABS,Absorption,nEntries_ABS);
  }
  if( core_rIndex > 0 )
  {
    const G4int nEntries_RI = 2;
    G4double PhotonEnergy_RI[nEntries_RI] = { 1.*eV, 10.*eV };
    G4double RefractiveIndex[nEntries_RI] = { core_rIndex, core_rIndex };
    
    CoMaterial -> GetMaterialPropertiesTable() -> RemoveProperty("RINDEX");
    CoMaterial -> GetMaterialPropertiesTable() -> AddProperty("RINDEX",PhotonEnergy_RI,RefractiveIndex,nEntries_RI);
  }
  
  
  if( capillary_absLength > 0 )
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = { 1.*eV, 10.*eV };
    G4double Absorption[nEntries_ABS] = { capillary_absLength*mm, capillary_absLength*mm };
    
    CaMaterial -> GetMaterialPropertiesTable() -> RemoveProperty("ABSLENGTH");
    CaMaterial -> GetMaterialPropertiesTable() -> AddProperty("ABSLENGTH",PhotonEnergy_ABS,Absorption,nEntries_ABS);
  }
  if( capillary_rIndex > 0 )
  {
    const G4int nEntries_RI = 2;
    G4double PhotonEnergy_RI[nEntries_RI] = { 1.*eV, 10.*eV };
    G4double RefractiveIndex[nEntries_RI] = { capillary_rIndex, capillary_rIndex };
    
    CaMaterial -> GetMaterialPropertiesTable() -> RemoveProperty("RINDEX");
    CaMaterial -> GetMaterialPropertiesTable() -> AddProperty("RINDEX",PhotonEnergy_RI,RefractiveIndex,nEntries_RI);
  }
  
  
  if( cladding_absLength > 0 )
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = { 1.*eV, 10.*eV };
    G4double Absorption[nEntries_ABS] = { cladding_absLength*mm, cladding_absLength*mm };
    
    ClMaterial -> GetMaterialPropertiesTable() -> RemoveProperty("ABSLENGTH");
    ClMaterial -> GetMaterialPropertiesTable() -> AddProperty("ABSLENGTH",PhotonEnergy_ABS,Absorption,nEntries_ABS);
  }
  if( cladding_rIndex > 0 )
  {
    const G4int nEntries_RI = 2;
    G4double PhotonEnergy_RI[nEntries_RI] = { 1.*eV, 10.*eV };
    G4double RefractiveIndex[nEntries_RI] = { cladding_rIndex, cladding_rIndex };
    
    ClMaterial -> GetMaterialPropertiesTable() -> RemoveProperty("RINDEX");
    ClMaterial -> GetMaterialPropertiesTable() -> AddProperty("RINDEX",PhotonEnergy_RI,RefractiveIndex,nEntries_RI);
  }
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 



void DetectorConstruction::ConstructField () 
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl ;
  
  static G4TransportationManager * trMgr = G4TransportationManager::GetTransportationManager () ; 
  
  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager * globalFieldMgr = trMgr->GetFieldManager () ;
  
  if( !B_field_IsInitialized )
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector(0.0522*B_field_intensity,0.0522*B_field_intensity,0.9973*B_field_intensity);   
    
    B_field = new G4UniformMagField (fieldVector) ; 
    globalFieldMgr->SetDetectorField (B_field) ;
    globalFieldMgr->CreateChordFinder (B_field) ;
    globalFieldMgr->GetChordFinder ()->SetDeltaChord (0.005 * mm) ;
    B_field_IsInitialized = true ;
  }
  
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl ;
  return ;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
