#include "MyMaterials.hh"
#include "G4NistManager.hh"

using namespace CLHEP;



MyMaterials::MyMaterials()
{}



MyMaterials::~MyMaterials()
{}



G4Material* MyMaterials::Air()
{
  G4double a, z, density;
  G4int nelements;
  
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  
  G4Material* mat = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  mat->AddElement(N, 70.*perCent);
  mat->AddElement(O, 30.*perCent);
  /*
  const G4int nEntries_RI = 42;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV, 
      4.5085*eV, 4.9594*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003, 1.0003, 1.0003, 
      1.0003, 1.0003 };*/

  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.0003,   1.0003,   1.0003,    1.0003,    1.0003,    1.0003,   1.0003,   1.0003};
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX",PhotonEnergy_RI,RefractiveIndex,nEntries_RI);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}



G4Material* MyMaterials::Water()
{
  G4double a, z, density;
  G4int nelements;
  
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  
  G4Material* mat = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  mat->AddElement(H, 2);
  mat->AddElement(O, 1);
  
  const G4int nEntries = 33;
  G4double PhotonEnergy[nEntries] =
    { 0.100*eV, 2.034*eV, 2.068*eV, 2.103*eV,
      2.139*eV, 2.177*eV, 2.216*eV, 2.256*eV,
      2.298*eV, 2.341*eV, 2.386*eV, 2.433*eV,
      2.481*eV, 2.532*eV, 2.585*eV, 2.640*eV,
      2.697*eV, 2.757*eV, 2.820*eV, 2.885*eV,
      2.954*eV, 3.026*eV, 3.102*eV, 3.181*eV,
      3.265*eV, 3.353*eV, 3.446*eV, 3.545*eV,
      3.649*eV, 3.760*eV, 3.877*eV, 4.002*eV,
      4.136*eV };
  G4double RefractiveIndex[nEntries] =
    { 1.3435, 1.3435, 1.3440, 1.3445,
      1.3450, 1.3455, 1.3460, 1.3465,
      1.3470, 1.3475, 1.3480, 1.3485,
      1.3492, 1.3500, 1.3505, 1.3510,
      1.3518, 1.3522, 1.3530, 1.3535,
      1.3540, 1.3545, 1.3550, 1.3555,
      1.3560, 1.3568, 1.3572, 1.3580,
      1.3585, 1.3590, 1.3595, 1.3600,
      1.3608};
  G4double Absorption[nEntries] =
    {  3.448*m,  3.448*m,  4.082*m,  6.329*m,
       9.174*m, 12.346*m, 13.889*m, 15.152*m,
      17.241*m, 18.868*m, 20.000*m, 26.316*m,
      35.714*m, 45.455*m, 47.619*m, 52.632*m,
      52.632*m, 55.556*m, 52.632*m, 52.632*m,
      47.619*m, 45.455*m, 41.667*m, 37.037*m,
      33.333*m, 30.000*m, 28.500*m, 27.000*m,
      24.500*m, 22.000*m, 19.500*m, 17.500*m,
      14.500*m };
  G4double FastComponent[nEntries] =
    { 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00,
      1.00 };
  G4double SlowComponent[nEntries] =
    { 0.01, 0.01, 1.00, 2.00,
      3.00, 4.00, 5.00, 6.00,
      7.00, 8.00, 9.00, 8.00,
      7.00, 6.00, 4.00, 3.00,
      2.00, 1.00, 0.01, 1.00,
      2.00, 3.00, 4.00, 5.00,
      6.00, 7.00, 8.00, 9.00,
      8.00, 7.00, 6.00, 5.00,
      4.00 };
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX",        PhotonEnergy, RefractiveIndex, nEntries);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy, Absorption,      nEntries);
  myMPT->AddProperty("FASTCOMPONENT", PhotonEnergy, FastComponent,   nEntries);
  myMPT->AddProperty("SLOWCOMPONENT", PhotonEnergy, SlowComponent,   nEntries);
  
  myMPT->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  //myMPT->AddConstProperty("ELECTRONSCINTILLATIONYIELD",50./MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT->AddConstProperty("YIELDRATIO",0.8);
  
  mat->SetMaterialPropertiesTable(myMPT);

  return mat;
}



G4Material* MyMaterials::Vacuum()
{
  G4double a, z, density;
  G4int nelements;

  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  
  G4Material* mat = new G4Material("Vacuum", density=0.001*mg/cm3, nelements=2);
  mat->AddElement(N, 70.*perCent);
  mat->AddElement(O, 30.*perCent);
  
  const G4int nEntries = 3;
  G4double PhotonEnergy[nEntries] =
    { 0.0001*eV, 1.00*eV,100.00*eV };
  G4double RefractiveIndex[nEntries] =
    { 1.00, 1.00, 1.00 };
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex, nEntries);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}



G4Material* MyMaterials::Silicon()
{
  G4double a, z, density;
  
  G4Element* Si = new G4Element("Silicon", "Si", z=14., a=28.09*g/mole);
  
  G4Material* mat = new G4Material("Silicon", density=2.33*g/cm3,1);
  mat->AddElement(Si,1);
  /*
  const G4int nEntries = 4;
  G4double PhotonEnergy[nEntries] =
    { 0.0001*eV, 1.0*eV, 1.84*eV, 6.26*eV };
  G4double RefractiveIndex[nEntries] =
    { 4.0, 4.0, 4.0, 4.0 };
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 4.0,   4.0,   4.0,    4.0,    4.0,   4.0,   4.0,   4.0};

  G4double Absorption[nEntries_RI] =
    { 0.1*mm, 0.1*mm, 0.1*mm, 0.1*mm, 0.1*mm, 0.1*mm, 0.1*mm, 0.1*mm};
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  myMPT->AddProperty("ABSLENGTH", PhotonEnergy_RI,  Absorption,      nEntries_RI);
  
  mat->SetMaterialPropertiesTable(myMPT);

  return mat;
}

G4Material* MyMaterials::PlasticO2WLS()	//O2 - orange green-to-wls
{
  G4double a, z, density;
  
  G4Element* H = new G4Element("Hydrogen", "H", z=1., a= 1.01*g/mole);
  G4Element* C = new G4Element("Carbon"  , "C", z=6., a=12.01*g/mole);
  
  G4Material* mat = new G4Material ("PlasticO2WLS", density = 1.*g/cm3,2);
  mat->AddElement(H,4);
  mat->AddElement(C,2);
  
  const G4int nEntries_RI = 11;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0 * eV, 2.0 * eV, 2.5 * eV, 3.0 * eV,
      3.5 * eV, 4.0 * eV, 4.5 * eV, 5.0 * eV,
      5.5 * eV, 6.0 * eV, 6.26 * eV };
  G4double RefractiveIndex[nEntries_RI] =
    { 1.5, 1.5, 1.5, 1.5,
      1.5, 1.5, 1.5, 1.5,
      1.5, 1.5, 1.5 };
  
  const G4int nEntries_ABS = 4;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.0 * eV, 1.84 * eV, 4.08 * eV, 6.26 * eV };
  G4double Absorption[nEntries_ABS] =
    { 1000.*mm, 1000.*mm, 1000.*mm, 1000. *mm };
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  myMPT->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  
  mat->SetMaterialPropertiesTable (myMPT);
  
  return mat;
}



G4Material* MyMaterials::Quartz()
{
  G4double a, z, density;
  
  G4Element* Si = new G4Element("Silicon", "Si", z = 14., a = 28.09* g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O",  z =  8., a = 16.00* g/mole);
  
  G4Material* mat = new G4Material ("Quartz", density = 2.2*g/cm3,2);
  mat->AddElement(Si,1);
  mat->AddElement(O,2);
  /*
  const G4int nEntries_RI = 11;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0 * eV, 2.0 * eV, 2.5 * eV, 3.0 * eV,
      3.5 * eV, 4.0 * eV, 4.5 * eV, 5.0 * eV,
      5.5 * eV, 6.0 * eV, 6.26 * eV };
  G4double RefractiveIndex[nEntries_RI] =
    { 1.53, 1.54, 1.55, 1.56,
      1.56, 1.57, 1.59, 1.60,
      1.62, 1.64, 1.65 };
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.55,   1.55,   1.55,    1.55,    1.55,    1.55,    1.55,   1.55};

  
  const G4int nEntries_ABS = 4;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.0 * eV, 1.84 * eV, 4.08 * eV, 6.26 * eV };
  G4double Absorption[nEntries_ABS] =
    { 138.*mm, 138.*mm, 138.*mm, 138. *mm };
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  myMPT->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  
  mat->SetMaterialPropertiesTable (myMPT);
  
  return mat;
}


G4Material* MyMaterials::ZnO_PS()
{

  G4double a, z, density, fractionmass;

  // ZnO
  G4Element *Zn = new G4Element ("Zinc",  "Zn",  z = 30., a = 65.409 * g / mole);
  G4Element* O  = new G4Element("Oxygen",  "O",  z =  8., a = 16.00* g/mole);

  G4Material *ZnO = new G4Material ("ZnO", density = 5.6 * g / cm3, 2, kStateSolid);
  ZnO->AddElement (Zn, 1);
  ZnO->AddElement (O, 1);

  // PS   
  G4Element* H = new G4Element("Hydrogen", "H", z=1., a= 1.01*g/mole);
  G4Element* C = new G4Element("Carbon"  , "C", z=6., a=12.01*g/mole);

  G4Material *PS = new G4Material ("PS", density = 1.04 * g / cm3, 2);
  PS->AddElement (C, 8);
  PS->AddElement (H, 8);

  G4Material* mat = new G4Material("ZnO_PS", 1.5*g/cm3, 2);
  mat->AddMaterial(ZnO, fractionmass=10*perCent);
  mat->AddMaterial(PS , fractionmass=90*perCent);
/*
  const G4int nEntries_RI = 11;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0 * eV, 2.0 * eV, 2.5 * eV, 3.0 * eV,
      3.5 * eV, 4.0 * eV, 4.5 * eV, 5.0 * eV,
      5.5 * eV, 6.0 * eV, 6.26 * eV };
  G4double RefractiveIndex[nEntries_RI] =
    { 1.55, 1.55, 1.55, 1.55,
      1.55, 1.55, 1.55, 1.55,
      1.55, 1.55, 1.55 };
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.55,   1.55,   1.55,    1.55,    1.55,    1.55,    1.55,   1.55};
  
  const G4int nEntries_ABS = 4;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.0 * eV, 1.84 * eV, 4.08 * eV, 6.26 * eV };
  G4double Absorption[nEntries_ABS] =
    { 1.*mm, 1.*mm, 1.*mm, 1. *mm };
  G4double Rayleigh[nEntries_ABS]       = { 0.1*mm, 0.1*mm, 0.1*mm, 0.1*mm};
  
  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  mt->AddProperty("RAYLEIGH",   PhotonEnergy_ABS, Rayleigh,        nEntries_ABS);

  mat->SetMaterialPropertiesTable(mt);
  return mat;
}

G4Material* MyMaterials::SiO2()
{
  G4double a, z, density;

  G4Element* Si = new G4Element("Silicon", "Si", z = 14., a = 28.09* g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O",  z =  8., a = 16.00* g/mole);

  G4Material* mat = new G4Material ("SiO2Ce", density = 2.65*g/cm3,2);
  mat->AddElement(Si,1);
  mat->AddElement(O,2);
  
  const G4int nEntries_RI = 42;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV, 
      4.5085*eV, 4.9594*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.4476, 1.4476, 1.4478, 1.4481, 
      1.4483, 1.4486, 1.4489, 1.4492, 
      1.4495, 1.4498, 1.4501, 1.4504, 
      1.4507, 1.4511, 1.4514, 1.4518, 
      1.4521, 1.4525, 1.4529, 1.4533, 
      1.4538, 1.4542, 1.4547, 1.4553, 
      1.4559, 1.4565, 1.4572, 1.4580, 
      1.4589, 1.4599, 1.4610, 1.4623, 
      1.4638, 1.4656, 1.4676, 1.4701, 
      1.4731, 1.4769, 1.4816, 1.4878, 
      1.4960, 1.5074 };
  
  const G4int nEntries_ABS = 4;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.0 * eV, 1.84 * eV, 4.08 * eV, 6.26 * eV };
  G4double Absorption[nEntries_ABS] =
    { 138.*mm, 138.*mm, 138.*mm, 138. *mm };
  
  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  
  mat->SetMaterialPropertiesTable(mt);
  return mat;
}

G4Material* MyMaterials::SiO2_Ce()
{
  G4double a, z, density;

  G4Element* Si = new G4Element("Silicon", "Si", z = 14., a = 28.09* g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O",  z =  8., a = 16.00* g/mole);

  G4Material* mat = new G4Material ("SiO2Ce", density = 2.65*g/cm3,2);
  mat->AddElement(Si,1);
  mat->AddElement(O,2);
  
  const G4int nEntries_RI = 42;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV, 
      4.5085*eV, 4.9594*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.4476, 1.4476, 1.4478, 1.4481, 
      1.4483, 1.4486, 1.4489, 1.4492, 
      1.4495, 1.4498, 1.4501, 1.4504, 
      1.4507, 1.4511, 1.4514, 1.4518, 
      1.4521, 1.4525, 1.4529, 1.4533, 
      1.4538, 1.4542, 1.4547, 1.4553, 
      1.4559, 1.4565, 1.4572, 1.4580, 
      1.4589, 1.4599, 1.4610, 1.4623, 
      1.4638, 1.4656, 1.4676, 1.4701, 
      1.4731, 1.4769, 1.4816, 1.4878, 
      1.4960, 1.5074 };
  
  const G4int nEntries_ABS = 4;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.0 * eV, 1.84 * eV, 4.08 * eV, 6.26 * eV };
  G4double Absorption[nEntries_ABS] =
    { 138.*mm, 138.*mm, 138.*mm, 138. *mm };

  const G4int NUMENTRIES_1 = 5;
  G4double FAST_Energy[NUMENTRIES_1]    = {1.8*eV,1.90*eV,2.7*eV,2.88*eV,4.08*eV};
  G4double FAST_COMPONENT[NUMENTRIES_1] = {0.00,1.00,2.0,1.0,0.00};

  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);

  mt->AddConstProperty("SCINTILLATIONYIELD",40000/MeV);
  mt->AddConstProperty("RESOLUTIONSCALE",8.5);
  mt->AddConstProperty("FASTTIMECONSTANT",55.*ns);
  mt->AddConstProperty("YIELDRATIO",1.0);
  mt->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
  
  mat->SetMaterialPropertiesTable(mt);
  return mat;
}

G4Material* MyMaterials::AFO_Ce()
{
  G4double a, z, density;

  G4Element* Si = new G4Element("Silicon", "Si", z = 14., a = 28.09* g/mole);
  G4Element* Ba = new G4Element("Barium",  "Ba", z = 56., a = 137.327* g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O" , z =  8., a = 16.00* g/mole);

  G4Material* mat = new G4Material ("AFO_Ce", density = 4.1*g/cm3,3);
  mat->AddElement(Si,1);
  mat->AddElement(Ba,1);
  mat->AddElement(O,2);
  
  const G4int nEntries_RI = 43;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV, 
      4.5085*eV, 4.9594*eV, 6.4*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.4476, 1.4476, 1.4478, 1.4481, 
      1.4483, 1.4486, 1.4489, 1.4492, 
      1.4495, 1.4498, 1.4501, 1.4504, 
      1.4507, 1.4511, 1.4514, 1.4518, 
      1.4521, 1.4525, 1.4529, 1.4533, 
      1.4538, 1.4542, 1.4547, 1.4553, 
      1.4559, 1.4565, 1.4572, 1.4580, 
      1.4589, 1.4599, 1.4610, 1.4623, 
      1.4638, 1.4656, 1.4676, 1.4701, 
      1.4731, 1.4769, 1.4816, 1.4878, 
      1.4960, 1.5074, 1.51 };
  
  const G4int nEntries_ABS = 143;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    {1.37778*eV, 1.38547*eV, 1.39326*eV, 1.40113*eV, 1.40909*eV, 1.41714*eV, 1.42529*eV, 1.43353*eV, 1.44186*eV, 1.45029*eV, 1.45882*eV, 1.46746*eV, 1.47619*eV, 1.48503*eV, 1.49398*eV, 1.50303*eV, 1.5122*eV, 1.52147*eV, 1.53086*eV, 1.54037*eV, 1.55*eV, 1.55975*eV, 1.56962*eV, 1.57962*eV, 1.58974*eV, 1.6*eV, 1.61039*eV, 1.62092*eV, 1.63158*eV, 1.64238*eV, 1.65333*eV, 1.66443*eV, 1.67568*eV, 1.68707*eV, 1.69863*eV, 1.71034*eV, 1.72222*eV, 1.73427*eV, 1.74648*eV, 1.75887*eV, 1.77143*eV, 1.78417*eV, 1.7971*eV, 1.81022*eV, 1.82353*eV, 1.83704*eV, 1.85075*eV, 1.86466*eV, 1.87879*eV, 1.89313*eV, 1.90769*eV, 1.92248*eV, 1.9375*eV, 1.95276*eV, 1.96825*eV, 1.984*eV, 2*eV, 2.01626*eV, 2.03279*eV, 2.04959*eV, 2.06667*eV, 2.08403*eV, 2.10169*eV, 2.11966*eV, 2.13793*eV, 2.15652*eV, 2.17544*eV, 2.19469*eV, 2.21429*eV, 2.23423*eV, 2.25455*eV, 2.27523*eV, 2.2963*eV, 2.31776*eV, 2.33962*eV, 2.3619*eV, 2.38462*eV, 2.40777*eV, 2.43137*eV, 2.45545*eV, 2.48*eV, 2.50505*eV, 2.53061*eV, 2.5567*eV, 2.58333*eV, 2.61053*eV, 2.6383*eV, 2.66667*eV, 2.69565*eV, 2.72527*eV, 2.75556*eV, 2.78652*eV, 2.81818*eV, 2.85057*eV, 2.88372*eV, 2.91765*eV, 2.95238*eV, 2.98795*eV, 3.02439*eV, 3.06173*eV, 3.1*eV, 3.13924*eV, 3.17949*eV, 3.22078*eV, 3.26316*eV, 3.30667*eV, 3.35135*eV, 3.39726*eV, 3.44444*eV, 3.49296*eV, 3.54286*eV, 3.5942*eV, 3.64706*eV, 3.70149*eV, 3.75758*eV, 3.81538*eV, 3.875*eV, 3.93651*eV, 4*eV, 4.06557*eV, 4.14716*eV, 4.21769*eV, 4.29066*eV, 4.3662*eV, 4.44444*eV, 4.52555*eV, 4.60967*eV, 4.69697*eV, 4.78764*eV, 4.88189*eV, 4.97992*eV, 5.08197*eV, 5.18828*eV, 5.29915*eV, 5.41485*eV, 5.53571*eV, 5.6621*eV, 5.8216*eV, 5.96154*eV, 6.10837*eV, 6.26263*eV, 6.42487*eV};

  G4double Absorption[nEntries_ABS] =
    {799.821*mm, 799.821*mm, 799.821*mm, 799.821*mm, 799.821*mm, 799.821*mm, 799.821*mm, 799.821*mm, 799.821*mm, 799.821*mm, 672.103*mm, 648.835*mm, 799.821*mm, 799.821*mm, 629.971*mm, 844.864*mm, 613.74*mm, 905.347*mm, 799.821*mm, 799.821*mm, 999.877*mm, 648.86*mm, 787.055*mm, 669.728*mm, 712.126*mm, 634.441*mm, 755.491*mm, 619.314*mm, 704.613*mm, 765.401*mm, 947.116*mm, 652.787*mm, 862.252*mm, 616.461*mm, 570.426*mm, 643.712*mm, 688.632*mm, 711.359*mm, 654.807*mm, 809.302*mm, 625.79*mm, 593.457*mm, 764.57*mm, 791.575*mm, 658.785*mm, 677.396*mm, 611.824*mm, 616.835*mm, 687.348*mm, 594.259*mm, 641.453*mm, 630.082*mm, 659.429*mm, 565.723*mm, 785.166*mm, 594.032*mm, 766.44*mm, 553.251*mm, 748.675*mm, 578.02*mm, 721.11*mm, 723.68*mm, 614.537*mm, 665.457*mm, 595.855*mm, 586.003*mm, 547.263*mm, 580.18*mm, 606.99*mm, 569.882*mm, 541.281*mm, 538.117*mm, 479.382*mm, 519.495*mm, 520.828*mm, 530.076*mm, 508.279*mm, 433.86*mm, 470.744*mm, 476.178*mm, 463.675*mm, 423.752*mm, 418.55*mm, 479.866*mm, 431.962*mm, 440.813*mm, 429.812*mm, 435.138*mm, 399.863*mm, 390.322*mm, 408.896*mm, 421.817*mm, 358.89*mm, 391.001*mm, 377.763*mm, 361.11*mm, 356.551*mm, 365.89*mm, 361.076*mm, 335.163*mm, 340.372*mm, 324.448*mm, 338.888*mm, 327.038*mm, 318.322*mm, 285.459*mm, 271.714*mm, 265.171*mm, 241.793*mm, 181.132*mm, 103.204*mm, 43.6429*mm, 16.3309*mm, 6.13723*mm, 2.39139*mm, 0.5*mm, 0.5*mm, 0.5*mm, 1.4214*mm, 1.33911*mm, 0.5*mm, 1.46426*mm, 1.46522*mm, 0.5*mm, 0.5*mm, 1.28414*mm, 1.30373*mm, 1.28356*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 1.3184*mm, 0.5*mm, 0.5*mm, 1.42229*mm, 0.5*mm, 0.5*mm, 1.08768*mm, 0.5*mm};
  
  //emission spectrum
  const G4int NUMENTRIES_1 = 93;
  G4double FAST_Energy[NUMENTRIES_1]    = 
  {3.875*eV, 3.84496*eV, 3.81538*eV, 3.78626*eV, 3.75758*eV, 3.72932*eV, 3.70149*eV, 3.67407*eV, 3.64706*eV, 3.62044*eV, 3.5942*eV, 3.56835*eV, 3.54286*eV, 3.51773*eV, 3.49296*eV, 3.46853*eV, 3.44444*eV, 3.42069*eV, 3.39726*eV, 3.37415*eV, 3.35135*eV, 3.32886*eV, 3.30667*eV, 3.28477*eV, 3.26316*eV, 3.24183*eV, 3.22078*eV, 3.2*eV, 3.17949*eV, 3.15924*eV, 3.13924*eV, 3.1195*eV, 3.1*eV, 3.08075*eV, 3.06173*eV, 3.04294*eV, 3.02439*eV, 3.00606*eV, 2.98795*eV, 2.97006*eV, 2.95238*eV, 2.93491*eV, 2.91765*eV, 2.90058*eV, 2.88372*eV, 2.86705*eV, 2.85057*eV, 2.83429*eV, 2.81818*eV, 2.80226*eV, 2.78652*eV, 2.77095*eV, 2.75556*eV, 2.74033*eV, 2.72527*eV, 2.71038*eV, 2.69565*eV, 2.68108*eV, 2.66667*eV, 2.65241*eV, 2.6383*eV, 2.62434*eV, 2.61053*eV, 2.59686*eV, 2.58333*eV, 2.56995*eV, 2.5567*eV, 2.54359*eV, 2.53061*eV, 2.51777*eV, 2.50505*eV, 2.49246*eV, 2.48*eV, 2.46766*eV, 2.45545*eV, 2.44335*eV, 2.43137*eV, 2.41951*eV, 2.40777*eV, 2.39614*eV, 2.38462*eV, 2.37321*eV, 2.3619*eV, 2.35071*eV, 2.33962*eV, 2.32864*eV, 2.31776*eV, 2.30698*eV, 2.2963*eV, 2.28571*eV, 2.27523*eV, 2.26484*eV, 2.25455*eV};
    
  G4double FAST_COMPONENT[NUMENTRIES_1] = 
  {0.00208044, 0.00333397, 0.00752478, 0.0234724, 0.0532366, 0.0854065, 0.114836, 0.163069, 0.234962, 0.318021, 0.405736, 0.500799, 0.589805, 0.668213, 0.732568, 0.788987, 0.834531, 0.864747, 0.89069, 0.903711, 0.907175, 0.899195, 0.882101, 0.855825, 0.823148, 0.782556, 0.737222, 0.689739, 0.638437, 0.586504, 0.536268, 0.4859, 0.437051, 0.39146, 0.348727, 0.312368, 0.278995, 0.247668, 0.220126, 0.191708, 0.168009, 0.145147, 0.126585, 0.109526, 0.0933741, 0.0803073, 0.06865, 0.0585006, 0.0500656, 0.0424515, 0.0359967, 0.0306027, 0.0261705, 0.0223074, 0.0186536, 0.0159328, 0.013442, 0.0116061, 0.00991835, 0.00854082, 0.00764762, 0.00673076, 0.00568263, 0.0052557, 0.00473967, 0.00425339, 0.00396469, 0.00349449, 0.00340015, 0.00295026, 0.00287799, 0.00269137, 0.00239846, 0.00223671, 0.00231405, 0.00205452, 0.00210516, 0.00198273, 0.00191159, 0.00195633, 0.00204159, 0.00186753, 0.00196798, 0.00194337, 0.00188266, 0.00193113, 0.00180067, 0.00173014, 0.00172164, 0.00176757, 0.00155514, 0.00175762, 0.00174964};

  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);

  mt->AddConstProperty("SCINTILLATIONYIELD", 400/MeV);
  mt->AddConstProperty("RESOLUTIONSCALE",1);
  mt->AddConstProperty("FASTTIMECONSTANT", 50.*ns);
//  mt->AddConstProperty("YIELDRATIO",1.0);
  mt->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
    
  mat->SetMaterialPropertiesTable(mt);
  return mat;
}

G4Material* MyMaterials::AFO_undoped()
{
  G4double a, z, density;

  G4Element* Si = new G4Element("Silicon", "Si", z = 14., a = 28.09* g/mole);
  G4Element* Ba = new G4Element("Barium",  "Ba", z = 56., a = 137.327* g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O" , z =  8., a = 16.00* g/mole);

  G4Material* mat = new G4Material ("AFO_undoped", density = 4.1*g/cm3,3);
  mat->AddElement(Si,1);
  mat->AddElement(Ba,1);
  mat->AddElement(O,2);
  
  const G4int nEntries_RI = 43;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV, 
      4.5085*eV, 4.9594*eV, 6.4*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.4476, 1.4476, 1.4478, 1.4481, 
      1.4483, 1.4486, 1.4489, 1.4492, 
      1.4495, 1.4498, 1.4501, 1.4504, 
      1.4507, 1.4511, 1.4514, 1.4518, 
      1.4521, 1.4525, 1.4529, 1.4533, 
      1.4538, 1.4542, 1.4547, 1.4553, 
      1.4559, 1.4565, 1.4572, 1.4580, 
      1.4589, 1.4599, 1.4610, 1.4623, 
      1.4638, 1.4656, 1.4676, 1.4701, 
      1.4731, 1.4769, 1.4816, 1.4878, 
      1.4960, 1.5074, 1.51 };
  
  const G4int nEntries_ABS = 143;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    {1.37778*eV, 1.38547*eV, 1.39326*eV, 1.40113*eV, 1.40909*eV, 1.41714*eV, 1.42529*eV, 1.43353*eV, 1.44186*eV, 1.45029*eV, 1.45882*eV, 1.46746*eV, 1.47619*eV, 1.48503*eV, 1.49398*eV, 1.50303*eV, 1.5122*eV, 1.52147*eV, 1.53086*eV, 1.54037*eV, 1.55*eV, 1.55975*eV, 1.56962*eV, 1.57962*eV, 1.58974*eV, 1.6*eV, 1.61039*eV, 1.62092*eV, 1.63158*eV, 1.64238*eV, 1.65333*eV, 1.66443*eV, 1.67568*eV, 1.68707*eV, 1.69863*eV, 1.71034*eV, 1.72222*eV, 1.73427*eV, 1.74648*eV, 1.75887*eV, 1.77143*eV, 1.78417*eV, 1.7971*eV, 1.81022*eV, 1.82353*eV, 1.83704*eV, 1.85075*eV, 1.86466*eV, 1.87879*eV, 1.89313*eV, 1.90769*eV, 1.92248*eV, 1.9375*eV, 1.95276*eV, 1.96825*eV, 1.984*eV, 2*eV, 2.01626*eV, 2.03279*eV, 2.04959*eV, 2.06667*eV, 2.08403*eV, 2.10169*eV, 2.11966*eV, 2.13793*eV, 2.15652*eV, 2.17544*eV, 2.19469*eV, 2.21429*eV, 2.23423*eV, 2.25455*eV, 2.27523*eV, 2.2963*eV, 2.31776*eV, 2.33962*eV, 2.3619*eV, 2.38462*eV, 2.40777*eV, 2.43137*eV, 2.45545*eV, 2.48*eV, 2.50505*eV, 2.53061*eV, 2.5567*eV, 2.58333*eV, 2.61053*eV, 2.6383*eV, 2.66667*eV, 2.69565*eV, 2.72527*eV, 2.75556*eV, 2.78652*eV, 2.81818*eV, 2.85057*eV, 2.88372*eV, 2.91765*eV, 2.95238*eV, 2.98795*eV, 3.02439*eV, 3.06173*eV, 3.1*eV, 3.13924*eV, 3.17949*eV, 3.22078*eV, 3.26316*eV, 3.30667*eV, 3.35135*eV, 3.39726*eV, 3.44444*eV, 3.49296*eV, 3.54286*eV, 3.5942*eV, 3.64706*eV, 3.70149*eV, 3.75758*eV, 3.81538*eV, 3.875*eV, 3.93651*eV, 4*eV, 4.06557*eV, 4.13333*eV, 4.20339*eV, 4.27586*eV, 4.35088*eV, 4.42857*eV, 4.50909*eV, 4.59259*eV, 4.67925*eV, 4.76923*eV, 4.86275*eV, 4.96*eV, 5.06122*eV, 5.16667*eV, 5.2766*eV, 5.3913*eV, 5.51111*eV, 5.63636*eV, 5.76744*eV, 5.90476*eV, 6.04878*eV, 6.2*eV, 6.35897*eV, 6.52632*eV};

  G4double Absorption[nEntries_ABS] =
    {727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 799.821*mm, 799.821*mm, 863.016*mm, 7137.52*mm, 752.307*mm, 645.178*mm, 788.825*mm, 554.496*mm, 799.821*mm, 938.506*mm, 799.821*mm, 463.043*mm, 799.821*mm, 595.76*mm, 672.641*mm, 727.29*mm, 732.717*mm, 637.631*mm, 542.757*mm, 730.837*mm, 799.821*mm, 642.004*mm, 732.557*mm, 572.838*mm, 668.199*mm, 633.246*mm, 579.953*mm, 724.578*mm, 679.925*mm, 809.302*mm, 697.23*mm, 574.274*mm, 667.979*mm, 996.469*mm, 778.907*mm, 589.811*mm, 682.036*mm, 614.903*mm, 612.425*mm, 558.8*mm, 600.388*mm, 573.264*mm, 682.24*mm, 620.87*mm, 807.574*mm, 513.77*mm, 810.301*mm, 605.708*mm, 663.732*mm, 569.737*mm, 759.523*mm, 618.363*mm, 681.626*mm, 596.102*mm, 633.495*mm, 659.046*mm, 541.37*mm, 586.942*mm, 596.202*mm, 574.758*mm, 545.677*mm, 518.743*mm, 528.427*mm, 527.69*mm, 520.828*mm, 551.011*mm, 526.645*mm, 477.234*mm, 483.831*mm, 494.08*mm, 467.447*mm, 426.827*mm, 421.471*mm, 463.217*mm, 466.91*mm, 425.735*mm, 429.812*mm, 422.051*mm, 413.414*mm, 379.887*mm, 389.232*mm, 383.074*mm, 375.402*mm, 349.485*mm, 373.998*mm, 360.416*mm, 396.534*mm, 336.949*mm, 353.571*mm, 327.523*mm, 322.734*mm, 288.729*mm, 334.029*mm, 288.477*mm, 302.728*mm, 280.321*mm, 283.31*mm, 268.642*mm, 254.899*mm, 294.522*mm, 300.912*mm, 281.016*mm, 236.868*mm, 232.818*mm, 221.732*mm, 246.612*mm, 168.529*mm, 151.129*mm, 129.74*mm, 108.552*mm, 90.7699*mm, 77.7944*mm, 68.669*mm, 58.5175*mm, 51.2097*mm, 43.4176*mm, 37.2326*mm, 31.6875*mm, 27.0823*mm, 23.8538*mm, 20.9881*mm, 19.217*mm, 17.7166*mm, 16.6901*mm, 15.8046*mm, 15.1983*mm, 14.825*mm, 14.2654*mm, 13.2542*mm, 11.9379*mm, 10.1158*mm, 8.19016*mm, 6.31106*mm};

  //emission spectrum
  const G4int NUMENTRIES_1 = 91;
  G4double FAST_Energy[NUMENTRIES_1]    = 
  {4.27586*eV, 4.23208*eV, 4.19628*eV, 4.16107*eV, 4.12646*eV, 4.09241*eV, 4.05892*eV, 4.02597*eV, 3.99356*eV, 3.96166*eV, 3.93027*eV, 3.89937*eV, 3.86895*eV, 3.83901*eV, 3.80952*eV, 3.78049*eV, 3.75189*eV, 3.72372*eV, 3.69598*eV, 3.66864*eV, 3.6417*eV, 3.61516*eV, 3.589*eV, 3.56322*eV, 3.5378*eV, 3.51275*eV, 3.48805*eV, 3.46369*eV, 3.43967*eV, 3.41598*eV, 3.39261*eV, 3.36957*eV, 3.34683*eV, 3.3244*eV, 3.30226*eV, 3.28042*eV, 3.25887*eV, 3.2376*eV, 3.2166*eV, 3.19588*eV, 3.17542*eV, 3.15522*eV, 3.13527*eV, 3.11558*eV, 3.09613*eV, 3.07692*eV, 3.05795*eV, 3.03922*eV, 3.02071*eV, 3.00242*eV, 2.98436*eV, 2.96651*eV, 2.94887*eV, 2.93144*eV, 2.91422*eV, 2.8972*eV, 2.88037*eV, 2.86374*eV, 2.8473*eV, 2.83105*eV, 2.81498*eV, 2.7991*eV, 2.78339*eV, 2.76786*eV, 2.7525*eV, 2.73731*eV, 2.72228*eV, 2.70742*eV, 2.69273*eV, 2.67819*eV, 2.6638*eV, 2.64957*eV, 2.63549*eV, 2.62156*eV, 2.60778*eV, 2.59414*eV, 2.58065*eV, 2.56729*eV, 2.55407*eV, 2.54098*eV, 2.52803*eV, 2.51521*eV, 2.50252*eV, 2.48996*eV, 2.47752*eV, 2.46521*eV, 2.45302*eV, 2.44094*eV, 2.42899*eV, 2.41715*eV, 2.40543*eV};
    
  G4double FAST_COMPONENT[NUMENTRIES_1] = 
  {0.00222164, 0.0947739, 0.133648, 0.18667, 0.249395, 0.31931, 0.391404, 0.462693, 0.534068, 0.597178, 0.661313, 0.718209, 0.767109, 0.815039, 0.855814, 0.882384, 0.901889, 0.909292, 0.903367, 0.890746, 0.869911, 0.835148, 0.795036, 0.756399, 0.705759, 0.657028, 0.605187, 0.552639, 0.506292, 0.457564, 0.410797, 0.370532, 0.326186, 0.286707, 0.252722, 0.2199, 0.195966, 0.168503, 0.146915, 0.126209, 0.111001, 0.0974553, 0.0850033, 0.0742505, 0.0655606, 0.0567462, 0.0516565, 0.0465925, 0.0418254, 0.0373131, 0.0345874, 0.0315096, 0.0295812, 0.0278235, 0.0267798, 0.0254359, 0.0242347, 0.0225283, 0.021659, 0.0202679, 0.0192306, 0.0184972, 0.0183235, 0.016838, 0.01585, 0.0161146, 0.0150103, 0.0144915, 0.0145025, 0.0138657, 0.0134015, 0.0125441, 0.012509, 0.0122657, 0.0121878, 0.0118191, 0.0116185, 0.0109438, 0.0105743, 0.0103989, 0.0098789, 0.0100983, 0.00930592, 0.00923501, 0.009329, 0.00916632, 0.00919149, 0.00894258, 0.00847391, 0.00884528, 0.00826064};

  
  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);
  mt->AddConstProperty("SCINTILLATIONYIELD", 100/MeV);
  mt->AddConstProperty("RESOLUTIONSCALE",1);
  mt->AddConstProperty("FASTTIMECONSTANT", 5.*ns);
//  mt->AddConstProperty("YIELDRATIO",1.0);
  mt->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
  
  mat->SetMaterialPropertiesTable(mt);
  return mat;
}

G4Material* MyMaterials::PbF2()
{
  G4double a, z, density;

  G4Element* Pb = new G4Element("Lead", "Pb", z = 82., a = 207.21* g/mole);
  G4Element* F = new G4Element("Fluorine",  "F", z = 9., a = 18.9984* g/mole);

  G4Material* mat = new G4Material ("PbF2", density = 7.77*g/cm3,2);
  mat->AddElement(Pb,1);
  mat->AddElement(F,2);
  
  const G4int nEntries_RI = 13;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    {0.135*eV, 0.177*eV, 0.248*eV, 0.413*eV, 1.24*eV, 1.377*eV, 1.55*eV, 1.77*eV, 2.066*eV, 2.48*eV, 3.1*eV, 4.*eV, 4.9*eV};
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.655, 1.685, 1.708, 1.724, 1.742, 1.745, 1.749, 1.755, 1.765, 1.782, 1.818, 1.89, 1.90};
  
  const G4int nEntries_ABS = 120;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    {1.37778*eV, 1.38547*eV, 1.39326*eV, 1.40113*eV, 1.40909*eV, 1.41714*eV, 1.42529*eV, 1.43353*eV, 1.44186*eV, 1.45029*eV, 1.45882*eV, 1.46746*eV, 1.47619*eV, 1.48503*eV, 1.49398*eV, 1.50303*eV, 1.5122*eV, 1.52147*eV, 1.53086*eV, 1.54037*eV, 1.55*eV, 1.55975*eV, 1.56962*eV, 1.57962*eV, 1.58974*eV, 1.6*eV, 1.61039*eV, 1.62092*eV, 1.63158*eV, 1.64238*eV, 1.65333*eV, 1.66443*eV, 1.67568*eV, 1.68707*eV, 1.69863*eV, 1.71034*eV, 1.72222*eV, 1.73427*eV, 1.74648*eV, 1.75887*eV, 1.77143*eV, 1.78417*eV, 1.7971*eV, 1.81022*eV, 1.82353*eV, 1.83704*eV, 1.85075*eV, 1.86466*eV, 1.87879*eV, 1.89313*eV, 1.90769*eV, 1.92248*eV, 1.9375*eV, 1.95276*eV, 1.96825*eV, 1.984*eV, 2*eV, 2.01626*eV, 2.03279*eV, 2.04959*eV, 2.06667*eV, 2.08403*eV, 2.10169*eV, 2.11966*eV, 2.13793*eV, 2.15652*eV, 2.17544*eV, 2.19469*eV, 2.21429*eV, 2.23423*eV, 2.25455*eV, 2.27523*eV, 2.2963*eV, 2.31776*eV, 2.33962*eV, 2.3619*eV, 2.38462*eV, 2.40777*eV, 2.43137*eV, 2.45545*eV, 2.48*eV, 2.50505*eV, 2.53061*eV, 2.5567*eV, 2.58333*eV, 2.61053*eV, 2.6383*eV, 2.66667*eV, 2.69565*eV, 2.72527*eV, 2.75556*eV, 2.78652*eV, 2.81818*eV, 2.85057*eV, 2.88372*eV, 2.91765*eV, 2.95238*eV, 2.98795*eV, 3.02439*eV, 3.06173*eV, 3.1*eV, 3.13924*eV, 3.17949*eV, 3.22078*eV, 3.26316*eV, 3.30667*eV, 3.35135*eV, 3.39726*eV, 3.44444*eV, 3.49296*eV, 3.54286*eV, 3.5942*eV, 3.64706*eV, 3.70149*eV, 3.75758*eV, 3.81538*eV, 3.875*eV, 3.93651*eV, 4*eV, 4.06557*eV};

  G4double Absorption[nEntries_ABS] =
    {727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 727.29*mm, 799.821*mm, 799.821*mm, 863.016*mm, 7137.52*mm, 752.307*mm, 645.178*mm, 788.825*mm, 554.496*mm, 799.821*mm, 938.506*mm, 799.821*mm, 463.043*mm, 799.821*mm, 595.76*mm, 672.641*mm, 727.29*mm, 732.717*mm, 637.631*mm, 542.757*mm, 730.837*mm, 799.821*mm, 642.004*mm, 732.557*mm, 572.838*mm, 668.199*mm, 633.246*mm, 579.953*mm, 724.578*mm, 679.925*mm, 809.302*mm, 697.23*mm, 574.274*mm, 667.979*mm, 996.469*mm, 778.907*mm, 589.811*mm, 682.036*mm, 614.903*mm, 612.425*mm, 558.8*mm, 600.388*mm, 573.264*mm, 682.24*mm, 620.87*mm, 607.574*mm, 513.77*mm, 610.301*mm, 605.708*mm, 663.732*mm, 569.737*mm, 759.523*mm, 618.363*mm, 681.626*mm, 596.102*mm, 633.495*mm, 659.046*mm, 541.37*mm, 586.942*mm, 596.202*mm, 574.758*mm, 545.677*mm, 518.743*mm, 528.427*mm, 527.69*mm, 520.828*mm, 551.011*mm, 526.645*mm, 477.234*mm, 483.831*mm, 494.08*mm, 467.447*mm, 426.827*mm, 421.471*mm, 463.217*mm, 466.91*mm, 425.735*mm, 429.812*mm, 422.051*mm, 413.414*mm, 379.887*mm, 389.232*mm, 383.074*mm, 375.402*mm, 349.485*mm, 373.998*mm, 360.416*mm, 396.534*mm, 336.949*mm, 353.571*mm, 327.523*mm, 322.734*mm, 288.729*mm, 334.029*mm, 288.477*mm, 200*mm, 100*mm, 90*mm, 80*mm, 70*mm, 60*mm, 50*mm, 40*mm, 30*mm, 20*mm, 10*mm, 5*mm, 1*mm, 0.*mm, 0*mm, 0.*mm};

/*
  //emission spectrum
  const G4int NUMENTRIES_1 = 91;
  G4double FAST_Energy[NUMENTRIES_1]    = 
  {4.27586*eV, 4.23208*eV, 4.19628*eV, 4.16107*eV, 4.12646*eV, 4.09241*eV, 4.05892*eV, 4.02597*eV, 3.99356*eV, 3.96166*eV, 3.93027*eV, 3.89937*eV, 3.86895*eV, 3.83901*eV, 3.80952*eV, 3.78049*eV, 3.75189*eV, 3.72372*eV, 3.69598*eV, 3.66864*eV, 3.6417*eV, 3.61516*eV, 3.589*eV, 3.56322*eV, 3.5378*eV, 3.51275*eV, 3.48805*eV, 3.46369*eV, 3.43967*eV, 3.41598*eV, 3.39261*eV, 3.36957*eV, 3.34683*eV, 3.3244*eV, 3.30226*eV, 3.28042*eV, 3.25887*eV, 3.2376*eV, 3.2166*eV, 3.19588*eV, 3.17542*eV, 3.15522*eV, 3.13527*eV, 3.11558*eV, 3.09613*eV, 3.07692*eV, 3.05795*eV, 3.03922*eV, 3.02071*eV, 3.00242*eV, 2.98436*eV, 2.96651*eV, 2.94887*eV, 2.93144*eV, 2.91422*eV, 2.8972*eV, 2.88037*eV, 2.86374*eV, 2.8473*eV, 2.83105*eV, 2.81498*eV, 2.7991*eV, 2.78339*eV, 2.76786*eV, 2.7525*eV, 2.73731*eV, 2.72228*eV, 2.70742*eV, 2.69273*eV, 2.67819*eV, 2.6638*eV, 2.64957*eV, 2.63549*eV, 2.62156*eV, 2.60778*eV, 2.59414*eV, 2.58065*eV, 2.56729*eV, 2.55407*eV, 2.54098*eV, 2.52803*eV, 2.51521*eV, 2.50252*eV, 2.48996*eV, 2.47752*eV, 2.46521*eV, 2.45302*eV, 2.44094*eV, 2.42899*eV, 2.41715*eV, 2.40543*eV};
    
  G4double FAST_COMPONENT[NUMENTRIES_1] = 
  {0.00222164, 0.0947739, 0.133648, 0.18667, 0.249395, 0.31931, 0.391404, 0.462693, 0.534068, 0.597178, 0.661313, 0.718209, 0.767109, 0.815039, 0.855814, 0.882384, 0.901889, 0.909292, 0.903367, 0.890746, 0.869911, 0.835148, 0.795036, 0.756399, 0.705759, 0.657028, 0.605187, 0.552639, 0.506292, 0.457564, 0.410797, 0.370532, 0.326186, 0.286707, 0.252722, 0.2199, 0.195966, 0.168503, 0.146915, 0.126209, 0.111001, 0.0974553, 0.0850033, 0.0742505, 0.0655606, 0.0567462, 0.0516565, 0.0465925, 0.0418254, 0.0373131, 0.0345874, 0.0315096, 0.0295812, 0.0278235, 0.0267798, 0.0254359, 0.0242347, 0.0225283, 0.021659, 0.0202679, 0.0192306, 0.0184972, 0.0183235, 0.016838, 0.01585, 0.0161146, 0.0150103, 0.0144915, 0.0145025, 0.0138657, 0.0134015, 0.0125441, 0.012509, 0.0122657, 0.0121878, 0.0118191, 0.0116185, 0.0109438, 0.0105743, 0.0103989, 0.0098789, 0.0100983, 0.00930592, 0.00923501, 0.009329, 0.00916632, 0.00919149, 0.00894258, 0.00847391, 0.00884528, 0.00826064};*/

  
  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
//  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);
//  mt->AddConstProperty("SCINTILLATIONYIELD", 100/MeV);
//  mt->AddConstProperty("RESOLUTIONSCALE",1);
//  mt->AddConstProperty("FASTTIMECONSTANT", 5.*ns);
//  mt->AddConstProperty("YIELDRATIO",1.0);
//  mt->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
  
  mat->SetMaterialPropertiesTable(mt);
  return mat;
}



G4Material* MyMaterials::YAG_Ce()
{
  G4double a, z, density;
  
  G4Element* Y  = new G4Element("Silicon",   "Y",  z = 39., a = 88.01* g/mole);
  G4Element* Al = new G4Element("Aluminium", "Al", z = 13., a = 28.09* g/mole);
  G4Element* O  = new G4Element("Oxygen",    "O",  z =  8., a = 16.00* g/mole);
  
  G4Material* mat = new G4Material ("YAG_Ce", density = 4.6*g/cm3,3);
  mat->AddElement(Y,3);
  mat->AddElement(Al,5);
  mat->AddElement(O,12);
  /*
  const G4int nEntries_RI = 40;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV}; //{,      4.5085*eV, 4.9594*eV, 6.5000*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.8115, 1.8116, 1.8118, 1.8122, 
      1.8126, 1.8131, 1.8135, 1.8140, 
      1.8144, 1.8149, 1.8154, 1.8160, 
      1.8165, 1.8171, 1.8177, 1.8184, 
      1.8191, 1.8198, 1.8206, 1.8214, 
      1.8223, 1.8233, 1.8244, 1.8256, 
      1.8269, 1.8284, 1.8300, 1.8318, 
      1.8338, 1.8362, 1.8388, 1.8419, 
      1.8455, 1.8497, 1.8548, 1.8608, 
      1.8683, 1.8775, 1.8892, 1.9045}; //, 1.9249, 1.9532, 1.9600 };
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.84,   1.84,   1.84,    1.84,    1.84,    1.84,    1.84,   1.84};

  //G4double Rayleigh[nEntries_RI] =
  //  { 138.*mm, 138.*mm, 138.*mm };
  
   //intrinsic absorption spectrum
  const G4int nEntries_ABS = 84;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    {4.42857*eV, 4.35088*eV, 4.27586*eV, 4.20339*eV, 4.13333*eV, 4.06557*eV, 4*eV, 3.93651*eV, 3.875*eV, 
    3.81538*eV, 3.75758*eV, 3.70149*eV, 3.64706*eV, 3.5942*eV, 3.54286*eV, 3.49296*eV, 3.44444*eV, 
    3.39726*eV, 3.35135*eV, 3.30667*eV, 3.26316*eV, 3.22078*eV, 3.17949*eV, 3.13924*eV, 3.1*eV, 3.06173*eV, 
    3.02439*eV, 2.98795*eV, 2.95238*eV, 2.91765*eV, 2.88372*eV, 2.85057*eV, 2.81818*eV, 2.78652*eV, 
    2.75556*eV, 2.72527*eV, 2.69565*eV, 2.66667*eV, 2.6383*eV, 2.61053*eV, 2.58333*eV, 2.5567*eV, 
    2.53061*eV, 2.50505*eV, 2.48*eV, 2.45545*eV, 2.43137*eV, 2.40777*eV, 2.38462*eV, 2.3619*eV, 
    2.33962*eV, 2.31776*eV, 2.2963*eV, 2.27523*eV, 2.25455*eV, 2.23423*eV, 2.21429*eV, 2.19469*eV, 
    2.17544*eV, 2.15652*eV, 2.13793*eV, 2.11966*eV, 2.10169*eV, 2.08403*eV, 2.06667*eV, 2.04959*eV, 
    2.03279*eV, 2.01626*eV, 2*eV, 1.984*eV, 1.96825*eV, 1.95276*eV, 1.9375*eV, 1.92248*eV, 1.90769*eV, 
    1.89313*eV, 1.87879*eV, 1.86466*eV, 1.85075*eV, 1.83704*eV, 1.82353*eV, 1.81022*eV, 1.7971*eV, 1.78417*eV};
    
  G4double Absorption[nEntries_ABS] =
    {12.7389*mm, 7.96854*mm, 9.79472*mm, 11.3642*mm, 13.9544*mm, 14.078*mm, 9.48293*mm, 4.60862*mm, 
    0.45205*mm, 0.14906*mm, 0.09215*mm, 0.05898*mm, 0.09349*mm, 0.09448*mm, 4.00614*mm, 14.055*mm, 
    59.4712*mm, 177.765*mm, 322.779*mm, 407.883*mm, 423.384*mm, 250.626*mm, 107.986*mm, 40.2593*mm, 
    16.1867*mm, 7.08721*mm, 3.3827*mm, 1.85802*mm, 0.42111*mm, 0.348*mm, 0.36044*mm, 0.37582*mm, 
    0.39199*mm, 0.33268*mm, 0.27998*mm, 0.33348*mm, 0.21815*mm, 0.28273*mm, 0.37124*mm, 0.4064*mm, 
    1.09597*mm, 1.45447*mm, 2.03141*mm, 3.62354*mm, 6.92578*mm, 13.9067*mm, 29.0131*mm, 60.2374*mm, 
    129.476*mm, 226.564*mm, 553.51*mm, 516.632*mm, 600.*mm, 700.*mm, 800.*mm, 900.*mm, 
    1000.*mm, 1000.*mm, 1000.*mm, 1000.2*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 
    1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 
    1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 
    1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm};
  
  //emission spectrum
  const G4int NUMENTRIES_1 = 137;
  G4double FAST_Energy[NUMENTRIES_1]    = 
  {2.72827*eV, 2.71335*eV, 2.69859*eV, 2.68398*eV, 2.66954*eV, 2.65525*eV, 2.64111*eV, 2.62712*eV, 2.61328*eV, 2.59958*eV, 
    2.58603*eV, 2.57261*eV, 2.55934*eV, 2.5462*eV, 2.5332*eV, 2.52033*eV, 2.50758*eV, 2.49497*eV, 2.48248*eV, 2.47012*eV, 
    2.45788*eV, 2.44576*eV, 2.43376*eV, 2.42188*eV, 2.41011*eV, 2.39845*eV, 2.38691*eV, 2.37548*eV, 2.36416*eV, 2.35294*eV, 
    2.34183*eV, 2.33083*eV, 2.31993*eV, 2.30912*eV, 2.29842*eV, 2.28782*eV, 2.27732*eV, 2.26691*eV, 2.2566*eV, 2.24638*eV, 
    2.23625*eV, 2.22621*eV, 2.21626*eV, 2.20641*eV, 2.19663*eV, 2.18695*eV, 2.17735*eV, 2.16783*eV, 2.1584*eV, 2.14905*eV, 
    2.13978*eV, 2.13058*eV, 2.12147*eV, 2.11244*eV, 2.10348*eV, 2.09459*eV, 2.08579*eV, 2.07705*eV, 2.06839*eV, 2.0598*eV, 
    2.05128*eV, 2.04283*eV, 2.03445*eV, 2.02614*eV, 2.0179*eV, 2.00972*eV, 2.00161*eV, 1.99357*eV, 1.98559*eV, 1.97767*eV, 
    1.96982*eV, 1.96203*eV, 1.95429*eV, 1.94662*eV, 1.93901*eV, 1.93146*eV, 1.92397*eV, 1.91654*eV, 1.90916*eV, 1.90184*eV, 
    1.89458*eV, 1.88737*eV, 1.88021*eV, 1.87311*eV, 1.86606*eV, 1.85907*eV, 1.85213*eV, 1.84524*eV, 1.8384*eV, 1.83161*eV, 
    1.82487*eV, 1.81818*eV, 1.81154*eV, 1.80495*eV, 1.7984*eV, 1.79191*eV, 1.78546*eV, 1.77905*eV, 1.77269*eV, 1.76638*eV, 
    1.76011*eV, 1.75389*eV, 1.74771*eV, 1.74157*eV, 1.73548*eV, 1.72943*eV, 1.72342*eV, 1.71745*eV, 1.71153*eV, 1.70564*eV, 
    1.69979*eV, 1.69399*eV, 1.68822*eV, 1.6825*eV, 1.67681*eV, 1.67116*eV, 1.66555*eV, 1.65997*eV, 1.65444*eV, 1.64894*eV, 
    1.64347*eV, 1.63804*eV, 1.63265*eV, 1.6273*eV, 1.62198*eV, 1.61669*eV, 1.61144*eV, 1.60622*eV, 1.60103*eV, 1.59588*eV, 
    1.59076*eV, 1.58568*eV, 1.58062*eV, 1.5756*eV, 1.57061*eV, 1.56566*eV, 1.56073*eV};
    
  G4double FAST_COMPONENT[NUMENTRIES_1] = 
  {0.017312, 0.0096453, 0.00790803, 0.00810917, 0.0088267, 0.0100443, 0.0119759, 0.013618, 0.0158952, 0.0205216, 0.0276587, 
    0.0355978, 0.047533, 0.063484, 0.0846373, 0.110553, 0.142739, 0.179262, 0.220629, 0.265465, 0.316428, 0.373613, 0.435236,
    0.503049, 0.573143, 0.645777, 0.716553, 0.781489, 0.830571, 0.868259, 0.892522, 0.905547, 0.90921, 0.904013, 0.894359, 
    0.879165, 0.863961, 0.848755, 0.832279, 0.814811, 0.796789, 0.781315, 0.763563, 0.741355, 0.724224, 0.705468, 0.687985, 
    0.668306, 0.65034, 0.62812, 0.608418, 0.588743, 0.569351, 0.545681, 0.521421, 0.500633, 0.47682, 0.455015, 0.432658, 
    0.408202, 0.385464, 0.363256, 0.341688, 0.321938, 0.302229, 0.283179, 0.264294, 0.244777, 0.227309, 0.208308, 0.193561, 
    0.178659, 0.163398, 0.14938, 0.137934, 0.126634, 0.116962, 0.106936, 0.0974455, 0.0876783, 0.0803322, 0.073087, 0.0664133, 
    0.0601629, 0.0542353, 0.0497304, 0.0440716, 0.0397871, 0.035986, 0.0335943, 0.0305447, 0.027748, 0.0251957, 0.0232941, 
    0.0210163, 0.0181465, 0.0167443, 0.014946, 0.0138611, 0.0129008, 0.0123556, 0.0118624, 0.011176, 0.0102985, 0.00944158, 
    0.00946978, 0.00862405, 0.00771896, 0.00694318, 0.00643815, 0.00590943, 0.00537023, 0.00446022, 0.00436113, 0.0043437, 
    0.00408437, 0.00385623, 0.0032281, 0.00318978, 0.00319752, 0.00311753, 0.00256044, 0.00262745, 0.00253195, 0.00259158, 
    0.00233708, 0.00213198, 0.00234155, 0.0024483, 0.00209324, 0.0018311, 0.00175891, 0.00183815, 0.00192849, 0.00163606, 
    0.00136759, 0.00170809};
  
  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);
  
  mt->AddConstProperty("SCINTILLATIONYIELD", 33200/MeV);
  mt->AddConstProperty("RESOLUTIONSCALE",1);
  mt->AddConstProperty("FASTTIMECONSTANT", 90.*ns);
  mt->AddConstProperty("YIELDRATIO",1.0);
  mt->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
  
  mat->SetMaterialPropertiesTable(mt);
  return mat;
}

G4Material* MyMaterials::YAG_Ce_Mg()
{
  G4double a, z, density;
  
  G4Element* Y  = new G4Element("Silicon",   "Y",  z = 39., a = 88.01* g/mole);
  G4Element* Al = new G4Element("Aluminium", "Al", z = 13., a = 28.09* g/mole);
  G4Element* O  = new G4Element("Oxygen",    "O",  z =  8., a = 16.00* g/mole);
  
  G4Material* mat = new G4Material ("YAG_Ce_Mg", density = 4.6*g/cm3,3);
  mat->AddElement(Y,3);
  mat->AddElement(Al,5);
  mat->AddElement(O,12);
/*  
  const G4int nEntries_RI = 40;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV}; //{,      4.5085*eV, 4.9594*eV, 6.5000*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.8115, 1.8116, 1.8118, 1.8122, 
      1.8126, 1.8131, 1.8135, 1.8140, 
      1.8144, 1.8149, 1.8154, 1.8160, 
      1.8165, 1.8171, 1.8177, 1.8184, 
      1.8191, 1.8198, 1.8206, 1.8214, 
      1.8223, 1.8233, 1.8244, 1.8256, 
      1.8269, 1.8284, 1.8300, 1.8318, 
      1.8338, 1.8362, 1.8388, 1.8419, 
      1.8455, 1.8497, 1.8548, 1.8608, 
      1.8683, 1.8775, 1.8892, 1.9045}; //, 1.9249, 1.9532, 1.9600 };
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.84,   1.84,   1.84,    1.84,    1.84,    1.84,    1.84,   1.84};

  //G4double Rayleigh[nEntries_RI] =
  //  { 138.*mm, 138.*mm, 138.*mm };
  
   //intrinsic absorption spectrum
  const G4int nEntries_ABS = 84;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    {4.42857*eV, 4.35088*eV, 4.27586*eV, 4.20339*eV, 4.13333*eV, 4.06557*eV, 4*eV, 3.93651*eV, 3.875*eV, 
    3.81538*eV, 3.75758*eV, 3.70149*eV, 3.64706*eV, 3.5942*eV, 3.54286*eV, 3.49296*eV, 3.44444*eV, 
    3.39726*eV, 3.35135*eV, 3.30667*eV, 3.26316*eV, 3.22078*eV, 3.17949*eV, 3.13924*eV, 3.1*eV, 3.06173*eV, 
    3.02439*eV, 2.98795*eV, 2.95238*eV, 2.91765*eV, 2.88372*eV, 2.85057*eV, 2.81818*eV, 2.78652*eV, 
    2.75556*eV, 2.72527*eV, 2.69565*eV, 2.66667*eV, 2.6383*eV, 2.61053*eV, 2.58333*eV, 2.5567*eV, 
    2.53061*eV, 2.50505*eV, 2.48*eV, 2.45545*eV, 2.43137*eV, 2.40777*eV, 2.38462*eV, 2.3619*eV, 
    2.33962*eV, 2.31776*eV, 2.2963*eV, 2.27523*eV, 2.25455*eV, 2.23423*eV, 2.21429*eV, 2.19469*eV, 
    2.17544*eV, 2.15652*eV, 2.13793*eV, 2.11966*eV, 2.10169*eV, 2.08403*eV, 2.06667*eV, 2.04959*eV, 
    2.03279*eV, 2.01626*eV, 2*eV, 1.984*eV, 1.96825*eV, 1.95276*eV, 1.9375*eV, 1.92248*eV, 1.90769*eV, 
    1.89313*eV, 1.87879*eV, 1.86466*eV, 1.85075*eV, 1.83704*eV, 1.82353*eV, 1.81022*eV, 1.7971*eV, 1.78417*eV};
    
  G4double Absorption[nEntries_ABS] =
    {12.7389*mm, 7.96854*mm, 9.79472*mm, 11.3642*mm, 13.9544*mm, 14.078*mm, 9.48293*mm, 4.60862*mm, 
    2.45205*mm, 2.14906*mm, 2.09215*mm, 2.05898*mm, 2.09349*mm, 2.09448*mm, 4.00614*mm, 14.055*mm, 
    59.4712*mm, 177.765*mm, 322.779*mm, 407.883*mm, 423.384*mm, 250.626*mm, 107.986*mm, 40.2593*mm, 
    16.1867*mm, 7.08721*mm, 3.3827*mm, 1.85802*mm, 1.42111*mm, 1.348*mm, 1.36044*mm, 1.37582*mm, 
    1.39199*mm, 1.33268*mm, 1.27998*mm, 1.33348*mm, 1.21815*mm, 1.28273*mm, 1.37124*mm, 1.4064*mm, 
    1.39597*mm, 1.45447*mm, 2.03141*mm, 3.62354*mm, 6.92578*mm, 13.9067*mm, 29.0131*mm, 60.2374*mm, 
    129.476*mm, 226.564*mm, 553.51*mm, 516.632*mm, 600.*mm, 700.*mm, 800.*mm, 900.*mm, 
    1000.*mm, 1000.*mm, 1000.*mm, 1000.2*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 
    1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 
    1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 
    1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm};
  
  //emission spectrum
  const G4int NUMENTRIES_1 = 137;
  G4double FAST_Energy[NUMENTRIES_1]    = 
  {2.72827*eV, 2.71335*eV, 2.69859*eV, 2.68398*eV, 2.66954*eV, 2.65525*eV, 2.64111*eV, 2.62712*eV, 2.61328*eV, 2.59958*eV, 
    2.58603*eV, 2.57261*eV, 2.55934*eV, 2.5462*eV, 2.5332*eV, 2.52033*eV, 2.50758*eV, 2.49497*eV, 2.48248*eV, 2.47012*eV, 
    2.45788*eV, 2.44576*eV, 2.43376*eV, 2.42188*eV, 2.41011*eV, 2.39845*eV, 2.38691*eV, 2.37548*eV, 2.36416*eV, 2.35294*eV, 
    2.34183*eV, 2.33083*eV, 2.31993*eV, 2.30912*eV, 2.29842*eV, 2.28782*eV, 2.27732*eV, 2.26691*eV, 2.2566*eV, 2.24638*eV, 
    2.23625*eV, 2.22621*eV, 2.21626*eV, 2.20641*eV, 2.19663*eV, 2.18695*eV, 2.17735*eV, 2.16783*eV, 2.1584*eV, 2.14905*eV, 
    2.13978*eV, 2.13058*eV, 2.12147*eV, 2.11244*eV, 2.10348*eV, 2.09459*eV, 2.08579*eV, 2.07705*eV, 2.06839*eV, 2.0598*eV, 
    2.05128*eV, 2.04283*eV, 2.03445*eV, 2.02614*eV, 2.0179*eV, 2.00972*eV, 2.00161*eV, 1.99357*eV, 1.98559*eV, 1.97767*eV, 
    1.96982*eV, 1.96203*eV, 1.95429*eV, 1.94662*eV, 1.93901*eV, 1.93146*eV, 1.92397*eV, 1.91654*eV, 1.90916*eV, 1.90184*eV, 
    1.89458*eV, 1.88737*eV, 1.88021*eV, 1.87311*eV, 1.86606*eV, 1.85907*eV, 1.85213*eV, 1.84524*eV, 1.8384*eV, 1.83161*eV, 
    1.82487*eV, 1.81818*eV, 1.81154*eV, 1.80495*eV, 1.7984*eV, 1.79191*eV, 1.78546*eV, 1.77905*eV, 1.77269*eV, 1.76638*eV, 
    1.76011*eV, 1.75389*eV, 1.74771*eV, 1.74157*eV, 1.73548*eV, 1.72943*eV, 1.72342*eV, 1.71745*eV, 1.71153*eV, 1.70564*eV, 
    1.69979*eV, 1.69399*eV, 1.68822*eV, 1.6825*eV, 1.67681*eV, 1.67116*eV, 1.66555*eV, 1.65997*eV, 1.65444*eV, 1.64894*eV, 
    1.64347*eV, 1.63804*eV, 1.63265*eV, 1.6273*eV, 1.62198*eV, 1.61669*eV, 1.61144*eV, 1.60622*eV, 1.60103*eV, 1.59588*eV, 
    1.59076*eV, 1.58568*eV, 1.58062*eV, 1.5756*eV, 1.57061*eV, 1.56566*eV, 1.56073*eV};
    
  G4double FAST_COMPONENT[NUMENTRIES_1] = 
  {0.017312, 0.0096453, 0.00790803, 0.00810917, 0.0088267, 0.0100443, 0.0119759, 0.013618, 0.0158952, 0.0205216, 0.0276587, 
    0.0355978, 0.047533, 0.063484, 0.0846373, 0.110553, 0.142739, 0.179262, 0.220629, 0.265465, 0.316428, 0.373613, 0.435236,
    0.503049, 0.573143, 0.645777, 0.716553, 0.781489, 0.830571, 0.868259, 0.892522, 0.905547, 0.90921, 0.904013, 0.894359, 
    0.879165, 0.863961, 0.848755, 0.832279, 0.814811, 0.796789, 0.781315, 0.763563, 0.741355, 0.724224, 0.705468, 0.687985, 
    0.668306, 0.65034, 0.62812, 0.608418, 0.588743, 0.569351, 0.545681, 0.521421, 0.500633, 0.47682, 0.455015, 0.432658, 
    0.408202, 0.385464, 0.363256, 0.341688, 0.321938, 0.302229, 0.283179, 0.264294, 0.244777, 0.227309, 0.208308, 0.193561, 
    0.178659, 0.163398, 0.14938, 0.137934, 0.126634, 0.116962, 0.106936, 0.0974455, 0.0876783, 0.0803322, 0.073087, 0.0664133, 
    0.0601629, 0.0542353, 0.0497304, 0.0440716, 0.0397871, 0.035986, 0.0335943, 0.0305447, 0.027748, 0.0251957, 0.0232941, 
    0.0210163, 0.0181465, 0.0167443, 0.014946, 0.0138611, 0.0129008, 0.0123556, 0.0118624, 0.011176, 0.0102985, 0.00944158, 
    0.00946978, 0.00862405, 0.00771896, 0.00694318, 0.00643815, 0.00590943, 0.00537023, 0.00446022, 0.00436113, 0.0043437, 
    0.00408437, 0.00385623, 0.0032281, 0.00318978, 0.00319752, 0.00311753, 0.00256044, 0.00262745, 0.00253195, 0.00259158, 
    0.00233708, 0.00213198, 0.00234155, 0.0024483, 0.00209324, 0.0018311, 0.00175891, 0.00183815, 0.00192849, 0.00163606, 
    0.00136759, 0.00170809};
  
  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);
  
  mt->AddConstProperty("SCINTILLATIONYIELD", 17000/MeV);
  mt->AddConstProperty("RESOLUTIONSCALE",5);
  mt->AddConstProperty("FASTTIMECONSTANT", 60.*ns);
  mt->AddConstProperty("YIELDRATIO",1.0);
  mt->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
  
  mat->SetMaterialPropertiesTable(mt);
  return mat;
}



G4Material* MyMaterials::GAGG_Ce() // Gadolinium Gallium Aluminum Garnet - undoped
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z = 8.,  a = 16.00  *g/mole);
  G4Element* Ga = new G4Element("Gallium",  "Ga", z = 31., a = 69.723 *g/mole);
  G4Element* Gd = new G4Element("Gadolinio","Gd", z = 64., a = 157.25 *g/mole);
  G4Element* Al = new G4Element("Aluminum", "Al", z = 13., a = 28.09  *g/mole);
  
  G4Material* mat = new G4Material("GAG_Ce", density = 6.63 *g/cm3, 4);
  mat->AddElement(Ga,3);
  mat->AddElement(Gd,3);
  mat->AddElement(Al,2);
  mat->AddElement(O,12);
/*
  const G4int nEntries_RI = 9;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    {0.1*eV, 1.0*eV, 1.550*eV, 1.771*eV, 1.907*eV, 2.105*eV, 2.917*eV, 3.1*eV, 4.1*eV}; //6.5*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    {1.890, 1.892, 1.897, 1.898, 1.899, 1.903, 1.937, 1.950, 1.983}; //2.400};
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.92,   1.92,   1.92,    1.92,    1.92,    1.92,    1.92,   1.92};
  
 //intrinsic absorption spectrum
  const G4int nEntries_ABS = 122;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    {1.55*eV, 1.55975*eV, 1.56962*eV, 1.57962*eV, 1.58974*eV, 1.6*eV, 1.61039*eV, 1.62092*eV, 1.63158*eV, 1.64238*eV, 1.65333*eV, 1.66443*eV, 1.67568*eV, 1.68707*eV, 1.69863*eV, 1.71034*eV, 1.72222*eV, 1.73427*eV, 1.74648*eV, 1.75887*eV, 1.77143*eV, 1.78417*eV, 1.7971*eV, 1.81022*eV, 1.82353*eV, 1.83704*eV, 1.85075*eV, 1.86466*eV, 1.87879*eV, 1.89313*eV, 1.90769*eV, 1.92248*eV, 1.9375*eV, 1.95276*eV, 1.96825*eV, 1.984*eV, 2*eV, 2.01626*eV, 2.03279*eV, 2.04959*eV, 2.06667*eV, 2.08403*eV, 2.10169*eV, 2.11966*eV, 2.13793*eV, 2.15652*eV, 2.17544*eV, 2.19469*eV, 2.21429*eV, 2.23423*eV, 2.25455*eV, 2.27523*eV, 2.2963*eV, 2.31776*eV, 2.33962*eV, 2.3619*eV, 2.38462*eV, 2.40777*eV, 2.43137*eV, 2.45545*eV, 2.48*eV, 2.50505*eV, 2.53061*eV, 2.5567*eV, 2.58333*eV, 2.61053*eV, 2.6383*eV, 2.66667*eV, 2.69565*eV, 2.72527*eV, 2.75556*eV, 2.78652*eV, 2.81818*eV, 2.85057*eV, 2.88372*eV, 2.92453*eV, 2.95943*eV, 3.00971*eV, 3.04668*eV, 3.08458*eV, 3.12343*eV, 3.16327*eV, 3.20413*eV, 3.24607*eV, 3.28912*eV, 3.33333*eV, 3.37875*eV, 3.42541*eV, 3.47339*eV, 3.52273*eV, 3.57349*eV, 3.63636*eV, 3.69048*eV, 3.74622*eV, 3.80368*eV, 3.86293*eV, 3.92405*eV, 3.98714*eV, 4.05229*eV, 4.1196*eV, 4.18919*eV, 4.26117*eV, 4.33566*eV, 4.41281*eV, 4.49275*eV, 4.57565*eV, 4.66165*eV, 4.75096*eV, 4.84375*eV, 4.94024*eV, 5.04065*eV, 5.18828*eV, 5.32189*eV, 5.4386*eV, 5.61086*eV, 5.76744*eV, 5.93301*eV, 6.07843*eV, 6.23116*eV, 6.39175*eV, 6.52632*eV};
    
  G4double Absorption[nEntries_ABS] =
    {1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 942.751*mm, 845.85*mm, 821.483*mm, 691.245*mm, 627.478*mm, 671.166*mm, 660.321*mm, 652.537*mm, 633.443*mm, 617.952*mm, 595.034*mm, 570.832*mm, 557.456*mm, 518.822*mm, 517.287*mm, 503.613*mm, 500.73*mm, 471.285*mm, 458.997*mm, 438.515*mm, 425.546*mm, 420.399*mm, 418.651*mm, 417.589*mm, 416.005*mm, 410.057*mm, 408.373*mm, 406.603*mm, 405.734*mm, 398.214*mm, 390.34*mm, 385.343*mm, 380.658*mm, 371.93*mm, 370.279*mm, 367.872*mm, 333.091*mm, 313.972*mm, 275.303*mm, 233.395*mm, 178.615*mm, 120.95*mm, 73.2142*mm, 41.3708*mm, 22.3063*mm, 12.0744*mm, 6.61098*mm, 3.74066*mm, 2.22006*mm, 1.43869*mm, 1.11985*mm, 1.31942*mm, 0.5*mm, 0.05*mm, 0.05*mm, 0.05*mm, 0.05*mm, 0.05*mm, 0.05*mm, 0.05*mm, 0.5*mm, 1.06781*mm, 1.78258*mm, 3.36962*mm, 6.81599*mm, 13.7903*mm, 23.6781*mm, 27.9876*mm, 23.5133*mm, 15.59*mm, 8.50471*mm, 3.964*mm, 1.86267*mm, 0.5*mm, 1.28132*mm, 1.03618*mm, 1.13582*mm, 0.5*mm, 1.36591*mm, 2.0196*mm, 2.45141*mm, 2.52958*mm, 2.53377*mm, 2.43753*mm, 2.262*mm, 1.99166*mm, 1.69886*mm, 1.43058*mm, 1.31941*mm, 1.25531*mm, 0.5*mm, 1.25407*mm, 1.25137*mm, 0.5*mm, 0.05*mm, 0.5*mm, 0.922364*mm, 1.28007*mm, 1.33808*mm, 1.40961*mm, 1.31044*mm, 1.48699*mm, 1.25399*mm, 1.2*mm, 0.1*mm};
  
  //emission spectrum
  const G4int NUMENTRIES_1 = 132;
  G4double FAST_Energy[NUMENTRIES_1]    = 
  {2.6327*eV, 2.6188*eV, 2.60504*eV, 2.59143*eV, 2.57796*eV, 2.56463*eV, 2.55144*eV, 2.53838*eV, 2.52546*eV, 2.51266*eV, 2.5*eV, 2.48746*eV, 2.47505*eV, 2.46276*eV, 2.45059*eV, 2.43854*eV, 2.42661*eV, 2.4148*eV, 2.4031*eV, 2.39151*eV, 2.38004*eV, 2.36867*eV, 2.35741*eV, 2.34626*eV, 2.33522*eV, 2.32427*eV, 2.31343*eV, 2.30269*eV, 2.29205*eV, 2.28151*eV, 2.27106*eV, 2.26071*eV, 2.25045*eV, 2.24029*eV, 2.23022*eV, 2.22023*eV, 2.21034*eV, 2.20053*eV, 2.19081*eV, 2.18118*eV, 2.17163*eV, 2.16216*eV, 2.15278*eV, 2.14347*eV, 2.13425*eV, 2.12511*eV, 2.11604*eV, 2.10705*eV, 2.09814*eV, 2.0893*eV, 2.08054*eV, 2.07185*eV, 2.06323*eV, 2.05468*eV, 2.0462*eV, 2.0378*eV, 2.02946*eV, 2.02119*eV, 2.01299*eV, 2.00485*eV, 1.99678*eV, 1.98877*eV, 1.98083*eV, 1.97295*eV, 1.96513*eV, 1.95738*eV, 1.94969*eV, 1.94205*eV, 1.93448*eV, 1.92696*eV, 1.9195*eV, 1.9121*eV, 1.90476*eV, 1.89748*eV, 1.89024*eV, 1.88307*eV, 1.87595*eV, 1.86888*eV, 1.86186*eV, 1.8549*eV, 1.84799*eV, 1.84113*eV, 1.83432*eV, 1.82756*eV, 1.82085*eV, 1.81419*eV, 1.80758*eV, 1.80102*eV, 1.7945*eV, 1.78803*eV, 1.78161*eV, 1.77523*eV, 1.7689*eV, 1.76262*eV, 1.75637*eV, 1.75018*eV, 1.74402*eV, 1.73791*eV, 1.73184*eV, 1.72582*eV, 1.71983*eV, 1.71389*eV, 1.70799*eV, 1.70213*eV, 1.69631*eV, 1.69052*eV, 1.68478*eV, 1.67908*eV, 1.67341*eV, 1.66779*eV, 1.6622*eV, 1.65665*eV, 1.65113*eV, 1.64565*eV, 1.64021*eV, 1.63481*eV, 1.62943*eV, 1.6241*eV, 1.6188*eV, 1.61353*eV, 1.6083*eV, 1.6031*eV, 1.59794*eV, 1.59281*eV, 1.58771*eV, 1.58264*eV, 1.57761*eV, 1.57261*eV, 1.56764*eV, 1.5627*eV, 1.55779*eV, 1.55291*eV};

    
  G4double FAST_COMPONENT[NUMENTRIES_1] = 
  {0.0348132, 0.041346, 0.0462213, 0.055994, 0.0670797, 0.0809222, 0.0962054, 0.119476, 0.147153, 0.182414, 0.217745, 0.256009, 0.299913, 0.34846, 0.391284, 0.436086, 0.487627, 0.538087, 0.591007, 0.64409, 0.68982, 0.740803, 0.785007, 0.825553, 0.855418, 0.87689, 0.893165, 0.909735, 0.904133, 0.904835, 0.904515, 0.896044, 0.885449, 0.872733, 0.857353, 0.837832, 0.829267, 0.810443, 0.78659, 0.767397, 0.745877, 0.730203, 0.705547, 0.682982, 0.665159, 0.645723, 0.621589, 0.595471, 0.576994, 0.556964, 0.533997, 0.507302, 0.484761, 0.463904, 0.439846, 0.418253, 0.404643, 0.384661, 0.360237, 0.339777, 0.324817, 0.304963, 0.287939, 0.271009, 0.252333, 0.235159, 0.226375, 0.212185, 0.196952, 0.184789, 0.1707, 0.159992, 0.147724, 0.140217, 0.130929, 0.121987, 0.113251, 0.102937, 0.0964418, 0.0913118, 0.0863073, 0.0811948, 0.0732201, 0.0693917, 0.0642312, 0.0575242, 0.0533638, 0.0488039, 0.0474854, 0.0466519, 0.0428456, 0.0431906, 0.0402339, 0.0363488, 0.0340234, 0.0328062, 0.0313949, 0.0292326, 0.0259991, 0.0247748, 0.0240735, 0.0210389, 0.021114, 0.0211829, 0.0192813, 0.0187306, 0.0185142, 0.0171277, 0.0167855, 0.0169199, 0.016618, 0.0150835, 0.0137522, 0.0137759, 0.0131568, 0.0134337, 0.0135361, 0.0132448, 0.0122409, 0.0119099, 0.011421, 0.0115286, 0.0103885, 0.0105544, 0.0118189, 0.0109937, 0.0102842, 0.0103254, 0.0101991, 0.00931661, 0.00833321, 0.00805133};

  //** scintillation and optical properties **//
  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);
  
  mt->AddConstProperty("SCINTILLATIONYIELD", 53000/MeV);
  mt->AddConstProperty("RESOLUTIONSCALE", 1);
  mt->AddConstProperty("FASTTIMECONSTANT", 100.*ns);
  mt->AddConstProperty("SLOWTIMECONSTANT", 300.*ns);
  mt->AddConstProperty("YIELDRATIO",0.7);
  mt->AddConstProperty("FASTSCINTILLATIONRISETIME", 1.8*ns);

  mat->SetMaterialPropertiesTable(mt);
  
  return mat;
}

G4Material* MyMaterials::GAGG_Ce_Mg() // Gadolinium Gallium Aluminum Garnet - undoped
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z = 8.,  a = 16.00  *g/mole);
  G4Element* Ga = new G4Element("Gallium",  "Ga", z = 31., a = 69.723 *g/mole);
  G4Element* Gd = new G4Element("Gadolinio","Gd", z = 64., a = 157.25 *g/mole);
  G4Element* Al = new G4Element("Aluminum", "Al", z = 13., a = 28.09  *g/mole);
  
  G4Material* mat = new G4Material("GAG_Ce_Mg", density = 6.63 *g/cm3, 4);
  mat->AddElement(Ga,3);
  mat->AddElement(Gd,3);
  mat->AddElement(Al,2);
  mat->AddElement(O,12);
/*
  const G4int nEntries_RI = 9;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    {0.1*eV, 1.0*eV, 1.550*eV, 1.771*eV, 1.907*eV, 2.105*eV, 2.917*eV, 3.1*eV, 4.1*eV}; //6.5*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    {1.890, 1.892, 1.897, 1.898, 1.899, 1.903, 1.937, 1.950, 1.983}; //2.400};
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.92,   1.92,   1.92,    1.92,    1.92,    1.92,    1.92,   1.92};
  
 //intrinsic absorption spectrum
  const G4int nEntries_ABS = 122;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    {1.55*eV, 1.55975*eV, 1.56962*eV, 1.57962*eV, 1.58974*eV, 1.6*eV, 1.61039*eV, 1.62092*eV, 1.63158*eV, 1.64238*eV, 1.65333*eV, 1.66443*eV, 1.67568*eV, 1.68707*eV, 1.69863*eV, 1.71034*eV, 1.72222*eV, 1.73427*eV, 1.74648*eV, 1.75887*eV, 1.77143*eV, 1.78417*eV, 1.7971*eV, 1.81022*eV, 1.82353*eV, 1.83704*eV, 1.85075*eV, 1.86466*eV, 1.87879*eV, 1.89313*eV, 1.90769*eV, 1.92248*eV, 1.9375*eV, 1.95276*eV, 1.96825*eV, 1.984*eV, 2*eV, 2.01626*eV, 2.03279*eV, 2.04959*eV, 2.06667*eV, 2.08403*eV, 2.10169*eV, 2.11966*eV, 2.13793*eV, 2.15652*eV, 2.17544*eV, 2.19469*eV, 2.21429*eV, 2.23423*eV, 2.25455*eV, 2.27523*eV, 2.2963*eV, 2.31776*eV, 2.33962*eV, 2.3619*eV, 2.38462*eV, 2.40777*eV, 2.43137*eV, 2.45545*eV, 2.48*eV, 2.50505*eV, 2.53061*eV, 2.5567*eV, 2.58333*eV, 2.61053*eV, 2.6383*eV, 2.66667*eV, 2.69565*eV, 2.72527*eV, 2.75556*eV, 2.78652*eV, 2.81818*eV, 2.85057*eV, 2.88372*eV, 2.91765*eV, 2.95238*eV, 2.98795*eV, 3.02439*eV, 3.06173*eV, 3.1*eV, 3.13924*eV, 3.17949*eV, 3.22078*eV, 3.26316*eV, 3.30667*eV, 3.35135*eV, 3.39726*eV, 3.44444*eV, 3.49296*eV, 3.54286*eV, 3.60465*eV, 3.65782*eV, 3.72372*eV, 3.78049*eV, 3.83901*eV, 3.89937*eV, 3.96166*eV, 4.02597*eV, 4.09241*eV, 4.16107*eV, 4.23208*eV, 4.30556*eV, 4.38163*eV, 4.46043*eV, 4.55882*eV, 4.64419*eV, 4.75096*eV, 4.84375*eV, 4.94024*eV, 5.04065*eV, 5.14523*eV, 5.25424*eV, 5.36797*eV, 5.48673*eV, 5.61086*eV, 5.74074*eV, 5.87678*eV, 6.01942*eV, 6.16915*eV, 6.32653*eV, 6.49215*eV};
    
  G4double Absorption[nEntries_ABS] =
    {1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 1000.*mm, 950.*mm, 900.*mm, 850.*mm, 800.*mm, 750.*mm, 700.*mm, 650.*mm, 600.*mm, 600.*mm, 550.*mm, 525.*mm, 500.*mm, 495.113*mm, 493.925*mm, 490.185*mm, 485.183*mm, 480.656*mm, 475.745*mm, 474.191*mm, 470.37*mm, 465.699*mm, 460.696*mm, 455.925*mm, 450.141*mm, 445.712*mm, 440.798*mm, 435.2*mm, 430.659*mm, 425.727*mm, 420.328*mm, 415.261*mm, 410.212*mm, 405.257*mm, 402.131*mm, 401.244*mm, 400.65*mm, 395.193*mm, 390.745*mm, 385.259*mm, 379.265*mm, 373.136*mm, 374.806*mm, 326.031*mm, 311.609*mm, 281.734*mm, 217.074*mm, 156.221*mm, 101.25*mm, 59.9635*mm, 32.9009*mm, 17.6943*mm, 9.62468*mm, 5.33133*mm, 3.06042*mm, 1.86737*mm, 0.931539*mm, 1.17313*mm, 1.33421*mm, 0.5*mm, 0.5*mm, 1.0692*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.994237*mm, 1.11628*mm, 0.5*mm, 1.38984*mm, 2.55123*mm, 5.06262*mm, 10.5534*mm, 20.0187*mm, 24.5102*mm, 17.6148*mm, 10.163*mm, 5.66707*mm, 3.08815*mm, 1.63268*mm, 0.5*mm, 1.0365*mm, 0.5*mm, 0.5*mm, 1.06797*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.926825*mm, 0.5*mm, 0.5*mm, 0.926322*mm, 0.925969*mm, 0.5*mm, 0.5*mm, 1.19071*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 0.5*mm, 1.25215*mm, 1.07403*mm, 1.16123*mm, 1.22498*mm};
  
  //emission spectrum
  const G4int NUMENTRIES_1 = 132;
  G4double FAST_Energy[NUMENTRIES_1]    = 
  {2.6327*eV, 2.6188*eV, 2.60504*eV, 2.59143*eV, 2.57796*eV, 2.56463*eV, 2.55144*eV, 2.53838*eV, 2.52546*eV, 2.51266*eV, 2.5*eV, 2.48746*eV, 2.47505*eV, 2.46276*eV, 2.45059*eV, 2.43854*eV, 2.42661*eV, 2.4148*eV, 2.4031*eV, 2.39151*eV, 2.38004*eV, 2.36867*eV, 2.35741*eV, 2.34626*eV, 2.33522*eV, 2.32427*eV, 2.31343*eV, 2.30269*eV, 2.29205*eV, 2.28151*eV, 2.27106*eV, 2.26071*eV, 2.25045*eV, 2.24029*eV, 2.23022*eV, 2.22023*eV, 2.21034*eV, 2.20053*eV, 2.19081*eV, 2.18118*eV, 2.17163*eV, 2.16216*eV, 2.15278*eV, 2.14347*eV, 2.13425*eV, 2.12511*eV, 2.11604*eV, 2.10705*eV, 2.09814*eV, 2.0893*eV, 2.08054*eV, 2.07185*eV, 2.06323*eV, 2.05468*eV, 2.0462*eV, 2.0378*eV, 2.02946*eV, 2.02119*eV, 2.01299*eV, 2.00485*eV, 1.99678*eV, 1.98877*eV, 1.98083*eV, 1.97295*eV, 1.96513*eV, 1.95738*eV, 1.94969*eV, 1.94205*eV, 1.93448*eV, 1.92696*eV, 1.9195*eV, 1.9121*eV, 1.90476*eV, 1.89748*eV, 1.89024*eV, 1.88307*eV, 1.87595*eV, 1.86888*eV, 1.86186*eV, 1.8549*eV, 1.84799*eV, 1.84113*eV, 1.83432*eV, 1.82756*eV, 1.82085*eV, 1.81419*eV, 1.80758*eV, 1.80102*eV, 1.7945*eV, 1.78803*eV, 1.78161*eV, 1.77523*eV, 1.7689*eV, 1.76262*eV, 1.75637*eV, 1.75018*eV, 1.74402*eV, 1.73791*eV, 1.73184*eV, 1.72582*eV, 1.71983*eV, 1.71389*eV, 1.70799*eV, 1.70213*eV, 1.69631*eV, 1.69052*eV, 1.68478*eV, 1.67908*eV, 1.67341*eV, 1.66779*eV, 1.6622*eV, 1.65665*eV, 1.65113*eV, 1.64565*eV, 1.64021*eV, 1.63481*eV, 1.62943*eV, 1.6241*eV, 1.6188*eV, 1.61353*eV, 1.6083*eV, 1.6031*eV, 1.59794*eV, 1.59281*eV, 1.58771*eV, 1.58264*eV, 1.57761*eV, 1.57261*eV, 1.56764*eV, 1.5627*eV, 1.55779*eV, 1.55291*eV};

    
  G4double FAST_COMPONENT[NUMENTRIES_1] = 
  {0.0348132, 0.041346, 0.0462213, 0.055994, 0.0670797, 0.0809222, 0.0962054, 0.119476, 0.147153, 0.182414, 0.217745, 0.256009, 0.299913, 0.34846, 0.391284, 0.436086, 0.487627, 0.538087, 0.591007, 0.64409, 0.68982, 0.740803, 0.785007, 0.825553, 0.855418, 0.87689, 0.893165, 0.909735, 0.904133, 0.904835, 0.904515, 0.896044, 0.885449, 0.872733, 0.857353, 0.837832, 0.829267, 0.810443, 0.78659, 0.767397, 0.745877, 0.730203, 0.705547, 0.682982, 0.665159, 0.645723, 0.621589, 0.595471, 0.576994, 0.556964, 0.533997, 0.507302, 0.484761, 0.463904, 0.439846, 0.418253, 0.404643, 0.384661, 0.360237, 0.339777, 0.324817, 0.304963, 0.287939, 0.271009, 0.252333, 0.235159, 0.226375, 0.212185, 0.196952, 0.184789, 0.1707, 0.159992, 0.147724, 0.140217, 0.130929, 0.121987, 0.113251, 0.102937, 0.0964418, 0.0913118, 0.0863073, 0.0811948, 0.0732201, 0.0693917, 0.0642312, 0.0575242, 0.0533638, 0.0488039, 0.0474854, 0.0466519, 0.0428456, 0.0431906, 0.0402339, 0.0363488, 0.0340234, 0.0328062, 0.0313949, 0.0292326, 0.0259991, 0.0247748, 0.0240735, 0.0210389, 0.021114, 0.0211829, 0.0192813, 0.0187306, 0.0185142, 0.0171277, 0.0167855, 0.0169199, 0.016618, 0.0150835, 0.0137522, 0.0137759, 0.0131568, 0.0134337, 0.0135361, 0.0132448, 0.0122409, 0.0119099, 0.011421, 0.0115286, 0.0103885, 0.0105544, 0.0118189, 0.0109937, 0.0102842, 0.0103254, 0.0101991, 0.00931661, 0.00833321, 0.00805133};

  //** scintillation and optical properties **//
  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty ("RINDEX",    PhotonEnergy_RI,  RefractiveIndex, nEntries_RI);
  mt->AddProperty ("ABSLENGTH", PhotonEnergy_ABS, Absorption,      nEntries_ABS);
  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);
  
  mt->AddConstProperty("SCINTILLATIONYIELD", 42000/MeV);
  mt->AddConstProperty("RESOLUTIONSCALE", 1);
  mt->AddConstProperty("FASTTIMECONSTANT", 50.*ns);
  mt->AddConstProperty("SLOWTIMECONSTANT", 200.*ns);
  mt->AddConstProperty("YIELDRATIO",0.6);
  mt->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.04*ns);

  mat->SetMaterialPropertiesTable(mt);
  
  return mat;
}






G4Material* MyMaterials::Brass()
{
  G4double a, z, density;

  G4Element* Cu = new G4Element("Copper", "Cu", z=29., a=63.546*g/mole);
  G4Element* Zn = new G4Element("Zinc",   "Zn", z=30., a=65.409*g/mole);
  
  G4Material* mat = new G4Material("Brass", density=8.73*g/cm3,2);
  mat->AddElement(Cu,0.75);
  mat->AddElement(Zn,0.25);
  
  return mat;
}



G4Material* MyMaterials::Aluminium()
{
  G4NistManager* man = G4NistManager::Instance();
  G4Element* Al = man->FindOrBuildElement("Al");
  
  G4Material* mat = new G4Material("Aluminium",2.700*g/cm3,1);
  mat->AddElement(Al,100.*perCent);
  
  return mat;
}



G4Material* MyMaterials::Iron()
{
  G4NistManager* man = G4NistManager::Instance();
  G4Element* Fe = man->FindOrBuildElement("Fe");
  
  G4Material* mat = new G4Material("Iron",7.874*g/cm3,1);
  mat->AddElement(Fe,100.*perCent);
  
  return mat;
}



G4Material* MyMaterials::Lead()
{
  G4NistManager* man = G4NistManager::Instance();
  G4Element* Pb = man->FindOrBuildElement("Pb");
  
  G4Material* mat = new G4Material("Lead",11.342*g/cm3,1);
  mat->AddElement(Pb,100.*perCent);
  
  return mat;
}



G4Material* MyMaterials::Tungsten()
{
  G4NistManager* man = G4NistManager::Instance();
  G4Element* W = man->FindOrBuildElement("W");
  G4Element* Ni = man->FindOrBuildElement("Ni");
  G4Element* Cu = man->FindOrBuildElement("Cu");
  
  G4Material* mat = new G4Material("Tungsten",17.*g/cm3,3);
  mat->AddElement(W,90.*perCent);
  mat->AddElement(Ni,5.*perCent);
  mat->AddElement(Cu,5.*perCent);
  
  return mat;
}



G4Material* MyMaterials::CopperTungstenAlloy(const G4double& WFrac)
{
  G4NistManager* man = G4NistManager::Instance();
  G4Element* W = man->FindOrBuildElement("W");
  G4Element* Cu = man->FindOrBuildElement("Cu");
  
  G4double rho_Cu = 8.96;
  G4double rho_W = 19.25;
  G4double rho = (1.-WFrac)*rho_Cu + WFrac*rho_W;
  G4Material* mat = new G4Material("CopperTungstenAlloy",rho*g/cm3,2);
  mat->AddElement(Cu,1.-WFrac);
  mat->AddElement(W,WFrac);
  
  return mat;
}



G4Material* MyMaterials::OpticalGrease()
{
  G4double a, z, density;
  G4Element* H = new G4Element("Hydrogen", "H", z=1., a= 1.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a=16.00*g/mole);
  G4Element* C = new G4Element("Carbon"  , "C", z=6., a=12.01*g/mole);

  G4Material* mat = new G4Material("Grease", density=1.0*g/cm3,3);
  mat->AddElement(C,1);
  mat->AddElement(H,1);
  mat->AddElement(O,1);
/*
  const G4int nEntries_RI = 11;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0 * eV, 2.0 * eV, 2.5 * eV, 3.0 * eV,
      3.5 * eV, 4.0 * eV, 4.5 * eV, 5.0 * eV,
      5.5 * eV, 6.0 * eV, 6.26 * eV };
  G4double RefractiveIndex[nEntries_RI] =
    { 1.48, 1.49, 1.50, 1.50,
      1.51, 1.52, 1.54, 1.55,
      1.57, 1.59, 1.60 }; 
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.50,   1.50,   1.50,    1.50,    1.50,    1.50,    1.50,   1.50};

  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX",    PhotonEnergy_RI, RefractiveIndex, nEntries_RI);
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}

G4Material* MyMaterials::OpticalGrease155()
{
  G4double a, z, density;
  G4Element* H = new G4Element("Hydrogen", "H", z=1., a= 1.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a=16.00*g/mole);
  G4Element* C = new G4Element("Carbon"  , "C", z=6., a=12.01*g/mole);

  G4Material* mat = new G4Material("Grease155", density=1.0*g/cm3,3);
  mat->AddElement(C,1);
  mat->AddElement(H,1);
  mat->AddElement(O,1);
/*
  const G4int nEntries_RI = 11;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0 * eV, 2.0 * eV, 2.5 * eV, 3.0 * eV,
      3.5 * eV, 4.0 * eV, 4.5 * eV, 5.0 * eV,
      5.5 * eV, 6.0 * eV, 6.26 * eV };
  G4double RefractiveIndex[nEntries_RI] =
    { 1.53, 1.54, 1.55, 1.56,
      1.56, 1.57, 1.59, 1.60,
      1.62, 1.64, 1.65 };
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.55,   1.55,   1.55,    1.55,    1.55,    1.55,    1.55,   1.55};

  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX",    PhotonEnergy_RI, RefractiveIndex, nEntries_RI);  
  mat->SetMaterialPropertiesTable(myMPT);
  return mat;
}


G4Material* MyMaterials::MeltMount168()
{
  G4double a, z, density;
  G4Element* H = new G4Element("Hydrogen", "H", z=1., a= 1.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a=16.00*g/mole);
  G4Element* C = new G4Element("Carbon"  , "C", z=6., a=12.01*g/mole);

  G4Material* mat = new G4Material("MeltMount168", density=1.0*g/cm3,3);
  mat->AddElement(C,1);
  mat->AddElement(H,1);
  mat->AddElement(O,1);

/*
  const G4int nEntries_RI = 11;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0 * eV, 2.0 * eV, 2.5 * eV, 3.0 * eV,
      3.5 * eV, 4.0 * eV, 4.5 * eV, 5.0 * eV,
      5.5 * eV, 6.0 * eV, 6.26 * eV };
  G4double RefractiveIndex[nEntries_RI] =
    { 1.66, 1.67, 1.68, 1.69,
      1.69, 1.70, 1.72, 1.73,
      1.75, 1.76, 1.78 };*/

  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.68,   1.68,   1.68,    1.68,    1.68,    1.68,    1.68,   1.68};


  const G4int NUMENTRIES_2 = 5;
  G4double ABS_Energy[NUMENTRIES_2] = { 1.0*eV, 1.84*eV, 2.48*eV, 2.75*eV,   3.02*eV};
  G4double ABS_LENGTH[NUMENTRIES_2] = { 1000.*mm, 500.*mm, 100.*mm, 10.*mm,  0.*mm}; //cut of around 430 nm so zero transmission at 430

  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX",    PhotonEnergy_RI, RefractiveIndex, nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     ABS_Energy,  ABS_LENGTH,     NUMENTRIES_2);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}

G4Material* MyMaterials::DSB_Ce()  // Nanostructured glass ceramics scintillator DSB:Ce
{
  G4double a, z, density;

  G4Element* Si = new G4Element("Silicon", "Si", z = 14., a = 28.09* g/mole);
  G4Element* Ba = new G4Element("Barium",  "Ba", z = 56., a = 137.327* g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O" , z =  8., a = 16.00* g/mole);

  G4Material* mat = new G4Material("DSB_Ce", density=4*g/cm3, 3);
  mat->AddElement(Si,1);
  mat->AddElement(Ba,1);
  mat->AddElement(O,2);

  // large band between 470 (2.64 eV) and 630 nm (1.97 eV) (mean 535 nm, 2.32)
  const G4int NUMENTRIES_1 = 5;
  G4double FAST_Energy[NUMENTRIES_1]    = {1.8*eV,1.90*eV,2.7*eV,2.88*eV,4.08*eV};
  G4double FAST_COMPONENT[NUMENTRIES_1] = {0.00,1.00,2.0,1.0,0.00};
  
  const G4int nEntries_RI = 42;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV, 
      4.5085*eV, 4.9594*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.4476, 1.4476, 1.4478, 1.4481, 
      1.4483, 1.4486, 1.4489, 1.4492, 
      1.4495, 1.4498, 1.4501, 1.4504, 
      1.4507, 1.4511, 1.4514, 1.4518, 
      1.4521, 1.4525, 1.4529, 1.4533, 
      1.4538, 1.4542, 1.4547, 1.4553, 
      1.4559, 1.4565, 1.4572, 1.4580, 
      1.4589, 1.4599, 1.4610, 1.4623, 
      1.4638, 1.4656, 1.4676, 1.4701, 
      1.4731, 1.4769, 1.4816, 1.4878, 
      1.4960, 1.5074 };
  
  const G4int NUMENTRIES_2 = 4;
  G4double ABS_Energy[NUMENTRIES_2] = { 1.0*eV, 1.84*eV, 4.08*eV, 6.26*eV };
  G4double ABS_LENGTH[NUMENTRIES_2] = { 500.*mm, 500.*mm, 500.*mm, 500.*mm }; //138 original
  //G4double Rayleigh[NUMENTRIES_2]       = { 138.*mm, 138.*mm, 138.*mm};

  G4MaterialPropertiesTable* mt = new G4MaterialPropertiesTable();
  mt->AddProperty("FASTCOMPONENT", FAST_Energy, FAST_COMPONENT, NUMENTRIES_1);
  mt->AddProperty("RINDEX",        PhotonEnergy_RI, RefractiveIndex, nEntries_RI);
  mt->AddProperty("ABSLENGTH",     ABS_Energy,  ABS_LENGTH,     NUMENTRIES_2);
  //mt->AddProperty("RAYLEIGH",      ABS_Energy,  Rayleigh,     NUMENTRIES_2);

  mt->AddConstProperty("SCINTILLATIONYIELD",1800/MeV);
  mt->AddConstProperty("RESOLUTIONSCALE",8.5);
  mt->AddConstProperty("FASTTIMECONSTANT",50.*ns);
  mt->AddConstProperty("YIELDRATIO",1.0);
  mt->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);

  mat->SetMaterialPropertiesTable(mt);
  return mat;
}


G4Material* MyMaterials::LuAG_undoped() // Lutetium Aluminum Garnet - undoped
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z=8.,  a= 16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71., a=174.97*g/mole);
  G4Element* Al = new G4Element("Aluminum", "Al", z=13., a= 28.09*g/mole);
  
  G4Material* mat = new G4Material("LuAG_undoped", density=6.7*g/cm3,3);
  mat->AddElement(Lu,3);
  mat->AddElement(Al,5);
  mat->AddElement(O,12);
  /*  
  const G4int nEntries_RI = 40;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV};//, 4.5085*eV, 4.9594*eV, 5.5*eV   , 6.5*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.8210, 1.8212, 1.8215, 1.8219, 
      1.8223, 1.8227, 1.8231, 1.8236, 
      1.8240, 1.8245, 1.8250, 1.8255, 
      1.8261, 1.8266, 1.8272, 1.8279, 
      1.8285, 1.8293, 1.8300, 1.8308, 
      1.8317, 1.8327, 1.8338, 1.8349, 
      1.8362, 1.8376, 1.8392, 1.8410, 
      1.8430, 1.8453, 1.8479, 1.8509, 
      1.8545, 1.8587, 1.8637, 1.8699, 
      1.8774, 1.8869, 1.8991, 1.9152}; //, 1.9374, 1.9694, 1.98  , 2. };
*/  
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.86,   1.86,   1.86,    1.86,    1.86,    1.86,    1.86,   1.86};

  //G4double Rayleigh[nEntries_RI] =
  //  { 138.*mm, 138.*mm, 138.*mm };
  
  const G4int nEntries_ABS = 203;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.55388*eV, 1.55975*eV, 1.56566*eV, 1.57161*eV, 1.57761*eV, 1.58365*eV, 1.58974*eV, 1.59588*eV, 1.60207*eV, 1.6083*eV, 1.61458*eV, 1.62092*eV, 1.6273*eV, 1.63373*eV, 1.64021*eV, 1.64675*eV, 1.65333*eV, 1.65997*eV, 1.66667*eV, 1.67341*eV, 1.68022*eV, 1.68707*eV, 1.69399*eV, 1.70096*eV, 1.70799*eV, 1.71508*eV, 1.72222*eV, 1.72943*eV, 1.73669*eV, 1.74402*eV, 1.75141*eV, 1.75887*eV, 1.76638*eV, 1.77396*eV, 1.78161*eV, 1.78932*eV, 1.7971*eV, 1.80495*eV, 1.81287*eV, 1.82085*eV, 1.82891*eV, 1.83704*eV, 1.84524*eV, 1.85351*eV, 1.86186*eV, 1.87029*eV, 1.87879*eV, 1.88737*eV, 1.89602*eV, 1.90476*eV, 1.91358*eV, 1.92248*eV, 1.93146*eV, 1.94053*eV, 1.94969*eV, 1.95893*eV, 1.96825*eV, 1.97767*eV, 1.98718*eV, 1.99678*eV, 2.00647*eV, 2.01626*eV, 2.02614*eV, 2.03612*eV, 2.0462*eV, 2.05638*eV, 2.06667*eV, 2.07705*eV, 2.08754*eV, 2.09814*eV, 2.10884*eV, 2.11966*eV, 2.13058*eV, 2.14162*eV, 2.15278*eV, 2.16405*eV, 2.17544*eV, 2.18695*eV, 2.19858*eV, 2.21034*eV, 2.22222*eV, 2.23423*eV, 2.24638*eV, 2.25865*eV, 2.27106*eV, 2.28361*eV, 2.2963*eV, 2.30912*eV, 2.3221*eV, 2.33522*eV, 2.34848*eV, 2.3619*eV, 2.37548*eV, 2.38921*eV, 2.4031*eV, 2.41715*eV, 2.43137*eV, 2.44576*eV, 2.46032*eV, 2.47505*eV, 2.48996*eV, 2.50505*eV, 2.52033*eV, 2.53579*eV, 2.55144*eV, 2.56729*eV, 2.58333*eV, 2.59958*eV, 2.61603*eV, 2.6327*eV, 2.64957*eV, 2.66667*eV, 2.68398*eV, 2.70153*eV, 2.7193*eV, 2.73731*eV, 2.75556*eV, 2.77405*eV, 2.79279*eV, 2.81179*eV, 2.83105*eV, 2.85057*eV, 2.87037*eV, 2.89044*eV, 2.9108*eV, 2.93144*eV, 2.95238*eV, 2.97362*eV, 2.99517*eV, 3.01703*eV, 3.03922*eV, 3.06173*eV, 3.08458*eV, 3.10777*eV, 3.13131*eV, 3.15522*eV, 3.17949*eV, 3.20413*eV, 3.22917*eV, 3.25459*eV, 3.28042*eV, 3.30667*eV, 3.33333*eV, 3.36043*eV, 3.38798*eV, 3.41598*eV, 3.44444*eV, 3.47339*eV, 3.50282*eV, 3.53276*eV, 3.56322*eV, 3.5942*eV, 3.62573*eV, 3.65782*eV, 3.69048*eV, 3.72372*eV, 3.75758*eV, 3.79205*eV, 3.82716*eV, 3.86293*eV, 3.89937*eV, 3.93651*eV, 3.97436*eV, 4.01294*eV, 4.05229*eV, 4.09241*eV, 4.13333*eV, 4.17508*eV, 4.21769*eV, 4.26117*eV, 4.30556*eV, 4.35088*eV, 4.39716*eV, 4.44444*eV, 4.49275*eV, 4.54212*eV, 4.59259*eV, 4.64419*eV, 4.69697*eV, 4.75096*eV, 4.8062*eV, 4.86275*eV, 4.92063*eV, 4.97992*eV, 5.04065*eV, 5.10288*eV, 5.16667*eV, 5.23207*eV, 5.29915*eV, 5.36797*eV, 5.4386*eV, 5.51111*eV, 5.58559*eV, 5.6621*eV, 5.74074*eV, 5.8216*eV, 5.90476*eV, 5.99034*eV, 6.07843*eV, 6.16915*eV, 6.26263*eV, 6.35897*eV, 6.45833*eV};


  G4double Absorption[nEntries_ABS] =
    { 1.000*m, 1.000*m, 01.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 0.95*m, 0.90*m, 0.85*m, 0.80*m, 0.75*m, 0.70*m, 0.65*m, 0.60*m, 0.55*m, 0.55*m, 0.54*m, 0.53653*m, 0.52467*m, 0.513653*m, 0.509889*m, 0.508994*m, 0.508051*m, 0.505419*m, 0.50765*m, 0.504605*m, 0.495435*m, 0.49046*m, 0.49047*m, 0.49327*m, 0.485816*m, 0.480906*m, 0.485679*m, 0.482268*m, 0.485718*m, 0.486185*m, 0.488607*m, 0.480324*m, 0.484148*m, 0.484234*m, 0.475095*m, 0.470936*m, 0.47703*m, 0.47404*m, 0.47822*m, 0.47574*m, 0.466905*m, 0.468739*m, 0.460301*m, 0.46061*m, 0.462486*m, 0.468237*m, 0.456888*m, 0.457251*m, 0.450251*m, 0.454983*m, 0.458911*m, 0.449242*m, 0.443889*m, 0.430268*m, 0.428066*m, 0.422729*m, 0.427122*m, 0.425943*m, 0.419851*m, 0.414304*m, 0.416655*m, 0.41294*m, 0.417711*m, 0.418304*m, 0.412927*m, 0.413118*m, 0.414686*m, 0.415758*m, 0.410793*m, 0.414454*m, 0.417264*m, 0.416738*m, 0.40941*m, 0.403909*m, 0.406116*m, 0.405126*m, 0.406951*m, 0.401353*m, 0.40883*m, 0.402408*m, 0.408617*m, 0.40139*m, 0.4021*m, 0.403939*m, 0.394537*m, 0.39998*m, 0.396536*m, 0.387493*m, 0.38305*m, 0.389145*m, 0.382008*m, 0.386039*m, 0.38835*m, 0.374416*m, 0.37622*m, 0.372372*m, 0.355058*m, 0.344739*m, 0.336844*m, 0.320373*m, 0.314119*m, 0.308814*m, 0.309788*m, 0.308062*m, 0.300838*m, 0.303496*m, 0.309025*m, 0.303451*m, 0.304254*m, 0.309881*m, 0.3053*m, 0.307437*m, 0.3006*m, 0.30426*m, 0.30292*m, 0.309979*m, 0.3009*m, 0.301757*m, 0.306717*m, 0.295276*m, 0.284917*m, 0.2726*m, 0.24962*m, 0.237705*m, 0.228777*m, 0.228844*m, 0.212415*m, 0.217188*m, 0.216965*m, 0.216849*m, 0.219055*m, 0.213229*m, 0.21006*m, 0.21604*m, 0.21683*m, 0.200704*m, 0.186524*m, 0.173764*m, 0.162931*m, 0.152322*m, 0.132755*m, 0.128307*m, 0.113216*m, 0.1077*m, 0.0993899*m, 0.0898436*m, 0.083731*m, 0.0787301*m, 0.0735095*m, 0.0691203*m, 0.0664335*m, 0.0634703*m, 0.0617119*m, 0.0595375*m, 0.0574144*m, 0.0550354*m, 0.0517719*m, 0.0488715*m, 0.0454854*m, 0.0425021*m, 0.0381949*m, 0.0340453*m, 0.0302877*m, 0.0276678*m, 0.0263617*m, 0.0248874*m, 0.0231268*m, 0.0212017*m, 0.0190633*m, 0.01683*m, 0.0145831*m, 0.0125656*m, 0.0107907*m, 0.00928679*m, 0.00804274*m, 0.00706246*m, 0.00634986*m, 0.00583677*m, 0.00545425*m, 0.00515469*m, 0.00493076*m, 0.00474148*m, 0.00457494*m, 0.00441682*m, 0.00417158*m, 0.00384893*m, 0.00344919*m, 0.00290348*m, 0.00219676*m};


  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  //myMPT->AddProperty("RAYLEIGH",      PhotonEnergy_RI,   Rayleigh,        nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_ABS);
//  myMPT->AddConstProperty("SCINTILLATIONYIELD",0/MeV);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}


G4Material* MyMaterials::LuAG_Ce() // Lutetium Aluminum Garnet - Ce-doped
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z=8.,  a= 16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71., a=174.97*g/mole);
  G4Element* Al = new G4Element("Aluminum", "Al", z=13., a= 28.09*g/mole);
  
  G4Material* mat = new G4Material("LuAG_Ce", density=6.7*g/cm3,3);
  mat->AddElement(Lu,3);
  mat->AddElement(Al,5);
  mat->AddElement(O,12);
  
  const G4int nEntries_FAST = 458;
  G4double PhotonEnergy_FAST[nEntries_FAST] =
    { 1.82487*eV, 1.82622*eV, 1.82756*eV, 1.82891*eV, 1.83026*eV, 1.83161*eV, 1.83296*eV, 1.83432*eV, 1.83568*eV, 1.83704*eV,
      1.8384*eV, 1.83976*eV, 1.84113*eV, 1.8425*eV, 1.84387*eV, 1.84524*eV, 1.84661*eV, 1.84799*eV, 1.84937*eV, 1.85075*eV,
      1.85213*eV, 1.85351*eV, 1.8549*eV, 1.85629*eV, 1.85768*eV, 1.85907*eV, 1.86047*eV, 1.86186*eV, 1.86326*eV, 1.86466*eV,
      1.86606*eV, 1.86747*eV, 1.86888*eV, 1.87029*eV, 1.8717*eV, 1.87311*eV, 1.87453*eV, 1.87595*eV, 1.87737*eV, 1.87879*eV,
      1.88021*eV, 1.88164*eV, 1.88307*eV, 1.8845*eV, 1.88593*eV, 1.88737*eV, 1.8888*eV, 1.89024*eV, 1.89169*eV, 1.89313*eV,
      1.89458*eV, 1.89602*eV, 1.89748*eV, 1.89893*eV, 1.90038*eV, 1.90184*eV, 1.9033*eV, 1.90476*eV, 1.90623*eV, 1.90769*eV,
      1.90916*eV, 1.91063*eV, 1.9121*eV, 1.91358*eV, 1.91506*eV, 1.91654*eV, 1.91802*eV, 1.9195*eV, 1.92099*eV, 1.92248*eV,
      1.92397*eV, 1.92547*eV, 1.92696*eV, 1.92846*eV, 1.92996*eV, 1.93146*eV, 1.93297*eV, 1.93448*eV, 1.93599*eV, 1.9375*eV,
      1.93901*eV, 1.94053*eV, 1.94205*eV, 1.94357*eV, 1.9451*eV, 1.94662*eV, 1.94815*eV, 1.94969*eV, 1.95122*eV, 1.95276*eV,
      1.95429*eV, 1.95584*eV, 1.95738*eV, 1.95893*eV, 1.96047*eV, 1.96203*eV, 1.96358*eV, 1.96513*eV, 1.96669*eV, 1.96825*eV,
      1.96982*eV, 1.97138*eV, 1.97295*eV, 1.97452*eV, 1.9761*eV, 1.97767*eV, 1.97925*eV, 1.98083*eV, 1.98241*eV, 1.984*eV,
      1.98559*eV, 1.98718*eV, 1.98877*eV, 1.99037*eV, 1.99197*eV, 1.99357*eV, 1.99517*eV, 1.99678*eV, 1.99839*eV, 2*eV,
      2.00161*eV, 2.00323*eV, 2.00485*eV, 2.00647*eV, 2.0081*eV, 2.00972*eV, 2.01135*eV, 2.01299*eV, 2.01462*eV, 2.01626*eV,
      2.0179*eV, 2.01954*eV, 2.02119*eV, 2.02284*eV, 2.02449*eV, 2.02614*eV, 2.0278*eV, 2.02946*eV, 2.03112*eV, 2.03279*eV,
      2.03445*eV, 2.03612*eV, 2.0378*eV, 2.03947*eV, 2.04115*eV, 2.04283*eV, 2.04452*eV, 2.0462*eV, 2.04789*eV, 2.04959*eV,
      2.05128*eV, 2.05298*eV, 2.05468*eV, 2.05638*eV, 2.05809*eV, 2.0598*eV, 2.06151*eV, 2.06323*eV, 2.06495*eV, 2.06667*eV,
      2.06839*eV, 2.07012*eV, 2.07185*eV, 2.07358*eV, 2.07531*eV, 2.07705*eV, 2.07879*eV, 2.08054*eV, 2.08228*eV, 2.08403*eV,
      2.08579*eV, 2.08754*eV, 2.0893*eV, 2.09106*eV, 2.09283*eV, 2.09459*eV, 2.09637*eV, 2.09814*eV, 2.09992*eV, 2.10169*eV,
      2.10348*eV, 2.10526*eV, 2.10705*eV, 2.10884*eV, 2.11064*eV, 2.11244*eV, 2.11424*eV, 2.11604*eV, 2.11785*eV, 2.11966*eV,
      2.12147*eV, 2.12329*eV, 2.12511*eV, 2.12693*eV, 2.12876*eV, 2.13058*eV, 2.13242*eV, 2.13425*eV, 2.13609*eV, 2.13793*eV,
      2.13978*eV, 2.14162*eV, 2.14347*eV, 2.14533*eV, 2.14719*eV, 2.14905*eV, 2.15091*eV, 2.15278*eV, 2.15465*eV, 2.15652*eV,
      2.1584*eV, 2.16028*eV, 2.16216*eV, 2.16405*eV, 2.16594*eV, 2.16783*eV, 2.16973*eV, 2.17163*eV, 2.17353*eV, 2.17544*eV,
      2.17735*eV, 2.17926*eV, 2.18118*eV, 2.1831*eV, 2.18502*eV, 2.18695*eV, 2.18888*eV, 2.19081*eV, 2.19275*eV, 2.19469*eV,
      2.19663*eV, 2.19858*eV, 2.20053*eV, 2.20249*eV, 2.20444*eV, 2.20641*eV, 2.20837*eV, 2.21034*eV, 2.21231*eV, 2.21429*eV,
      2.21626*eV, 2.21825*eV, 2.22023*eV, 2.22222*eV, 2.22422*eV, 2.22621*eV, 2.22821*eV, 2.23022*eV, 2.23222*eV, 2.23423*eV,
      2.23625*eV, 2.23827*eV, 2.24029*eV, 2.24231*eV, 2.24434*eV, 2.24638*eV, 2.24841*eV, 2.25045*eV, 2.2525*eV, 2.25455*eV,
      2.2566*eV, 2.25865*eV, 2.26071*eV, 2.26277*eV, 2.26484*eV, 2.26691*eV, 2.26898*eV, 2.27106*eV, 2.27314*eV, 2.27523*eV,
      2.27732*eV, 2.27941*eV, 2.28151*eV, 2.28361*eV, 2.28571*eV, 2.28782*eV, 2.28994*eV, 2.29205*eV, 2.29417*eV, 2.2963*eV,
      2.29842*eV, 2.30056*eV, 2.30269*eV, 2.30483*eV, 2.30698*eV, 2.30912*eV, 2.31128*eV, 2.31343*eV, 2.31559*eV, 2.31776*eV,
      2.31993*eV, 2.3221*eV, 2.32427*eV, 2.32645*eV, 2.32864*eV, 2.33083*eV, 2.33302*eV, 2.33522*eV, 2.33742*eV, 2.33962*eV,
      2.34183*eV, 2.34405*eV, 2.34626*eV, 2.34848*eV, 2.35071*eV, 2.35294*eV, 2.35518*eV, 2.35741*eV, 2.35966*eV, 2.3619*eV,
      2.36416*eV, 2.36641*eV, 2.36867*eV, 2.37094*eV, 2.37321*eV, 2.37548*eV, 2.37776*eV, 2.38004*eV, 2.38232*eV, 2.38462*eV,
      2.38691*eV, 2.38921*eV, 2.39151*eV, 2.39382*eV, 2.39614*eV, 2.39845*eV, 2.40077*eV, 2.4031*eV, 2.40543*eV, 2.40777*eV,
      2.41011*eV, 2.41245*eV, 2.4148*eV, 2.41715*eV, 2.41951*eV, 2.42188*eV, 2.42424*eV, 2.42661*eV, 2.42899*eV, 2.43137*eV,
      2.43376*eV, 2.43615*eV, 2.43854*eV, 2.44094*eV, 2.44335*eV, 2.44576*eV, 2.44817*eV, 2.45059*eV, 2.45302*eV, 2.45545*eV,
      2.45788*eV, 2.46032*eV, 2.46276*eV, 2.46521*eV, 2.46766*eV, 2.47012*eV, 2.47258*eV, 2.47505*eV, 2.47752*eV, 2.48*eV,
      2.48248*eV, 2.48497*eV, 2.48746*eV, 2.48996*eV, 2.49246*eV, 2.49497*eV, 2.49748*eV, 2.5*eV, 2.50252*eV, 2.50505*eV,
      2.50758*eV, 2.51012*eV, 2.51266*eV, 2.51521*eV, 2.51777*eV, 2.52033*eV, 2.52289*eV, 2.52546*eV, 2.52803*eV, 2.53061*eV,
      2.5332*eV, 2.53579*eV, 2.53838*eV, 2.54098*eV, 2.54359*eV, 2.5462*eV, 2.54882*eV, 2.55144*eV, 2.55407*eV, 2.5567*eV,
      2.55934*eV, 2.56198*eV, 2.56463*eV, 2.56729*eV, 2.56995*eV, 2.57261*eV, 2.57529*eV, 2.57796*eV, 2.58065*eV, 2.58333*eV,
      2.58603*eV, 2.58873*eV, 2.59143*eV, 2.59414*eV, 2.59686*eV, 2.59958*eV, 2.60231*eV, 2.60504*eV, 2.60778*eV, 2.61053*eV,
      2.61328*eV, 2.61603*eV, 2.6188*eV, 2.62156*eV, 2.62434*eV, 2.62712*eV, 2.6299*eV, 2.6327*eV, 2.63549*eV, 2.6383*eV,
      2.64111*eV, 2.64392*eV, 2.64674*eV, 2.64957*eV, 2.65241*eV, 2.65525*eV, 2.65809*eV, 2.66094*eV, 2.6638*eV, 2.66667*eV,
      2.66954*eV, 2.67241*eV, 2.6753*eV, 2.67819*eV, 2.68108*eV, 2.68398*eV, 2.68689*eV, 2.6898*eV, 2.69273*eV, 2.69565*eV,
      2.69859*eV, 2.70153*eV, 2.70447*eV, 2.70742*eV, 2.71038*eV, 2.71335*eV, 2.71632*eV, 2.7193*eV, 2.72228*eV, 2.72527*eV,
      2.72827*eV, 2.73128*eV, 2.73429*eV, 2.73731*eV, 2.74033*eV, 2.74336*eV, 2.7464*eV, 2.74945*eV };
  G4double FastComponent[nEntries_FAST] =
    { 5.81332e-05, 6.44431e-05, 5.14981e-05, 5.53578e-05, 7.63256e-05, 7.53282e-05, 7.58269e-05, 8.97693e-05, 7.76917e-05, 7.38103e-05,
      7.78435e-05, 7.09481e-05, 7.49162e-05, 8.77528e-05, 8.86852e-05, 9.01596e-05, 7.3355e-05, 8.61916e-05, 8.31125e-05, 9.63177e-05,
      9.64045e-05, 8.96609e-05, 0.000118934, 0.000122446, 0.000112017, 8.10092e-05, 9.10487e-05, 9.54287e-05, 0.000102975, 0.000102996,
      0.00010833, 9.44529e-05, 9.82259e-05, 0.000117372, 0.000121601, 0.00011206, 0.000123183, 0.000126371, 0.000114987, 0.000121687,
      0.00011065, 0.000131879, 0.000124766, 0.000119606, 0.000146146, 0.000145279, 0.000141441, 0.000148553, 0.000156012, 0.000149746,
      0.000163168, 0.000161043, 0.000174898, 0.000182661, 0.000175918, 0.000175939, 0.000171169, 0.000159807, 0.0001726, 0.000178866,
      0.000175028, 0.000190836, 0.000208768, 0.000179539, 0.000198165, 0.000197644, 0.000199509, 0.000202545, 0.000218005, 0.000208031,
      0.000212584, 0.000219848, 0.000234961, 0.000249945, 0.000232078, 0.000224814, 0.000229476, 0.000248232, 0.000290948, 0.000269568,
      0.000289062, 0.000288346, 0.000276572, 0.000287999, 0.000306235, 0.000306452, 0.000293724, 0.000325186, 0.000335356, 0.000332277,
      0.000323885, 0.000335464, 0.000335724, 0.000372413, 0.000366406, 0.000352919, 0.000341297, 0.000400384, 0.00038796, 0.000373518,
      0.000396785, 0.000419986, 0.00042712, 0.000413654, 0.000429158, 0.000443621, 0.000482629, 0.000489546, 0.00047736, 0.000457151,
      0.000534973, 0.000505028, 0.000518277, 0.000507934, 0.000520879, 0.000552884, 0.00054859, 0.00057396, 0.000561015, 0.000557329,
      0.000603645, 0.000608111, 0.000628646, 0.000624894, 0.000634196, 0.000647055, 0.000670148, 0.000680122, 0.000679211, 0.000696731,
      0.000703887, 0.000728888, 0.000708722, 0.000749509, 0.000780863, 0.000772819, 0.000771908, 0.000782294, 0.000797256, 0.000792767,
      0.000837197, 0.00086081, 0.000856668, 0.000894267, 0.000905629, 0.000906323, 0.00097068, 0.00095175, 0.000955696, 0.000969986,
      0.000980069, 0.00103508, 0.00108445, 0.00109566, 0.00108985, 0.00108554, 0.00116115, 0.00116232, 0.00123272, 0.00118337,
      0.0012059, 0.00121104, 0.00122034, 0.00134806, 0.00132759, 0.00131655, 0.00134255, 0.0014004, 0.00139422, 0.00140613,
      0.00148002, 0.00146782, 0.00151511, 0.00152601, 0.00156979, 0.00156765, 0.00161945, 0.00161797, 0.00164276, 0.00167585,
      0.00163298, 0.00169692, 0.00173064, 0.00185434, 0.0018602, 0.00183116, 0.00184591, 0.00187262, 0.00185005, 0.00187863,
      0.00193908, 0.00196593, 0.00204381, 0.00207335, 0.00214098, 0.00216439, 0.00214946, 0.00215717, 0.00218378, 0.00220102,
      0.00226733, 0.00225952, 0.00232741, 0.0023407, 0.00237913, 0.0023986, 0.00246035, 0.00246658, 0.00256504, 0.00255051,
      0.00260585, 0.00261381, 0.00263094, 0.00263575, 0.00272168, 0.00272411, 0.00270891, 0.00276706, 0.00281004, 0.00290874,
      0.00298958, 0.00287925, 0.00292917, 0.00294574, 0.00308388, 0.00300632, 0.00300664, 0.00304351, 0.00310511, 0.00314704,
      0.00307382, 0.0031967, 0.00324223, 0.0032804, 0.00328582, 0.00328493, 0.00322873, 0.00335108, 0.00344131, 0.00348589,
      0.00347806, 0.00350257, 0.00354257, 0.00361777, 0.00364609, 0.00357236, 0.00361374, 0.0036504, 0.00367545, 0.00370214,
      0.00372157, 0.00380769, 0.00386073, 0.00378243, 0.0038225, 0.00388172, 0.00388896, 0.0039206, 0.00382433, 0.0039355,
      0.00394595, 0.00403936, 0.00412466, 0.00398097, 0.00412995, 0.00410556, 0.00412607, 0.00409762, 0.00417612, 0.00419223,
      0.00418902, 0.00420858, 0.00418737, 0.00438864, 0.00428859, 0.00424058, 0.00428341, 0.00434349, 0.00426255, 0.00425739,
      0.00426099, 0.00435654, 0.00430615, 0.0043434, 0.00436442, 0.00443317, 0.00453702, 0.00458465, 0.00451718, 0.00454829,
      0.004526, 0.00444941, 0.00459784, 0.00461657, 0.00464116, 0.00463936, 0.00462128, 0.00464869, 0.00472074, 0.00464255,
      0.00463531, 0.00464357, 0.00472312, 0.00471226, 0.00472876, 0.00475929, 0.00477558, 0.00477493, 0.00476745, 0.00480513,
      0.00488634, 0.00489984, 0.00491695, 0.0049675, 0.00488809, 0.00492643, 0.0048836, 0.00497446, 0.00506449, 0.00503294,
      0.00507216, 0.00511015, 0.00528854, 0.00508509, 0.00508214, 0.00515293, 0.00521213, 0.00535207, 0.00521807, 0.00530294,
      0.00523004, 0.00531701, 0.00543113, 0.00540544, 0.00534221, 0.00529414, 0.00536786, 0.00530663, 0.00540221, 0.0054366,
      0.00546841, 0.00534307, 0.0053802, 0.00543647, 0.00542493, 0.00540019, 0.0054354, 0.00542636, 0.00540218, 0.00539761,
      0.00546834, 0.00538748, 0.0054119, 0.00524333, 0.0052661, 0.00528475, 0.00527178, 0.00527688, 0.00527451, 0.0051803,
      0.00532525, 0.00516377, 0.00502179, 0.00498588, 0.00493792, 0.00504124, 0.00497812, 0.00480997, 0.00484698, 0.00475996,
      0.00467631, 0.0046375, 0.00452982, 0.00445893, 0.00443842, 0.00431203, 0.004365, 0.00422527, 0.00416564, 0.00407752,
      0.00394289, 0.00401572, 0.00385061, 0.00368295, 0.00359867, 0.00356532, 0.00350603, 0.00343337, 0.00328738, 0.00314308,
      0.00308744, 0.00307742, 0.00298919, 0.00291293, 0.00276594, 0.00274861, 0.00263172, 0.00244029, 0.00243841, 0.0023837,
      0.00222602, 0.00220015, 0.00206933, 0.0019725, 0.00198421, 0.00188412, 0.00176243, 0.00169384, 0.00163604, 0.00154739,
      0.00148061, 0.00135909, 0.00128767, 0.00121056, 0.00116466, 0.00113142, 0.00102363, 0.000933341, 0.000903136, 0.00086764,
      0.000834096, 0.000722275, 0.000704733, 0.000665919, 0.000615007, 0.00057151, 0.000541261, 0.000484668, 0.000465868, 0.000435923,
      0.000406347, 0.000362091, 0.000334272, 0.000302766, 0.000280562, 0.000268874, 0.000248102, 0.000231557, 0.0002133, 0.0001981,
      0.000193611, 0.000166442, 0.000154017, 0.000139056, 0.0001301, 0.00011928, 0.000122511, 0.000106704, 8.5411e-05, 8.49339e-05,
      8.20717e-05, 6.96905e-05, 6.09304e-05, 5.20402e-05, 5.76345e-05, 5.77646e-05, 4.39089e-05, 6.16243e-05 };
  /*
  const G4int nEntries_RI = 40;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV};//, 4.5085*eV, 4.9594*eV, 5.5*eV   , 6.5*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.8210, 1.8212, 1.8215, 1.8219, 
      1.8223, 1.8227, 1.8231, 1.8236, 
      1.8240, 1.8245, 1.8250, 1.8255, 
      1.8261, 1.8266, 1.8272, 1.8279, 
      1.8285, 1.8293, 1.8300, 1.8308, 
      1.8317, 1.8327, 1.8338, 1.8349, 
      1.8362, 1.8376, 1.8392, 1.8410, 
      1.8430, 1.8453, 1.8479, 1.8509, 
      1.8545, 1.8587, 1.8637, 1.8699, 
      1.8774, 1.8869, 1.8991, 1.9152}; //, 1.9374, 1.9694, 1.98  , 2. };
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.86,   1.86,   1.86,    1.86,    1.86,    1.86,    1.86,   1.86};
  
  const G4int nEntries_ABS = 89;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.78417*eV, 1.7971*eV, 1.81022*eV, 1.82353*eV, 1.83704*eV, 1.85075*eV, 1.86466*eV, 1.87879*eV, 1.89313*eV, 1.90769*eV,
      1.92248*eV, 1.9375*eV, 1.95276*eV, 1.96825*eV, 1.984*eV, 2*eV, 2.01626*eV, 2.03279*eV, 2.04959*eV, 2.06667*eV,
      2.08403*eV, 2.10169*eV, 2.11966*eV, 2.13793*eV, 2.15652*eV, 2.17544*eV, 2.19469*eV, 2.21429*eV, 2.23423*eV, 2.25455*eV,
      2.27523*eV, 2.2963*eV, 2.31776*eV, 2.33962*eV, 2.3619*eV, 2.38462*eV, 2.40777*eV, 2.43137*eV, 2.45545*eV, 2.48*eV,
      2.50505*eV, 2.53061*eV, 2.5567*eV, 2.58333*eV, 2.61053*eV, 2.6383*eV, 2.66667*eV, 2.69565*eV, 2.72527*eV, 2.75556*eV,
      2.78652*eV, 2.81818*eV, 2.85057*eV, 2.88372*eV, 2.91765*eV, 2.95238*eV, 2.98795*eV, 3.02439*eV, 3.06173*eV, 3.1*eV,
      3.13924*eV, 3.17949*eV, 3.22078*eV, 3.26316*eV, 3.30667*eV, 3.35135*eV, 3.39726*eV, 3.44444*eV, 3.49296*eV, 3.54286*eV,
      3.5942*eV, 3.64706*eV, 3.70149*eV, 3.75758*eV, 3.81538*eV, 3.875*eV, 3.93651*eV, 4*eV, 4.06557*eV, 4.13333*eV,
      4.20339*eV, 4.27586*eV, 4.35088*eV, 4.42857*eV, 4.50909*eV, 4.59259*eV, 4.67925*eV, 4.76923*eV, 4.86275*eV };
  G4double Absorption[nEntries_ABS] =
    { 1.000*m, 1.000*m, 1.0000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m,
      1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 0.90*m, 0.85*m, 0.80*m, 0.75*m,
      0.70*m, 0.65*m, 0.60*m, 0.58*m, 0.57*m, 0.56*m, 0.55*m, 0.54*m, 0.54*m, 0.54*m,
      0.531547*m, 0.53*m, 0.5251*m, 0.510*m, 0.493902*m, 0.422165*m, 0.354962*m, 0.255139*m, 0.015762*m, 0.00827965*m,
      0.00409174*m, 0.001007*m, 0.000100463*m, 0.000188694*m, 0.000636339*m, 0.0084497*m, 0.0224574*m, 0.0509883*m, 0.262914*m, 0.0571499*m,
      0.0830375*m, 0.378696*m, 0.0528428*m, 0.0661874*m, 0.0930821*m, 0.0672707*m, 0.0015385*m, 0.00067752*m, 0.0005306*m, 0.00079596*m,
      0.0177025*m, 0.0411282*m, 0.0919861*m, 0.149875*m, 0.132761*m, 0.006419*m, 0.0026548*m, 0.00092619*m, 0.00092168*m, 0.0024256*m,
      0.0089517*m, 0.0076384*m, 0.0052649*m, 0.0019203*m, 0.00082616*m, 0.00074327*m, 0.0015009*m, 0.0029903*m, 0.0040526*m, 0.0037371*m,
      0.00032887*m, 0.00021734*m, 0.00014992*m, 0.00015645*m, 0.00011908*m, 0.00010775*m, 0.00011081*m, 0.00015907*m, 0.0001793*m };
  

  const G4int nEntries_SCY = 12;
  G4double ElectronEnergy_SCY[nEntries_SCY] =
    { 0.000*MeV, 0.015*MeV, 0.020*MeV, 0.030*MeV,
      0.040*MeV, 0.060*MeV, 0.080*MeV, 0.090*MeV,
      0.105*MeV, 0.300*MeV, 0.500*MeV, 1.000*MeV };
  G4double ScintilYield[nEntries_SCY] =
    { 0.10, 0.46, 0.60, 0.68,
      0.74, 0.80, 0.82, 0.84,
      0.87, 0.96, 0.98, 1.00 };
  for(int i = 0; i < nEntries_SCY; i++)
    ScintilYield[i] = 15000.0*MeV*ScintilYield[i]*ElectronEnergy_SCY[i];
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("FASTCOMPONENT", PhotonEnergy_FAST, FastComponent,   nEntries_FAST);
  myMPT->AddProperty("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  //myMPT->AddProperty("RAYLEIGH",      PhotonEnergy_RI,   Rayleigh,        nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_ABS);
  //myMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", ElectronEnergy_SCY, ScintilYield, nEntries_SCY);
  myMPT->AddConstProperty("SCINTILLATIONYIELD",23000/MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE",8.5);
  myMPT->AddConstProperty("FASTTIMECONSTANT",55.*ns);
  myMPT->AddConstProperty("SLOWTIMECONSTANT",300.*ns);
  myMPT->AddConstProperty("YIELDRATIO", 0.5);
  myMPT->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}

G4Material* MyMaterials::LuAG_Ce_Mg() // Lutetium Aluminum Garnet - Ce-doped
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z=8.,  a= 16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71., a=174.97*g/mole);
  G4Element* Al = new G4Element("Aluminum", "Al", z=13., a= 28.09*g/mole);
  
  G4Material* mat = new G4Material("LuAG_Ce_Mg", density=6.7*g/cm3,3);
  mat->AddElement(Lu,3);
  mat->AddElement(Al,5);
  mat->AddElement(O,12);
  
  const G4int nEntries_FAST = 458;
  G4double PhotonEnergy_FAST[nEntries_FAST] =
    { 1.82487*eV, 1.82622*eV, 1.82756*eV, 1.82891*eV, 1.83026*eV, 1.83161*eV, 1.83296*eV, 1.83432*eV, 1.83568*eV, 1.83704*eV,
      1.8384*eV, 1.83976*eV, 1.84113*eV, 1.8425*eV, 1.84387*eV, 1.84524*eV, 1.84661*eV, 1.84799*eV, 1.84937*eV, 1.85075*eV,
      1.85213*eV, 1.85351*eV, 1.8549*eV, 1.85629*eV, 1.85768*eV, 1.85907*eV, 1.86047*eV, 1.86186*eV, 1.86326*eV, 1.86466*eV,
      1.86606*eV, 1.86747*eV, 1.86888*eV, 1.87029*eV, 1.8717*eV, 1.87311*eV, 1.87453*eV, 1.87595*eV, 1.87737*eV, 1.87879*eV,
      1.88021*eV, 1.88164*eV, 1.88307*eV, 1.8845*eV, 1.88593*eV, 1.88737*eV, 1.8888*eV, 1.89024*eV, 1.89169*eV, 1.89313*eV,
      1.89458*eV, 1.89602*eV, 1.89748*eV, 1.89893*eV, 1.90038*eV, 1.90184*eV, 1.9033*eV, 1.90476*eV, 1.90623*eV, 1.90769*eV,
      1.90916*eV, 1.91063*eV, 1.9121*eV, 1.91358*eV, 1.91506*eV, 1.91654*eV, 1.91802*eV, 1.9195*eV, 1.92099*eV, 1.92248*eV,
      1.92397*eV, 1.92547*eV, 1.92696*eV, 1.92846*eV, 1.92996*eV, 1.93146*eV, 1.93297*eV, 1.93448*eV, 1.93599*eV, 1.9375*eV,
      1.93901*eV, 1.94053*eV, 1.94205*eV, 1.94357*eV, 1.9451*eV, 1.94662*eV, 1.94815*eV, 1.94969*eV, 1.95122*eV, 1.95276*eV,
      1.95429*eV, 1.95584*eV, 1.95738*eV, 1.95893*eV, 1.96047*eV, 1.96203*eV, 1.96358*eV, 1.96513*eV, 1.96669*eV, 1.96825*eV,
      1.96982*eV, 1.97138*eV, 1.97295*eV, 1.97452*eV, 1.9761*eV, 1.97767*eV, 1.97925*eV, 1.98083*eV, 1.98241*eV, 1.984*eV,
      1.98559*eV, 1.98718*eV, 1.98877*eV, 1.99037*eV, 1.99197*eV, 1.99357*eV, 1.99517*eV, 1.99678*eV, 1.99839*eV, 2*eV,
      2.00161*eV, 2.00323*eV, 2.00485*eV, 2.00647*eV, 2.0081*eV, 2.00972*eV, 2.01135*eV, 2.01299*eV, 2.01462*eV, 2.01626*eV,
      2.0179*eV, 2.01954*eV, 2.02119*eV, 2.02284*eV, 2.02449*eV, 2.02614*eV, 2.0278*eV, 2.02946*eV, 2.03112*eV, 2.03279*eV,
      2.03445*eV, 2.03612*eV, 2.0378*eV, 2.03947*eV, 2.04115*eV, 2.04283*eV, 2.04452*eV, 2.0462*eV, 2.04789*eV, 2.04959*eV,
      2.05128*eV, 2.05298*eV, 2.05468*eV, 2.05638*eV, 2.05809*eV, 2.0598*eV, 2.06151*eV, 2.06323*eV, 2.06495*eV, 2.06667*eV,
      2.06839*eV, 2.07012*eV, 2.07185*eV, 2.07358*eV, 2.07531*eV, 2.07705*eV, 2.07879*eV, 2.08054*eV, 2.08228*eV, 2.08403*eV,
      2.08579*eV, 2.08754*eV, 2.0893*eV, 2.09106*eV, 2.09283*eV, 2.09459*eV, 2.09637*eV, 2.09814*eV, 2.09992*eV, 2.10169*eV,
      2.10348*eV, 2.10526*eV, 2.10705*eV, 2.10884*eV, 2.11064*eV, 2.11244*eV, 2.11424*eV, 2.11604*eV, 2.11785*eV, 2.11966*eV,
      2.12147*eV, 2.12329*eV, 2.12511*eV, 2.12693*eV, 2.12876*eV, 2.13058*eV, 2.13242*eV, 2.13425*eV, 2.13609*eV, 2.13793*eV,
      2.13978*eV, 2.14162*eV, 2.14347*eV, 2.14533*eV, 2.14719*eV, 2.14905*eV, 2.15091*eV, 2.15278*eV, 2.15465*eV, 2.15652*eV,
      2.1584*eV, 2.16028*eV, 2.16216*eV, 2.16405*eV, 2.16594*eV, 2.16783*eV, 2.16973*eV, 2.17163*eV, 2.17353*eV, 2.17544*eV,
      2.17735*eV, 2.17926*eV, 2.18118*eV, 2.1831*eV, 2.18502*eV, 2.18695*eV, 2.18888*eV, 2.19081*eV, 2.19275*eV, 2.19469*eV,
      2.19663*eV, 2.19858*eV, 2.20053*eV, 2.20249*eV, 2.20444*eV, 2.20641*eV, 2.20837*eV, 2.21034*eV, 2.21231*eV, 2.21429*eV,
      2.21626*eV, 2.21825*eV, 2.22023*eV, 2.22222*eV, 2.22422*eV, 2.22621*eV, 2.22821*eV, 2.23022*eV, 2.23222*eV, 2.23423*eV,
      2.23625*eV, 2.23827*eV, 2.24029*eV, 2.24231*eV, 2.24434*eV, 2.24638*eV, 2.24841*eV, 2.25045*eV, 2.2525*eV, 2.25455*eV,
      2.2566*eV, 2.25865*eV, 2.26071*eV, 2.26277*eV, 2.26484*eV, 2.26691*eV, 2.26898*eV, 2.27106*eV, 2.27314*eV, 2.27523*eV,
      2.27732*eV, 2.27941*eV, 2.28151*eV, 2.28361*eV, 2.28571*eV, 2.28782*eV, 2.28994*eV, 2.29205*eV, 2.29417*eV, 2.2963*eV,
      2.29842*eV, 2.30056*eV, 2.30269*eV, 2.30483*eV, 2.30698*eV, 2.30912*eV, 2.31128*eV, 2.31343*eV, 2.31559*eV, 2.31776*eV,
      2.31993*eV, 2.3221*eV, 2.32427*eV, 2.32645*eV, 2.32864*eV, 2.33083*eV, 2.33302*eV, 2.33522*eV, 2.33742*eV, 2.33962*eV,
      2.34183*eV, 2.34405*eV, 2.34626*eV, 2.34848*eV, 2.35071*eV, 2.35294*eV, 2.35518*eV, 2.35741*eV, 2.35966*eV, 2.3619*eV,
      2.36416*eV, 2.36641*eV, 2.36867*eV, 2.37094*eV, 2.37321*eV, 2.37548*eV, 2.37776*eV, 2.38004*eV, 2.38232*eV, 2.38462*eV,
      2.38691*eV, 2.38921*eV, 2.39151*eV, 2.39382*eV, 2.39614*eV, 2.39845*eV, 2.40077*eV, 2.4031*eV, 2.40543*eV, 2.40777*eV,
      2.41011*eV, 2.41245*eV, 2.4148*eV, 2.41715*eV, 2.41951*eV, 2.42188*eV, 2.42424*eV, 2.42661*eV, 2.42899*eV, 2.43137*eV,
      2.43376*eV, 2.43615*eV, 2.43854*eV, 2.44094*eV, 2.44335*eV, 2.44576*eV, 2.44817*eV, 2.45059*eV, 2.45302*eV, 2.45545*eV,
      2.45788*eV, 2.46032*eV, 2.46276*eV, 2.46521*eV, 2.46766*eV, 2.47012*eV, 2.47258*eV, 2.47505*eV, 2.47752*eV, 2.48*eV,
      2.48248*eV, 2.48497*eV, 2.48746*eV, 2.48996*eV, 2.49246*eV, 2.49497*eV, 2.49748*eV, 2.5*eV, 2.50252*eV, 2.50505*eV,
      2.50758*eV, 2.51012*eV, 2.51266*eV, 2.51521*eV, 2.51777*eV, 2.52033*eV, 2.52289*eV, 2.52546*eV, 2.52803*eV, 2.53061*eV,
      2.5332*eV, 2.53579*eV, 2.53838*eV, 2.54098*eV, 2.54359*eV, 2.5462*eV, 2.54882*eV, 2.55144*eV, 2.55407*eV, 2.5567*eV,
      2.55934*eV, 2.56198*eV, 2.56463*eV, 2.56729*eV, 2.56995*eV, 2.57261*eV, 2.57529*eV, 2.57796*eV, 2.58065*eV, 2.58333*eV,
      2.58603*eV, 2.58873*eV, 2.59143*eV, 2.59414*eV, 2.59686*eV, 2.59958*eV, 2.60231*eV, 2.60504*eV, 2.60778*eV, 2.61053*eV,
      2.61328*eV, 2.61603*eV, 2.6188*eV, 2.62156*eV, 2.62434*eV, 2.62712*eV, 2.6299*eV, 2.6327*eV, 2.63549*eV, 2.6383*eV,
      2.64111*eV, 2.64392*eV, 2.64674*eV, 2.64957*eV, 2.65241*eV, 2.65525*eV, 2.65809*eV, 2.66094*eV, 2.6638*eV, 2.66667*eV,
      2.66954*eV, 2.67241*eV, 2.6753*eV, 2.67819*eV, 2.68108*eV, 2.68398*eV, 2.68689*eV, 2.6898*eV, 2.69273*eV, 2.69565*eV,
      2.69859*eV, 2.70153*eV, 2.70447*eV, 2.70742*eV, 2.71038*eV, 2.71335*eV, 2.71632*eV, 2.7193*eV, 2.72228*eV, 2.72527*eV,
      2.72827*eV, 2.73128*eV, 2.73429*eV, 2.73731*eV, 2.74033*eV, 2.74336*eV, 2.7464*eV, 2.74945*eV };
  G4double FastComponent[nEntries_FAST] =
    { 5.81332e-05, 6.44431e-05, 5.14981e-05, 5.53578e-05, 7.63256e-05, 7.53282e-05, 7.58269e-05, 8.97693e-05, 7.76917e-05, 7.38103e-05,
      7.78435e-05, 7.09481e-05, 7.49162e-05, 8.77528e-05, 8.86852e-05, 9.01596e-05, 7.3355e-05, 8.61916e-05, 8.31125e-05, 9.63177e-05,
      9.64045e-05, 8.96609e-05, 0.000118934, 0.000122446, 0.000112017, 8.10092e-05, 9.10487e-05, 9.54287e-05, 0.000102975, 0.000102996,
      0.00010833, 9.44529e-05, 9.82259e-05, 0.000117372, 0.000121601, 0.00011206, 0.000123183, 0.000126371, 0.000114987, 0.000121687,
      0.00011065, 0.000131879, 0.000124766, 0.000119606, 0.000146146, 0.000145279, 0.000141441, 0.000148553, 0.000156012, 0.000149746,
      0.000163168, 0.000161043, 0.000174898, 0.000182661, 0.000175918, 0.000175939, 0.000171169, 0.000159807, 0.0001726, 0.000178866,
      0.000175028, 0.000190836, 0.000208768, 0.000179539, 0.000198165, 0.000197644, 0.000199509, 0.000202545, 0.000218005, 0.000208031,
      0.000212584, 0.000219848, 0.000234961, 0.000249945, 0.000232078, 0.000224814, 0.000229476, 0.000248232, 0.000290948, 0.000269568,
      0.000289062, 0.000288346, 0.000276572, 0.000287999, 0.000306235, 0.000306452, 0.000293724, 0.000325186, 0.000335356, 0.000332277,
      0.000323885, 0.000335464, 0.000335724, 0.000372413, 0.000366406, 0.000352919, 0.000341297, 0.000400384, 0.00038796, 0.000373518,
      0.000396785, 0.000419986, 0.00042712, 0.000413654, 0.000429158, 0.000443621, 0.000482629, 0.000489546, 0.00047736, 0.000457151,
      0.000534973, 0.000505028, 0.000518277, 0.000507934, 0.000520879, 0.000552884, 0.00054859, 0.00057396, 0.000561015, 0.000557329,
      0.000603645, 0.000608111, 0.000628646, 0.000624894, 0.000634196, 0.000647055, 0.000670148, 0.000680122, 0.000679211, 0.000696731,
      0.000703887, 0.000728888, 0.000708722, 0.000749509, 0.000780863, 0.000772819, 0.000771908, 0.000782294, 0.000797256, 0.000792767,
      0.000837197, 0.00086081, 0.000856668, 0.000894267, 0.000905629, 0.000906323, 0.00097068, 0.00095175, 0.000955696, 0.000969986,
      0.000980069, 0.00103508, 0.00108445, 0.00109566, 0.00108985, 0.00108554, 0.00116115, 0.00116232, 0.00123272, 0.00118337,
      0.0012059, 0.00121104, 0.00122034, 0.00134806, 0.00132759, 0.00131655, 0.00134255, 0.0014004, 0.00139422, 0.00140613,
      0.00148002, 0.00146782, 0.00151511, 0.00152601, 0.00156979, 0.00156765, 0.00161945, 0.00161797, 0.00164276, 0.00167585,
      0.00163298, 0.00169692, 0.00173064, 0.00185434, 0.0018602, 0.00183116, 0.00184591, 0.00187262, 0.00185005, 0.00187863,
      0.00193908, 0.00196593, 0.00204381, 0.00207335, 0.00214098, 0.00216439, 0.00214946, 0.00215717, 0.00218378, 0.00220102,
      0.00226733, 0.00225952, 0.00232741, 0.0023407, 0.00237913, 0.0023986, 0.00246035, 0.00246658, 0.00256504, 0.00255051,
      0.00260585, 0.00261381, 0.00263094, 0.00263575, 0.00272168, 0.00272411, 0.00270891, 0.00276706, 0.00281004, 0.00290874,
      0.00298958, 0.00287925, 0.00292917, 0.00294574, 0.00308388, 0.00300632, 0.00300664, 0.00304351, 0.00310511, 0.00314704,
      0.00307382, 0.0031967, 0.00324223, 0.0032804, 0.00328582, 0.00328493, 0.00322873, 0.00335108, 0.00344131, 0.00348589,
      0.00347806, 0.00350257, 0.00354257, 0.00361777, 0.00364609, 0.00357236, 0.00361374, 0.0036504, 0.00367545, 0.00370214,
      0.00372157, 0.00380769, 0.00386073, 0.00378243, 0.0038225, 0.00388172, 0.00388896, 0.0039206, 0.00382433, 0.0039355,
      0.00394595, 0.00403936, 0.00412466, 0.00398097, 0.00412995, 0.00410556, 0.00412607, 0.00409762, 0.00417612, 0.00419223,
      0.00418902, 0.00420858, 0.00418737, 0.00438864, 0.00428859, 0.00424058, 0.00428341, 0.00434349, 0.00426255, 0.00425739,
      0.00426099, 0.00435654, 0.00430615, 0.0043434, 0.00436442, 0.00443317, 0.00453702, 0.00458465, 0.00451718, 0.00454829,
      0.004526, 0.00444941, 0.00459784, 0.00461657, 0.00464116, 0.00463936, 0.00462128, 0.00464869, 0.00472074, 0.00464255,
      0.00463531, 0.00464357, 0.00472312, 0.00471226, 0.00472876, 0.00475929, 0.00477558, 0.00477493, 0.00476745, 0.00480513,
      0.00488634, 0.00489984, 0.00491695, 0.0049675, 0.00488809, 0.00492643, 0.0048836, 0.00497446, 0.00506449, 0.00503294,
      0.00507216, 0.00511015, 0.00528854, 0.00508509, 0.00508214, 0.00515293, 0.00521213, 0.00535207, 0.00521807, 0.00530294,
      0.00523004, 0.00531701, 0.00543113, 0.00540544, 0.00534221, 0.00529414, 0.00536786, 0.00530663, 0.00540221, 0.0054366,
      0.00546841, 0.00534307, 0.0053802, 0.00543647, 0.00542493, 0.00540019, 0.0054354, 0.00542636, 0.00540218, 0.00539761,
      0.00546834, 0.00538748, 0.0054119, 0.00524333, 0.0052661, 0.00528475, 0.00527178, 0.00527688, 0.00527451, 0.0051803,
      0.00532525, 0.00516377, 0.00502179, 0.00498588, 0.00493792, 0.00504124, 0.00497812, 0.00480997, 0.00484698, 0.00475996,
      0.00467631, 0.0046375, 0.00452982, 0.00445893, 0.00443842, 0.00431203, 0.004365, 0.00422527, 0.00416564, 0.00407752,
      0.00394289, 0.00401572, 0.00385061, 0.00368295, 0.00359867, 0.00356532, 0.00350603, 0.00343337, 0.00328738, 0.00314308,
      0.00308744, 0.00307742, 0.00298919, 0.00291293, 0.00276594, 0.00274861, 0.00263172, 0.00244029, 0.00243841, 0.0023837,
      0.00222602, 0.00220015, 0.00206933, 0.0019725, 0.00198421, 0.00188412, 0.00176243, 0.00169384, 0.00163604, 0.00154739,
      0.00148061, 0.00135909, 0.00128767, 0.00121056, 0.00116466, 0.00113142, 0.00102363, 0.000933341, 0.000903136, 0.00086764,
      0.000834096, 0.000722275, 0.000704733, 0.000665919, 0.000615007, 0.00057151, 0.000541261, 0.000484668, 0.000465868, 0.000435923,
      0.000406347, 0.000362091, 0.000334272, 0.000302766, 0.000280562, 0.000268874, 0.000248102, 0.000231557, 0.0002133, 0.0001981,
      0.000193611, 0.000166442, 0.000154017, 0.000139056, 0.0001301, 0.00011928, 0.000122511, 0.000106704, 8.5411e-05, 8.49339e-05,
      8.20717e-05, 6.96905e-05, 6.09304e-05, 5.20402e-05, 5.76345e-05, 5.77646e-05, 4.39089e-05, 6.16243e-05 };
  /*
  const G4int nEntries_RI = 44;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV, 
      4.5085*eV, 4.9594*eV, 5.5*eV   , 6.5*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.8210, 1.8212, 1.8215, 1.8219, 
      1.8223, 1.8227, 1.8231, 1.8236, 
      1.8240, 1.8245, 1.8250, 1.8255, 
      1.8261, 1.8266, 1.8272, 1.8279, 
      1.8285, 1.8293, 1.8300, 1.8308, 
      1.8317, 1.8327, 1.8338, 1.8349, 
      1.8362, 1.8376, 1.8392, 1.8410, 
      1.8430, 1.8453, 1.8479, 1.8509, 
      1.8545, 1.8587, 1.8637, 1.8699, 
      1.8774, 1.8869, 1.8991, 1.9152, 
      1.9374, 1.9694, 1.98  , 2. };
*/

  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.86,   1.86,   1.86,    1.86,    1.86,    1.86,    1.86,   1.86};
  
  const G4int nEntries_ABS = 89;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.78417*eV, 1.7971*eV, 1.81022*eV, 1.82353*eV, 1.83704*eV, 1.85075*eV, 1.86466*eV, 1.87879*eV, 1.89313*eV, 1.90769*eV,
      1.92248*eV, 1.9375*eV, 1.95276*eV, 1.96825*eV, 1.984*eV, 2*eV, 2.01626*eV, 2.03279*eV, 2.04959*eV, 2.06667*eV,
      2.08403*eV, 2.10169*eV, 2.11966*eV, 2.13793*eV, 2.15652*eV, 2.17544*eV, 2.19469*eV, 2.21429*eV, 2.23423*eV, 2.25455*eV,
      2.27523*eV, 2.2963*eV, 2.31776*eV, 2.33962*eV, 2.3619*eV, 2.38462*eV, 2.40777*eV, 2.43137*eV, 2.45545*eV, 2.48*eV,
      2.50505*eV, 2.53061*eV, 2.5567*eV, 2.58333*eV, 2.61053*eV, 2.6383*eV, 2.66667*eV, 2.69565*eV, 2.72527*eV, 2.75556*eV,
      2.78652*eV, 2.81818*eV, 2.85057*eV, 2.88372*eV, 2.91765*eV, 2.95238*eV, 2.98795*eV, 3.02439*eV, 3.06173*eV, 3.1*eV,
      3.13924*eV, 3.17949*eV, 3.22078*eV, 3.26316*eV, 3.30667*eV, 3.35135*eV, 3.39726*eV, 3.44444*eV, 3.49296*eV, 3.54286*eV,
      3.5942*eV, 3.64706*eV, 3.70149*eV, 3.75758*eV, 3.81538*eV, 3.875*eV, 3.93651*eV, 4*eV, 4.06557*eV, 4.13333*eV,
      4.20339*eV, 4.27586*eV, 4.35088*eV, 4.42857*eV, 4.50909*eV, 4.59259*eV, 4.67925*eV, 4.76923*eV, 4.86275*eV };
  G4double Absorption[nEntries_ABS] =
    { 1.000*m, 1.000*m, 1.0000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m,
      1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 1.000*m, 0.90*m, 0.85*m, 0.80*m, 0.75*m,
      0.70*m, 0.65*m, 0.60*m, 0.58*m, 0.57*m, 0.56*m, 0.55*m, 0.54*m, 0.54*m, 0.54*m,
      0.531547*m, 0.53*m, 0.5251*m, 0.510*m, 0.493902*m, 0.422165*m, 0.354962*m, 0.255139*m, 0.151762*m, 0.0827965*m,
      0.0409174*m, 0.02007*m, 0.0100463*m, 0.00588694*m, 0.00636339*m, 0.0084497*m, 0.0224574*m, 0.0509883*m, 0.262914*m, 0.0571499*m,
      0.0830375*m, 0.378696*m, 0.0528428*m, 0.0661874*m, 0.0930821*m, 0.0672707*m, 0.0152385*m, 0.00676752*m, 0.00538106*m, 0.00799596*m,
      0.0177025*m, 0.0411282*m, 0.0919861*m, 0.149875*m, 0.132761*m, 0.068419*m, 0.0246548*m, 0.00922619*m, 0.00902168*m, 0.0264256*m,
      0.0839517*m, 0.0796384*m, 0.0552649*m, 0.0197203*m, 0.00872616*m, 0.00764327*m, 0.0153009*m, 0.0299903*m, 0.0403526*m, 0.0377371*m,
      0.0322887*m, 0.0251734*m, 0.0194992*m, 0.0145645*m, 0.0112908*m, 0.0100775*m, 0.0112081*m, 0.0158907*m, 0.019793*m };
  
  const G4int nEntries_SCY = 12;
  G4double ElectronEnergy_SCY[nEntries_SCY] =
    { 0.000*MeV, 0.015*MeV, 0.020*MeV, 0.030*MeV,
      0.040*MeV, 0.060*MeV, 0.080*MeV, 0.090*MeV,
      0.105*MeV, 0.300*MeV, 0.500*MeV, 1.000*MeV };
  G4double ScintilYield[nEntries_SCY] =
    { 0.10, 0.46, 0.60, 0.68,
      0.74, 0.80, 0.82, 0.84,
      0.87, 0.96, 0.98, 1.00 };
  for(int i = 0; i < nEntries_SCY; i++)
    ScintilYield[i] = 15000.0*MeV*ScintilYield[i]*ElectronEnergy_SCY[i];
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("FASTCOMPONENT", PhotonEnergy_FAST, FastComponent,   nEntries_FAST);
  myMPT->AddProperty("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  //myMPT->AddProperty("RAYLEIGH",      PhotonEnergy_RI,   Rayleigh,        nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_ABS);
  //myMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", ElectronEnergy_SCY, ScintilYield, nEntries_SCY);
  myMPT->AddConstProperty("SCINTILLATIONYIELD",16000/MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE",8.5);
  myMPT->AddConstProperty("FASTTIMECONSTANT",55.*ns);
//  myMPT->AddConstProperty("SLOWTIMECONSTANT",0.*ns);
  myMPT->AddConstProperty("YIELDRATIO",1.0);
  myMPT->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}



G4Material* MyMaterials::LuAG_Pr() // Lutetium Aluminum Garnet - 
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z= 8., a= 16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71., a=174.97*g/mole);
  G4Element* Al = new G4Element("Aluminum", "Al", z=13., a= 28.09*g/mole);
  
  G4Material* mat = new G4Material("LuAG_Pr", density=6.7*g/cm3,3);
  mat->AddElement(Lu,3);
  mat->AddElement(Al,5);
  mat->AddElement(O,12);
  
  //  2 bands at 290nm (4.28eV) and 350nm (3.54eV) about 50% of the light in each.
  const G4int nEntries_FAST = 152;
  G4double PhotonEnergy_FAST[nEntries_FAST] = 
    {4.27586*eV, 4.26117*eV, 4.24658*eV, 4.23208*eV, 4.21769*eV, 4.20339*eV, 4.18919*eV, 4.17508*eV, 4.16107*eV, 4.14716*eV, 4.13333*eV, 4.1196*eV, 4.10596*eV, 4.09241*eV, 4.07895*eV, 4.06557*eV, 4.05229*eV, 4.03909*eV, 4.02597*eV, 4.01294*eV, 4*eV, 3.98714*eV, 3.97436*eV, 3.96166*eV, 3.94904*eV, 3.93651*eV, 3.92405*eV, 3.91167*eV, 3.89937*eV, 3.88715*eV, 3.875*eV, 3.86293*eV, 3.85093*eV, 3.83901*eV, 3.82716*eV, 3.81538*eV, 3.80368*eV, 3.79205*eV, 3.78049*eV, 3.769*eV, 3.75758*eV, 3.74622*eV, 3.73494*eV, 3.72372*eV, 3.71257*eV, 3.70149*eV, 3.69048*eV, 3.67953*eV, 3.66864*eV, 3.65782*eV, 3.64706*eV, 3.63636*eV, 3.62573*eV, 3.61516*eV, 3.60465*eV, 3.5942*eV, 3.58382*eV, 3.57349*eV, 3.56322*eV, 3.55301*eV, 3.54286*eV, 3.53276*eV, 3.52273*eV, 3.51275*eV, 3.50282*eV, 3.49296*eV, 3.48315*eV, 3.47339*eV, 3.46369*eV, 3.45404*eV, 3.44444*eV, 3.4349*eV, 3.42541*eV, 3.41598*eV, 3.40659*eV, 3.39726*eV, 3.38798*eV, 3.37875*eV, 3.36957*eV, 3.36043*eV, 3.35135*eV, 3.34232*eV, 3.33333*eV, 3.3244*eV, 3.31551*eV, 3.30667*eV, 3.29787*eV, 3.28912*eV, 3.28042*eV, 3.27177*eV, 3.26316*eV, 3.25459*eV, 3.24607*eV, 3.2376*eV, 3.22917*eV, 3.22078*eV, 3.21244*eV, 3.20413*eV, 3.19588*eV, 3.18766*eV, 3.17949*eV, 3.17136*eV, 3.16327*eV, 3.15522*eV, 3.14721*eV, 3.13924*eV, 3.13131*eV, 3.12343*eV, 3.11558*eV, 3.10777*eV, 3.1*eV, 3.09227*eV, 3.08458*eV, 3.07692*eV, 3.06931*eV, 3.06173*eV, 3.05419*eV, 3.04668*eV, 3.03922*eV, 3.03178*eV, 3.02439*eV, 3.01703*eV, 3.00971*eV, 3.00242*eV, 2.99517*eV, 2.98795*eV, 2.98077*eV, 2.97362*eV, 2.96651*eV, 2.95943*eV, 2.95238*eV, 2.94537*eV, 2.93839*eV, 2.93144*eV, 2.92453*eV, 2.91765*eV, 2.9108*eV, 2.90398*eV, 2.8972*eV, 2.89044*eV, 2.88372*eV, 2.87703*eV, 2.87037*eV, 2.86374*eV, 2.85714*eV, 2.85057*eV, 2.84404*eV, 2.83753*eV, 2.83105*eV, 2.8246*eV, 2.81818*eV, 2.81498*eV};

  G4double FastComponent[nEntries_FAST] = 
    {0.0432883, 0.0544879, 0.0728316, 0.0984743, 0.131255, 0.168052, 0.209737, 0.258707, 0.319857, 0.388254, 0.457261, 0.521868, 0.585468, 0.647954, 0.708267, 0.761854, 0.806207, 0.846097, 0.875271, 0.888461, 0.89361, 0.901929, 0.909677, 0.909872, 0.905667, 0.903999, 0.90391, 0.901108, 0.888578, 0.868688, 0.854241, 0.850029, 0.846092, 0.833223, 0.816599, 0.806819, 0.799456, 0.787796, 0.766995, 0.745956, 0.723528, 0.701548, 0.680262, 0.658418, 0.635332, 0.616613, 0.601253, 0.581025, 0.563156, 0.549948, 0.542195, 0.529153, 0.519084, 0.513169, 0.515516, 0.520102, 0.523, 0.527527, 0.532145, 0.538273, 0.552397, 0.564101, 0.571792, 0.577502, 0.592028, 0.604143, 0.612318, 0.623396, 0.636049, 0.644146, 0.649778, 0.651031, 0.653231, 0.659969, 0.673727, 0.68391, 0.68999, 0.697213, 0.703817, 0.706272, 0.704779, 0.707771, 0.710137, 0.706503, 0.698464, 0.693413, 0.692696, 0.689932, 0.684025, 0.673953, 0.663181, 0.65151, 0.639958, 0.624073, 0.610787, 0.594851, 0.575486, 0.552245, 0.535806, 0.521441, 0.504155, 0.481418, 0.456575, 0.431041, 0.413058, 0.398687, 0.382835, 0.360953, 0.337387, 0.317093, 0.303603, 0.289973, 0.270853, 0.250078, 0.235426, 0.224264, 0.212303, 0.199736, 0.186369, 0.173554, 0.161876, 0.15471, 0.148464, 0.140449, 0.131234, 0.12331, 0.115828, 0.107644, 0.100514, 0.0965623, 0.0943724, 0.0904375, 0.0846862, 0.0784364, 0.0720687, 0.0648027, 0.0592781, 0.0577879, 0.0574301, 0.0547233, 0.0515322, 0.0501597, 0.0493648, 0.0503557, 0.0510614, 0.0495324, 0.044173, 0.0387929, 0.0347393, 0.0334513, 0.0326892, 0.0321444};
  
/*
  const G4int nEntries_RI = 44;
  G4double PhotonEnergy_RI[nEntries_RI] = 
    { 0.1000*eV, 1.0000*eV, 1.0121*eV, 1.0332*eV, 
      1.0552*eV, 1.0781*eV, 1.1021*eV, 1.1271*eV, 
      1.1533*eV, 1.1808*eV, 1.2096*eV, 1.2398*eV, 
      1.2716*eV, 1.3051*eV, 1.3404*eV, 1.3776*eV, 
      1.4170*eV, 1.4586*eV, 1.5028*eV, 1.5498*eV, 
      1.5998*eV, 1.6531*eV, 1.7101*eV, 1.7712*eV, 
      1.8368*eV, 1.9074*eV, 1.9837*eV, 2.0664*eV, 
      2.1562*eV, 2.2543*eV, 2.3616*eV, 2.4797*eV, 
      2.6102*eV, 2.7552*eV, 2.9173*eV, 3.0996*eV, 
      3.3062*eV, 3.5424*eV, 3.8149*eV, 4.1328*eV, 
      4.5085*eV, 4.9594*eV, 5.5*eV   , 6.53*eV };
  
  G4double RefractiveIndex[nEntries_RI] = 
    { 1.8210, 1.8212, 1.8215, 1.8219, 
      1.8223, 1.8227, 1.8231, 1.8236, 
      1.8240, 1.8245, 1.8250, 1.8255, 
      1.8261, 1.8266, 1.8272, 1.8279, 
      1.8285, 1.8293, 1.8300, 1.8308, 
      1.8317, 1.8327, 1.8338, 1.8349, 
      1.8362, 1.8376, 1.8392, 1.8410, 
      1.8430, 1.8453, 1.8479, 1.8509, 
      1.8545, 1.8587, 1.8637, 1.8699, 
      1.8774, 1.8869, 1.8991, 1.9152, 
      1.9374, 1.9694, 1.98  , 2. };
*/
  //fixed index of refraction
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.86,   1.86,   1.86,    1.86,    1.86,    1.86,    1.86,   1.86};

  const G4int nEntries_ABS = 86;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 4.42857*eV, 4.35088*eV, 4.27586*eV, 4.20339*eV, 4.13333*eV, 4.06557*eV, 4*eV, 3.93651*eV, 3.875*eV, 3.81538*eV, 3.75758*eV, 3.70149*eV, 3.64706*eV, 3.5942*eV, 3.54286*eV, 3.49296*eV, 3.44444*eV, 3.39726*eV, 3.35135*eV, 3.30667*eV, 3.26316*eV, 3.22078*eV, 3.17949*eV, 3.13924*eV, 3.1*eV, 3.06173*eV, 3.02439*eV, 2.98795*eV, 2.95238*eV, 2.91765*eV, 2.88372*eV, 2.85057*eV, 2.81818*eV, 2.78652*eV, 2.75556*eV, 2.72527*eV, 2.69565*eV, 2.66667*eV, 2.6383*eV, 2.61053*eV, 2.58333*eV, 2.5567*eV, 2.53061*eV, 2.50505*eV, 2.48*eV, 2.45545*eV, 2.43137*eV, 2.40777*eV, 2.38462*eV, 2.3619*eV, 2.33962*eV, 2.31776*eV, 2.2963*eV, 2.27523*eV, 2.25455*eV, 2.23423*eV, 2.21429*eV, 2.19469*eV, 2.17544*eV, 2.15652*eV, 2.13793*eV, 2.11966*eV, 2.10169*eV, 2.08403*eV, 2.06667*eV, 2.04959*eV, 2.03279*eV, 2.01626*eV, 2*eV, 1.984*eV, 1.96825*eV, 1.95276*eV, 1.9375*eV, 1.92248*eV, 1.90769*eV, 1.89313*eV, 1.87879*eV, 1.86466*eV, 1.85075*eV, 1.83704*eV, 1.82353*eV, 1.81022*eV, 1.7971*eV, 1.78417*eV, 1.77143*eV, 1.*eV};

  G4double Absorption[nEntries_ABS] =
    { 1.0*mm, 3.56338*mm, 3.24351*mm, 4.00004*mm, 14.0716*mm, 38.9298*mm, 56.2426*mm, 64.8461*mm, 84.1463*mm, 83.1221*mm, 86.9579*mm, 94.4561*mm, 98.4977*mm, 104.277*mm, 139.249*mm, 133.6*mm, 182.28*mm, 163.112*mm, 204.257*mm, 185.78*mm, 241.729*mm, 295.022*mm, 334.904*mm, 440.211*mm, 303.927*mm, 409.643*mm, 291.351*mm, 249.281*mm, 220.063*mm, 140.16*mm, 85.426*mm, 37.4797*mm, 22.7967*mm, 9.87634*mm, 19.3781*mm, 23.7056*mm, 71.9216*mm, 34.3367*mm, 47.3344*mm, 122.49*mm, 25.9357*mm, 218.548*mm, 414.841*mm, 600.652*mm, 885.534*mm, 776.05*mm, 810.135*mm, 1108.04*mm, 782.709*mm, 906.707*mm, 964.158*mm, 1463*mm, 939.49*mm, 1118.37*mm, 839.318*mm, 1592.08*mm, 762.109*mm, 823.407*mm, 441.257*mm, 189.253*mm, 124.454*mm, 63.6358*mm, 229.853*mm, 429.362*mm, 39.1411*mm, 112.491*mm, 566.704*mm, 935.834*mm, 793.498*mm, 1031.76*mm, 1032.27*mm, 1065.16*mm, 1160*mm, 1019.89*mm, 1079.36*mm, 1072.12*mm, 1027.8*mm, 1073.86*mm, 967.219*mm, 1017.6*mm, 1038.24*mm, 1064.4*mm, 1065.45*mm, 1098.97*mm, 982.938*mm, 1100*mm};

  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("FASTCOMPONENT", PhotonEnergy_FAST, FastComponent,   nEntries_FAST);
  myMPT->AddProperty("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  //myMPT->AddProperty("RAYLEIGH",      PhotonEnergy_RI,   Rayleigh,        nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_RI);
  myMPT->AddConstProperty("SCINTILLATIONYIELD",7000/MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE",6.4);
  myMPT->AddConstProperty("FASTTIMECONSTANT",20.*ns);
  myMPT->AddConstProperty("YIELDRATIO",1.0);
  myMPT->AddConstProperty("FASTSCINTILLATIONRISETIME",0.5*ns);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}



G4Material* MyMaterials::LYSO()
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z= 8., a= 16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71., a=174.97*g/mole);
  G4Element* Si = new G4Element("Silicon",  "Si", z=14., a= 28.09*g/mole);
  
  G4Material* mat = new G4Material("LYSO", density=7.1*g/cm3,3,kStateSolid);
  mat->AddElement(Lu,2);
  mat->AddElement(Si,1);
  mat->AddElement(O,5);
  
  const G4int nEntries_FAST = 261;
  G4double PhotonEnergy_FAST[nEntries_FAST] = 
    { 1.77169*eV, 1.77266*eV, 1.77558*eV, 1.77851*eV, 1.78145*eV, 1.78539*eV, 1.79033*eV, 1.7963*eV, 1.80231*eV, 1.80836*eV,
      1.81445*eV, 1.82058*eV, 1.82882*eV, 1.83401*eV, 1.84553*eV, 1.85293*eV, 1.86147*eV, 1.869*eV, 1.87769*eV, 1.89308*eV,
      1.90536*eV, 1.92007*eV, 1.93039*eV, 1.94901*eV, 1.95846*eV, 1.9668*eV, 1.97884*eV, 1.99102*eV, 2.00088*eV, 2.01209*eV,
      2.02596*eV, 2.03617*eV, 2.04519*eV, 2.0569*eV, 2.06611*eV, 2.0794*eV, 2.09151*eV, 2.10239*eV, 2.112*eV, 2.1231*eV,
      2.13431*eV, 2.14565*eV, 2.15566*eV, 2.16868*eV, 2.18038*eV, 2.19519*eV, 2.21171*eV, 2.2193*eV, 2.23619*eV, 2.23464*eV,
      2.24395*eV, 2.25806*eV, 2.27234*eV, 2.28358*eV, 2.29493*eV, 2.30475*eV, 2.31631*eV, 2.32463*eV, 2.33134*eV, 2.33809*eV,
      2.34487*eV, 2.35856*eV, 2.36719*eV, 2.37939*eV, 2.38642*eV, 2.40238*eV, 2.41134*eV, 2.424*eV, 2.43312*eV, 2.44047*eV,
      2.44786*eV, 2.46278*eV, 2.47788*eV, 2.48741*eV, 2.49317*eV, 2.49702*eV, 2.50282*eV, 2.50865*eV, 2.5145*eV, 2.52038*eV,
      2.52432*eV, 2.53223*eV, 2.5362*eV, 2.54619*eV, 2.55424*eV, 2.56031*eV, 2.56437*eV, 2.57049*eV, 2.57663*eV, 2.58487*eV,
      2.59317*eV, 2.59734*eV, 2.60571*eV, 2.61414*eV, 2.61414*eV, 2.61837*eV, 2.62262*eV, 2.62475*eV, 2.62902*eV, 2.63331*eV,
      2.63545*eV, 2.63976*eV, 2.64191*eV, 2.64841*eV, 2.65493*eV, 2.6593*eV, 2.66149*eV, 2.66588*eV, 2.67914*eV, 2.67914*eV,
      2.68136*eV, 2.68136*eV, 2.68359*eV, 2.68805*eV, 2.68805*eV, 2.68805*eV, 2.69477*eV, 2.69477*eV, 2.69702*eV, 2.70153*eV,
      2.70605*eV, 2.71286*eV, 2.71742*eV, 2.71971*eV, 2.722*eV, 2.722*eV, 2.72429*eV, 2.72889*eV, 2.72889*eV, 2.73351*eV,
      2.73814*eV, 2.74279*eV, 2.74512*eV, 2.74979*eV, 2.75213*eV, 2.75447*eV, 2.75917*eV, 2.75682*eV, 2.76389*eV, 2.76626*eV,
      2.76389*eV, 2.76626*eV, 2.77338*eV, 2.77576*eV, 2.78533*eV, 2.79255*eV, 2.79738*eV, 2.80223*eV, 2.80466*eV, 2.80709*eV,
      2.80953*eV, 2.80953*eV, 2.81934*eV, 2.8218*eV, 2.82673*eV, 2.83168*eV, 2.84164*eV, 2.84916*eV, 2.85419*eV, 2.8643*eV,
      2.86684*eV, 2.87449*eV, 2.87705*eV, 2.87961*eV, 2.88475*eV, 2.88733*eV, 2.8925*eV, 2.89509*eV, 2.90028*eV, 2.90549*eV,
      2.90811*eV, 2.91073*eV, 2.91335*eV, 2.91335*eV, 2.91335*eV, 2.91861*eV, 2.92125*eV, 2.92125*eV, 2.92389*eV, 2.92654*eV,
      2.92654*eV, 2.92919*eV, 2.92919*eV, 2.93185*eV, 2.93451*eV, 2.93717*eV, 2.93985*eV, 2.94252*eV, 2.9452*eV, 2.94789*eV,
      2.94789*eV, 2.94789*eV, 2.95058*eV, 2.95868*eV, 2.96411*eV, 2.96955*eV, 2.97228*eV, 2.97228*eV, 2.96955*eV, 2.97228*eV,
      2.97502*eV, 2.97776*eV, 2.97502*eV, 2.9805*eV, 2.9805*eV, 2.9805*eV, 2.98601*eV, 2.99154*eV, 2.99431*eV, 2.99431*eV,
      2.99708*eV, 2.99431*eV, 2.99708*eV, 3.00544*eV, 3.00824*eV, 3.00824*eV, 3.00824*eV, 3.00824*eV, 3.01385*eV, 3.0223*eV,
      3.02797*eV, 3.03081*eV, 3.02797*eV, 3.03365*eV, 3.03081*eV, 3.03081*eV, 3.0365*eV, 3.03935*eV, 3.04221*eV, 3.04795*eV,
      3.04795*eV, 3.05083*eV, 3.05371*eV, 3.05949*eV, 3.06239*eV, 3.06529*eV, 3.0682*eV, 3.06529*eV, 3.07112*eV, 3.0682*eV,
      3.07696*eV, 3.08283*eV, 3.0976*eV, 3.09464*eV, 3.09464*eV, 3.10653*eV, 3.11252*eV, 3.11852*eV, 3.12757*eV, 3.13668*eV,
      3.14583*eV, 3.15813*eV, 3.16741*eV, 3.17675*eV, 3.20828*eV, 3.23719*eV, 3.26664*eV, 3.28656*eV, 3.31351*eV, 3.34783*eV,
      3.38287*eV };
  G4double FastComponent[nEntries_FAST] = 
    { 0.011691, 0.011691, 0.011691, 0.0146138, 0.0146138, 0.0146138, 0.011691, 0.011691, 0.00876827, 0.00876827,
      0.00584551, 0.00584551, 0.00584551, 0.00292276, 0.00876827, 0.0146138, 0.0146138, 0.0146138, 0.0204593, 0.023382,
      0.0263048, 0.0204593, 0.0204593, 0.023382, 0.0292276, 0.0321503, 0.0350731, 0.0379958, 0.0379958, 0.0379958,
      0.0350731, 0.0379958, 0.0409186, 0.0438413, 0.0526096, 0.0584551, 0.0643006, 0.0730689, 0.0730689, 0.0818372,
      0.0906054, 0.0964509, 0.0993737, 0.105219, 0.111065, 0.122756, 0.125678, 0.146138, 0.146138, 0.160752,
      0.157829, 0.163674, 0.184134, 0.192902, 0.20167, 0.219207, 0.230898, 0.242589, 0.25428, 0.265971,
      0.274739, 0.292276, 0.306889, 0.315658, 0.321503, 0.350731, 0.368267, 0.385804, 0.397495, 0.415031,
      0.432568, 0.458873, 0.482255, 0.496868, 0.514405, 0.529019, 0.549478, 0.564092, 0.581628, 0.593319,
      0.602088, 0.616701, 0.637161, 0.660543, 0.681002, 0.71023, 0.736534, 0.756994, 0.777453, 0.806681,
      0.844676, 0.868058, 0.891441, 0.9119, 0.938205, 0.955741, 0.984969, 1.0142, 1.03173, 1.05511,
      1.07557, 1.11649, 1.13695, 1.15741, 1.17495, 1.19248, 1.21002, 1.22756, 1.27432, 1.2977,
      1.31524, 1.32985, 1.36785, 1.40292, 1.39415, 1.4, 1.41754, 1.44092, 1.47015, 1.48476,
      1.50814, 1.5286, 1.54906, 1.56952, 1.58998, 1.61921, 1.63967, 1.66597, 1.68935, 1.71566,
      1.73904, 1.76242, 1.77996, 1.80042, 1.8238, 1.83549, 1.85303, 1.8618, 1.87933, 1.89979,
      1.91733, 1.92902, 1.95825, 1.98163, 2.01378, 2.03424, 2.0547, 2.07808, 2.09562, 2.11023,
      2.12484, 2.13361, 2.15407, 2.15699, 2.15992, 2.16576, 2.16868, 2.16868, 2.16284, 2.15699,
      2.14823, 2.13946, 2.12484, 2.11023, 2.08977, 2.06639, 2.04593, 2.02839, 2.01086, 1.98455,
      1.96409, 1.94948, 1.93194, 1.91733, 1.90271, 1.87641, 1.86472, 1.8501, 1.83841, 1.82088,
      1.79749, 1.77119, 1.75073, 1.73027, 1.70689, 1.68058, 1.65428, 1.6309, 1.60167, 1.57244,
      1.55491, 1.53152, 1.50522, 1.47891, 1.45261, 1.43215, 1.40877, 1.38831, 1.362, 1.33862,
      1.31232, 1.28601, 1.27432, 1.25678, 1.21587, 1.19541, 1.17203, 1.14864, 1.12234, 1.10772,
      1.08434, 1.06096, 1.0142, 0.987891, 0.967432, 0.938205, 0.9119, 0.879749, 0.853445, 0.82714,
      0.786221, 0.765762, 0.739457, 0.716075, 0.681002, 0.660543, 0.637161, 0.60501, 0.581628, 0.552401,
      0.531942, 0.505637, 0.485177, 0.458873, 0.435491, 0.412109, 0.379958, 0.356576, 0.336117, 0.309812,
      0.280585, 0.25428, 0.207516, 0.175365, 0.157829, 0.13737, 0.119833, 0.0993737, 0.0759916, 0.0613779,
      0.0526096, 0.0350731, 0.0263048, 0.011691, 0.00876827, 0.00876827, 0.011691, 0.011691, 0.011691, 0.00876827,
      0.011691 };
  
  const G4int nEntries_RI = 3;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0*eV, 1.84*eV, 4.5*eV};
  G4double RefractiveIndex[nEntries_RI] =
    //{ 1.75, 1.82, 1.88};
    { 1.82, 1.82, 1.82};
  //G4double Rayleigh[nEntries_RI] =
  //  { 138.*mm, 138.*mm, 138.*mm};
  
  const G4int nEntries_ABS = 9;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.00*eV , 2.82*eV , 2.88*eV , 2.95*eV , 3.02*eV  , 3.10*eV  , 3.18*eV  , 3.26*eV , 4.08*eV };
  G4double Absorption[nEntries_ABS] =
    { 438.*mm , 438.*mm , 413.*mm , 375.*mm , 263.*mm  , 87.5*mm  , 11.5*mm  , 1.0*mm  , 1.0*mm  };
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("FASTCOMPONENT", PhotonEnergy_FAST, FastComponent,  nEntries_FAST);
  myMPT->AddProperty("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  //myMPT->AddProperty("RAYLEIGH",      PhotonEnergy_ABS,  Rayleigh,        nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_ABS);
  myMPT->AddConstProperty("SCINTILLATIONYIELD",40000./MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE",3.4);
  myMPT->AddConstProperty("FASTTIMECONSTANT",40.*ns);
  myMPT->AddConstProperty("YIELDRATIO",1.0);
  myMPT->AddConstProperty("FASTSCINTILLATIONRISETIME",0.1*ns);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}






G4Material* MyMaterials::LYSO_lowLY()
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z=8.,  a= 16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71., a=174.97*g/mole);
  G4Element* Si = new G4Element("Silicon",  "Si", z=14., a= 28.09*g/mole);
  
  G4Material* mat = new G4Material("LYSO_lowLY", density=7.1*g/cm3,3,kStateSolid);
  mat->AddElement(Lu,2);
  mat->AddElement(Si,1);
  mat->AddElement(O,5);
  
  const G4int nEntries_FAST = 261;
  G4double PhotonEnergy_FAST[nEntries_FAST] = 
    { 1.77169*eV, 1.77266*eV, 1.77558*eV, 1.77851*eV, 1.78145*eV, 1.78539*eV, 1.79033*eV, 1.7963*eV, 1.80231*eV, 1.80836*eV,
      1.81445*eV, 1.82058*eV, 1.82882*eV, 1.83401*eV, 1.84553*eV, 1.85293*eV, 1.86147*eV, 1.869*eV, 1.87769*eV, 1.89308*eV,
      1.90536*eV, 1.92007*eV, 1.93039*eV, 1.94901*eV, 1.95846*eV, 1.9668*eV, 1.97884*eV, 1.99102*eV, 2.00088*eV, 2.01209*eV,
      2.02596*eV, 2.03617*eV, 2.04519*eV, 2.0569*eV, 2.06611*eV, 2.0794*eV, 2.09151*eV, 2.10239*eV, 2.112*eV, 2.1231*eV,
      2.13431*eV, 2.14565*eV, 2.15566*eV, 2.16868*eV, 2.18038*eV, 2.19519*eV, 2.21171*eV, 2.2193*eV, 2.23619*eV, 2.23464*eV,
      2.24395*eV, 2.25806*eV, 2.27234*eV, 2.28358*eV, 2.29493*eV, 2.30475*eV, 2.31631*eV, 2.32463*eV, 2.33134*eV, 2.33809*eV,
      2.34487*eV, 2.35856*eV, 2.36719*eV, 2.37939*eV, 2.38642*eV, 2.40238*eV, 2.41134*eV, 2.424*eV, 2.43312*eV, 2.44047*eV,
      2.44786*eV, 2.46278*eV, 2.47788*eV, 2.48741*eV, 2.49317*eV, 2.49702*eV, 2.50282*eV, 2.50865*eV, 2.5145*eV, 2.52038*eV,
      2.52432*eV, 2.53223*eV, 2.5362*eV, 2.54619*eV, 2.55424*eV, 2.56031*eV, 2.56437*eV, 2.57049*eV, 2.57663*eV, 2.58487*eV,
      2.59317*eV, 2.59734*eV, 2.60571*eV, 2.61414*eV, 2.61414*eV, 2.61837*eV, 2.62262*eV, 2.62475*eV, 2.62902*eV, 2.63331*eV,
      2.63545*eV, 2.63976*eV, 2.64191*eV, 2.64841*eV, 2.65493*eV, 2.6593*eV, 2.66149*eV, 2.66588*eV, 2.67914*eV, 2.67914*eV,
      2.68136*eV, 2.68136*eV, 2.68359*eV, 2.68805*eV, 2.68805*eV, 2.68805*eV, 2.69477*eV, 2.69477*eV, 2.69702*eV, 2.70153*eV,
      2.70605*eV, 2.71286*eV, 2.71742*eV, 2.71971*eV, 2.722*eV, 2.722*eV, 2.72429*eV, 2.72889*eV, 2.72889*eV, 2.73351*eV,
      2.73814*eV, 2.74279*eV, 2.74512*eV, 2.74979*eV, 2.75213*eV, 2.75447*eV, 2.75917*eV, 2.75682*eV, 2.76389*eV, 2.76626*eV,
      2.76389*eV, 2.76626*eV, 2.77338*eV, 2.77576*eV, 2.78533*eV, 2.79255*eV, 2.79738*eV, 2.80223*eV, 2.80466*eV, 2.80709*eV,
      2.80953*eV, 2.80953*eV, 2.81934*eV, 2.8218*eV, 2.82673*eV, 2.83168*eV, 2.84164*eV, 2.84916*eV, 2.85419*eV, 2.8643*eV,
      2.86684*eV, 2.87449*eV, 2.87705*eV, 2.87961*eV, 2.88475*eV, 2.88733*eV, 2.8925*eV, 2.89509*eV, 2.90028*eV, 2.90549*eV,
      2.90811*eV, 2.91073*eV, 2.91335*eV, 2.91335*eV, 2.91335*eV, 2.91861*eV, 2.92125*eV, 2.92125*eV, 2.92389*eV, 2.92654*eV,
      2.92654*eV, 2.92919*eV, 2.92919*eV, 2.93185*eV, 2.93451*eV, 2.93717*eV, 2.93985*eV, 2.94252*eV, 2.9452*eV, 2.94789*eV,
      2.94789*eV, 2.94789*eV, 2.95058*eV, 2.95868*eV, 2.96411*eV, 2.96955*eV, 2.97228*eV, 2.97228*eV, 2.96955*eV, 2.97228*eV,
      2.97502*eV, 2.97776*eV, 2.97502*eV, 2.9805*eV, 2.9805*eV, 2.9805*eV, 2.98601*eV, 2.99154*eV, 2.99431*eV, 2.99431*eV,
      2.99708*eV, 2.99431*eV, 2.99708*eV, 3.00544*eV, 3.00824*eV, 3.00824*eV, 3.00824*eV, 3.00824*eV, 3.01385*eV, 3.0223*eV,
      3.02797*eV, 3.03081*eV, 3.02797*eV, 3.03365*eV, 3.03081*eV, 3.03081*eV, 3.0365*eV, 3.03935*eV, 3.04221*eV, 3.04795*eV,
      3.04795*eV, 3.05083*eV, 3.05371*eV, 3.05949*eV, 3.06239*eV, 3.06529*eV, 3.0682*eV, 3.06529*eV, 3.07112*eV, 3.0682*eV,
      3.07696*eV, 3.08283*eV, 3.0976*eV, 3.09464*eV, 3.09464*eV, 3.10653*eV, 3.11252*eV, 3.11852*eV, 3.12757*eV, 3.13668*eV,
      3.14583*eV, 3.15813*eV, 3.16741*eV, 3.17675*eV, 3.20828*eV, 3.23719*eV, 3.26664*eV, 3.28656*eV, 3.31351*eV, 3.34783*eV,
      3.38287*eV };
  G4double FastComponent[nEntries_FAST] = 
    { 0.011691, 0.011691, 0.011691, 0.0146138, 0.0146138, 0.0146138, 0.011691, 0.011691, 0.00876827, 0.00876827,
      0.00584551, 0.00584551, 0.00584551, 0.00292276, 0.00876827, 0.0146138, 0.0146138, 0.0146138, 0.0204593, 0.023382,
      0.0263048, 0.0204593, 0.0204593, 0.023382, 0.0292276, 0.0321503, 0.0350731, 0.0379958, 0.0379958, 0.0379958,
      0.0350731, 0.0379958, 0.0409186, 0.0438413, 0.0526096, 0.0584551, 0.0643006, 0.0730689, 0.0730689, 0.0818372,
      0.0906054, 0.0964509, 0.0993737, 0.105219, 0.111065, 0.122756, 0.125678, 0.146138, 0.146138, 0.160752,
      0.157829, 0.163674, 0.184134, 0.192902, 0.20167, 0.219207, 0.230898, 0.242589, 0.25428, 0.265971,
      0.274739, 0.292276, 0.306889, 0.315658, 0.321503, 0.350731, 0.368267, 0.385804, 0.397495, 0.415031,
      0.432568, 0.458873, 0.482255, 0.496868, 0.514405, 0.529019, 0.549478, 0.564092, 0.581628, 0.593319,
      0.602088, 0.616701, 0.637161, 0.660543, 0.681002, 0.71023, 0.736534, 0.756994, 0.777453, 0.806681,
      0.844676, 0.868058, 0.891441, 0.9119, 0.938205, 0.955741, 0.984969, 1.0142, 1.03173, 1.05511,
      1.07557, 1.11649, 1.13695, 1.15741, 1.17495, 1.19248, 1.21002, 1.22756, 1.27432, 1.2977,
      1.31524, 1.32985, 1.36785, 1.40292, 1.39415, 1.4, 1.41754, 1.44092, 1.47015, 1.48476,
      1.50814, 1.5286, 1.54906, 1.56952, 1.58998, 1.61921, 1.63967, 1.66597, 1.68935, 1.71566,
      1.73904, 1.76242, 1.77996, 1.80042, 1.8238, 1.83549, 1.85303, 1.8618, 1.87933, 1.89979,
      1.91733, 1.92902, 1.95825, 1.98163, 2.01378, 2.03424, 2.0547, 2.07808, 2.09562, 2.11023,
      2.12484, 2.13361, 2.15407, 2.15699, 2.15992, 2.16576, 2.16868, 2.16868, 2.16284, 2.15699,
      2.14823, 2.13946, 2.12484, 2.11023, 2.08977, 2.06639, 2.04593, 2.02839, 2.01086, 1.98455,
      1.96409, 1.94948, 1.93194, 1.91733, 1.90271, 1.87641, 1.86472, 1.8501, 1.83841, 1.82088,
      1.79749, 1.77119, 1.75073, 1.73027, 1.70689, 1.68058, 1.65428, 1.6309, 1.60167, 1.57244,
      1.55491, 1.53152, 1.50522, 1.47891, 1.45261, 1.43215, 1.40877, 1.38831, 1.362, 1.33862,
      1.31232, 1.28601, 1.27432, 1.25678, 1.21587, 1.19541, 1.17203, 1.14864, 1.12234, 1.10772,
      1.08434, 1.06096, 1.0142, 0.987891, 0.967432, 0.938205, 0.9119, 0.879749, 0.853445, 0.82714,
      0.786221, 0.765762, 0.739457, 0.716075, 0.681002, 0.660543, 0.637161, 0.60501, 0.581628, 0.552401,
      0.531942, 0.505637, 0.485177, 0.458873, 0.435491, 0.412109, 0.379958, 0.356576, 0.336117, 0.309812,
      0.280585, 0.25428, 0.207516, 0.175365, 0.157829, 0.13737, 0.119833, 0.0993737, 0.0759916, 0.0613779,
      0.0526096, 0.0350731, 0.0263048, 0.011691, 0.00876827, 0.00876827, 0.011691, 0.011691, 0.011691, 0.00876827,
      0.011691 };
  
  const G4int nEntries_RI = 3;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0*eV, 1.84*eV, 4.5*eV};
  G4double RefractiveIndex[nEntries_RI] =
    //{ 1.75, 1.82, 1.88};
    { 1.82, 1.82, 1.82};
  //G4double Rayleigh[nEntries_RI]
  //  = { 138.*mm, 138.*mm, 138.*mm};
  
  const G4int nEntries_ABS = 9;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.00*eV , 2.82*eV , 2.88*eV , 2.95*eV , 3.02*eV  , 3.10*eV  , 3.18*eV  , 3.26*eV , 4.08*eV };
  G4double Absorption[nEntries_ABS] =
    { 438.*mm , 438.*mm , 413.*mm , 375.*mm , 263.*mm  , 87.5*mm  , 11.5*mm  , 1.0*mm  , 1.0*mm  };
  
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("FASTCOMPONENT", PhotonEnergy_FAST, FastComponent,  nEntries_FAST);
  myMPT->AddProperty("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  //myMPT->AddProperty("RAYLEIGH",      PhotonEnergy_RI,   Rayleigh,        nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_ABS);
  myMPT->AddConstProperty("SCINTILLATIONYIELD",10./MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE",3.4);
  myMPT->AddConstProperty("FASTTIMECONSTANT",40.*ns);
  myMPT->AddConstProperty("YIELDRATIO",1.0);
  myMPT->AddConstProperty("FASTSCINTILLATIONRISETIME",0.1*ns);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}



G4Material* MyMaterials::LSO()
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z= 8., a= 16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71., a=174.97*g/mole);
  G4Element* Si = new G4Element("Silicon",  "Si", z=14., a= 28.09*g/mole);
  
  G4Material* mat = new G4Material("LSO", density=7.4*g/cm3,3);
  mat->AddElement(Lu,2);
  mat->AddElement(Si,1);
  mat->AddElement(O,5);
  
  const G4int nEntries_FAST = 192;
  G4double PhotonEnergy_FAST[nEntries_FAST] = 
    { 1.83966*eV, 1.84948*eV, 1.86274*eV, 1.87507*eV, 1.88413*eV, 1.90369*eV, 1.91187*eV, 1.92129*eV, 1.92962*eV, 1.93922*eV,
      1.95258*eV, 1.96365*eV, 1.97986*eV, 1.99124*eV, 2.00533*eV, 2.02618*eV, 2.04747*eV, 2.06101*eV, 2.07472*eV, 2.09424*eV,
      2.11269*eV, 2.12565*eV, 2.14466*eV, 2.16251*eV, 2.17914*eV, 2.19602*eV, 2.21317*eV, 2.22422*eV, 2.24021*eV, 2.25479*eV,
      2.26462*eV, 2.27785*eV, 2.29462*eV, 2.30821*eV, 2.32024*eV, 2.33588*eV, 2.34643*eV, 2.35529*eV, 2.37322*eV, 2.38594*eV,
      2.3896*eV, 2.39879*eV, 2.40805*eV, 2.41365*eV, 2.4268*eV, 2.44009*eV, 2.45161*eV, 2.46518*eV, 2.47693*eV, 2.48483*eV,
      2.49477*eV, 2.50479*eV, 2.51692*eV, 2.53123*eV, 2.5457*eV, 2.54986*eV, 2.55613*eV, 2.56033*eV, 2.56665*eV, 2.58796*eV,
      2.59658*eV, 2.60091*eV, 2.60309*eV, 2.60744*eV, 2.614*eV, 2.62059*eV, 2.62943*eV, 2.6361*eV, 2.64057*eV, 2.64729*eV,
      2.65632*eV, 2.66085*eV, 2.6654*eV, 2.66997*eV, 2.67684*eV, 2.67684*eV, 2.68839*eV, 2.69303*eV, 2.70237*eV, 2.70471*eV,
      2.71177*eV, 2.72124*eV, 2.72362*eV, 2.73077*eV, 2.73077*eV, 2.73317*eV, 2.73797*eV, 2.74279*eV, 2.74762*eV, 2.7549*eV,
      2.7549*eV, 2.75978*eV, 2.75978*eV, 2.76468*eV, 2.76713*eV, 2.77205*eV, 2.77699*eV, 2.77699*eV, 2.77947*eV, 2.78941*eV,
      2.79692*eV, 2.80195*eV, 2.80699*eV, 2.8146*eV, 2.81714*eV, 2.8248*eV, 2.8325*eV, 2.83507*eV, 2.85063*eV, 2.85847*eV,
      2.86635*eV, 2.86899*eV, 2.87428*eV, 2.87959*eV, 2.88225*eV, 2.89027*eV, 2.89295*eV, 2.89833*eV, 2.90103*eV, 2.90915*eV,
      2.91186*eV, 2.91731*eV, 2.92278*eV, 2.92278*eV, 2.92553*eV, 2.93103*eV, 2.93103*eV, 2.93103*eV, 2.94487*eV, 2.94487*eV,
      2.94766*eV, 2.95324*eV, 2.95604*eV, 2.95885*eV, 2.95604*eV, 2.96166*eV, 2.96447*eV, 2.97012*eV, 2.96166*eV, 2.97295*eV,
      2.98434*eV, 2.98434*eV, 2.98148*eV, 2.98434*eV, 2.99006*eV, 2.9872*eV, 2.99006*eV, 2.9872*eV, 2.99006*eV, 2.99869*eV,
      3.00447*eV, 3.00737*eV, 3.0161*eV, 3.01902*eV, 3.0161*eV, 3.0161*eV, 3.01318*eV, 3.01318*eV, 3.02194*eV, 3.02781*eV,
      3.03666*eV, 3.03666*eV, 3.03666*eV, 3.04556*eV, 3.05152*eV, 3.05152*eV, 3.05451*eV, 3.05451*eV, 3.05451*eV, 3.06051*eV,
      3.05751*eV, 3.07258*eV, 3.07258*eV, 3.07561*eV, 3.08169*eV, 3.09085*eV, 3.08779*eV, 3.09085*eV, 3.09699*eV, 3.10935*eV,
      3.10625*eV, 3.1218*eV, 3.12807*eV, 3.13121*eV, 3.14067*eV, 3.15657*eV, 3.16941*eV, 3.19213*eV, 3.21849*eV, 3.24529*eV,
      3.27255*eV, 3.28981*eV };
  G4double FastComponent[nEntries_FAST] = 
    { 0.0121475, 0.0121475, 0.0151844, 0.0151844, 0.0151844, 0.0182213, 0.0182213, 0.0182213, 0.024295, 0.024295,
      0.0212581, 0.0212581, 0.0303688, 0.0303688, 0.0303688, 0.0425163, 0.0516269, 0.0607375, 0.0698482, 0.072885,
      0.0850325, 0.0941432, 0.106291, 0.127549, 0.130586, 0.142733, 0.163991, 0.179176, 0.19436, 0.212581,
      0.224729, 0.239913, 0.252061, 0.273319, 0.297614, 0.318872, 0.34013, 0.355315, 0.376573, 0.38872,
      0.413015, 0.4282, 0.440347, 0.458568, 0.47679, 0.507158, 0.531453, 0.567896, 0.595228, 0.628633,
      0.652928, 0.68026, 0.71974, 0.759219, 0.77744, 0.813883, 0.835141, 0.859436, 0.886768, 0.920174,
      0.956616, 0.990022, 1.00521, 1.01735, 1.04165, 1.06898, 1.09328, 1.11757, 1.15098, 1.17223,
      1.2026, 1.23297, 1.26334, 1.29067, 1.32104, 1.37874, 1.40304, 1.43341, 1.46074, 1.49414,
      1.52451, 1.56095, 1.60043, 1.63991, 1.67028, 1.69761, 1.72191, 1.7462, 1.77354, 1.81605,
      1.84946, 1.88286, 1.88286, 1.88894, 1.9102, 1.94056, 1.98308, 2.00434, 2.03167, 2.07419,
      2.10759, 2.13189, 2.15315, 2.16833, 2.17744, 2.19566, 2.20781, 2.20781, 2.21996, 2.21692,
      2.20477, 2.18959, 2.16833, 2.14403, 2.11367, 2.08026, 2.04685, 2.01649, 1.98308, 1.94056,
      1.90716, 1.87679, 1.84642, 1.80998, 1.77354, 1.73406, 1.70369, 1.66421, 1.60651, 1.53362,
      1.5154, 1.49111, 1.46985, 1.44252, 1.4243, 1.39696, 1.36356, 1.318, 1.26941, 1.21171,
      1.16616, 1.13275, 1.09935, 1.12972, 1.11453, 1.08416, 1.05683, 1.02343, 0.993059, 0.956616,
      0.929284, 0.895879, 0.87462, 0.835141, 0.801735, 0.77744, 0.747072, 0.704555, 0.67115, 0.640781,
      0.595228, 0.570933, 0.540564, 0.510195, 0.473753, 0.443384, 0.419089, 0.394794, 0.373536, 0.34013,
      0.318872, 0.276356, 0.252061, 0.203471, 0.185249, 0.163991, 0.142733, 0.127549, 0.112364, 0.0911063,
      0.072885, 0.0577007, 0.0425163, 0.0303688, 0.024295, 0.00911063, 0.00607375, 0.00607375, 0.00303688, 0.00303688,
      0.00911063, 0.00911063 };

/*  
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.844*eV, 3.06*eV, 3.4*eV, 4.1*eV};// 6.5*eV};
  G4double RefractiveIndex[nEntries_RI] =
    { 1.79, 1.80, 1.806, 1.813, 1.822, 1.833, 1.84, 1.86}; //1.96};
*/
  const G4int nEntries_RI = 8;
  G4double PhotonEnergy_RI[nEntries_RI] =  { 0.1*eV, 1.0*eV, 2.26*eV, 2.55*eV, 2.84*eV, 3.06*eV, 3.4*eV, 4.1*eV};
  G4double RefractiveIndex[nEntries_RI] =  { 1.81,   1.81,   1.81,    1.81,    1.81,    1.81,    1.81,   1.81};


  //G4double Rayleigh[nEntries_RI] =
  //  { 138.*mm, 138.*mm, 138.*mm};

  
  const G4int nEntries_ABS = 85;  
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 4.42857*eV, 4.35088*eV, 4.27586*eV, 4.20339*eV, 4.13333*eV, 4.06557*eV, 4*eV, 3.93651*eV, 3.875*eV, 3.81538*eV, 
      3.75758*eV, 3.70149*eV, 3.64706*eV, 3.5942*eV, 3.54286*eV, 3.49296*eV, 3.44444*eV, 3.39726*eV, 3.35135*eV, 
      3.30667*eV, 3.26316*eV, 3.22078*eV, 3.17949*eV, 3.13924*eV, 3.1*eV, 3.06173*eV, 3.02439*eV, 2.98795*eV, 
      2.95238*eV, 2.91765*eV, 2.88372*eV, 2.85057*eV, 2.81818*eV, 2.78652*eV, 2.75556*eV, 2.72527*eV, 2.69565*eV, 
      2.66667*eV, 2.6383*eV, 2.61053*eV, 2.58333*eV, 2.5567*eV, 2.53061*eV, 2.50505*eV, 2.48*eV, 2.45545*eV, 
      2.43137*eV, 2.40777*eV, 2.38462*eV, 2.3619*eV, 2.33962*eV, 2.31776*eV, 2.2963*eV, 2.27523*eV, 2.25455*eV, 
      2.23423*eV, 2.21429*eV, 2.19469*eV, 2.17544*eV, 2.15652*eV, 2.13793*eV, 2.11966*eV, 2.10169*eV, 2.08403*eV, 
      2.06667*eV, 2.04959*eV, 2.03279*eV, 2.01626*eV, 2*eV, 1.984*eV, 1.96825*eV, 1.95276*eV, 1.9375*eV, 1.92248*eV, 
      1.90769*eV, 1.89313*eV, 1.87879*eV, 1.86466*eV, 1.85075*eV, 1.83704*eV, 1.82353*eV, 1.81022*eV, 1.7971*eV, 1.78417*eV, 1.0*eV
    };
      
  G4double Absorption[nEntries_ABS] =
    { 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 
      0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm,
      0.5*mm, 0.93828*mm, 1.92292*mm, 3.19852*mm, 7.84337*mm, 19.4627*mm, 43.2508*mm, 90.9185*mm, 161.725*mm, 
      267.35*mm, 325.724*mm, 369.759*mm, 527.641*mm, 433.773*mm, 451.487*mm, 510.549*mm, 352.191*mm, 432.319*mm, 
      499.728*mm, 416.34*mm, 433.836*mm, 580.102*mm, 434.027*mm, 469.628*mm, 458.186*mm, 500.534*mm, 487.4*mm, 
      426.389*mm, 486.997*mm, 451.24*mm, 500.008*mm, 500.425*mm, 500.712*mm, 470.476*mm, 475.17*mm, 500.267*mm, 
      493.917*mm, 500.192*mm, 507.001*mm, 544.182*mm, 500.366*mm, 500.553*mm, 500.793*mm, 500.08*mm, 500.38*mm, 
      451.895*mm, 500.089*mm, 500.602*mm, 500.547*mm, 500.33*mm, 500.171*mm, 500.461*mm, 500.85*mm, 500.493*mm, 
      500.18*mm, 500.011*mm, 500.956*mm, 500.05*mm, 500.495*mm, 500.982*mm, 500.846*mm, 500.64*mm, 500.212*mm, 
      500.801*mm, 500.457*mm, 500.729*mm, 500.8*mm
    };
  
  const G4int nEntries_SCY = 12;
  G4double ElectronEnergy_SCY[nEntries_SCY] =
    { 0.000*MeV, 0.015*MeV, 0.020*MeV, 0.030*MeV,
      0.040*MeV, 0.060*MeV, 0.080*MeV, 0.090*MeV,
      0.105*MeV, 0.300*MeV, 0.500*MeV, 1.000*MeV };
  G4double ScintilYield[nEntries_SCY] =
    { 0.10, 0.46, 0.60, 0.68,
      0.74, 0.80, 0.82, 0.84,
      0.87,  0.96,  0.98,  1.00 };
  for(int i=0; i < nEntries_SCY; i++)    ScintilYield[i] = 40000.*MeV*ScintilYield[i]*ElectronEnergy_SCY[i];
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("FASTCOMPONENT", PhotonEnergy_FAST, FastComponent,   nEntries_FAST);
  myMPT->AddProperty("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_ABS);
  myMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", ElectronEnergy_SCY, ScintilYield, nEntries_SCY);
  myMPT->AddConstProperty("SCINTILLATIONYIELD", 40000./MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE", 3.2);
  myMPT->AddConstProperty("FASTTIMECONSTANT", 40.*ns);
  myMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.06*ns);
  myMPT->AddConstProperty("YIELDRATIO", 1.0);
  
  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}



G4Material* MyMaterials::PWO()
{
  G4double a, z, density;
  G4Element* Pb = new G4Element("Lead",     "Pb", z = 82., a = 207.21*g/mole);
  G4Element* W  = new G4Element("Tungsten", "W",  z = 74., a = 183.85*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O",  z =  8., a =  16.00*g/mole);
  
  G4Material* mat = new G4Material ("PWO", density = 8.28*g/cm3,3);
  mat->AddElement (Pb, 1);
  mat->AddElement (W, 1);
  mat->AddElement (O, 4);
  
  //my latest lab measurement
      const int nEntries_FAST = 507;
      G4double PhotonEnergy_FAST[nEntries_FAST] = { 1.90769*eV, 2.06667*eV, 2.16783*eV, 2.16973*eV, 2.17163*eV, 2.17353*eV, 2.17544*eV, 2.17735*eV, 2.17926*eV, 2.18118*eV, 2.1831*eV, 2.18502*eV,
       2.18695*eV, 2.18888*eV, 2.19081*eV, 2.19275*eV, 2.19469*eV, 2.19663*eV, 2.19858*eV, 2.20053*eV, 2.20249*eV, 2.20444*eV, 2.20641*eV, 2.20837*eV, 2.21034*eV, 2.21231*eV, 2.21429*eV,
        2.21626*eV, 2.21825*eV, 2.22023*eV, 2.22222*eV, 2.22422*eV, 2.22621*eV, 2.22821*eV, 2.23022*eV, 2.23222*eV, 2.23423*eV, 2.23625*eV, 2.23827*eV, 2.24029*eV, 2.24231*eV, 2.24434*eV,
         2.24638*eV, 2.24841*eV, 2.25045*eV, 2.2525*eV, 2.25455*eV, 2.2566*eV, 2.25865*eV, 2.26071*eV, 2.26277*eV, 2.26484*eV, 2.26691*eV, 2.26898*eV, 2.27106*eV, 2.27314*eV, 2.27523*eV,
          2.27732*eV, 2.27941*eV, 2.28151*eV, 2.28361*eV, 2.28571*eV, 2.28782*eV, 2.28994*eV, 2.29205*eV, 2.29417*eV, 2.2963*eV, 2.29842*eV, 2.30056*eV, 2.30269*eV, 2.30483*eV, 2.30698*eV,
           2.30912*eV, 2.31128*eV, 2.31343*eV, 2.31559*eV, 2.31776*eV, 2.31993*eV, 2.3221*eV, 2.32427*eV, 2.32645*eV, 2.32864*eV, 2.33083*eV, 2.33302*eV, 2.33522*eV, 2.33742*eV, 2.33962*eV,
            2.34183*eV, 2.34405*eV, 2.34626*eV, 2.34848*eV, 2.35071*eV, 2.35294*eV, 2.35518*eV, 2.35741*eV, 2.35966*eV, 2.3619*eV, 2.36416*eV, 2.36641*eV, 2.36867*eV, 2.37094*eV, 2.37321*eV,
             2.37548*eV, 2.37776*eV, 2.38004*eV, 2.38232*eV, 2.38462*eV, 2.38691*eV, 2.38921*eV, 2.39151*eV, 2.39382*eV, 2.39614*eV, 2.39845*eV, 2.40077*eV, 2.4031*eV, 2.40543*eV, 2.40777*eV,
              2.41011*eV, 2.41245*eV, 2.4148*eV, 2.41715*eV, 2.41951*eV, 2.42188*eV, 2.42424*eV, 2.42661*eV, 2.42899*eV, 2.43137*eV, 2.43376*eV, 2.43615*eV, 2.43854*eV, 2.44094*eV, 2.44335*eV,
               2.44576*eV, 2.44817*eV, 2.45059*eV, 2.45302*eV, 2.45545*eV, 2.45788*eV, 2.46032*eV, 2.46276*eV, 2.46521*eV, 2.46766*eV, 2.47012*eV, 2.47258*eV, 2.47505*eV, 2.47752*eV, 2.48*eV,
                2.48248*eV, 2.48497*eV, 2.48746*eV, 2.48996*eV, 2.49246*eV, 2.49497*eV, 2.49748*eV, 2.5*eV, 2.50252*eV, 2.50505*eV, 2.50758*eV, 2.51012*eV, 2.51266*eV, 2.51521*eV, 2.51777*eV,
                 2.52033*eV, 2.52289*eV, 2.52546*eV, 2.52803*eV, 2.53061*eV, 2.5332*eV, 2.53579*eV, 2.53838*eV, 2.54098*eV, 2.54359*eV, 2.5462*eV, 2.54882*eV, 2.55144*eV, 2.55407*eV, 2.5567*eV,
                  2.55934*eV, 2.56198*eV, 2.56463*eV, 2.56729*eV, 2.56995*eV, 2.57261*eV, 2.57529*eV, 2.57796*eV, 2.58065*eV, 2.58333*eV, 2.58603*eV, 2.58873*eV, 2.59143*eV, 2.59414*eV,
                   2.59686*eV, 2.59958*eV, 2.60231*eV, 2.60504*eV, 2.60778*eV, 2.61053*eV, 2.61328*eV, 2.61603*eV, 2.6188*eV, 2.62156*eV, 2.62434*eV, 2.62712*eV, 2.6299*eV, 2.6327*eV, 2.63549*eV, 2.6383*eV, 2.64111*eV, 2.64392*eV,                   2.64674*eV, 2.64957*eV, 2.65241*eV, 2.65525*eV, 2.65809*eV, 2.66094*eV, 2.6638*eV, 2.66667*eV, 2.66954*eV, 2.67241*eV, 2.6753*eV, 2.67819*eV, 2.68108*eV, 2.68398*eV, 2.68689*eV, 2.6898*eV, 2.69273*eV, 2.69565*eV, 2.69859*eV, 2.70153*eV, 2.70447*eV, 2.70742*eV, 2.71038*eV, 2.71335*eV, 2.71632*eV, 2.7193*eV, 2.72228*eV, 2.72527*eV, 2.72827*eV, 2.73128*eV, 2.73429*eV, 2.73731*eV, 2.74033*eV, 2.74336*eV, 2.7464*eV, 2.74945*eV, 2.7525*eV, 2.75556*eV, 2.75862*eV, 2.76169*eV, 2.76477*eV, 2.76786*eV, 2.77095*eV, 2.77405*eV, 2.77716*eV, 2.78027*eV, 2.78339*eV, 2.78652*eV, 2.78965*eV, 2.79279*eV, 2.79594*eV, 2.7991*eV, 2.80226*eV, 2.80543*eV, 2.80861*eV, 2.81179*eV, 2.81498*eV, 2.81818*eV, 2.82139*eV, 2.8246*eV, 2.82782*eV, 2.83105*eV, 2.83429*eV, 2.83753*eV, 2.84078*eV, 2.84404*eV, 2.8473*eV, 2.85057*eV, 2.85386*eV, 2.85714*eV, 2.86044*eV, 2.86374*eV, 2.86705*eV, 2.87037*eV, 2.8737*eV, 2.87703*eV, 2.88037*eV, 2.88372*eV, 2.88708*eV, 2.89044*eV, 2.89382*eV, 2.8972*eV, 2.90058*eV, 2.90398*eV, 2.90739*eV, 2.9108*eV, 2.91422*eV, 2.91765*eV, 2.92108*eV, 2.92453*eV, 2.92798*eV, 2.93144*eV, 2.93491*eV, 2.93839*eV, 2.94187*eV, 2.94537*eV, 2.94887*eV, 2.95238*eV, 2.9559*eV, 2.95943*eV, 2.96296*eV, 2.96651*eV, 2.97006*eV, 2.97362*eV, 2.97719*eV, 2.98077*eV, 2.98436*eV, 2.98795*eV, 2.99156*eV, 2.99517*eV, 2.99879*eV, 3.00242*eV, 3.00606*eV, 3.00971*eV, 3.01337*eV, 3.01703*eV, 3.02071*eV, 3.02439*eV, 3.02808*eV, 3.03178*eV, 3.0355*eV, 3.03922*eV, 3.04294*eV, 3.04668*eV, 3.05043*eV, 3.05419*eV, 3.05795*eV, 3.06173*eV, 3.06551*eV, 3.06931*eV, 3.07311*eV, 3.07692*eV, 3.08075*eV, 3.08458*eV, 3.08842*eV, 3.09227*eV, 3.09613*eV, 3.1*eV, 3.10388*eV, 3.10777*eV, 3.11167*eV, 3.11558*eV, 3.1195*eV, 3.12343*eV, 3.12736*eV, 3.13131*eV, 3.13527*eV, 3.13924*eV, 3.14322*eV, 3.14721*eV, 3.15121*eV, 3.15522*eV, 3.15924*eV, 3.16327*eV, 3.16731*eV, 3.17136*eV, 3.17542*eV, 3.17949*eV, 3.18357*eV, 3.18766*eV, 3.19176*eV, 3.19588*eV, 3.2*eV, 3.20413*eV, 3.20828*eV, 3.21244*eV, 3.2166*eV, 3.22078*eV, 3.22497*eV, 3.22917*eV, 3.23338*eV, 3.2376*eV, 3.24183*eV, 3.24607*eV, 3.25033*eV, 3.25459*eV, 3.25887*eV, 3.26316*eV, 3.26746*eV, 3.27177*eV, 3.27609*eV, 3.28042*eV, 3.28477*eV, 3.28912*eV, 3.29349*eV, 3.29787*eV, 3.30226*eV, 3.30667*eV, 3.31108*eV, 3.31551*eV, 3.31995*eV, 3.3244*eV, 3.32886*eV, 3.33333*eV, 3.33782*eV, 3.34232*eV, 3.34683*eV, 3.35135*eV, 3.35589*eV, 3.36043*eV, 3.36499*eV, 3.36957*eV, 3.37415*eV, 3.37875*eV, 3.38336*eV, 3.38798*eV, 3.39261*eV, 3.39726*eV, 3.40192*eV, 3.40659*eV, 3.41128*eV, 3.41598*eV, 3.42069*eV, 3.42541*eV, 3.43015*eV, 3.4349*eV, 3.43967*eV, 3.44444*eV, 3.44924*eV, 3.45404*eV, 3.45886*eV, 3.46369*eV, 3.46853*eV, 3.47339*eV, 3.47826*eV, 3.48315*eV, 3.48805*eV, 3.49296*eV, 3.49788*eV, 3.50282*eV, 3.50778*eV, 3.51275*eV, 3.51773*eV, 3.52273*eV, 3.52774*eV, 3.53276*eV, 3.5378*eV, 3.54286*eV, 3.54793*eV, 3.55301*eV, 3.55811*eV, 3.56322*eV, 3.56835*eV, 3.57349*eV, 3.57864*eV, 3.58382*eV, 3.589*eV, 3.5942*eV, 3.59942*eV, 3.60465*eV, 3.6099*eV, 3.61516*eV, 3.62044*eV, 3.62573*eV, 3.63104*eV, 3.63636*eV, 3.6417*eV, 3.64706*eV, 3.65243*eV, 3.65782*eV, 3.66322*eV, 3.66864*eV, 3.67407*eV, 3.67953*eV, 3.68499*eV, 3.69048*eV, 3.69598*eV, 3.70149*eV, 3.70703*eV, 3.71257*eV, 3.71814*eV, 3.72372*eV, 3.72932*eV, 3.73494*eV, 3.74057*eV, 3.74622*eV, 3.75189*eV, 3.75758*eV, 3.76328*eV, 3.769*eV, 3.77473*eV, 3.78049*eV, 3.78626*eV, 3.79205*eV, 3.79786*eV, 3.80368*eV, 3.80952*eV, 3.81538*eV, 3.82126*eV, 3.82716*eV, 3.83308*eV, 3.83901*eV, 3.84496*eV, 3.85093*eV, 3.85692*eV, 3.86293*eV, 3.86895*eV};
      
      G4double FastComponent[nEntries_FAST] ={ 0.00052381, 0.0152381, 0.0293924, 0.0269467, 0.0256848, 0.0254981, 0.0260076, 0.0268133, 0.0273648, 0.0280457, 0.0287067, 0.0289924, 0.029, 0.0288914, 0.0289476, 0.0291505, 0.0292914, 0.0290752, 0.0287514, 0.0286714, 0.0288343, 0.0293133, 0.0302143, 0.031561, 0.0330762, 0.034459, 0.0351705, 0.03478, 0.0340324, 0.0329552, 0.0314514, 0.0300962, 0.0292152, 0.028841, 0.0289524, 0.0294419, 0.0297876, 0.0303133, 0.0310962, 0.0317743, 0.0322238, 0.0325276, 0.0329733, 0.03336, 0.0337895, 0.0342257, 0.0347362, 0.0354695, 0.0363, 0.0369238, 0.036979, 0.036579, 0.0356076, 0.0343095, 0.032999, 0.0316876, 0.0306171, 0.0302533, 0.0308752, 0.0322067, 0.0339943, 0.0360067, 0.0384286, 0.0411333, 0.04334, 0.0450924, 0.0462695, 0.0472533, 0.0481019, 0.0483352, 0.0482181, 0.04864, 0.0499019, 0.0517543, 0.0543505, 0.0575267, 0.0607876, 0.0641314, 0.0667838, 0.0683514, 0.0693419, 0.0702543, 0.0710981, 0.0720552, 0.0736676, 0.0752762, 0.0773286, 0.0791752, 0.0807333, 0.082079, 0.0833629, 0.0845933, 0.0859524, 0.0877581, 0.0892943, 0.0910914, 0.0929019, 0.0952905, 0.0978371, 0.100682, 0.103224, 0.105718, 0.107852, 0.109241, 0.109795, 0.109354, 0.109503, 0.110292, 0.112409, 0.115256, 0.11961, 0.124909, 0.130732, 0.13613, 0.140314, 0.144119, 0.14665, 0.148469, 0.149407, 0.150265, 0.151147, 0.152253, 0.154354, 0.15737, 0.161138, 0.164838, 0.168873, 0.172971, 0.177416, 0.181659, 0.185039, 0.18823, 0.191871, 0.195512, 0.198159, 0.20033, 0.202124, 0.203936, 0.205858, 0.207766, 0.20994, 0.212648, 0.216484, 0.22057, 0.22497, 0.229337, 0.232884, 0.23596, 0.239326, 0.242921, 0.245595, 0.248666, 0.2521, 0.255899, 0.260253, 0.264172, 0.268089, 0.27221, 0.27688, 0.280961, 0.284639, 0.28845, 0.291873, 0.294659, 0.29668, 0.298346, 0.300048, 0.302431, 0.305583, 0.309162, 0.313943, 0.319651, 0.324986, 0.329565, 0.333289, 0.336565, 0.33924, 0.341681, 0.343442, 0.345074, 0.347293, 0.349824, 0.352279, 0.354626, 0.357285, 0.360091, 0.363057, 0.366477, 0.37011, 0.373876, 0.377648, 0.380938, 0.384106, 0.387288, 0.390137, 0.39199, 0.393832, 0.39634, 0.399255, 0.402285, 0.405265, 0.408943, 0.412963, 0.417683, 0.42203, 0.425902, 0.430032, 0.434444, 0.438648, 0.442674, 0.445648, 0.44771, 0.450105, 0.452818, 0.455475, 0.459257, 0.46497, 0.471776, 0.4804, 0.489469, 0.497992, 0.507131, 0.516506, 0.525018, 0.532744, 0.540899, 0.548129, 0.555484, 0.562471, 0.568386, 0.574721, 0.581799, 0.58892, 0.595488, 0.603247, 0.611145, 0.619414, 0.628492, 0.637634, 0.646924, 0.656813, 0.66681, 0.675362, 0.683059, 0.689791, 0.69601, 0.702039, 0.708148, 0.714768, 0.722988, 0.732344, 0.74167, 0.75073, 0.759785, 0.768961, 0.777719, 0.78556, 0.792373, 0.79914, 0.80549, 0.81113, 0.815989, 0.821, 0.826681, 0.832696, 0.838765, 0.844266, 0.850464, 0.85645, 0.862316, 0.86791, 0.873631, 0.879935, 0.887078, 0.895343, 0.903276, 0.910896, 0.917349, 0.922744, 0.927535, 0.931457, 0.93481, 0.938199, 0.942202, 0.946601, 0.950905, 0.954545, 0.957467, 0.96089, 0.964401, 0.967275, 0.968981, 0.97016, 0.971976, 0.973287, 0.973908, 0.973841, 0.97475, 0.977297, 0.981698, 0.986047, 0.990052, 0.994456, 0.998224, 1.0007, 1.00106, 1.00068, 1.00017, 1.00024, 1.00012, 0.99895, 0.997525, 0.995587, 0.994011, 0.991989, 0.990154, 0.988663, 0.987812, 0.988681, 0.990095, 0.990958, 0.990571, 0.990366, 0.989786, 0.988225, 0.985411, 0.981308, 0.977582, 0.973715, 0.96889, 0.963269, 0.958067, 0.954727, 0.952326, 0.95059, 0.949261, 0.949456, 0.949965, 0.948949, 0.946274, 0.94231, 0.937687, 0.93214, 0.926558, 0.920441, 0.915683, 0.912037, 0.908574, 0.904968, 0.901723, 0.898805, 0.895016, 0.891226, 0.886231, 0.880648, 0.874508, 0.867607, 0.859656, 0.851617, 0.844196, 0.83623, 0.828706, 0.822149, 0.817083, 0.812778, 0.808321, 0.803222, 0.798333, 0.793735, 0.787804, 0.780487, 0.772463, 0.764901, 0.75783, 0.750741, 0.743238, 0.737445, 0.73221, 0.725617, 0.717075, 0.707011, 0.696076, 0.684175, 0.670404, 0.65462, 0.640174, 0.627405, 0.6152, 0.603558, 0.592237, 0.58155, 0.57139, 0.559997, 0.546431, 0.532181, 0.517833, 0.503294, 0.488553, 0.474083, 0.460749, 0.449591, 0.439908, 0.431058, 0.42282, 0.414699, 0.406633, 0.398634, 0.39069, 0.382162, 0.373201, 0.364355, 0.355435, 0.346777, 0.337376, 0.32759, 0.31762, 0.307741, 0.29736, 0.286301, 0.274514, 0.262031, 0.249232, 0.235434, 0.220187, 0.204732, 0.18955, 0.174084, 0.158599, 0.143255, 0.128051, 0.113468, 0.0996657, 0.0862448, 0.0741762, 0.0637238, 0.054501, 0.0470933, 0.0412562, 0.0365495, 0.0324981, 0.0291943, 0.0259467, 0.0229876, 0.0201476, 0.0172495, 0.0144133, 0.0121181, 0.010861, 0.0100343, 0.00974476, 0.0103733, 0.0119886, 0.01364, 0.0151286, 0.0161257, 0.0168276, 0.0176267, 0.0177667, 0.0169867, 0.01598, 0.015241, 0.0144143, 0.0135886, 0.0125457, 0.0115524, 0.0113305, 0.0114295, 0.0114038, 0.0114352, 0.01208, 0.0132114, 0.0141905, 0.0147667, 0.0149648, 0.0148695, 0.0140505, 0.0127952, 0.0109514, 0.00864667, 0.00670762, 0.00527143, 0.0046019, 0.00473524, 0.00552476, 0.0065, 0.00768667, 0.0084381, 0.00831333, 0.00752286, 0.0062181, 0.00454952, 0.00287905, 0.00136476, 0.000487619, 0.000487619, 0.000514286, 0.000467619, 0.000337143, 0.00047619, 0.00104, 0.00124, 0.000652381, 0.0015, 0.00581905, 0.0120495, 0.0200286};
      


  const G4int nEntries_RI = 32;
  G4double PhotonEnergy_RI[nEntries_RI] =    
{ 0.0001*eV , 0.5*eV, 1.0*eV, 1.1*eV, 1.2*eV, 1.3*eV, 1.4*eV, 1.5*eV, 1.6*eV, 1.7*eV, 1.8*eV, 1.9*eV, 2.0*eV, 2.1*eV, 2.2*eV, 2.3*eV, 2.4*eV, 2.5*eV, 2.6*eV, 2.7*eV, 2.8*eV, 2.9*eV, 3.0*eV, 3.1*eV, 3.2*eV, 3.3*eV, 3.4*eV, 3.5*eV, 3.6*eV, 3.7*eV, 3.8*eV, 3.9*eV};
  G4double RefractiveIndex[nEntries_RI] =    
{ 2.17701, 2.1804, 2.19103, 2.19414, 2.19762, 2.20149, 2.20577, 2.2105, 2.21571, 2.22144, 2.22773, 2.23463, 2.24221, 2.25053, 2.25968, 2.26974, 2.28084, 2.2931, 2.30668, 2.32177, 2.3386, 2.35745, 2.37866, 2.40268, 2.43006, 2.4615, 2.49794, 2.54063, 2.59128, 2.6523, 2.72722, 2.82141};


  const G4int nEntries_ABS = 50;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    {       
      1.56942 * eV, 1.58954 * eV, 1.61018 * eV, 1.63137 * eV, 1.65312 * eV, 1.67546 * eV, 1.69841 * eV, 1.722 * eV, 1.74626 * eV, 1.7712 * eV,
      1.79687 * eV, 1.8233 * eV, 1.85051 * eV, 1.87855 * eV, 1.90745 * eV, 1.93725 * eV, 1.968 * eV, 1.99974 * eV, 2.03253 * eV, 2.0664 * eV,
      2.10143 * eV, 2.13766 * eV, 2.17516 * eV, 2.214 * eV, 2.25426 * eV, 2.296 * eV, 2.33932 * eV, 2.38431 * eV, 2.43106 * eV, 2.47968 * eV,
      2.53029 * eV, 2.583 * eV, 2.63796 * eV, 2.69531 * eV, 2.7552 * eV, 2.81782 * eV, 2.88335 * eV, 2.952 * eV, 3.024 * eV, 3.0996 * eV,
      3.17908 * eV, 3.26274 * eV, 3.35092 * eV, 3.44401 * eV, 3.54241 * eV, 3.64659 * eV, 3.7571 * eV, 3.87451 * eV, 3.99949 * eV, 4.13281 * eV };

      double att0 = 1;
      G4double Absorption[nEntries_ABS] =
    { 
      390.8 *att0*mm, 390.9 *att0*mm, 390.7 *att0*mm, 390.2 *att0*mm, 390.7 *att0*mm, 390.5 *att0*mm, 390.6 *att0*mm, 390.7 *att0*mm, 390.3 *att0*mm, 390.2 *att0*mm,
      390.8 *att0*mm, 390.9 *att0*mm, 390.7 *att0*mm, 390.2 *att0*mm, 390.7 *att0*mm, 390.5 *att0*mm, 390.6 *att0*mm, 390.7 *att0*mm, 390.3 *att0*mm, 390.2 *att0*mm,
      390.8 *att0*mm, 390.5 *att0*mm, 390.3 *att0*mm, 390.4 *att0*mm, 390.3 *att0*mm, 390.8 *att0*mm, 390.9 *att0*mm, 390.9 *att0*mm, 390.4 *att0*mm, 390.9 *att0*mm,
      390 *att0*mm, 390.2 *att0*mm, 390.1 *att0*mm, 345.3 *att0*mm, 298.9 *att0*mm, 256.7 *att0*mm, 219.8 *att0*mm, 185.4 *att0*mm, 150.9 *att0*mm, 116.4 *att0*mm,
      84.8 *att0*mm, 59.4 *att0*mm, 41.1 *att0*mm, 0 *att0*mm, 0 *att0*mm, 0 *att0*mm, 0 *att0*mm, 0 *att0*mm, 0 *att0*mm, 0 *att0*mm };
  
  const G4int nEntries_SCY = 12;
  G4double ElectronEnergy_SCY[nEntries_SCY] =
    { 0.000 * MeV, 0.015 * MeV, 0.020 * MeV, 0.030 * MeV,
      0.040 * MeV, 0.060 * MeV, 0.080 * MeV, 0.090 * MeV,
      0.105 * MeV, 0.300 * MeV, 0.500 * MeV, 1.000 * MeV };
      
  G4double ScintilYield[nEntries_SCY] =
    { 0.10, 0.46, 0.60, 0.68,
      0.74, 0.80, 0.82, 0.84,
      0.87,  0.96,  0.98,  1.00 };
      
  for(int i = 0; i < nEntries_SCY; i++)    ScintilYield[i] = 0.3 * MeV * ScintilYield[i] * ElectronEnergy_SCY[i];
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty ("FASTCOMPONENT", PhotonEnergy_FAST, FastComponent,   nEntries_FAST);
  myMPT->AddProperty ("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  myMPT->AddProperty ("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_ABS);
  myMPT->AddProperty ("ELECTRONSCINTILLATIONYIELD", ElectronEnergy_SCY, ScintilYield, nEntries_SCY);
//  myMPT->AddConstProperty ("SCINTILLATIONYIELD", 1000/MeV );//for 10% of detected cherenkov in APDs
  myMPT->AddConstProperty ("SCINTILLATIONYIELD", 450/MeV );//for 10% of detected cherenkov in PMTs
  myMPT->AddConstProperty ("RESOLUTIONSCALE", 1.0); //3.2 default value
  myMPT->AddConstProperty ("FASTTIMECONSTANT", 5.*ns);
  myMPT->AddConstProperty ("SLOWTIMECONSTANT", 15.*ns);
  myMPT->AddConstProperty ("YIELDRATIO", 0.3);
  myMPT->AddConstProperty ("FASTSCINTILLATIONRISETIME", 0.01 * ns);	//careful on rise time
  
  mat->SetMaterialPropertiesTable (myMPT);
  
  return mat;
}

G4Material* MyMaterials::ToyLSO(float toy_ly, float toy_decay, float toy_rise)
{
  G4double a, z, density;
  G4Element*  O = new G4Element("Oxygen",   "O",  z= 8., a= 16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71., a=174.97*g/mole);
  G4Element* Si = new G4Element("Silicon",  "Si", z=14., a= 28.09*g/mole);
  
  G4Material* mat = new G4Material("ToyLSO", density=7.4*g/cm3,3,kStateSolid);
  mat->AddElement(Lu,2);
  mat->AddElement(Si,1);
  mat->AddElement(O,5);
  
  const G4int nEntries_FAST = 261;
  G4double PhotonEnergy_FAST[nEntries_FAST] = 
    { 1.77169*eV, 1.77266*eV, 1.77558*eV, 1.77851*eV, 1.78145*eV, 1.78539*eV, 1.79033*eV, 1.7963*eV, 1.80231*eV, 1.80836*eV,
      1.81445*eV, 1.82058*eV, 1.82882*eV, 1.83401*eV, 1.84553*eV, 1.85293*eV, 1.86147*eV, 1.869*eV, 1.87769*eV, 1.89308*eV,
      1.90536*eV, 1.92007*eV, 1.93039*eV, 1.94901*eV, 1.95846*eV, 1.9668*eV, 1.97884*eV, 1.99102*eV, 2.00088*eV, 2.01209*eV,
      2.02596*eV, 2.03617*eV, 2.04519*eV, 2.0569*eV, 2.06611*eV, 2.0794*eV, 2.09151*eV, 2.10239*eV, 2.112*eV, 2.1231*eV,
      2.13431*eV, 2.14565*eV, 2.15566*eV, 2.16868*eV, 2.18038*eV, 2.19519*eV, 2.21171*eV, 2.2193*eV, 2.23619*eV, 2.23464*eV,
      2.24395*eV, 2.25806*eV, 2.27234*eV, 2.28358*eV, 2.29493*eV, 2.30475*eV, 2.31631*eV, 2.32463*eV, 2.33134*eV, 2.33809*eV,
      2.34487*eV, 2.35856*eV, 2.36719*eV, 2.37939*eV, 2.38642*eV, 2.40238*eV, 2.41134*eV, 2.424*eV, 2.43312*eV, 2.44047*eV,
      2.44786*eV, 2.46278*eV, 2.47788*eV, 2.48741*eV, 2.49317*eV, 2.49702*eV, 2.50282*eV, 2.50865*eV, 2.5145*eV, 2.52038*eV,
      2.52432*eV, 2.53223*eV, 2.5362*eV, 2.54619*eV, 2.55424*eV, 2.56031*eV, 2.56437*eV, 2.57049*eV, 2.57663*eV, 2.58487*eV,
      2.59317*eV, 2.59734*eV, 2.60571*eV, 2.61414*eV, 2.61414*eV, 2.61837*eV, 2.62262*eV, 2.62475*eV, 2.62902*eV, 2.63331*eV,
      2.63545*eV, 2.63976*eV, 2.64191*eV, 2.64841*eV, 2.65493*eV, 2.6593*eV, 2.66149*eV, 2.66588*eV, 2.67914*eV, 2.67914*eV,
      2.68136*eV, 2.68136*eV, 2.68359*eV, 2.68805*eV, 2.68805*eV, 2.68805*eV, 2.69477*eV, 2.69477*eV, 2.69702*eV, 2.70153*eV,
      2.70605*eV, 2.71286*eV, 2.71742*eV, 2.71971*eV, 2.722*eV, 2.722*eV, 2.72429*eV, 2.72889*eV, 2.72889*eV, 2.73351*eV,
      2.73814*eV, 2.74279*eV, 2.74512*eV, 2.74979*eV, 2.75213*eV, 2.75447*eV, 2.75917*eV, 2.75682*eV, 2.76389*eV, 2.76626*eV,
      2.76389*eV, 2.76626*eV, 2.77338*eV, 2.77576*eV, 2.78533*eV, 2.79255*eV, 2.79738*eV, 2.80223*eV, 2.80466*eV, 2.80709*eV,
      2.80953*eV, 2.80953*eV, 2.81934*eV, 2.8218*eV, 2.82673*eV, 2.83168*eV, 2.84164*eV, 2.84916*eV, 2.85419*eV, 2.8643*eV,
      2.86684*eV, 2.87449*eV, 2.87705*eV, 2.87961*eV, 2.88475*eV, 2.88733*eV, 2.8925*eV, 2.89509*eV, 2.90028*eV, 2.90549*eV,
      2.90811*eV, 2.91073*eV, 2.91335*eV, 2.91335*eV, 2.91335*eV, 2.91861*eV, 2.92125*eV, 2.92125*eV, 2.92389*eV, 2.92654*eV,
      2.92654*eV, 2.92919*eV, 2.92919*eV, 2.93185*eV, 2.93451*eV, 2.93717*eV, 2.93985*eV, 2.94252*eV, 2.9452*eV, 2.94789*eV,
      2.94789*eV, 2.94789*eV, 2.95058*eV, 2.95868*eV, 2.96411*eV, 2.96955*eV, 2.97228*eV, 2.97228*eV, 2.96955*eV, 2.97228*eV,
      2.97502*eV, 2.97776*eV, 2.97502*eV, 2.9805*eV, 2.9805*eV, 2.9805*eV, 2.98601*eV, 2.99154*eV, 2.99431*eV, 2.99431*eV,
      2.99708*eV, 2.99431*eV, 2.99708*eV, 3.00544*eV, 3.00824*eV, 3.00824*eV, 3.00824*eV, 3.00824*eV, 3.01385*eV, 3.0223*eV,
      3.02797*eV, 3.03081*eV, 3.02797*eV, 3.03365*eV, 3.03081*eV, 3.03081*eV, 3.0365*eV, 3.03935*eV, 3.04221*eV, 3.04795*eV,
      3.04795*eV, 3.05083*eV, 3.05371*eV, 3.05949*eV, 3.06239*eV, 3.06529*eV, 3.0682*eV, 3.06529*eV, 3.07112*eV, 3.0682*eV,
      3.07696*eV, 3.08283*eV, 3.0976*eV, 3.09464*eV, 3.09464*eV, 3.10653*eV, 3.11252*eV, 3.11852*eV, 3.12757*eV, 3.13668*eV,
      3.14583*eV, 3.15813*eV, 3.16741*eV, 3.17675*eV, 3.20828*eV, 3.23719*eV, 3.26664*eV, 3.28656*eV, 3.31351*eV, 3.34783*eV,
      3.38287*eV };
  G4double FastComponent[nEntries_FAST] = 
    { 0.011691, 0.011691, 0.011691, 0.0146138, 0.0146138, 0.0146138, 0.011691, 0.011691, 0.00876827, 0.00876827,
      0.00584551, 0.00584551, 0.00584551, 0.00292276, 0.00876827, 0.0146138, 0.0146138, 0.0146138, 0.0204593, 0.023382,
      0.0263048, 0.0204593, 0.0204593, 0.023382, 0.0292276, 0.0321503, 0.0350731, 0.0379958, 0.0379958, 0.0379958,
      0.0350731, 0.0379958, 0.0409186, 0.0438413, 0.0526096, 0.0584551, 0.0643006, 0.0730689, 0.0730689, 0.0818372,
      0.0906054, 0.0964509, 0.0993737, 0.105219, 0.111065, 0.122756, 0.125678, 0.146138, 0.146138, 0.160752,
      0.157829, 0.163674, 0.184134, 0.192902, 0.20167, 0.219207, 0.230898, 0.242589, 0.25428, 0.265971,
      0.274739, 0.292276, 0.306889, 0.315658, 0.321503, 0.350731, 0.368267, 0.385804, 0.397495, 0.415031,
      0.432568, 0.458873, 0.482255, 0.496868, 0.514405, 0.529019, 0.549478, 0.564092, 0.581628, 0.593319,
      0.602088, 0.616701, 0.637161, 0.660543, 0.681002, 0.71023, 0.736534, 0.756994, 0.777453, 0.806681,
      0.844676, 0.868058, 0.891441, 0.9119, 0.938205, 0.955741, 0.984969, 1.0142, 1.03173, 1.05511,
      1.07557, 1.11649, 1.13695, 1.15741, 1.17495, 1.19248, 1.21002, 1.22756, 1.27432, 1.2977,
      1.31524, 1.32985, 1.36785, 1.40292, 1.39415, 1.4, 1.41754, 1.44092, 1.47015, 1.48476,
      1.50814, 1.5286, 1.54906, 1.56952, 1.58998, 1.61921, 1.63967, 1.66597, 1.68935, 1.71566,
      1.73904, 1.76242, 1.77996, 1.80042, 1.8238, 1.83549, 1.85303, 1.8618, 1.87933, 1.89979,
      1.91733, 1.92902, 1.95825, 1.98163, 2.01378, 2.03424, 2.0547, 2.07808, 2.09562, 2.11023,
      2.12484, 2.13361, 2.15407, 2.15699, 2.15992, 2.16576, 2.16868, 2.16868, 2.16284, 2.15699,
      2.14823, 2.13946, 2.12484, 2.11023, 2.08977, 2.06639, 2.04593, 2.02839, 2.01086, 1.98455,
      1.96409, 1.94948, 1.93194, 1.91733, 1.90271, 1.87641, 1.86472, 1.8501, 1.83841, 1.82088,
      1.79749, 1.77119, 1.75073, 1.73027, 1.70689, 1.68058, 1.65428, 1.6309, 1.60167, 1.57244,
      1.55491, 1.53152, 1.50522, 1.47891, 1.45261, 1.43215, 1.40877, 1.38831, 1.362, 1.33862,
      1.31232, 1.28601, 1.27432, 1.25678, 1.21587, 1.19541, 1.17203, 1.14864, 1.12234, 1.10772,
      1.08434, 1.06096, 1.0142, 0.987891, 0.967432, 0.938205, 0.9119, 0.879749, 0.853445, 0.82714,
      0.786221, 0.765762, 0.739457, 0.716075, 0.681002, 0.660543, 0.637161, 0.60501, 0.581628, 0.552401,
      0.531942, 0.505637, 0.485177, 0.458873, 0.435491, 0.412109, 0.379958, 0.356576, 0.336117, 0.309812,
      0.280585, 0.25428, 0.207516, 0.175365, 0.157829, 0.13737, 0.119833, 0.0993737, 0.0759916, 0.0613779,
      0.0526096, 0.0350731, 0.0263048, 0.011691, 0.00876827, 0.00876827, 0.011691, 0.011691, 0.011691, 0.00876827,
      0.011691 };
  
  
  const G4int nEntries_RI = 6;
  G4double PhotonEnergy_RI[nEntries_RI] =
    { 1.0*eV, 2.26*eV, 2.55*eV, 2.844*eV, 3.06*eV, 3.4*eV};
  G4double RefractiveIndex[nEntries_RI] =
 //   { 1.75, 1.82, 1.88};
    { 1.82, 1.806, 1.813, 1.822, 1.833, 1.82};
  
  const G4int nEntries_ABS = 9;
  G4double PhotonEnergy_ABS[nEntries_ABS] =
    { 1.00*eV , 2.82*eV , 2.88*eV , 2.95*eV , 3.02*eV  , 3.10*eV  , 3.18*eV  , 3.26*eV , 4.08*eV };
  G4double Absorption[nEntries_ABS] =
    { 438.*mm , 438.*mm , 413.*mm , 375.*mm , 263.*mm  , 87.5*mm  , 11.5*mm  , 1.0*mm  , 1.0*mm  };
  
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("FASTCOMPONENT", PhotonEnergy_FAST, FastComponent,  nEntries_FAST);
  myMPT->AddProperty("RINDEX",        PhotonEnergy_RI,   RefractiveIndex, nEntries_RI);
  myMPT->AddProperty("ABSLENGTH",     PhotonEnergy_ABS,  Absorption,      nEntries_ABS);
  myMPT->AddConstProperty ("RESOLUTIONSCALE", 1.0);
  myMPT->AddConstProperty ("YIELDRATIO", 1);
  myMPT->AddConstProperty("SCINTILLATIONYIELD",toy_ly/MeV);
  myMPT->AddConstProperty("FASTTIMECONSTANT",toy_decay*ns);
  myMPT->AddConstProperty("FASTSCINTILLATIONRISETIME",toy_rise*ns);

  mat->SetMaterialPropertiesTable(myMPT);
  
  return mat;
}



G4double MyMaterials::CalculateSellmeier (int size, G4double indexZero, G4double *nVec, G4double *lVec, G4double wavelength)
{
  /*------http://gentitfx.fr/SLitrani/code/SLitraniCode/TLitSellmeier.html----*/
  
  float partial = indexZero * indexZero;
  float sum = 0;
  for (int i = 0; i < size; i++)
  {
    sum += nVec[i] * nVec[i] / (1 - lVec[i] * lVec[i] / (wavelength*wavelength));
  }
  
  partial += sum;
  partial += 1;
  
  //G4cout << "Wavelength: " << wavelength << " -> rifr. index: " << sqrt(partial) << G4endl;
  
  return sqrt(partial);
}



G4double MyMaterials::fromEvToNm (G4double energy)
{
  return 1239.84187 / energy;
}

G4double MyMaterials::fromNmToEv (G4double wavelength)
{
  return 1239.84187 / wavelength;
}
