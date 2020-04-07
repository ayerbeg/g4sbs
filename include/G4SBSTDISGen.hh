#ifndef G4SBSTDISGen_HH
#define G4SBSTDISGen_HH

// I will try to make it as clean as possible

//#include "TFile.h"
//#include "TTree.h"
//#include "TChain.h"

#include "G4SystemOfUnits.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "sbstypes.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4SBSIO.hh"
#include "G4SBSEventGen.hh"

using namespace CLHEP;

class G4SBSTDISGen //should I use inheritance (CA)
{
public: 
  G4SBSTDISGen();
  ~G4SBSTDISGen();

  G4SBSIO *fIO;
 
   void Test(Kine_t, Nucl_t, G4LorentzVector, G4LorentzVector);

  void Generate(Kine_t, Nucl_t, G4LorentzVector , G4LorentzVector);

  G4double PhotoD_XS(G4double E_photon); //PhotoDisintegration Cross Section
  
  G4double VXPhoton_flux(double photon_energy, double beam_energy);
 
  G4double QuasiElasticXS(G4double beam_energy, G4int z1, G4int n1, G4int partID, G4double momentum, G4double angle);
  G4double FluxCorrection(G4LorentzVector nf_lab, G4LorentzVector ni, G4LorentzVector ei, G4LorentzVector ei_Nrest);

  G4double ElasticXS(G4double beam_energy, G4double scatter_e_energy, G4double theta, Nucl_t nu);

  G4double MottXS(G4double, G4double);
  G4double DipoleFF();
  G4double GE( Nucl_t);
  G4double GM( Nucl_t);
  G4double tau();
  void qe_epc();
  TGraph *gD;

  G4double alpha()
  {return fine_structure_const;} // fine structure, to shorten the definition (useless?)

  void setFinalNucleon(Nucl_t ft) //ft is a temporary variable
  {tFinalNucl = ft;}

  Nucl_t tGetFinalNucleon(){ return tFinalNucl; }

  // momentum-energy of the scattered electron to be sent to PrimaryGeneratorAction
  G4ThreeVector tGetElectronP(){ return tElectronP; }
  G4double tGetElectronE(){ return tElectronE; }

  // momentum-energy of the scattered nucleon to be sent to PrimaryGeneratorAction
  G4ThreeVector tGetNucleonP(){ return tNucleonP; }
  G4double tGetNucleonE(){ return tNucleonE; }

  G4double tGetnu(){ return nu; }
  G4double tGetQ2(){ return Q2; }
  G4double tGetW2(){ return W2; }
  G4double tGetxbj(){ return xbj; }

  G4LorentzVector tGetelef_lab(){ return tElectronf_lab; }
  G4LorentzVector tGetnucf_lab(){ return tNucleonf_lab; }





private:

  //maybe we need to deleta most of these definitions

  G4double PI;       // 3.1415...
  G4double e;        // electron charge
  G4double m_deu;    // Deuterom mass (MeV)
  G4double m_pro;    // Deuterom mass (MeV)
  G4double deu_bind; // deuterium binding energy (MeV)
  G4double h_bar;    // reduced Planck Constant (MeV s)
  G4double c;        // speed of light
  G4double m_e;      // mass of the electron
  //maybe we need to deleta most of these definitions


  G4double Mp;



  G4double iGE; // internal calculus GE
  G4double iGM; // internal calculus GM

  G4double Q2; // momentum transfer

  // Nucl_t tNuclType, tFinalNucl;
 
  G4double FluxC; //flux correction from Nucleon Rest to Lab frames
  
  Nucl_t tNuclType, tFinalNucl;
 
  G4double W2, xbj;
  G4ThreeVector tElectronP, tNucleonP;
  G4double tElectronE, tNucleonE;
  G4double  sigma;
  G4double nu;
  G4LorentzVector tElectronf_lab, tNucleonf_lab;

};

extern G4SBSTDISGen *tdishandler;

#endif//G4SBSTDISGen_HH
