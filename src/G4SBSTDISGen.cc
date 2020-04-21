// TDIS Generators class 
// Carlos Ayerbe (gayoso@jlab.org)


#include "TMath.h"
#include <vector>


#include "CLHEP/Random/RandGauss.h"

#include "globals.hh"
#include "CLHEP/Random/RandFlat.h"
#include "G4SBSTDISGen.hh"
#include "G4SBSIO.hh"
#include "G4SystemOfUnits.hh"
#include "G4SBSEventGen.hh"


// Should be in the header? (CA)
extern"C"{
  double epc_func_(double *, int*, int*, int*, double *, double *);
}

// WHY SHOULD BE DECLARED GLOBAL?

using namespace CLHEP;


G4SBSTDISGen::G4SBSTDISGen()
{
  G4cout<<"<G4SBSTDISGen::G4SBSTDISGen()>: initializing"<<G4endl;
  // A bunch of constants (some of them should be replaced by CLHEP)

  PI = TMath::Pi();  //3.14159

  e = 1.602e-19;     // electron charge
  m_deu = 1875.6;    // Deuterom mass (MeV)
  m_pro = 938.3;     // Deuterom mass (MeV)
  deu_bind = 2.224;  // deuterium binding energy (MeV)
  h_bar = 6.582e-22; // reduced Planck Constant (MeV s)
  c = 299792458;     // speed of light
  m_e = 0.511;       // mass of the electron (MeV)
 

  // these are useless, we can get the mass of the particle
  // with ---> ni_Nrest.m()
  Mp =  proton_mass_c2;      // proton mass 
  Mn = neutron_mass_c2; 
  //  iQ2 = 0;

  Q2 = 0;

  // proposal numbers
  tThMin = 5.0*deg;
  tThMax = 45.0*deg;
  tPhMin = -12.0*deg;
  tPhMax = 12.0*deg;

  tEeMin = 0.05 *GeV;
  tEeMax = 5.0 *GeV;

}

G4SBSTDISGen::~G4SBSTDISGen()
{
 G4cout<<"G4SBSTDISGen::~G4SBSTDISGen()"<<G4endl;
}

G4SBSTDISGen *tdishandler=NULL;


void G4SBSTDISGen::Test(Kine_t KinType, Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni)
{

  G4cout<<"Welcome to the generators"<<G4endl;
  G4cout << "GENERATETDIS" << G4endl;


  if(nucl == 0)
    {
      G4cout << "nucl: "<<  nucl << " A PROTON" << G4endl;
      Mt = Mp;
    }
  else
    {
      G4cout << "nucl: "<<  nucl << " A NEUTRON" << G4endl;
      Mt = Mn;
    }
  G4cout << "ni: "<<  ni << G4endl;
  G4cout << "ei: "<<  ei << G4endl;

  G4cout<<"kine: "<<KinType<<G4endl; 

  Generate(KinType, nucl, ei, ni);
  //  qe_epc(); 



  return;
}




void G4SBSTDISGen::Generate(Kine_t tKineType, Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni)
{
  G4cout<<"Entering  <G4SBSTDISGen::Generate>"<<G4endl;
  G4cout<<"kine: "<< tKineType <<G4endl; 
  G4cout<<"nucl: "<< nucl <<G4endl; 
  G4cout<<"ei: "<< ei <<G4endl; 
  G4cout<<"ni: "<< ni <<G4endl; 
    
  
   


  // this method receive the 4-vector beam and 4-vector nucleon
  // WHICH could be at rest (H likewise) or Fermi smeared (CHECK!)

  // So, we have two frames, Lab system and Nucleon rest
  // both frames coincides IF the nucleon is NOT Fermi smeared

 

  // This creates a boost to a rest Nucleon in case it is Fermi smeared
  G4ThreeVector boost_Nrest = ni.boostVector();
  // boostVector() return the spatial coordinates divided by the time component
  //if the nucleon is on rest, it will return 0

  // Total initial 4-momentum in lab
  G4LorentzVector Pisum_lab = ei + ni;

  // 4-momentum vectors in Nucleon rest frame
  // defined FIRST from the lab frame, them boosted to Nucleon rest
  G4LorentzVector ei_Nrest = ei;
  G4LorentzVector ni_Nrest = ni;

  ei_Nrest.boost( -boost_Nrest );
  ni_Nrest.boost( -boost_Nrest );

  G4cout<<"ei_Nrest: "<< ei_Nrest <<G4endl; 
  G4cout<<"ni_Nrest: "<< ni_Nrest <<G4endl; 


  // Note that, if the Nucleon is not smeared, 
  // ei_Nrest = ei; ni_Nrest = ni 


  //electron beam energy in Nucleon rest frame
  G4double Ebeam_Nrest = ei_Nrest.e();

  //electron beam energy in Lab frame
  G4double Ebeam_lab = ei.e();

  G4cout<<"ei_Nrest.e(): "<< ei_Nrest.e() <<G4endl; 
  //  G4cout<<"ei.m2(): "<< ei.m2() <<G4endl; 

  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  //GENERATE KINEMATICS (copy from G4SBSEventGen)
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


  //Generate both theta and phi angles in the LAB frame:
  G4double th = acos( CLHEP::RandFlat::shoot(cos(tThMax), cos(tThMin)) );
  G4double ph = CLHEP::RandFlat::shoot(tPhMin, tPhMax); 
  G4cout<<"th: "<<th/deg<<" ph: "<<ph/deg<< G4endl;


  //unit vector in the direction of the scattered electron in the LAB frame:
  G4ThreeVector kfhat_lab( sin(th)*cos(ph), sin(th)*sin(ph),  cos(th) );
  G4cout<<"kfhat_lab: "<< kfhat_lab <<G4endl; 

  // Outgoing energy of the scattered electron in the LAB frame accounting for 
  // the initial nucleon motion (no off-shell or binding energy corrections, just Fermi momentum)

  if ( tKineType != tElastic && tKineType != tQuasiElastic ) // QE is a special case of Elastic
    {
      Eeprime_lab = CLHEP::RandFlat::shoot(tEeMin, tEeMax );
    }
  else
    {
      G4cout<< " Ela or QE energy"<< G4endl;
      Eeprime_lab = (ei.e()*(ni.e()-ni.pz()))/(ei.e()*(1.-cos(th))+ni.e()-ni.vect().dot(kfhat_lab));
    }

  // DELETE
  // G4ThreeVector kfhat_labMAX( sin(45*deg)*cos(ph), sin(45*deg)*sin(ph),  cos(45*deg) );
  // G4ThreeVector kfhat_labMIN ( sin(  5*deg)*cos(ph), sin(  5*deg)*sin(ph),  cos(  5*deg) );
  // G4double    Eeprime_lab_MAX = (ei.e()*(ni.e()-ni.pz()))/(ei.e()*(1.-cos(45*deg))+ni.e()-ni.vect().dot(kfhat_labMAX));
  // G4double    Eeprime_lab_MIN  = (ei.e()*(ni.e()-ni.pz()))/(ei.e()*(1.-cos(  5*deg))+ni.e()-ni.vect().dot(kfhat_labMIN));

  // G4cout<<"Eeprime_lab_MAX: "<< Eeprime_lab_MAX << " Eeprime_lab_MIN: "<<  Eeprime_lab_MIN<<G4endl;

  Peprime_lab = sqrt(pow(Eeprime_lab, 2) - ei.m2());

  //  G4cout<<"Eeprime_lab: "<< Eeprime_lab <<G4endl; 
  // G4cout<<"Peprime_lab: "<< Peprime_lab <<G4endl; 
  
  // 3-momentum scatter electron in lab frame:
  G4ThreeVector kf_lab = Peprime_lab*kfhat_lab;

  //  G4cout<<"kf_lab: "<< kf_lab <<G4endl; 

  //Four-momentum of scattered electron in the LAB frame:
  G4LorentzVector ef_lab( kf_lab, Eeprime_lab );
  tElectronf_lab = ef_lab; //to rootfile

  //  G4cout<<"ef_lab: "<< ef_lab <<G4endl; 

  //q vector in the LAB frame (virtual photon momentum):
  G4LorentzVector q_lab = ei - ef_lab;
  //  G4cout<<"q_lab: "<<q_lab<<G4endl;


  // Q2 definition
  Q2 = -q_lab.m2(); // set the variable to be used along the whole class

  nu =  q_lab.e(); // return Energy component
  //  G4cout<<"nu: "<<nu<<G4endl;

  // Calculate four-momentum of scattered electron boosted to the 
  // nucleon REST frame for cross section calculation:
 
  G4LorentzVector ef_Nrest = ef_lab; // first define the 4-vector...
  ef_Nrest.boost( -boost_Nrest ); // then boost to the Nucleon rest frame

  //  G4cout<<"ef_Nrest: "<< ef_Nrest <<G4endl; 

  //scattered electron energy in Nucleon Rest frame
  G4double Eprime_Nrest = ef_Nrest.e();  

  //  G4cout<<"ef_Nrest.e(): "<< ef_Nrest.e() <<G4endl; 

  // this is the angle used to calculate Cross-sections 
  // at least Elastic and Mott (for the moment)
  G4double th_Nrest = acos( ei_Nrest.vect().unit().dot( ef_Nrest.vect().unit()) );

  G4LorentzVector q_Nrest = ei_Nrest - ef_Nrest;
  //  G4cout<<"q_Nrest: "<< q_Nrest <<G4endl; 


  G4LorentzVector nf_Nrest = ni_Nrest + q_Nrest;

  //  G4cout<<"nf_Nrest: "<<nf_Nrest<<G4endl;  

  G4double x_Nrest = -q_Nrest.m2()/(2.0*ni_Nrest.dot( q_Nrest ) );
  
  //boost back to LAB frame
  
  G4LorentzVector nf_lab = nf_Nrest;
  nf_lab.boost( boost_Nrest );
  
  //  tNucleonf_lab = nf_lab; // to rootfile

  tNucleonf_lab = nf_Nrest ;// temporary, just to test. to rootfile
  


  // Note: we need to check the units of this correction
  // and take care when is applicable. Like this, it is only useful in the elastic case.
  FluxC = FluxCorrection(nf_lab,  ni, ei, ei_Nrest); // calculate the flux correction for this kinematics

  W2 = (q_lab+ni_Nrest).m2();


  G4cout<<"nf 3-mon: "<< nf_Nrest .vect().mag() <<" nf ene: "<< nf_Nrest .e()<<G4endl;
  //****** more kinematic values*****

  // 

  // double xa =  Q2/(2*Mt*nu);
  // double ya =  nu/EBeam;

  //********************************


  //These values are set fixed now for Deuterium target with the epc code
  // later on, the target used should be carried here

  G4int z1 = 1; //atomic number (number of protons)
  G4int n1 = 1; //THIS IS NUMBER OF NEUTRONS!!! 
  G4int partID = (Nucl_t) nucl; 

  G4cout<<"kine: "<<tKineType<<G4endl; 


  // In principle, each kinematic case just need to carry the 
  // values calculated until here. Then each case will be 
  // treated differently, maybe just one line, maybe a whole 
  // method (as the pion structure generator will need)

  switch(tKineType){
  case tElastic:
    xbj = 1.0;
    ELAsigma = ElasticXS(Ebeam_Nrest, Eprime_Nrest, th_Nrest, nucl);
    G4cout<<"enter in TDIS elastic: "<<ELAsigma/barn<<G4endl;
    break;

  case tQuasiElastic:
    //always in Nucleon Rest Frame
    thN_Nrest = nf_Nrest.theta() ; //nucleon theta angle in N rest frame
    xbj = Q2 / (2.0*nu*ni_Nrest.m() ); //x Bjorken, CHECK!!! it is in Lab frame. 

    //    G4cout<<"standard xbj: "<< Q2/(2.*(ni_Nrest.dot(q_Nrest) ))<<G4endl;

    QEsigma = QuasiElasticXS(Ebeam_Nrest, z1, n1, partID, nf_Nrest .vect().mag(), thN_Nrest);

    G4cout<<"[TDIS] quasielastic sigma (uB/(MeV sr): "<<QEsigma<<G4endl;


    // This is from definition BUT it is NOT clear for me. Nevertheless, it coincides with the epc code
    // value 'tp' which I believe is the kinetic energy  (CA)
    G4cout<<" KE from definition:  "<< sqrt( pow(nf_Nrest.vect().mag(),2) + ni_Nrest.m2()) - ni_Nrest.m()<<G4endl; 

    KE =  sqrt( pow(nf_Nrest.vect().mag(),2) + ni_Nrest.m2()) - ni_Nrest.m();

    break;

  case tInelastic:
    //Calculate the boosted value of Bjorken x:
    xbj = -q_Nrest.m2()/(2.0*ni_Nrest.dot( q_Nrest ) );//not sure
    break;

  default:
    break;
  }

  // momentum-energy of the scattered electron to be sent to PrimaryGeneratorAction
  tElectronP = ef_lab.vect();
  tElectronE = ef_lab.e();

  // momentum-energy of the scattered nucleon to be sent to PrimaryGeneratorAction
  tNucleonP = nf_lab.vect();
  tNucleonE = nf_lab.e();

  tFinalNucl = nucl;
  G4cout<<"nucleon sent: "<<tFinalNucl<<G4endl;


  return;
}



G4double G4SBSTDISGen::FluxCorrection( G4LorentzVector nf_lab, G4LorentzVector ni, G4LorentzVector ei, G4LorentzVector ei_Nrest)
{

  G4double beta            = nf_lab.vect().mag()/nf_lab.e(); //beta = p/e
  G4double costheta_eN_lab = (ei.vect().unit() ).dot( ni.vect().unit() );
  G4double betaN_lab       = ni.beta();
  G4double gammaN_lab      = ni.gamma();
  G4double flux_Nrest      = 4.0*ni.m()*ei_Nrest.e();
  G4double flux_lab        = 4.0*ei.e()*ni.e()*sqrt( 2.0*(1.0-betaN_lab*costheta_eN_lab) - pow(gammaN_lab,-2) );

//The lines above already converted the cross section to GEANT4 units. Now this has dimensions of area, is expressed in the lab frame, and is differential in solid angle only! (note from EventGen)

  return flux_Nrest/flux_lab;
}




G4double G4SBSTDISGen::ElasticXS(G4double beam_energy, G4double scatter_e_energy, G4double theta, Nucl_t nu)
{
 

  // I redefined the variables, just for facility notation
  // I prefer to keep large names at the beggining to clear know what are we talking about
  G4double E = beam_energy;
  G4double E_prime = scatter_e_energy;


  //Mott x Recoil fraction x form factors relation

  G4double eSigma = MottXS(theta, E)* (E_prime/E)*
    ( (pow(GE(nu),2)+tau()*pow(GM(nu),2)/(1.0+tau()) + 2.0*tau()*pow(GM(nu),2)*pow(tan(theta/2.0),2) )); 
  // Dimensions of area 
  
  return eSigma*FluxC;
}



G4double G4SBSTDISGen::PhotoD_XS(G4double E_photon)
{
  // E_photon; // virtual photon energy in MeV

  // Bethe-Peirls Deuterium photodesintegration cross-section
  // http://www1.jinr.ru/Pepan_letters/panl_2013_3/16_did.pdf

  if(E_photon < deu_bind)
    {
      return 0; //to extend the X axis to 0
    }
  else
    {

      // a is the constant factor of the Bethe-Peirls cross-section 
      // READ NOTE BELOW!!
      
      //  a = (8*PI)/(3*M);//h_bar,c, e = 1 
      //  a = (8*PI*pow(e,2)*h_bar)/(3*M*c);

      // b is the factor energy-dependent of the Bethe-Peirls cross-section 
      G4double b = (sqrt(deu_bind) * sqrt(pow(E_photon - deu_bind,3)))/pow(E_photon,3);
	
      return (2.4/0.056195)*b; // READ NOTE BELOW!!
    }

  // NOTE:
  // The paper indicates the maximum is 2.4mb @4.4MeV (which is true from  
  // real data base)
  // IN PRINCIPLE XS is in fm² (actually in MeV⁻² doing the
  // dimensional analysis with e=h=c=1) 
  // BUT considering the values of the constant units in many systems
  // there is no a consensus about their values.
  // Following a suggestion from Rey Cruz, I calculated the constant value
  // knowing the maximum. 2.4mb@4.4MeV. In other words /sigma = K * f(E)
  // where f(E) is the factor energy dependent.
  // thus K = 2.4 mb / 0.056195 MeV^-1
}


G4double G4SBSTDISGen::VXPhoton_flux(double E_photon, double E_beam)
{
  // virtual photon flux, check with Raffo's CLAS12 J/psi photoproduction proposal

  G4double Q2max = 0.3; // (MeV/c)^2 Cut-off value. The formula is not so sensitive 
  // to this value. It was tested with 0.3, 3, 30 and 300 (MeV/c)²
  // and the change in flux is neglegible
  
  G4double x = E_photon / E_beam;

  if (x >= 1) return 0;
  
  G4double Q2min = m_e *x*x/(1-x);

  G4double term1 = (1 - x + x*x/2);
  G4double term2 = log(Q2max/Q2min);
  G4double term3 = alpha() /(E_beam*x*TMath::Pi());
  
  return term3 * ( term1 * term2 - (1-x)) ;
}



G4double G4SBSTDISGen::QuasiElasticXS(G4double beam_energy, G4int z1, G4int n1, G4int partID, G4double momentum, G4double angle)
{
  // This function calls a fortran function based in the 
  // EPC code by Lightbody and O'Connell, 
  // J.S. O'Connell and J.W.Lightbody, Computers in Physics 2,57(1988).
  // the function used here was modified by O. Rondon 
  // and adapted as a function by C. Ayerbe

  // the arguments needed (in order):
  // beam energy (in MeV)
  // atomic number
  // number of nucleons
  // kind of particle // 1=proton, -1=neutron, 2=pi+, -2=pi-, 0=pi0
  // momentum of the emmited particle
  // angle of emission (in deg, it is converted in rad inside the function)

  angle = angle*(180./TMath::Pi());

  G4cout<<"QuasiElasticXS: beam_energy: "<<beam_energy<< " momentum: "<<momentum<<" angle: "<<angle<<G4endl;

  // EPC notation (just neutron/proton now)
  if (partID == 0)
    {
      partID = 1; //proton
    }
  else
    {
      partID = -1; // neutron
    }

  return  epc_func_(&beam_energy, &z1, &n1, &partID, &momentum, &angle);

}


G4double  G4SBSTDISGen::MottXS(G4double theta, G4double beam_energy)
{
  //the theta angle used is the theta in the NUCLEON REST frame!!

  G4double M1 = cos(theta/2)*alpha();

  G4double M2 = 2*beam_energy*pow(sin(theta/2),2);

  //MXS: Mott Cross Section. Keep in mind, you don't care units.
  G4double MXS = pow(hbarc, 2) * pow(M1/M2, 2);

  return MXS; 

  // According to Andrew, if we don't touch the units, Geant4 works internally with them
  // in other words, BECAUSE we give some units to certain numbers, GEANT4 knows
  // how to convert to the units WE WANT. 
  // For example, the final result here will be in mm^2, if we want barns
  // just MXS/barn and the final number is in such unit.
}


G4double G4SBSTDISGen::DipoleFF()
{
  G4double Q2t = Q2;

  return pow(1.0 + Q2t/(0.71*GeV*GeV), -2.0);
}


G4double G4SBSTDISGen::GE( Nucl_t nucleon)
{
   //I need to find where this comes from

  if (nucleon == 0) //a proton
    {
      iGE = (1.0-0.24*tau())/
	(1.0 + 10.98*tau() + 12.82*pow(tau(),2) + 21.97*pow(tau(),3));
    }
  else
    {
      iGE = (1.520*tau() + 2.629*pow(tau(),2) + 3.055*pow(tau(),2)*DipoleFF())/
	(1.0+5.222*tau()+0.040*pow(tau(),2)+11.438*pow(tau(),3));
    }
  
  return iGE;

}

G4double G4SBSTDISGen::GM( Nucl_t nucleon)
{
  if (nucleon == 0) //a proton
    {
      iGM = 2.79*(1.0+0.12*tau())/
	(1.0 + 10.97*tau() + 18.86*pow(tau(),2) + 6.55*pow(tau(),3) );
    }
  else
    {
      iGM = -1.913*(1.0+2.33*tau())/ 
	(1.0 + 14.72*tau() + 24.20*pow(tau(),2) + 84.1*pow(tau(),3));
    }
  
  return iGM;
  
}

G4double G4SBSTDISGen::tau()
{
  G4double Q2t = Q2;

  return Q2t/(4.0*Mp*Mp);
}

