#include "G4SBSRICHoutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSRICHoutput::G4SBSRICHoutput(){
  //Set sensible default values for threshold and time window:

  timewindow = 100.0*ns; 
  threshold  = 0.5; //single photo-electron threshold!
  

  Clear();
}

G4SBSRICHoutput::~G4SBSRICHoutput(){
  ;
}

void G4SBSRICHoutput::Clear(){
  nhits_RICH = 0;
  ntracks_RICH = 0;
  
  PMTnumber.clear();
  row.clear();
  col.clear();
  //NumPhotons.clear();
  NumPhotoelectrons.clear();
  Time_avg.clear();
  Time_rms.clear();
  mTrackNo.clear();

  xhit.clear();
  yhit.clear();
  zhit.clear();
  
  pxhit.clear();
  pyhit.clear();
  pzhit.clear();

  pvx.clear();
  pvy.clear();
  pvz.clear();

  ppx.clear();
  ppy.clear();
  ppz.clear();
  
  //  thetaC.clear();
  volume_flag.clear();
  
  mPID.clear();
  mvx.clear();
  mvy.clear();
  mvz.clear();
  mpx.clear();
  mpy.clear();
  mpz.clear();

}
