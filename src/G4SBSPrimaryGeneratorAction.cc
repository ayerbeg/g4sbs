#include "TBuffer.h"
#include "TMatrixTBase.h"
#include "TString.h"
#include "THashTable.h"

#include "G4SBSPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SBSEventGen.hh"

#include "G4SBSTDISGen.hh"

#include "G4SBSIO.hh"
#include "G4SBSRunAction.hh"
#include "sbstypes.hh"
#include "globals.hh"
#include "TVector3.h"
#include "Randomize.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

G4SBSPrimaryGeneratorAction::G4SBSPrimaryGeneratorAction()
{

 
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  GunParticleName = "e-";

  particle = particleTable->FindParticle(particleName="e-");

  particleGun->SetParticleDefinition(particle);

  particleGun->SetParticleMomentumDirection(G4ParticleMomentum(sin(-40.0*deg),0.0,cos(-40.0*deg)));
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  particleGun->SetParticlePolarization( G4ThreeVector(0,0,0) );
  
  GunParticleType = particle;
  GunPolarization = G4ThreeVector(0,0,0);
 
 
  G4cout << "*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*" << G4endl;
  G4cout<<"This is PrimaryGeneratorAction constructor"<< G4endl;
  sbsgen = new G4SBSEventGen();
  G4cout << "Initialize G4SBSEventGen" << G4endl;
  G4cout << "*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*" << G4endl;
   

  fUseGeantino = false;
}

G4SBSPrimaryGeneratorAction::~G4SBSPrimaryGeneratorAction()
{
  delete particleGun;
}

void G4SBSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4cout<<"Now we are in PrimaryGeneratorAction"<< G4endl;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  //ev_t evdata;

  //G4SBSPythiaOutput Primaries;
  
  // Several different types of scattering
  // Let's start with e'N elastic


  
  //  Roll up random values
  int ntries = 1;
  while( !sbsgen->GenerateEvent() ) //THIS IS THE CALL TO GENERATE EVENT!!
    { 
      ntries++; 
    }


//  G4cout << "Got event, ntries = " << ntries << G4endl;

  int ntries_run = RunAction->GetNtries();
  RunAction->SetNtries( ntries_run + ntries );

 //evdata = sbsgen->GetEventData();

  fIO->SetEventData(sbsgen->GetEventData());// THIS GUY SAVES THE DATA


  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  //PYTHIA6 event
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  if( sbsgen->GetKine() == kPYTHIA6 )
    { 
      G4SBSPythiaOutput Primaries = sbsgen->GetPythiaEvent();
      //cout << "G4SBSPrimaryGeneratorAction.cc" << endl;
      for( int ipart = 0; ipart<Primaries.Nprimaries; ipart++ )
	{
	  //cout << ipart << " " << Primaries.genflag[ipart] << " " << Primaries.PID[ipart] << " " 
	  //   << Primaries.Px[ipart] << " " << Primaries.Py[ipart] << " " << Primaries.Pz[ipart] << endl;
	  if( Primaries.genflag[ipart] != 0 )
	    {
	      particle = particleTable->FindParticle( Primaries.PID[ipart] );
	      particleGun->SetParticleDefinition(particle);
	      particleGun->SetNumberOfParticles(1);
	      particleGun->SetParticleEnergy( Primaries.E[ipart]-Primaries.M[ipart] );
	      particleGun->SetParticleMomentumDirection( G4ThreeVector( Primaries.Px[ipart], Primaries.Py[ipart], Primaries.Pz[ipart] ).unit() );
	      
	      //G4ThreeVector vertex( Primaries.vx[ipart], Primaries.vy[ipart], Primaries.vz[ipart] );
	      //vertex += sbsgen->GetV();
	      particleGun->SetParticlePosition( G4ThreeVector( Primaries.vx[ipart], Primaries.vy[ipart], Primaries.vz[ipart] ) );
	      particleGun->SetParticleTime( Primaries.t[ipart] );
	      particleGun->GeneratePrimaryVertex(anEvent);
	    } 
	}
      
      Primaries.ConvertToTreeUnits();
      fIO->SetPythiaOutput( Primaries );
      
      return;
    }
  
  
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // TDIS AcquMC
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  if( sbsgen->GetKine() == kAcquMC )
    { //AcquMC event:
      G4SBSAcquMCOutput Primaries = sbsgen->GetAcquMCEvent();
      
      particle = particleTable->FindParticle( Primaries.pid );
      particleGun->SetParticleDefinition(particle);
      particleGun->SetNumberOfParticles(1);
      particleGun->SetParticleEnergy( Primaries.E - particle->GetPDGMass() );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( Primaries.Px, Primaries.Py, Primaries.Pz ).unit() );
      particleGun->SetParticlePosition( G4ThreeVector( Primaries.Vx, Primaries.Vy, Primaries.Vz ) );
      particleGun->GeneratePrimaryVertex(anEvent);

      
      // Primaries.ConvertToTreeUnits();
      fIO->SetAcquMCOutput( Primaries );
      
      return;
    }


  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // electron or geantino
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  if( !fUseGeantino && sbsgen->GetKine() != kGun )
    { //first primary is an electron!
      particle = particleTable->FindParticle(particleName="e-");
    }
  else if( fUseGeantino )
    { //first primary is a geantino!
      particle = particleTable->FindParticle(particleName="chargedgeantino");
    } 
  else 
    { //first primary according to particle gun generator:
      particle = particleTable->FindParticle( GunParticleName );
      if( particle != 0 ) SetParticleType( particle );
      particle = GunParticleType;
    }
  
 //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // cosmics 
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  if( sbsgen->GetKine()==kCosmics )
    {
      particle = particleTable->FindParticle(particleName="mu-");
      if( fUseGeantino )
	{ //first primary is a geantino!
	  particle = particleTable->FindParticle(particleName="chargedgeantino");
	}
    }
  
  particleGun->SetParticleDefinition(particle);
  
  // I guess making a flag could make it easier? (CA)
  if ( sbsgen->GetKine() != tElastic && //new (CA)
       sbsgen->GetKine() != tInelastic && //new (CA)
       sbsgen->GetKine() != tQuasiElastic && //new (CA)
       sbsgen->GetKine() != tTDISKin //new (CA)
       )
    {
      particleGun->SetParticleMomentumDirection( sbsgen->GetElectronP().unit() );
    }
  else
    {
      particleGun->SetParticleMomentumDirection( tdishandler -> tGetElectronP().unit() );
    }
  
    

  if( sbsgen->GetKine() != kGun && sbsgen->GetKine() != kCosmics)
    { 
      
      if ( sbsgen->GetKine() != tElastic && //new (CA)
	   sbsgen->GetKine() != tInelastic && //new (CA)
	   sbsgen->GetKine() != tQuasiElastic && //new (CA)
	   sbsgen->GetKine() != tTDISKin //new (CA)
	   )
	{
	  particleGun->SetParticleEnergy(sbsgen->GetElectronE());
	}
      else
	{
	  particleGun->SetParticleEnergy(tdishandler->tGetElectronE());//(CA)
	}
      
    } 
  else 
    { //kGun!
      //SetParticleEnergy sets the ***kinetic energy*** of particles
      //GenerateGun() generates the ***momentum***; therefore we need:
      // T = E - M --> (T + M)^2 = p^2 + M^2 --> T^2 + 2MT = p^2 
      particleGun->SetParticleEnergy( sqrt( pow( sbsgen->GetElectronP().mag(), 2) + pow( GunParticleType->GetPDGMass(), 2 ) ) - GunParticleType->GetPDGMass() );
    }
  
  particleGun->SetParticlePosition(sbsgen->GetV());
  
  /*
    particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(-40.0*deg),0.0,cos(-40.0*deg)));
    particleGun->SetParticleEnergy(1.0*GeV);
    particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  */
  
  // Not necessarily kinematically allowed
  
  ev_t evdata = fIO->GetEventData();
  
  
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // Wiser
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  if( sbsgen->GetKine()!= kWiser )// is this correct?
    {
      particleGun->SetParticlePolarization( G4ThreeVector(0.0,0.0,0.0) );
      //G4cout << "Gun polarization for the primary electron: " << particleGun->GetParticlePolarization() << G4endl;
      if( sbsgen->GetKine() == kGun )
	{ //If a gun polarization is defined, transform polarization to TRANSPORT coordinates and set particle polarization:
	  gen_t gendata = fIO->GetGenData();
	  G4double sbsangle = gendata.thsbs;
	  
	  G4ThreeVector sbsaxis( -sin(sbsangle), 0.0, cos(sbsangle) );
	  G4ThreeVector yaxis(0,1,0);
	  G4ThreeVector xaxis = (yaxis.cross(sbsaxis)).unit();
	  
	  G4ThreeVector Pol_transport = GunPolarization.y() * xaxis - GunPolarization.x() * yaxis + GunPolarization.z() * sbsaxis;  
	  
	  evdata.Sx = GunPolarization.x();
	  evdata.Sy = GunPolarization.y();
	  evdata.Sz = GunPolarization.z();
	  fIO->SetEventData( evdata );
	  
	  particleGun->SetParticlePolarization( Pol_transport );
	  
	}
      particleGun->GeneratePrimaryVertex(anEvent);
    }
  
  
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // GMn Elastic Check
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  if( sbsgen->GetKine() == kGMnElasticCheck ) 
    { // In this case we
      // generate TWO nucleons, a proton and neutron simultaneously
      for(int i = 0; i < 2; i++ ) 
	{
	  if(i==0)
	    particle = particleTable->FindParticle(particleName="neutron");
	  else
	    particle = particleTable->FindParticle(particleName="proton");
	  
	  particleGun->SetParticleDefinition(particle);
	  
	  particleGun->SetParticleMomentumDirection(sbsgen->GetNucleonP().unit());
	  // This is KINETIC energy
	  particleGun->SetParticleEnergy(sbsgen->GetNucleonE()-particle->GetPDGMass());
	  particleGun->SetParticlePosition(sbsgen->GetV());
	  
	  // Generate the particle
	  particleGun->GeneratePrimaryVertex(anEvent);
	}
      return; // So that we don't generate any more particles by accident
    }
  
  
  
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // Any other kinematics
  // This is where TDISkin is generated
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  
 
  if( sbsgen->GetKine() != kSIDIS && 
      sbsgen->GetKine() != kWiser && 
      sbsgen->GetKine() != kGun   && 
      sbsgen->GetKine() != kBeam  && 
      sbsgen->GetKine() != kCosmics &&
      sbsgen->GetKine() != tElastic && //new (CA)
      sbsgen->GetKine() != tInelastic && //new (CA)
      sbsgen->GetKine() != tQuasiElastic && //new (CA)
      sbsgen->GetKine() != tTDISKin //new (CA)
      )
    { //Then we are generating a final nucleon
      
      G4cout<<"Creating kinematics"<< G4endl;
      switch( sbsgen->GetFinalNucleon() ){
      case kProton:
	G4cout<<"a Proton"<< G4endl;
	particle = particleTable->FindParticle(particleName="proton");
	break;
      case kNeutron:
	G4cout<<"a Neutron"<< G4endl;
	particle = particleTable->FindParticle(particleName="neutron");
	break;
      default:
	particle = particleTable->FindParticle(particleName="geantino");
	break;
      } 
      particleGun->SetParticleDefinition(particle);
      
      // Ensure we're doing something sensible for Geant4
      if( sbsgen->GetNucleonE()-particle->GetPDGMass() > 0.0 ) 
	{
	  particleGun->SetParticleMomentumDirection(sbsgen->GetNucleonP().unit());
	  // This is KINETIC energy
	  particleGun->SetParticleEnergy(sbsgen->GetNucleonE()-particle->GetPDGMass());
	  particleGun->SetParticlePosition(sbsgen->GetV());
	  
	  if( sbsgen->GetKine() == kElastic ) 
	    {
	      G4double muB = 0.5*eplus*hbar_Planck/(particle->GetPDGMass()/c_squared);
	      
	      // G4cout << "Particle magnetic moment = "
	      //        << particle->GetPDGMagneticMoment()/muB << G4endl;
	      
	      G4ThreeVector k_hat(0,0,1.0); // beam polarization unit vector
	      G4ThreeVector n_hat = ((sbsgen->GetNucleonP().unit()).cross(k_hat)).unit();
	      G4ThreeVector t_hat = n_hat.cross( (sbsgen->GetNucleonP().unit()) );
	      G4ThreeVector S_hat = (sbsgen->GetPl())*(sbsgen->GetNucleonP().unit()) + (sbsgen->GetPt())*t_hat;
	      
	      //	G4cout << "Initial polarization = " << S_hat << G4endl;
	      
	      particleGun->SetParticlePolarization( S_hat.unit() );
	      
	      gen_t gendata = fIO->GetGenData();
	      G4double sbsangle = gendata.thsbs;
	      
	      G4ThreeVector sbsaxis( -sin(sbsangle), 0, cos(sbsangle) );
	      G4ThreeVector yaxis(0,1,0);
	      G4ThreeVector xaxis = (yaxis.cross(sbsaxis)).unit();
	      
	      evdata.Sx = -S_hat.dot( yaxis );
	      evdata.Sy = S_hat.dot( xaxis );
	      evdata.Sz = S_hat.dot( sbsaxis );
	      
	      fIO->SetEventData( evdata );
	  
	    }
	}
    } 
  
  else 
    
    if( sbsgen->GetKine() == kSIDIS || sbsgen->GetKine() == kWiser )
      { //SIDIS case: generate a final hadron:
	switch( sbsgen->GetHadronType() ){
	case kPiPlus:
	  particle = particleTable->FindParticle(particleName="pi+");
	  break;
	case kPiMinus:
	  particle = particleTable->FindParticle(particleName="pi-");
	  break;
	case kPi0:
	  particle = particleTable->FindParticle(particleName="pi0");
	  break;
	case kKPlus:
	  particle = particleTable->FindParticle(particleName="kaon+");
	  break;
	case kKMinus:
	  particle = particleTable->FindParticle(particleName="kaon-");
	  break;
	case kP:
	  particle = particleTable->FindParticle(particleName="proton");
	break;
	case kPbar:
	  particle = particleTable->FindParticle(particleName="anti_proton");
	  break;
	default:
	  particle = particleTable->FindParticle(particleName="pi+");
	  break;
	}
	
	particleGun->SetParticleDefinition( particle );
	if( sbsgen->GetHadronE()-particle->GetPDGMass() > 0.0 ) 
	  {
	    particleGun->SetParticleMomentumDirection( sbsgen->GetHadronP().unit() );
	    particleGun->SetParticleEnergy( sbsgen->GetHadronE()-particle->GetPDGMass() );
	    particleGun->SetParticlePosition( sbsgen->GetV() );
	  }
      }
  
  
    else // new TDIS generators

      if( sbsgen->GetKine() == tElastic) 
      // sbsgen->GetKine() != kWiser && 
      // sbsgen->GetKine() != kGun   && 
      // sbsgen->GetKine() != kBeam  && 
      // sbsgen->GetKine() != kCosmics )
    { //Then we are generating a final nucleon

      G4cout<<"Creating TDIS kinematics"<< G4endl;
      G4cout<<"nucleon received: "<<tdishandler->tGetFinalNucleon()<<G4endl;
      switch( tdishandler->tGetFinalNucleon() )
	{
	case kProton:
	  G4cout<<"a tProton"<< G4endl;
	  particle = particleTable->FindParticle(particleName="proton");
	  break;
	case kNeutron:
	  G4cout<<"a tNeutron"<< G4endl;
	  particle = particleTable->FindParticle(particleName="neutron");
	  break;
	default:
	  particle = particleTable->FindParticle(particleName="geantino");
	  break;
	} 
      particleGun->SetParticleDefinition(particle);
      
      // Ensure we're doing something sensible for Geant4
      if( tdishandler -> tGetNucleonE()-particle->GetPDGMass() > 0.0 ) 
	{
	  particleGun->SetParticleMomentumDirection(tdishandler -> tGetNucleonP().unit());
	  // This is KINETIC energy
	  particleGun->SetParticleEnergy(tdishandler -> tGetNucleonE()-particle->GetPDGMass());
	  particleGun->SetParticlePosition(sbsgen->GetV()); //this is common to all kinematics
	}
      fIO->SetEventData( evdata ); 
    }









  // Only do final nucleon/hadron for generators other than
  // the generic beam generator
  if( sbsgen->GetKine() != kBeam && sbsgen->GetKine() != kGun && sbsgen->GetKine() != kDIS && sbsgen->GetKine() != kCosmics ){
    //G4cout << "Gun polarization = " << particleGun->GetParticlePolarization() << G4endl;
    particleGun->GeneratePrimaryVertex(anEvent);
  }
  
}

G4ParticleGun* G4SBSPrimaryGeneratorAction::GetParticleGun()
{
  return particleGun;
} 

