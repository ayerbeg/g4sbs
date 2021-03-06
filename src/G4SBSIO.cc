#include "TObjString.h"
#include "TString.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>

#include "G4SBSGlobalField.hh"
#include "G4SBSRun.hh"
#include "G4SBSIO.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSEArmBuilder.hh"
#include "G4SBSHArmBuilder.hh"
//#include "G4SDManager.hh"
#include <assert.h>
#include "sbstypes.hh"

G4SBSIO::G4SBSIO(){
    fTree = NULL;
    //InitializeTree(); //We want experiment-dependent ROOT tree! Don't invoke until after fdetcon->ConstructAll() has been invoked!
    // Default filename
    strcpy(fFilename, "g4sbsout.root");
    fFile = NULL;

    //Sometimes these default parameters don't match what was actually used in simulation if certain commands
    //aren't invoked. Moreover, these default values don't match the default values in
    //G4SBSDetectorConstruction, G4SBSEarmBuilder, G4SBSHArmBuilder, etc.
    //Therefore, if commands aren't invoked, these defaults can be wrong/misleading.
    //In principle the best way to make sure there is one consistent set of default values is to grab this
    //information from the G4SBSDetectorConstruction whenever it changes!
    //Easiest (but not necessarily best?) way is to store a pointer to G4SBSIO as a data member of G4SBSDetectorConstruction so they can
    //directly talk to each other?
    gendata.Ebeam = 2.2;
    gendata.thbb = 40.0*CLHEP::deg;
    gendata.dbb = 1.5;
    gendata.thsbs = 39.4*CLHEP::deg;
    gendata.dhcal = 17.0;
    gendata.voffhcal = 0.0;
    gendata.dsbs = 1.6;
    gendata.drich = 4.6;
    gendata.dsbstrkr = 4.3;

    KeepPartCALflags.clear();
    KeepHistoryflags.clear();

    Esum_histograms = NULL;
    PulseShape_histograms = NULL;
    
    fUsePythia = false;

    fNhistograms = 0;
}

G4SBSIO::~G4SBSIO(){
    if( fTree ){delete fTree;}
    fTree = NULL;

    if( fFile ){delete fFile;}
    fFile = NULL;

    if( Esum_histograms ){ delete Esum_histograms; }
    if( PulseShape_histograms ){ delete PulseShape_histograms; }
}

void G4SBSIO::SetGEMData( G4String SDname, G4SBSGEMoutput gd ){
  GEMdata[SDname] = gd;
}

void G4SBSIO::SetTrackData( G4String SDname, G4SBSTrackerOutput td ){
  trackdata[SDname] = td;
}

void G4SBSIO::SetCalData( G4String SDname, G4SBSCALoutput cd ){
  CALdata[SDname] = cd;
}

void G4SBSIO::SetRICHData( G4String SDname, G4SBSRICHoutput rd ){
  richdata[SDname] = rd;
}

void G4SBSIO::SetECalData( G4String SDname, G4SBSECaloutput ed ){
  ecaldata[SDname] = ed;
}

void G4SBSIO::InitializeTree(){
    if( fFile ){
	fFile->Close();
	delete fFile;
    }

    fFile = new TFile(fFilename, "RECREATE"); 

    if( fTree ){ delete fTree; }

    Esum_histograms = new TClonesArray("TH1F",10);
    PulseShape_histograms = new TClonesArray("TH1F",10);

    fNhistograms = 0;
    
    fTree = new TTree("T", "Geant4 SBS Simulation");
    fTree->Branch("ev", &evdata, "count/D:rate/D:solang/D:sigma/D:W2/D:xbj/D:Q2/D:th/D:ph/D:Aperp/D:Apar/D:Pt/D:Pl/D:vx/D:vy/D:vz/D:ep/D:np/D:epx/D:epy/D:epz/D:npx/D:npy/D:npz/D:nth/D:nph/D:pmperp/D:pmpar/D:pmparsm/D:z/D:phperp/D:phih/D:MX2/D:Sx/D:Sy/D:Sz/D:nucl/I:fnucl/I:hadr/I:earmaccept/I:harmaccept/I");
    //fTree->Branch("tr", &trdata, "x/D:y/D:xp/D:yp/D:tx/D:ty/D:txp/D:typ/D:hcal/I:bb/I:gemtr/I:hcx/D:hcy/D:bcx/D:bcy/D:hct/D:hctex/D:hclx/D:hcly/D:hclz/D:hcdang/D");
    //fTree->Branch("gen", &gendata, "thbb/D:thsbs/D:dbb/D:dsbs/D:dhcal/D:voffhcal/D:drich/D:dsbstrkr/D:Ebeam/D");

    //Instead of having the same tree structure as before, we want to dynamically generate tree branches depending on what kinds of detectors are present: Since we already require the ROOT libraries, we might as well use TStrings:
    
    //For all tree branches representing data in sensitive detectors, we want to grab the information from fdetcon->SDlist
    //Later, we will add other kinds of sensitive detectors:
    for( set<G4String>::iterator d = (fdetcon->SDlist).begin(); d != (fdetcon->SDlist).end(); d++ ){
      //for( G4int idet=0; idet<fdetcon->fSDman->G
      G4String SDname = *d;
      SDet_t SDtype = (fdetcon->SDtype)[SDname];

      G4cout << "Initializing tree branches for Sensitive Detector " << SDname.data() << G4endl;

      switch( SDtype ){
      case kGEM: //GEM: Add branches for the GEM AND "tracker" branches:
	//Create "GEM output" and "Tracker Output" data structures and associate them with this sensitive detector name:
	GEMdata[SDname] = G4SBSGEMoutput();
	trackdata[SDname] = G4SBSTrackerOutput();
	
	BranchGEM(SDname);
	
	break;
      case kCAL: //"CAL": Add appropriate branches:
	//Initialize "CAL output" data structure and associate with this sensitive detector:
	CALdata[SDname] = G4SBSCALoutput();
	BranchCAL(SDname);
	break;
      case kRICH: //"RICH"
	richdata[SDname] = G4SBSRICHoutput();
	BranchRICH(SDname);
	break;
      case kECAL: //"ECAL"
	ecaldata[SDname] = G4SBSECaloutput();
	BranchECAL(SDname);
	break;
      }
    }
    if( fUsePythia ){
      BranchPythia();
    }

    // // Tedious, but we want dynamically scaled
    // fTree->Branch("ht.ndata", &hitdata.ndata, "ht.ndata/I");
    // fTree->Branch("ht.gid", &hitdata.gid, "ht.gid[ht.ndata]/I");
    // fTree->Branch("ht.trkrid", &hitdata.trkrid, "ht.trkrid[ht.ndata]/I");
    // fTree->Branch("ht.x", &hitdata.x, "ht.x[ht.ndata]/D");
    // fTree->Branch("ht.y", &hitdata.y, "ht.y[ht.ndata]/D");
    // fTree->Branch("ht.z", &hitdata.z, "ht.z[ht.ndata]/D");
    // fTree->Branch("ht.t", &hitdata.t, "ht.t[ht.ndata]/D");
    // fTree->Branch("ht.vx", &hitdata.vx, "ht.vx[ht.ndata]/D");
    // fTree->Branch("ht.vy", &hitdata.vy, "ht.vy[ht.ndata]/D");
    // fTree->Branch("ht.vz", &hitdata.vz, "ht.vz[ht.ndata]/D");
    // fTree->Branch("ht.dx", &hitdata.dx, "ht.dx[ht.ndata]/D");
    // fTree->Branch("ht.dy", &hitdata.dy, "ht.dy[ht.ndata]/D");
    // fTree->Branch("ht.p", &hitdata.p, "ht.p[ht.ndata]/D");

    // fTree->Branch("ht.trid", &hitdata.trid, "ht.trid[ht.ndata]/I");
    // fTree->Branch("ht.pid", &hitdata.pid, "ht.pid[ht.ndata]/I");
    // fTree->Branch("ht.mid", &hitdata.mid, "ht.mid[ht.ndata]/I");
    // fTree->Branch("ht.edep", &hitdata.edep, "ht.edep[ht.ndata]/D");

    // fTree->Branch("ht.tx", &hitdata.tx, "ht.tx[ht.ndata]/D");
    // fTree->Branch("ht.ty", &hitdata.ty, "ht.ty[ht.ndata]/D");
    // fTree->Branch("ht.txp", &hitdata.txp, "ht.txp[ht.ndata]/D");
    // fTree->Branch("ht.typ", &hitdata.typ, "ht.typ[ht.ndata]/D");

    // fTree->Branch("hc.ndata", &caldata.hcndata, "hc.ndata/I");
    // fTree->Branch("hc.x", &caldata.hcx, "hc.x[hc.ndata]/D");
    // fTree->Branch("hc.y", &caldata.hcy, "hc.y[hc.ndata]/D");
    // fTree->Branch("hc.z", &caldata.hcz, "hc.z[hc.ndata]/D");
    // fTree->Branch("hc.e", &caldata.hce, "hc.e[hc.ndata]/D");
    // fTree->Branch("hc.t", &caldata.hct, "hc.t[hc.ndata]/D");
    // fTree->Branch("hc.vx", &caldata.hcvx, "hc.vx[hc.ndata]/D");
    // fTree->Branch("hc.vy", &caldata.hcvy, "hc.vy[hc.ndata]/D");
    // fTree->Branch("hc.vz", &caldata.hcvz, "hc.vz[hc.ndata]/D");
    // fTree->Branch("hc.trid", &caldata.hctrid, "hc.trid[hc.ndata]/I");
    // fTree->Branch("hc.mid", &caldata.hcmid, "hc.mid[hc.ndata]/I");
    // fTree->Branch("hc.pid", &caldata.hcpid, "hc.pid[hc.ndata]/I");

    // fTree->Branch("bc.ndata", &caldata.bcndata, "bc.ndata/I");
    // fTree->Branch("bc.x", &caldata.bcx, "bc.x[bc.ndata]/D");
    // fTree->Branch("bc.y", &caldata.bcy, "bc.y[bc.ndata]/D");
    // fTree->Branch("bc.z", &caldata.bcz, "bc.z[bc.ndata]/D");
    // fTree->Branch("bc.e", &caldata.bce, "bc.e[bc.ndata]/D");
    // fTree->Branch("bc.t", &caldata.bct, "bc.t[bc.ndata]/D");
    // fTree->Branch("bc.vx", &caldata.bcvx, "bc.vx[bc.ndata]/D");
    // fTree->Branch("bc.vy", &caldata.bcvy, "bc.vy[bc.ndata]/D");
    // fTree->Branch("bc.vz", &caldata.bcvz, "bc.vz[bc.ndata]/D");
    // fTree->Branch("bc.trid", &caldata.bctrid, "bc.trid[bc.ndata]/I");
    // fTree->Branch("bc.mid", &caldata.bcmid, "bc.mid[bc.ndata]/I");
    // fTree->Branch("bc.pid", &caldata.bcpid, "bc.pid[bc.ndata]/I");

    // //Declare new "vectorized" Tracking output branches:
    // fTree->Branch("ntracks", &(trackdata.ntracks), "ntracks/I");
    // fTree->Branch("trackerid", &(trackdata.TrackerID) );
    // fTree->Branch("trackid", &(trackdata.TrackTID) );
    // fTree->Branch("trackpid", &(trackdata.TrackPID) );
    // fTree->Branch("tracknhits", &(trackdata.NumHits) );
    // fTree->Branch("tracknplanes", &(trackdata.NumPlanes) );
    // fTree->Branch("trackndf", &(trackdata.NDF) );
    // fTree->Branch("trackchi2", &(trackdata.Chi2fit) );
    // fTree->Branch("trackchi2true", &(trackdata.Chi2true) );
    // fTree->Branch("trackx", &(trackdata.TrackX) );
    // fTree->Branch("tracky", &(trackdata.TrackY) );
    // fTree->Branch("trackxp", &(trackdata.TrackXp) );
    // fTree->Branch("trackyp", &(trackdata.TrackYp) );
    // fTree->Branch("trackt", &(trackdata.TrackT) );
    // fTree->Branch("trackp", &(trackdata.TrackP) );
    // fTree->Branch("trackxfit", &(trackdata.TrackXfit) );
    // fTree->Branch("trackyfit", &(trackdata.TrackYfit) );
    // fTree->Branch("trackxpfit", &(trackdata.TrackXpfit) );
    // fTree->Branch("trackypfit", &(trackdata.TrackYpfit) );

    // //Declare RICH-related branches of the tree:
    // //richdata are stored as STL vectors (basically dynamically sized arrays). Newer ROOT versions know how to handle this, older may not.
    // fTree->Branch("RICH_nhits", &(richdata.nhits_RICH), "nhits_RICH/I");
    // fTree->Branch("RICH_pmt", &(richdata.PMTnumber) );
    // fTree->Branch("RICH_row", &(richdata.row) );
    // fTree->Branch("RICH_col", &(richdata.col) );
    // fTree->Branch("RICH_nphe", &(richdata.NumPhotoelectrons) );
    // fTree->Branch("RICH_tavg", &(richdata.Time_avg) );
    // fTree->Branch("RICH_trms", &(richdata.Time_rms) );
    // fTree->Branch("RICH_mID", &(richdata.mTrackNo) );
    // fTree->Branch("RICH_vol", &(richdata.volume_flag) );
    // fTree->Branch("RICH_xhit", &(richdata.xhit) );
    // fTree->Branch("RICH_yhit", &(richdata.yhit) );
    // fTree->Branch("RICH_zhit", &(richdata.zhit) );
    // fTree->Branch("RICH_pxhit", &(richdata.pxhit) );
    // fTree->Branch("RICH_pyhit", &(richdata.pyhit) );
    // fTree->Branch("RICH_pzhit", &(richdata.pzhit) );
    // fTree->Branch("RICH_vxhit", &(richdata.pvx) );
    // fTree->Branch("RICH_vyhit", &(richdata.pvy) );
    // fTree->Branch("RICH_vzhit", &(richdata.pvz) );
    // fTree->Branch("RICH_vpxhit", &(richdata.ppx) );
    // fTree->Branch("RICH_vpyhit", &(richdata.ppy) );
    // fTree->Branch("RICH_vpzhit", &(richdata.ppz) );
    // //The RICH parent track data only gets filled if /tracking/storeTrajectory 1 has been issued by the user (because keeping track of this info is CPU and memory intensive).
    // fTree->Branch("RICH_ntracks", &(richdata.ntracks_RICH), "ntracks_RICH/I");
    // fTree->Branch("RICH_mPID", &(richdata.mPID) );
    // fTree->Branch("RICH_mTID", &(richdata.mTID) );
    // fTree->Branch("RICH_mMID", &(richdata.mMID) );
    // fTree->Branch("RICH_mvx", &(richdata.mvx) );
    // fTree->Branch("RICH_mvy", &(richdata.mvy) );
    // fTree->Branch("RICH_mvz", &(richdata.mvz) );
    // fTree->Branch("RICH_mpx", &(richdata.mpx) );
    // fTree->Branch("RICH_mpy", &(richdata.mpy) );
    // fTree->Branch("RICH_mpz", &(richdata.mpz) );

    // //ECal uses the same approach as RICH
    // fTree->Branch("ECal_nhits", &(ecaldata.nhits_ECal), "nhits_ECal/I");
    // fTree->Branch("ECal_pmt", &(ecaldata.PMTnumber) );
    // fTree->Branch("ECal_row", &(ecaldata.row) );
    // fTree->Branch("ECal_col", &(ecaldata.col) );
    // fTree->Branch("ECal_nphe", &(ecaldata.NumPhotoelectrons) );
    // fTree->Branch("ECal_tavg", &(ecaldata.Time_avg) );
    // fTree->Branch("ECal_trms", &(ecaldata.Time_rms) );
    //fTree->Print();
    return;
}

void G4SBSIO::FillTree(){
    if( !fTree ){ 
	fprintf(stderr, "Error %s: %s line %d - Trying to fill non-existant tree\n", __PRETTY_FUNCTION__, __FILE__, __LINE__ );
	return; 
    }

    fTree->Fill();
}

void G4SBSIO::WriteTree(){
    assert( fFile );
    assert( fTree );
    if( !fFile->IsOpen() ){
	G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << ": TFile not open" << G4endl;
	exit(1);
    }

    fFile->cd();
    fTree->Write("T", TObject::kOverwrite);

    Esum_histograms->Compress();
    Esum_histograms->Write();
    PulseShape_histograms->Compress();
    PulseShape_histograms->Write();
    
    G4SBSRun::GetRun()->GetData()->Write("run_data", TObject::kOverwrite);

    // Produce and write out field map graphics
    fGlobalField->DebugField( gendata.thbb, gendata.thsbs );
    for( vector<TH2F *>::iterator it = fGlobalField->fFieldPlots.begin(); it!= fGlobalField->fFieldPlots.end(); it++ ){
	(*it)->Write((*it)->GetName(), TObject::kOverwrite );
	delete (*it);
    }
    fGlobalField->fFieldPlots.clear();

    fTree->ResetBranchAddresses();
    delete fTree;
    fTree = NULL;

    fFile->Close();
    delete fFile;
    fFile = NULL;

    return;
}

void G4SBSIO::BranchGEM(G4String SDname="GEM"){
  TString branch_prefix = SDname.data();
  TString branch_name;
  
  branch_prefix.ReplaceAll("/",".");
 
  //Branches with raw GEM data:
  
  fTree->Branch( branch_name.Format( "%s.hit.nhits", branch_prefix.Data() ), &(GEMdata[SDname].nhits_GEM) );
  fTree->Branch( branch_name.Format( "%s.hit.plane", branch_prefix.Data() ), &(GEMdata[SDname].plane) );
  fTree->Branch( branch_name.Format( "%s.hit.strip", branch_prefix.Data() ), &(GEMdata[SDname].strip) );
  fTree->Branch( branch_name.Format( "%s.hit.x", branch_prefix.Data() ), &(GEMdata[SDname].x) );
  fTree->Branch( branch_name.Format( "%s.hit.y", branch_prefix.Data() ), &(GEMdata[SDname].y) );
  fTree->Branch( branch_name.Format( "%s.hit.z", branch_prefix.Data() ), &(GEMdata[SDname].z) );
  fTree->Branch( branch_name.Format( "%s.hit.polx", branch_prefix.Data() ), &(GEMdata[SDname].polx) );
  fTree->Branch( branch_name.Format( "%s.hit.poly", branch_prefix.Data() ), &(GEMdata[SDname].poly) );
  fTree->Branch( branch_name.Format( "%s.hit.polz", branch_prefix.Data() ), &(GEMdata[SDname].polz) );
  fTree->Branch( branch_name.Format( "%s.hit.t", branch_prefix.Data() ), &(GEMdata[SDname].t) );
  fTree->Branch( branch_name.Format( "%s.hit.trms", branch_prefix.Data() ), &(GEMdata[SDname].trms) );
  fTree->Branch( branch_name.Format( "%s.hit.tmin", branch_prefix.Data() ), &(GEMdata[SDname].tmin) );
  fTree->Branch( branch_name.Format( "%s.hit.tmax", branch_prefix.Data() ), &(GEMdata[SDname].tmax) );
  // fTree->Branch( branch_name.Format( "%s.hit.dx", branch_prefix.Data() ), &(GEMdata[SDname].dx) );
  // fTree->Branch( branch_name.Format( "%s.hit.dy", branch_prefix.Data() ), &(GEMdata[SDname].dy) );
  fTree->Branch( branch_name.Format( "%s.hit.tx", branch_prefix.Data() ), &(GEMdata[SDname].tx) );
  fTree->Branch( branch_name.Format( "%s.hit.ty", branch_prefix.Data() ), &(GEMdata[SDname].ty) );
  fTree->Branch( branch_name.Format( "%s.hit.xin", branch_prefix.Data() ), &(GEMdata[SDname].xin) );
  fTree->Branch( branch_name.Format( "%s.hit.yin", branch_prefix.Data() ), &(GEMdata[SDname].yin) );
  fTree->Branch( branch_name.Format( "%s.hit.zin", branch_prefix.Data() ), &(GEMdata[SDname].zin) );
  fTree->Branch( branch_name.Format( "%s.hit.xout", branch_prefix.Data() ), &(GEMdata[SDname].xout) );
  fTree->Branch( branch_name.Format( "%s.hit.yout", branch_prefix.Data() ), &(GEMdata[SDname].yout) );
  fTree->Branch( branch_name.Format( "%s.hit.zout", branch_prefix.Data() ), &(GEMdata[SDname].zout) );
  fTree->Branch( branch_name.Format( "%s.hit.txp", branch_prefix.Data() ), &(GEMdata[SDname].txp) );
  fTree->Branch( branch_name.Format( "%s.hit.typ", branch_prefix.Data() ), &(GEMdata[SDname].typ) );
  fTree->Branch( branch_name.Format( "%s.hit.xg", branch_prefix.Data() ), &(GEMdata[SDname].xg) );
  fTree->Branch( branch_name.Format( "%s.hit.yg", branch_prefix.Data() ), &(GEMdata[SDname].yg) );
  fTree->Branch( branch_name.Format( "%s.hit.zg", branch_prefix.Data() ), &(GEMdata[SDname].zg) );
  fTree->Branch( branch_name.Format( "%s.hit.trid", branch_prefix.Data() ), &(GEMdata[SDname].trid) );
  fTree->Branch( branch_name.Format( "%s.hit.mid", branch_prefix.Data() ), &(GEMdata[SDname].mid) );
  fTree->Branch( branch_name.Format( "%s.hit.pid", branch_prefix.Data() ), &(GEMdata[SDname].pid) );
  fTree->Branch( branch_name.Format( "%s.hit.vx", branch_prefix.Data() ), &(GEMdata[SDname].vx) );
  fTree->Branch( branch_name.Format( "%s.hit.vy", branch_prefix.Data() ), &(GEMdata[SDname].vy) );
  fTree->Branch( branch_name.Format( "%s.hit.vz", branch_prefix.Data() ), &(GEMdata[SDname].vz) );
  fTree->Branch( branch_name.Format( "%s.hit.p", branch_prefix.Data() ), &(GEMdata[SDname].p) );
  fTree->Branch( branch_name.Format( "%s.hit.edep", branch_prefix.Data() ), &(GEMdata[SDname].edep) );
  fTree->Branch( branch_name.Format( "%s.hit.beta", branch_prefix.Data() ), &(GEMdata[SDname].beta) );

  //Branches with "Tracker output" data:
  fTree->Branch( branch_name.Format("%s.Track.ntracks",branch_prefix.Data() ), &(trackdata[SDname].ntracks) );
  fTree->Branch( branch_name.Format("%s.Track.TID",branch_prefix.Data() ), &(trackdata[SDname].TrackTID) );
  fTree->Branch( branch_name.Format("%s.Track.PID",branch_prefix.Data() ), &(trackdata[SDname].TrackPID) );
  fTree->Branch( branch_name.Format("%s.Track.MID",branch_prefix.Data() ), &(trackdata[SDname].TrackMID) );
  fTree->Branch( branch_name.Format("%s.Track.NumHits",branch_prefix.Data() ), &(trackdata[SDname].NumHits) );
  fTree->Branch( branch_name.Format("%s.Track.NumPlanes",branch_prefix.Data() ), &(trackdata[SDname].NumPlanes) );
  fTree->Branch( branch_name.Format("%s.Track.NDF",branch_prefix.Data() ), &(trackdata[SDname].NDF) );
  fTree->Branch( branch_name.Format("%s.Track.Chi2fit",branch_prefix.Data() ), &(trackdata[SDname].Chi2fit) );
  fTree->Branch( branch_name.Format("%s.Track.Chi2true",branch_prefix.Data() ), &(trackdata[SDname].Chi2true) );
  fTree->Branch( branch_name.Format("%s.Track.X",branch_prefix.Data() ), &(trackdata[SDname].TrackX) );
  fTree->Branch( branch_name.Format("%s.Track.Y",branch_prefix.Data() ), &(trackdata[SDname].TrackY) );
  fTree->Branch( branch_name.Format("%s.Track.Xp",branch_prefix.Data() ), &(trackdata[SDname].TrackXp) );
  fTree->Branch( branch_name.Format("%s.Track.Yp",branch_prefix.Data() ), &(trackdata[SDname].TrackYp) );
  fTree->Branch( branch_name.Format("%s.Track.T",branch_prefix.Data() ), &(trackdata[SDname].TrackT) );
  fTree->Branch( branch_name.Format("%s.Track.P",branch_prefix.Data() ), &(trackdata[SDname].TrackP) );
  fTree->Branch( branch_name.Format("%s.Track.Sx",branch_prefix.Data() ), &(trackdata[SDname].TrackSx) );
  fTree->Branch( branch_name.Format("%s.Track.Sy",branch_prefix.Data() ), &(trackdata[SDname].TrackSy) );
  fTree->Branch( branch_name.Format("%s.Track.Sz",branch_prefix.Data() ), &(trackdata[SDname].TrackSz) );
  fTree->Branch( branch_name.Format("%s.Track.Xfit",branch_prefix.Data() ), &(trackdata[SDname].TrackXfit) );
  fTree->Branch( branch_name.Format("%s.Track.Yfit",branch_prefix.Data() ), &(trackdata[SDname].TrackYfit) );
  fTree->Branch( branch_name.Format("%s.Track.Xpfit",branch_prefix.Data() ), &(trackdata[SDname].TrackXpfit) );
  fTree->Branch( branch_name.Format("%s.Track.Ypfit",branch_prefix.Data() ), &(trackdata[SDname].TrackYpfit) );
  
  map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );

  if( it != KeepHistoryflags.end() && it->second ){
    //Branches with "Particle History" data:
    fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.npart) );
    fTree->Branch( branch_name.Format("%s.part.PID", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.PID) );
    fTree->Branch( branch_name.Format("%s.part.MID", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.MID) );
    fTree->Branch( branch_name.Format("%s.part.TID", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.TID) );
    fTree->Branch( branch_name.Format("%s.part.nbounce", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.nbounce) );
    fTree->Branch( branch_name.Format("%s.part.hitindex", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.hitindex) );
    fTree->Branch( branch_name.Format("%s.part.vx", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.vx) );
    fTree->Branch( branch_name.Format("%s.part.vy", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.vy) );
    fTree->Branch( branch_name.Format("%s.part.vz", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.vz) );
    fTree->Branch( branch_name.Format("%s.part.px", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.px) );
    fTree->Branch( branch_name.Format("%s.part.py", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.py) );
    fTree->Branch( branch_name.Format("%s.part.pz", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.pz) );
  }

  return;
}

void G4SBSIO::BranchCAL( G4String SDname="CAL" ){
  TString branch_prefix = SDname.data();
  TString branch_name;
  
  branch_prefix.ReplaceAll("/",".");

  TString histname;

  G4SBSCalSD *SDtemp = (G4SBSCalSD*) fdetcon->fSDman->FindSensitiveDetector( SDname );

  G4double threshtemp = SDtemp->GetEnergyThreshold();
  G4double gatewidthtemp = SDtemp->GetTimeWindow();
  G4int ntimebinstemp = SDtemp->GetNTimeBins();
  
  new( (*Esum_histograms)[fNhistograms] ) TH1F( histname.Format("%s.esum",branch_prefix.Data()), "",
						100, 0.0, 100.0*threshtemp );
  new( (*PulseShape_histograms)[fNhistograms] ) TH1F( histname.Format("%s.pulseshape",branch_prefix.Data()), "",
						      ntimebinstemp, 0.0,
						      gatewidthtemp );

  histogram_index[SDname] = fNhistograms;
  
  fNhistograms++;
			  
  //Define "hit" branches:
  fTree->Branch( branch_name.Format( "%s.det.esum", branch_prefix.Data() ), &(CALdata[SDname].Esum) );
  fTree->Branch( branch_name.Format( "%s.hit.nhits", branch_prefix.Data() ), &(CALdata[SDname].nhits_CAL) );
  fTree->Branch( branch_name.Format( "%s.hit.row", branch_prefix.Data() ), &(CALdata[SDname].row) );
  fTree->Branch( branch_name.Format( "%s.hit.col", branch_prefix.Data() ), &(CALdata[SDname].col) );
  fTree->Branch( branch_name.Format( "%s.hit.cell", branch_prefix.Data() ), &(CALdata[SDname].cell) );
  fTree->Branch( branch_name.Format( "%s.hit.plane", branch_prefix.Data() ), &(CALdata[SDname].plane) );
  fTree->Branch( branch_name.Format( "%s.hit.wire", branch_prefix.Data() ), &(CALdata[SDname].wire) );
  fTree->Branch( branch_name.Format( "%s.hit.xcell", branch_prefix.Data() ), &(CALdata[SDname].xcell) );
  fTree->Branch( branch_name.Format( "%s.hit.ycell", branch_prefix.Data() ), &(CALdata[SDname].ycell) );
  fTree->Branch( branch_name.Format( "%s.hit.zcell", branch_prefix.Data() ), &(CALdata[SDname].zcell) );
  fTree->Branch( branch_name.Format( "%s.hit.xcellg", branch_prefix.Data() ), &(CALdata[SDname].xcellg) );
  fTree->Branch( branch_name.Format( "%s.hit.ycellg", branch_prefix.Data() ), &(CALdata[SDname].ycellg) );
  fTree->Branch( branch_name.Format( "%s.hit.zcellg", branch_prefix.Data() ), &(CALdata[SDname].zcellg) );
  fTree->Branch( branch_name.Format( "%s.hit.xhit", branch_prefix.Data() ), &(CALdata[SDname].xhit) );
  fTree->Branch( branch_name.Format( "%s.hit.yhit", branch_prefix.Data() ), &(CALdata[SDname].yhit) );
  fTree->Branch( branch_name.Format( "%s.hit.zhit", branch_prefix.Data() ), &(CALdata[SDname].zhit) );
  fTree->Branch( branch_name.Format( "%s.hit.xhitg", branch_prefix.Data() ), &(CALdata[SDname].xhitg) );
  fTree->Branch( branch_name.Format( "%s.hit.yhitg", branch_prefix.Data() ), &(CALdata[SDname].yhitg) );
  fTree->Branch( branch_name.Format( "%s.hit.zhitg", branch_prefix.Data() ), &(CALdata[SDname].zhitg) );
  fTree->Branch( branch_name.Format( "%s.hit.sumedep", branch_prefix.Data() ), &(CALdata[SDname].sumedep) );
  fTree->Branch( branch_name.Format( "%s.hit.tavg", branch_prefix.Data() ), &(CALdata[SDname].tavg) );
  fTree->Branch( branch_name.Format( "%s.hit.trms", branch_prefix.Data() ), &(CALdata[SDname].trms) );
  fTree->Branch( branch_name.Format( "%s.hit.tmin", branch_prefix.Data() ), &(CALdata[SDname].tmin) );
  fTree->Branch( branch_name.Format( "%s.hit.tmax", branch_prefix.Data() ), &(CALdata[SDname].tmax) );

  map<G4String,G4bool>::iterator it = KeepPartCALflags.find( SDname );

  if( it != KeepPartCALflags.end() && it->second ){
    //Define "particle" branches:
    fTree->Branch( branch_name.Format( "%s.npart_CAL", branch_prefix.Data() ), &(CALdata[SDname].npart_CAL) );
    fTree->Branch( branch_name.Format( "%s.ihit", branch_prefix.Data() ), &(CALdata[SDname].ihit) );
    fTree->Branch( branch_name.Format( "%s.x", branch_prefix.Data() ), &(CALdata[SDname].x) );
    fTree->Branch( branch_name.Format( "%s.y", branch_prefix.Data() ), &(CALdata[SDname].y) );
    fTree->Branch( branch_name.Format( "%s.z", branch_prefix.Data() ), &(CALdata[SDname].z) );
    fTree->Branch( branch_name.Format( "%s.t", branch_prefix.Data() ), &(CALdata[SDname].t) );
    fTree->Branch( branch_name.Format( "%s.E", branch_prefix.Data() ), &(CALdata[SDname].E) );
    fTree->Branch( branch_name.Format( "%s.dt", branch_prefix.Data() ), &(CALdata[SDname].dt) );
    fTree->Branch( branch_name.Format( "%s.L", branch_prefix.Data() ), &(CALdata[SDname].L) );
    fTree->Branch( branch_name.Format( "%s.vx", branch_prefix.Data() ), &(CALdata[SDname].vx) );
    fTree->Branch( branch_name.Format( "%s.vy", branch_prefix.Data() ), &(CALdata[SDname].vy) );
    fTree->Branch( branch_name.Format( "%s.vz", branch_prefix.Data() ), &(CALdata[SDname].vz) );
    fTree->Branch( branch_name.Format( "%s.trid", branch_prefix.Data() ),  &(CALdata[SDname].trid) );
    fTree->Branch( branch_name.Format( "%s.mid", branch_prefix.Data() ), &(CALdata[SDname].mid) );
    fTree->Branch( branch_name.Format( "%s.pid", branch_prefix.Data() ), &(CALdata[SDname].pid) );
    fTree->Branch( branch_name.Format( "%s.p", branch_prefix.Data() ), &(CALdata[SDname].p) );
    fTree->Branch( branch_name.Format( "%s.px", branch_prefix.Data() ), &(CALdata[SDname].px) );
    fTree->Branch( branch_name.Format( "%s.py", branch_prefix.Data() ), &(CALdata[SDname].py) );
    fTree->Branch( branch_name.Format( "%s.pz", branch_prefix.Data() ), &(CALdata[SDname].pz) );
    fTree->Branch( branch_name.Format( "%s.edep", branch_prefix.Data() ), &(CALdata[SDname].edep) );
  }

  it = KeepHistoryflags.find( SDname );

  if( it != KeepHistoryflags.end() && it->second ){
    //Branches with "Particle History" data:
    fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.npart) );
    fTree->Branch( branch_name.Format("%s.part.PID", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.PID) );
    fTree->Branch( branch_name.Format("%s.part.MID", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.MID) );
    fTree->Branch( branch_name.Format("%s.part.TID", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.TID) );
    fTree->Branch( branch_name.Format("%s.part.nbounce", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.nbounce) );
    fTree->Branch( branch_name.Format("%s.part.hitindex", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.hitindex) );
    fTree->Branch( branch_name.Format("%s.part.vx", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.vx) );
    fTree->Branch( branch_name.Format("%s.part.vy", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.vy) );
    fTree->Branch( branch_name.Format("%s.part.vz", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.vz) );
    fTree->Branch( branch_name.Format("%s.part.px", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.px) );
    fTree->Branch( branch_name.Format("%s.part.py", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.py) );
    fTree->Branch( branch_name.Format("%s.part.pz", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.pz) );
  }

  return;
}

void G4SBSIO::BranchRICH(G4String SDname="RICH"){
  TString branch_prefix = SDname.data();
  TString branch_name;
  branch_prefix.ReplaceAll("/",".");
  
  //Branches for "hits": 
  
  fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(richdata[SDname].nhits_RICH) );
  fTree->Branch( branch_name.Format("%s.hit.PMT", branch_prefix.Data() ), &(richdata[SDname].PMTnumber) );
  fTree->Branch( branch_name.Format("%s.hit.row", branch_prefix.Data() ), &(richdata[SDname].row) );
  fTree->Branch( branch_name.Format("%s.hit.col", branch_prefix.Data() ), &(richdata[SDname].col) );
  fTree->Branch( branch_name.Format("%s.hit.xpmt", branch_prefix.Data() ), &(richdata[SDname].xpmt) );
  fTree->Branch( branch_name.Format("%s.hit.ypmt", branch_prefix.Data() ), &(richdata[SDname].ypmt) );
  fTree->Branch( branch_name.Format("%s.hit.zpmt", branch_prefix.Data() ), &(richdata[SDname].zpmt) );
  fTree->Branch( branch_name.Format( "%s.hit.xgpmt", branch_prefix.Data() ), &(richdata[SDname].xgpmt) );
  fTree->Branch( branch_name.Format( "%s.hit.ygpmt", branch_prefix.Data() ), &(richdata[SDname].ygpmt) );
  fTree->Branch( branch_name.Format( "%s.hit.zgpmt", branch_prefix.Data() ), &(richdata[SDname].zgpmt) );
  fTree->Branch( branch_name.Format("%s.hit.NumPhotoelectrons", branch_prefix.Data() ), &(richdata[SDname].NumPhotoelectrons) );
  fTree->Branch( branch_name.Format("%s.hit.Time_avg", branch_prefix.Data() ), &(richdata[SDname].Time_avg) );
  fTree->Branch( branch_name.Format("%s.hit.Time_rms", branch_prefix.Data() ), &(richdata[SDname].Time_rms) );
  fTree->Branch( branch_name.Format("%s.hit.Time_min", branch_prefix.Data() ), &(richdata[SDname].Time_min) );
  fTree->Branch( branch_name.Format("%s.hit.Time_max", branch_prefix.Data() ), &(richdata[SDname].Time_max) );
  fTree->Branch( branch_name.Format("%s.hit.mTrackNo", branch_prefix.Data() ), &(richdata[SDname].mTrackNo) );
  fTree->Branch( branch_name.Format("%s.hit.xhit", branch_prefix.Data() ), &(richdata[SDname].xhit) );
  fTree->Branch( branch_name.Format("%s.hit.yhit", branch_prefix.Data() ), &(richdata[SDname].yhit) );
  fTree->Branch( branch_name.Format("%s.hit.zhit", branch_prefix.Data() ), &(richdata[SDname].zhit) );
  fTree->Branch( branch_name.Format("%s.hit.pxhit", branch_prefix.Data() ), &(richdata[SDname].pxhit) );
  fTree->Branch( branch_name.Format("%s.hit.pyhit", branch_prefix.Data() ), &(richdata[SDname].pyhit) );
  fTree->Branch( branch_name.Format("%s.hit.pzhit", branch_prefix.Data() ), &(richdata[SDname].pzhit) );
  fTree->Branch( branch_name.Format("%s.hit.pvx", branch_prefix.Data() ), &(richdata[SDname].pvx) );
  fTree->Branch( branch_name.Format("%s.hit.pvy", branch_prefix.Data() ), &(richdata[SDname].pvy) );
  fTree->Branch( branch_name.Format("%s.hit.pvz", branch_prefix.Data() ), &(richdata[SDname].pvz) );
  fTree->Branch( branch_name.Format("%s.hit.ppx", branch_prefix.Data() ), &(richdata[SDname].ppx) );
  fTree->Branch( branch_name.Format("%s.hit.ppy", branch_prefix.Data() ), &(richdata[SDname].ppy) );
  fTree->Branch( branch_name.Format("%s.hit.ppz", branch_prefix.Data() ), &(richdata[SDname].ppz) );
  fTree->Branch( branch_name.Format("%s.hit.volume_flag", branch_prefix.Data() ), &(richdata[SDname].volume_flag) );
  
  //Branches for "tracks": This might be reorganized later:
  // branch_name.Format( "%s.ntracks_RICH", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].ntracks_RICH), "ntracks_RICH/I" );
  // branch_name.Format( "%s.mPID", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mPID) );
  // branch_name.Format( "%s.mTID", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mTID) );
  // branch_name.Format( "%s.mMID", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mMID) );
  // branch_name.Format( "%s.mvx", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mvx) );
  // branch_name.Format( "%s.mvy", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mvy) );
  // branch_name.Format( "%s.mvz", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mvz) );
  // branch_name.Format( "%s.mpx", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mpx) );
  // branch_name.Format( "%s.mpy", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mpy) );
  // branch_name.Format( "%s.mpz", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mpz) );

  map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );

  if( it != KeepHistoryflags.end() && it->second ){
    //Branches with "Particle History" data:
    fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.npart) );
    fTree->Branch( branch_name.Format("%s.part.PID", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.PID) );
    fTree->Branch( branch_name.Format("%s.part.MID", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.MID) );
    fTree->Branch( branch_name.Format("%s.part.TID", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.TID) );
    fTree->Branch( branch_name.Format("%s.part.nbounce", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.nbounce) );
    fTree->Branch( branch_name.Format("%s.part.hitindex", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.hitindex) );
    fTree->Branch( branch_name.Format("%s.part.vx", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.vx) );
    fTree->Branch( branch_name.Format("%s.part.vy", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.vy) );
    fTree->Branch( branch_name.Format("%s.part.vz", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.vz) );
    fTree->Branch( branch_name.Format("%s.part.px", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.px) );
    fTree->Branch( branch_name.Format("%s.part.py", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.py) );
    fTree->Branch( branch_name.Format("%s.part.pz", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.pz) );
    fTree->Branch( branch_name.Format("%s.part.Nphe_part", branch_prefix.Data() ), &(richdata[SDname].Nphe_part) );
  }
  return;
}

void G4SBSIO::BranchECAL(G4String SDname="ECAL"){
  TString branch_prefix = SDname.data();
  TString branch_name;
  branch_prefix.ReplaceAll("/",".");
  
  fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(ecaldata[SDname].nhits_ECal) );
  fTree->Branch( branch_name.Format("%s.hit.PMT", branch_prefix.Data() ), &(ecaldata[SDname].PMTnumber) );
  fTree->Branch( branch_name.Format("%s.hit.row", branch_prefix.Data() ), &(ecaldata[SDname].row) );
  fTree->Branch( branch_name.Format("%s.hit.col", branch_prefix.Data() ), &(ecaldata[SDname].col) );
  fTree->Branch( branch_name.Format("%s.hit.plane", branch_prefix.Data() ), &(ecaldata[SDname].plane) );
  fTree->Branch( branch_name.Format("%s.hit.xcell", branch_prefix.Data() ), &(ecaldata[SDname].xcell) );
  fTree->Branch( branch_name.Format("%s.hit.ycell", branch_prefix.Data() ), &(ecaldata[SDname].ycell) );
  fTree->Branch( branch_name.Format("%s.hit.zcell", branch_prefix.Data() ), &(ecaldata[SDname].zcell) );
  fTree->Branch( branch_name.Format("%s.hit.xgcell", branch_prefix.Data() ), &(ecaldata[SDname].xgcell) );
  fTree->Branch( branch_name.Format("%s.hit.ygcell", branch_prefix.Data() ), &(ecaldata[SDname].ygcell) );
  fTree->Branch( branch_name.Format("%s.hit.zgcell", branch_prefix.Data() ), &(ecaldata[SDname].zgcell) );
  fTree->Branch( branch_name.Format("%s.hit.NumPhotoelectrons", branch_prefix.Data() ), &(ecaldata[SDname].NumPhotoelectrons) );
  fTree->Branch( branch_name.Format("%s.hit.Time_avg", branch_prefix.Data() ), &(ecaldata[SDname].Time_avg) );
  fTree->Branch( branch_name.Format("%s.hit.Time_rms", branch_prefix.Data() ), &(ecaldata[SDname].Time_rms) );
  fTree->Branch( branch_name.Format("%s.hit.Time_min", branch_prefix.Data() ), &(ecaldata[SDname].Time_min) );
  fTree->Branch( branch_name.Format("%s.hit.Time_max", branch_prefix.Data() ), &(ecaldata[SDname].Time_max) );

}

void G4SBSIO::BranchPythia(){
  //Per-event variables:
  fTree->Branch("primaries.Sigma",&(Primaries.Sigma),"primaries.Sigma/D");
  fTree->Branch("primaries.Ebeam",&(Primaries.Ebeam),"primaries.Ebeam/D");
  fTree->Branch("primaries.Eprime",&(Primaries.Eprime),"primaries.Eprime/D");
  fTree->Branch("primaries.Q2",&(Primaries.Q2),"primaries.Q2/D");
  fTree->Branch("primaries.xbj",&(Primaries.xbj),"primaries.xbj/D");
  fTree->Branch("primaries.y",&(Primaries.y),"primaries.y/D");
  fTree->Branch("primaries.W2",&(Primaries.W2),"primaries.W2/D");
  fTree->Branch("primaries.theta_e",&(Primaries.theta_e),"primaries.theta_e/D");
  fTree->Branch("primaries.phi_e",&(Primaries.phi_e),"primaries.phi_e/D");
  fTree->Branch("primaries.px_e",&(Primaries.px_e),"primaries.px_e/D");
  fTree->Branch("primaries.py_e",&(Primaries.py_e),"primaries.py_e/D");
  fTree->Branch("primaries.pz_e",&(Primaries.pz_e),"primaries.pz_e/D");
  fTree->Branch("primaries.vx_e",&(Primaries.vx_e),"primaries.vx_e/D");
  fTree->Branch("primaries.vy_e",&(Primaries.vy_e),"primaries.vy_e/D");
  fTree->Branch("primaries.vz_e",&(Primaries.vz_e),"primaries.vz_e/D");
  fTree->Branch("primaries.Egamma",&(Primaries.Egamma),"primaries.Egamma/D");
  fTree->Branch("primaries.theta_gamma",&(Primaries.theta_gamma),"primaries.theta_gamma/D");
  fTree->Branch("primaries.phi_gamma",&(Primaries.phi_gamma),"primaries.phi_gamma/D");
  fTree->Branch("primaries.px_gamma",&(Primaries.px_gamma),"primaries.px_gamma/D");
  fTree->Branch("primaries.py_gamma",&(Primaries.py_gamma),"primaries.py_gamma/D");
  fTree->Branch("primaries.pz_gamma",&(Primaries.pz_gamma),"primaries.pz_gamma/D");
  fTree->Branch("primaries.vx_gamma",&(Primaries.vx_gamma),"primaries.vx_gamma/D");
  fTree->Branch("primaries.vy_gamma",&(Primaries.vy_gamma),"primaries.vy_gamma/D");
  fTree->Branch("primaries.vz_gamma",&(Primaries.vz_gamma),"primaries.vz_gamma/D");
  //Primary particle arrays:
  fTree->Branch("Primaries.Nprimaries",&(Primaries.Nprimaries),"Nprimaries/I");
  fTree->Branch("Primaries.PID",&(Primaries.PID));
  fTree->Branch("Primaries.genflag",&(Primaries.genflag));
  fTree->Branch("Primaries.Px",&(Primaries.Px));
  fTree->Branch("Primaries.Py",&(Primaries.Py));
  fTree->Branch("Primaries.Pz",&(Primaries.Pz));
  fTree->Branch("Primaries.vx",&(Primaries.vx));
  fTree->Branch("Primaries.vy",&(Primaries.vy));
  fTree->Branch("Primaries.vz",&(Primaries.vz));
  fTree->Branch("Primaries.M",&(Primaries.M));
  fTree->Branch("Primaries.E",&(Primaries.E));
  fTree->Branch("Primaries.P",&(Primaries.P));
  fTree->Branch("Primaries.t",&(Primaries.t));
  fTree->Branch("Primaries.theta",&(Primaries.theta));
  fTree->Branch("Primaries.phi",&(Primaries.phi));
}

void G4SBSIO::UpdateGenDataFromDetCon(){ //Go with whatever is in fdetcon as of run start for constant parameters of the run describing detector layout:
  gendata.thbb = fdetcon->fEArmBuilder->fBBang;
  gendata.dbb  = fdetcon->fEArmBuilder->fBBdist/CLHEP::m;
  gendata.thsbs = fdetcon->fHArmBuilder->f48D48ang;
  gendata.dsbs  = fdetcon->fHArmBuilder->f48D48dist/CLHEP::m;
  gendata.dhcal = fdetcon->fHArmBuilder->fHCALdist/CLHEP::m;
  gendata.voffhcal = fdetcon->fHArmBuilder->fHCALvertical_offset/CLHEP::m;
  gendata.hoffhcal = fdetcon->fHArmBuilder->fHCALhorizontal_offset/CLHEP::m;
  gendata.dlac = fdetcon->fHArmBuilder->fLACdist/CLHEP::m;
  gendata.vofflac = fdetcon->fHArmBuilder->fLACvertical_offset/CLHEP::m;
  gendata.hofflac = fdetcon->fHArmBuilder->fLAChorizontal_offset/CLHEP::m;
  gendata.drich = fdetcon->fHArmBuilder->fRICHdist/CLHEP::m;
  gendata.dsbstrkr = fdetcon->fHArmBuilder->fSBS_tracker_dist/CLHEP::m;
  gendata.sbstrkrpitch = fdetcon->fHArmBuilder->fSBS_tracker_pitch; //radians
}
