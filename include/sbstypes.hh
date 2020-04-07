#ifndef SBSTYPES_HH
#define SBSTYPES_HH

#include "TTimeStamp.h"

#define MAXTARG 5
#define NNUCL   2

#define __RUNSTR_LEN 255
#define __MAXFILE_LEN 1048576 // MB file size limit

// Include fictional neutron target
enum Targ_t { kH2, kD2, kLH2, kLD2, k3He, kNeutTarg };
enum Nucl_t { kProton, kNeutron };
enum Hadron_t { kPiPlus, kPiMinus, kPi0, kKPlus, kKMinus, kP, kPbar}; //Hadron types for SIDIS event generator
enum Kine_t { kElastic, // 0
	      kFlat, 
	      kInelastic,
	      kDIS, 
	      kBeam, 
	      kSIDIS, //5
	      kGun, 
	      kWiser, 
	      kPYTHIA6, 
	      kGMnElasticCheck, 
	      kCosmics, //10
	      kTDISKin, 
	      kAcquMC, 
	      kTDISGen, 
	      tElastic, 
	      tQuasiElastic, //15 
	      tInelastic, 
	      tTDISKin }; // TDIS, added kTDISGen, tElastic, tQuasiElastic, tInelastic, tTDISKin, used in my class (CA)



enum Exp_t { kGEp, kNeutronExp, kSIDISExp, kC16, kA1n, kTDIS, kNDVCS, kGEnRP};
enum Arm_t { kEarm, kHarm }; //Types for association of detector modules with spectrometer arms. Presently "E arm" and "H arm" are possible.
enum SDet_t { kGEM, kCAL, kRICH, kECAL, kmTPC }; //Types of sensitive detectors (others to be added later)

struct filedata_t {
    char filename[__RUNSTR_LEN];
    char hashsum[__RUNSTR_LEN];
    TTimeStamp timestamp;
};


#endif//SBSTYPES_HH
