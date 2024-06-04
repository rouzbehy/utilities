#include <vector>
#include <cmath>
#include "MARTINI.h"
#include "fastjet/ClusterSequence.hh"
using namespace std;
//using namespace std::chrono;
using namespace fastjet;
// constants:
#define PI 3.1415926535
//LHC:
#define PHOTON_ETA_CUT 0.8
#define CH_HAD_ETA_CUT 1.0
#define PARTON_ETA_CUT 2.0
//RHIC:
//#define PHOTON_ETA_CUT 0.35
//#define CH_HAD_ETA_CUT 0.5
//
#define JET_ETA_CUT 2.0
#define ETA_TRK_MAX 2.5

#define JET_PTMIN_FF 100
#define JET_PTMAX_FF 398
#define JET_ETA_FF 2.1

#define JET_PTMIN_SHAPE 100
#define JET_ETA_MIN_SHAPE 0.3
#define JET_ETA_MAX_SHAPE 2.0

#define PARTON_JET_PTMIN_SHAPE 20
#define PARTON_JET_ETA_MIN_SHAPE -1.0
#define PARTON_JET_ETA_MAX_SHAPE 1.0

#define trk_cut_shape 1
#define QUARK 0
#define GLUON 1

#define HARD 0
#define RECOIL 1
#define HOLE -1

// spectra before clustering
void bin_partons_by_z(vector<Parton>* plist, vector<double>& gluon_hist, vector<double>& fermion_hist, vector<double>& zbins, double pTRef, double tau);
void bin_parton(vector<Parton>* plist, vector<double>& hist, vector<double>& pTbins, int recoil_flag, int particle_type);
void bin_chgd_hads(Event& evt, vector<double>& hist, vector<double>& pTbins, double weight);
void bin_identified_had(Event& evt, vector<double>& hist, vector<double>& pTbins, double weight, int pid);

// preparation for clustering 
void cluster_prep_partonic(vector<Parton>* plist, vector<PseudoJet>& ps, int hole_flag);
void cluster_prep_inclusive(Event& evt, vector<PseudoJet>& hads);
void cluster_prep_charged(Event& evt, vector<PseudoJet>& hads);
void cluster_prep_inclusive_wth_trk_cut(Event& evt, vector<PseudoJet>& hads, double trk_cut);

// Subtract Holes from clustered jets
void subtract_holes_from_jets(vector<PseudoJet>& jets, vector<PseudoJet>& holes, double R);

// Binning Routines: Jets
void bin_jets(vector<PseudoJet>& jets, vector<double>& hist, vector<double>& pTbins);
int bin_jet_fragmentation(vector<PseudoJet>& jets, vector<double>& FFz, 
     vector<double>& zbins, vector<double>& FFpT, vector<double>& pTbins);
int bin_jet_shape(vector<PseudoJet>& jets, vector<double>& shape, vector<double>& rbins);

int bin_partonic_shape(vector<PseudoJet>& jets, vector<double>& shape, vector<double>& rbins, int particle_type);
// Photons
void bin_photons(vector<Parton>* plist, vector<double>& prompt, 
                 vector<double>& conv, vector<double>& brem, vector<double>& pTbins);
