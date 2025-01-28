#ifndef MyTauFinder_h
#define MyTauFinder_h 1

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TTree.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <string>

using namespace lcio;
using namespace marlin;

#define NTAU_MAX 10000

/** MyTauFinder processor for marlin.
 *
 * @author A. Muennich, CERN
 *
 */

struct MyParticle;

class MyTauFinder : public Processor {

public:
  virtual Processor *newProcessor() { return new MyTauFinder; }

  MyTauFinder();
  MyTauFinder(const MyTauFinder &) = delete;
  MyTauFinder &operator=(const MyTauFinder &) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

protected:
  /** Input collection name.
   */
  std::string _colNameMC{}, _colNameRECO{}, _incol{};
  std::string _colNameMCTruth{}, _colNameTauRecLink{};
  std::string _outcol{}, _outcolRest{};
  std::string _OutputFile_Signal{};
  int _nRun = -1;
  int _nEvt = -1;

  float _ptcut = 0.0, _ptseed = 0.0, _cosTcut = 0.0;
  float _coneAngle = 0.0, _isoAngle = 0.0, _isoE = 0.0;
  float _D0seedmin = 0.0, _D0seedmax = 0.0, _minv = 0.0;

  int _fail_pt_seed = 0, _fail_minv = 0, _fail_minv_neg = 0, _fail_Qtr = 0,
      _fail_isoE = 0, _mergeTries = 0;

  TFile *rootfile = NULL;
  TNtuple *failtuple = NULL;
  TTree *anatree = NULL;

  bool FindTau(std::vector<ReconstructedParticle *> &Qvec,
               std::vector<ReconstructedParticle *> &Nvec,
               std::vector<std::vector<ReconstructedParticle *>> &tauvec);

private:
  int _ntau{};
  int _ngood{};
  int _nfail_seed{};
  int _nfail_Qtrack{};
  int _nrej{};
  int _nrej_isoE{};
  int _nrej_minv{};
  int _nrej_nQ{};
  int _nrej_nQN{};

  float _tau_isoE[NTAU_MAX]{};
  float _tau_minv[NTAU_MAX]{};
  float _tau_pt[NTAU_MAX]{};
  float _tau_p[NTAU_MAX]{};
  float _tau_ene[NTAU_MAX]{};
  float _tau_phi[NTAU_MAX]{};
  float _tau_eta[NTAU_MAX]{};
  float _tau_nQ[NTAU_MAX]{};
  float _tau_nN[NTAU_MAX]{};
  float _tau_nQN[NTAU_MAX]{};
  int _tau_good[NTAU_MAX]{};
  // float _tau_nmerged[NTAU_MAX]{};

  /*float _rej_isoE[NTAU_MAX]{};
  float _rej_minv[NTAU_MAX]{};
  float _rej_pt[NTAU_MAX]{};
  float _rej_ene[NTAU_MAX]{};
  float _rej_nQ[NTAU_MAX]{};
  float _rej_nQN[NTAU_MAX]{};*/
  // float _rej_nmerged[NTAU_MAX]{};
};

#endif