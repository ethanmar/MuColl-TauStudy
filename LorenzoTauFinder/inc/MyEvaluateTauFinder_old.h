#ifndef MyEvaluateTauFinder_h
#define MyEvaluateTauFinder_h 1

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TTree.h"
#include "UTIL/LCRelationNavigator.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <string>

using namespace lcio;
using namespace marlin;
#define NTAU_MAX 100

struct MyParticle;

/**  My Evaluation processor for TauFinder
 *
 * @author A. Muennich, CERN, L. Valla, LIP
 */

class MyEvaluateTauFinder : public Processor {

public:
  virtual Processor *newProcessor() { return new MyEvaluateTauFinder; }

  MyEvaluateTauFinder();
  MyEvaluateTauFinder(const MyEvaluateTauFinder &) = delete;
  MyEvaluateTauFinder &operator=(const MyEvaluateTauFinder &) = delete;

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
  virtual void LoopDaughters(MCParticle *particle, double &Evis, double &pxvis,
                             double &pyvis, double &pzvis);
  virtual void LoopDaughtersRelation(MCParticle *particle,
                                     LCRelationNavigator *relationNavigatorTau,
                                     LCRelationNavigator *relationNavigatorMC,
                                     bool &ralToTau);
  virtual int FindDecayMode(MCParticle *tau);

protected:
  /** Input collection name.
   */
  std::string _colNameMC{}, _colNameTrack{};
  std::string _colNameMCTruth{}, _incol{}, _colNamePFORecLink{};
  std::string _colNameMCRecLink{}, _colNameTracksRecLink{},
      _colNameTauRecLink{};
  std::string _OutputFile_Signal{};
  float _bField = 0.0;
  int _nRun = -1;
  int _nEvt = -1;

  double _ntot_rec = 0.0;
  double _ntot_mc = 0.0;
  double _ntau_correct = 0.0;
  double _dEsum = 0.0;
  double _dEsumsq = 0.0;
  int _ndE = 0.0;

  float mcE[NTAU_MAX], mcPx[NTAU_MAX], mcPy[NTAU_MAX], mcPz[NTAU_MAX],
      mcPt[NTAU_MAX], mcPhi[NTAU_MAX], mcEta[NTAU_MAX], mcD0[NTAU_MAX];
  float mcE_vis[NTAU_MAX], mcPx_vis[NTAU_MAX], mcPy_vis[NTAU_MAX],
      mcPz_vis[NTAU_MAX], mcPt_vis[NTAU_MAX], mcPhi_vis[NTAU_MAX],
      mcEta_vis[NTAU_MAX];
  float recE[NTAU_MAX], recPx[NTAU_MAX], recPy[NTAU_MAX], recPz[NTAU_MAX],
      recPt[NTAU_MAX], recPhi[NTAU_MAX], recEta[NTAU_MAX], recD0[NTAU_MAX],
      recESeed[NTAU_MAX];
  int charge[NTAU_MAX], recNTracks[NTAU_MAX], recNQTracks[NTAU_MAX],
      mcDecayMode[NTAU_MAX];
  int matched[NTAU_MAX], missed[NTAU_MAX], fake[NTAU_MAX];
  int pfosPdg[5 * NTAU_MAX];
  float pfosDeltaR[5 * NTAU_MAX], pfosPt[5 * NTAU_MAX];

  int ntau_mc = 0, ntau_rec = 0, ntau_missed = 0, nfakes = 0, ntau_matched = 0,
      npfos = 0;

  TFile *rootfile = NULL;
  TTree *tree = NULL;
};

#endif