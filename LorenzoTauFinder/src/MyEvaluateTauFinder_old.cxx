#include "MyTauFinder/MyEvaluateTauFinder.h"
#include "MyTauFinder/HelixClass.h"
#include <algorithm>
#include <iomanip>
#include <iostream>

// NOTE: ALL THE getGeneratorStatus() WERE INCREASED BY ONE (THEY WERE PROBABLY
// WRONG)

using namespace std;

#ifdef MARLIN_USE_AIDA
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogramFactory.h>
#include <marlin/AIDAProcessor.h>
// #include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/LCObject.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/VertexImpl.h>
#include <UTIL/LCRelationNavigator.h>

#include <marlin/Global.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#define coutEv -1

#define coutUpToEv 0

using namespace lcio;
using namespace marlin;
using namespace UTIL;

MyEvaluateTauFinder aEvaluateTauFinder;

struct TAU { // Declare  struct type
  float E, px, py, pz, pt, phi, eta, D0, Evis, pxvis, pyvis, pzvis, ptvis,
      phivis, etavis, minv, Eseed;
  int NQ, N, decaymode, missed, matched, fake, charge; // Declare member types
};

bool My_AngleSort(TAU p1, TAU p2) { return fabs(p1.phi) > fabs(p2.phi); }

MyEvaluateTauFinder::MyEvaluateTauFinder() : Processor("MyEvaluateTauFinder") {
  // modify processor description
  _description = "MyEvaluateTauFinder checks performance of TauFinder and "
                 "writes output to root file.";

  // register steering parameters: name, description, class-variable, default
  // value

  registerProcessorParameter(
      "FileName_Signal", "Name of the Signal output file ", _OutputFile_Signal,
      std::string("EvalTauFinder_out.root"));

  registerProcessorParameter("B_Field",
                             "Value of Bz field (in T) at the origin ", _bField,
                             float(3.57));

  registerInputCollection(LCIO::MCPARTICLE, "MCCollectionName",
                          "Name of the MCParticle collection", _colNameMC,
                          std::string("MCParticlesSkimmed"));

  registerInputCollection(LCIO::LCRELATION, "RECOMCTRUTHCollectionName",
                          "Name of the MC Truth PFA collection",
                          _colNameMCTruth, std::string("RecoMCTruthLink"));

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "TauRecCollection",
                          "Collection of Tau Candidates", _incol,
                          std::string("TauRec_PFO"));

  registerInputCollection(LCIO::LCRELATION, "TauLinkCollectionName",
                          "Name of the link between tau and pfos",
                          _colNameTauRecLink, std::string("TauRecLink_PFO"));
}

void MyEvaluateTauFinder::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // usually a good idea to
  printParameters();

  _nRun = 0;
  _nEvt = 0;
  _ntot_rec = 0;
  _ntot_mc = 0;
  _ntau_correct = 0;
  _dEsum = 0;
  _dEsumsq = 0;
  _ndE = 0;

  ntau_mc = 0;
  ntau_rec = 0;
  ntau_missed = 0;
  nfakes = 0;
  ntau_matched = 0;

  rootfile = new TFile((_OutputFile_Signal).c_str(), "RECREATE");

  tree = new TTree("evtree", "evtree");

  tree->Branch("RunID", &_nRun, "RunID/I");
  tree->Branch("EvID", &_nEvt, "EvID/I");
  tree->Branch("nTausMC", &ntau_mc, "nTausMC/I");
  tree->Branch("nTausRec", &ntau_rec, "nTausRec/I");
  tree->Branch("nTausMiss", &ntau_missed, "nTausMiss/I");
  tree->Branch("nTausFake", &nfakes, "nTausFake/I");
  tree->Branch("Charge", charge, "Charge[nTausRec]/I");
  tree->Branch("mcDecayMode", mcDecayMode, "mcDecayMode[nTausMC]/I");
  tree->Branch("mcE", mcE, "mcE[nTausMC]/F");
  tree->Branch("mcPx", mcPx, "mcPx[nTausMC]/F");
  tree->Branch("mcPy", mcPy, "mcPy[nTausMC]/F");
  tree->Branch("mcPz", mcPz, "mcPz[nTausMC]/F");
  tree->Branch("mcPt", mcPt, "mcPt[nTausMC]/F");
  tree->Branch("mcPhi", mcPhi, "mcPhi[nTausMC]/F");
  tree->Branch("mcEta", mcEta, "mcEta[nTausMC]/F");
  tree->Branch("mcD0", mcD0, "mcD0[nTausMC]/F");
  tree->Branch("mcE_vis", mcE_vis, "mcE_vis[nTausMC]/F");
  tree->Branch("mcPx_vis", mcPx_vis, "mcPx_vis[nTausMC]/F");
  tree->Branch("mcPy_vis", mcPy_vis, "mcPy_vis[nTausMC]/F");
  tree->Branch("mcPz_vis", mcPz_vis, "mcPz_vis[nTausMC]/F");
  tree->Branch("mcPt_vis", mcPt_vis, "mcPt_vis[nTausMC]/F");
  tree->Branch("mcPhi_vis", mcPhi_vis, "mcPhi_vis[nTausMC]/F");
  tree->Branch("mcEta_vis", mcEta_vis, "mcEta_vis[nTausMC]/F");
  tree->Branch("recE", recE, "recE[nTausRec]/F");
  tree->Branch("recPx", recPx, "recPx[nTausRec]/F");
  tree->Branch("recPy", recPy, "recPy[nTausRec]/F");
  tree->Branch("recPz", recPz, "recPz[nTausRec]/F");
  tree->Branch("recPt", recPt, "recPt[nTausRec]/F");
  tree->Branch("recPhi", recPhi, "recPhi[nTausRec]/F");
  tree->Branch("recEta", recEta, "recEta[nTausRec]/F");
  tree->Branch("recD0", recD0, "recD0[nTausRec]/F");
  tree->Branch("recNTracks", recNTracks, "recNTracks[nTausRec]/I");
  tree->Branch("recNQTracks", recNQTracks, "recNQTracks[nTausRec]/I");
  tree->Branch("recESeed", recESeed, "recESeed[nTausRec]/F");

  tree->Branch("Matched", matched, "Matched[nTausRec]/I");
  tree->Branch("Missed", missed, "Missed[nTausMC]/I");
  tree->Branch("Fake", fake, "Matched[nTausRec]/I");
}

void MyEvaluateTauFinder::processRunHeader(LCRunHeader *) { _nRun++; }

void MyEvaluateTauFinder::processEvent(LCEvent *evt) {
  // this gets called for every event
  // usually the working horse ...

  LCCollection *colMC, *colMCTruth, *colTau;
  LCCollection *colTauRecLink;
  try {
    colMC = evt->getCollection(_colNameMC);
  } catch (Exception &e) {
    colMC = 0;
  }

  try {
    colTau = evt->getCollection(_incol);
  } catch (Exception &e) {
    colTau = 0;
  }

  try {
    colMCTruth = evt->getCollection(_colNameMCTruth);
  } catch (Exception &e) {
    colMCTruth = 0;
  }

  try {
    colTauRecLink = evt->getCollection(_colNameTauRecLink);
  } catch (Exception &e) {
    colTauRecLink = 0;
  }

  _nEvt = evt->getEventNumber();

  ntau_mc = 0;
  ntau_rec = 0;
  ntau_missed = 0;
  nfakes = 0;
  ntau_matched = 0;

  LCRelationNavigator *relationNavigatorTau = 0;
  LCRelationNavigator *relationNavigatorPFOMC = 0;
  if (colTauRecLink != 0) {
    relationNavigatorTau = new LCRelationNavigator(colTauRecLink);
  }

  if (colMCTruth != 0)
    relationNavigatorPFOMC = new LCRelationNavigator(colMCTruth);

  bool isfake = false;
  nfakes = 0;
  int LD_rec = 0;

  std::vector<TAU> rectauvec;

  if (colTau != 0) {
    int nT = colTau->getNumberOfElements();
    ntau_rec = nT;
    if (nT == 0) {
      // std::cout << _nEvt << '\n';
    }
    int LD = 0;
    if (_nEvt < coutUpToEv || _nEvt == coutEv)
      streamlog_out(DEBUG) << "EVENT " << _nEvt << " with " << nT << " taus"
                           << endl;
    HelixClass *helix = new HelixClass();
    HelixClass *mc_helix = new HelixClass();
    for (int k = 0; k < nT; k++) {
      ReconstructedParticle *tau =
          static_cast<ReconstructedParticle *>(colTau->getElementAt(k));
      const double *pvec = tau->getMomentum();
      double pt = sqrt(pvec[0] * pvec[0] + pvec[1] * pvec[1]);
      double p =
          sqrt(pvec[0] * pvec[0] + pvec[1] * pvec[1] + pvec[2] * pvec[2]);
      double px = pvec[0];
      double py = pvec[1];
      double pz = pvec[2];
      double phi = TMath::ATan2(pvec[1], pvec[0]);
      double theta = TMath::ATan2(pt, pvec[2]);
      double eta = -std::log(TMath::Tan(theta / 2.));

      // tauvec contains
      std::vector<ReconstructedParticle *> tauvec = tau->getParticles();
      double Eseed = 0., D0 = 0.;

      /*for (int daugh = 0; daugh < tauvec.size(); ++daugh)
      {
              std::cout << tauvec[daugh]->getType() << ", ";
      }
      std::cout << '\n';
      */

      int NQ = 0;
      for (unsigned int o = 0; o < tauvec.size(); o++) {
        // find seed track for D0
        if (tauvec[o]->getCharge() != 0) {
          NQ++;
          if (tauvec[o]->getEnergy() > Eseed &&
              tauvec[o]->getTracks().size() != 0) {
            D0 = (float)tauvec[o]->getTracks()[0]->getD0();
            Eseed = tauvec[o]->getEnergy();
          }
        }
        // check for leptonic decay
        if (abs(tauvec[o]->getType()) == 11 ||
            abs(tauvec[o]->getType()) == 13) {
          LD = 1;
          ++LD_rec;
        }
      }
      // float mom[3];
      // 	  float ver[3];

      // 	  for (int icomp=0; icomp<3; ++icomp) {
      // 	    mom[icomp]=(float)tau->getMomentum()[icomp];
      // 	    VertexImpl
      // *vtx=static_cast<VertexImpl*>(tau->getStartVertex()); 	    if(vtx)
      // 	      {
      // 		const float *vpos=vtx->getPosition();
      // 		ver[icomp]=vpos[icomp];
      // 	      }
      // 	    else
      // 	      ver[icomp]=0;
      // 	  }

      // 	  float charge = tau->getCharge();
      // 	  helix->Initialize_VP(ver,mom,charge,_bField);
      // double D0=fabs(helix->getD0());
      if (_nEvt < coutUpToEv || _nEvt == coutEv)
        streamlog_out(DEBUG) << tau->getEnergy() << " " << phi << " " << eta
                             << " " << D0 << " " << tauvec.size() << endl;

      TAU rtau;
      rtau.E = tau->getEnergy();
      rtau.px = px;
      rtau.py = py;
      rtau.pz = pz;
      rtau.pt = pt;
      rtau.phi = phi;
      rtau.eta = eta;
      rtau.D0 = D0;
      rtau.N = tauvec.size();
      rtau.NQ = NQ;
      rtau.charge = tau->getCharge();
      rtau.Eseed = Eseed;
      double mass_inv = 0.;
      if (tau->getEnergy() * tau->getEnergy() < p * p)
        mass_inv = tau->getEnergy() - sqrt(p * p);
      else
        mass_inv = sqrt(tau->getEnergy() * tau->getEnergy() - p * p);
      rtau.minv = mass_inv;
      // std::cout << "E * E = " << tau->getEnergy() * tau->getEnergy() << ", p
      // * p = " << p * p << ", mass_inv = " << mass_inv << '\n';

      _ntot_rec++;

      // follow the chain back to mc truth
      if (relationNavigatorTau) {
        bool istau = false;
        bool contaminated = false;
        MCParticle *mctau = NULL;
        EVENT::LCObjectVec relobjFROM =
            relationNavigatorTau->getRelatedToObjects(tau);
        // std::cout << "size = " << relobjFROM.size() << '\n';
        for (unsigned int o = 0; o < relobjFROM.size(); o++) {
          ReconstructedParticle *rec =
              static_cast<ReconstructedParticle *>(relobjFROM[o]);
          if (relationNavigatorPFOMC) {
            EVENT::LCObjectVec relobjMC =
                relationNavigatorPFOMC->getRelatedToObjects(rec);

            for (unsigned int m = 0; m < relobjMC.size(); m++) {
              MCParticle *mc = static_cast<MCParticle *>(relobjMC[m]);
              // check whether particles parent is really a tau:
              MCParticle *dummy = mc;
              MCParticle *parent = mc;
              // need to catch broken relations in DST file for background
              // particles
              if (mc == 0) {
                // std::cout<<"Broken Relation: "<<_nEvt<<" "<<rec->getType()<<"
                // "<<rec->getEnergy()<<std::endl;
                continue;
              }
              int size = mc->getParents().size();

              while (size != 0) {
                dummy = parent->getParents()[0];
                size = dummy->getParents().size();
                parent = dummy;
                if (abs(parent->getPDG()) == 15)
                  size = 0;
              }
              if (abs(parent->getPDG()) == 15) {
                istau = true;
                mctau = parent;
              } else
                contaminated = true;
            }
          }
        }

        // compare tau with mc truth
        if (mctau) {
          int mc_decaymode = FindDecayMode(mctau);

          float mc_mom[3];
          float mc_ver[3];
          const double *mc_pvec = mctau->getMomentum();
          double mc_pt =
              sqrt(mc_pvec[0] * mc_pvec[0] + mc_pvec[1] * mc_pvec[1]);
          double mc_phi = TMath::ATan2(mc_pvec[1], mc_pvec[0]);
          double mc_theta = TMath::ATan2(mc_pt, mc_pvec[2]);
          double mc_eta = -std::log(TMath::Tan(mc_theta / 2.));

          for (int icomp = 0; icomp < 3; ++icomp) {
            mc_mom[icomp] = (float)mctau->getMomentum()[icomp];
            mc_ver[icomp] = (float)mctau->getDaughters()[0]->getVertex()[icomp];
          }
          float mc_charge = mctau->getCharge();
          mc_helix->Initialize_VP(mc_ver, mc_mom, mc_charge, _bField);
          double mc_D0 = fabs(mc_helix->getD0());
          double Evis = 0., pxvis = 0., pyvis = 0., pzvis = 0.;

          LoopDaughters(mctau, Evis, pxvis, pyvis, pzvis);

          double ptvis = sqrt(pxvis * pxvis + pyvis * pyvis);

          // taumatchtuple->Fill();
          // if (!contaminated)
          // tauexacttuple->Fill();
          _dEsum += Evis - tau->getEnergy();
          _dEsumsq += (Evis - tau->getEnergy()) * (Evis - tau->getEnergy());
          _ndE++;
        }
        if (istau) {
          _ntau_correct++;
          ++ntau_matched;
        } else {
          int d1 = 0, d2 = 0, pdg = 0;
          for (unsigned int o = 0; o < relobjFROM.size(); o++) {
            ReconstructedParticle *rec =
                static_cast<ReconstructedParticle *>(relobjFROM[o]);
            if (relationNavigatorPFOMC) {
              EVENT::LCObjectVec relobj =
                  relationNavigatorPFOMC->getRelatedToObjects(rec);
              for (unsigned int m = 0; m < relobj.size(); m++) {
                MCParticle *mc = static_cast<MCParticle *>(relobj[m]);
                // need to catch broken relations in DST file for background
                // particles
                if (mc == 0)
                  continue;
                if (mc->getCharge() == 0)
                  continue;
                MCParticle *dummy = mc;
                MCParticle *parent = mc;
                int size = mc->getParents().size();
                while (size != 0) {
                  dummy = parent->getParents()[0];
                  size = dummy->getParents().size();
                  parent = dummy;
                  if (parent->getGeneratorStatus() == 2)
                    break;
                }
                pdg = parent->getPDG();
                if (parent->getDaughters().size())
                  d1 = parent->getDaughters()[0]->getPDG();
                if (parent->getDaughters().size() > 1)
                  d2 = parent->getDaughters()[1]->getPDG();
              }
            }
          }
          // faketuple->Fill();
          isfake = true;
          nfakes++;
        }

        rtau.fake = isfake;
        rtau.matched = istau;
      } // relNavTau
      rectauvec.push_back(rtau);
    }

    delete helix;
    delete mc_helix;
  }

  int D1 = 0, D2 = 0, D3 = 0, D4 = 0;
  std::vector<TAU> mctauvec;
  if (colMC != 0) {
    // std::cout << "ColMC if entered" << '\n';
    int nMCP = colMC->getNumberOfElements();
    if (_nEvt < coutUpToEv || _nEvt == coutEv)
      streamlog_out(DEBUG) << "MCTRUTH: " << endl;
    HelixClass *helix = new HelixClass();
    for (int k = 0; k < nMCP; k++) {
      MCParticle *particle = static_cast<MCParticle *>(colMC->getElementAt(k));

      if (particle->getGeneratorStatus() == 2 &&
          abs(particle->getPDG()) == 15 &&
          particle->getDaughters().size() > 1) {
        // std::cout << "getGenStatus if entered" << '\n';
        if (abs(particle->getDaughters()[0]->getPDG()) == 15)
          continue;
        ntau_mc++;
        _ntot_mc++;
        int mc_decaymode = FindDecayMode(particle);

        const double *pvec = particle->getMomentum();
        double pt = sqrt(pvec[0] * pvec[0] + pvec[1] * pvec[1]);
        double phi = TMath::ATan2(pvec[1], pvec[0]);
        // double theta = 180. / TMath::Pi() * TMath::ATan2(pt, pvec[2]);
        double theta = TMath::ATan2(pt, pvec[2]);
        double eta = -std::log(TMath::Tan(theta / 2.));

        float mom[3];
        float ver[3];
        for (int icomp = 0; icomp < 3; ++icomp) {
          mom[icomp] = (float)particle->getMomentum()[icomp];
          ver[icomp] = (float)particle->getDaughters()[0]->getVertex()[icomp];
        }

        float charge = particle->getCharge();
        helix->Initialize_VP(ver, mom, charge, _bField);
        double D0 = fabs(helix->getD0());
        double Evis = 0, pxvis = 0, pyvis = 0, pzvis = 0;

        LoopDaughters(particle, Evis, pxvis, pyvis, pzvis);

        double ptvis = sqrt(pxvis * pxvis + pyvis * pyvis);
        double phivis = TMath::ATan2(pyvis, pxvis);
        // double theta = 180. / TMath::Pi() * TMath::ATan2(pt, pvec[2]);
        double thetavis = TMath::ATan2(ptvis, pzvis);
        double etavis = -std::log(TMath::Tan(thetavis / 2.));

        TAU mctau;
        // std::cout << "mctau energy = " << Evis << '\n';
        mctau.E = particle->getEnergy();
        mctau.px = pvec[0];
        mctau.py = pvec[1];
        mctau.pz = pvec[2];
        mctau.pt = pt;
        mctau.phi = phi;
        mctau.eta = eta;
        mctau.D0 = D0;
        mctau.Evis = Evis;
        mctau.pxvis = pxvis;
        mctau.pyvis = pyvis;
        mctau.pzvis = pzvis;
        mctau.ptvis = ptvis;
        mctau.phivis = phivis;
        mctau.etavis = etavis;
        mctau.decaymode = mc_decaymode;

        if (_nEvt < coutUpToEv || _nEvt == coutEv)
          streamlog_out(DEBUG)
              << Evis << " " << phi << " " << eta << " " << D0 << endl;

        bool ismissed = false;

        // find out which mc taus do not have a link to the rec
        if (relationNavigatorPFOMC && relationNavigatorTau) {
          bool hasRel = false;
          LoopDaughtersRelation(particle, relationNavigatorTau,
                                relationNavigatorPFOMC, hasRel);
          if (!hasRel) {
            ntau_missed++;
            int d1 = 0, d2 = 0;
            if (particle->getDaughters().size() == 2) {
              d1 = particle->getDaughters()[0]->getPDG();
              d2 = particle->getDaughters()[1]->getPDG();
            }
            // mcmisstuple->Fill();
            ismissed = true;
            if (_nEvt < coutUpToEv || _nEvt == coutEv)
              streamlog_out(DEBUG) << "Missed: " << Evis << " " << D0 << " "
                                   << d1 << " " << d2 << endl;
          }
        }
        mctau.missed = ismissed;
        mctauvec.push_back(mctau);
      } // tau

      if (particle->getGeneratorStatus() < 3 &&
          abs(particle->getPDG()) == 24) {
        if (particle->getPDG() == 24 && particle->getDaughters().size() == 2) {
          D1 = particle->getDaughters()[0]->getPDG();
          D2 = particle->getDaughters()[1]->getPDG();
        }
        if (particle->getPDG() == -24 && particle->getDaughters().size() == 2) {
          D3 = particle->getDaughters()[0]->getPDG();
          D4 = particle->getDaughters()[1]->getPDG();
        }
      }
    }
    delete helix;
  }
  // if (isfake) topofaketuple->Fill();
  //  sort the mc t and rec taus for comparison
  std::sort(mctauvec.begin(), mctauvec.end(), My_AngleSort);
  std::sort(rectauvec.begin(), rectauvec.end(), My_AngleSort);

  for (unsigned int p = 0; p < mctauvec.size(); p++) {
    mcE[p] = mctauvec[p].E;
    mcPx[p] = mctauvec[p].px;
    mcPy[p] = mctauvec[p].py;
    mcPz[p] = mctauvec[p].pz;
    mcPt[p] = mctauvec[p].pt;
    mcPhi[p] = mctauvec[p].phi;
    mcEta[p] = mctauvec[p].eta;
    mcD0[p] = mctauvec[p].D0;
    mcE_vis[p] = mctauvec[p].Evis;
    mcPx_vis[p] = mctauvec[p].pxvis;
    mcPy_vis[p] = mctauvec[p].pyvis;
    mcPz_vis[p] = mctauvec[p].pzvis;
    mcPt_vis[p] = mctauvec[p].ptvis;
    mcPhi_vis[p] = mctauvec[p].phivis;
    mcEta_vis[p] = mctauvec[p].etavis;

    mcDecayMode[p] = mctauvec[p].decaymode;
    missed[p] = mctauvec[p].missed;
  }

  for (unsigned int p = 0; p < rectauvec.size(); p++) {
    recE[p] = rectauvec[p].E;
    recPx[p] = rectauvec[p].px;
    recPy[p] = rectauvec[p].py;
    recPz[p] = rectauvec[p].pz;
    recPt[p] = rectauvec[p].pt;
    recPhi[p] = rectauvec[p].phi;
    recEta[p] = rectauvec[p].eta;
    recD0[p] = rectauvec[p].D0;
    recESeed[p] = rectauvec[p].Eseed;
    recNTracks[p] = rectauvec[p].N;
    recNQTracks[p] = rectauvec[p].NQ;

    charge[p] = rectauvec[p].charge;
    matched[p] = rectauvec[p].matched;
    fake[p] = mctauvec[p].fake;
  }

  tree->Fill();

  _nEvt++;
  // cleanup
  delete relationNavigatorTau;
  delete relationNavigatorPFOMC;
}

void MyEvaluateTauFinder::LoopDaughters(MCParticle *particle, double &Evis,
                                        double &pxvis, double &pyvis,
                                        double &pzvis) {
  for (unsigned int d = 0; d < particle->getDaughters().size(); d++) {
    // std::cout << "particle->getDaughters().size() = " <<
    // particle->getDaughters().size() << '\n';
    MCParticle *daughter = particle->getDaughters()[d];
    // only particles visible in the detector
    if (daughter->hasLeftDetector() == false ||
        abs(daughter->getPDG()) == 13) {
      if (daughter->getGeneratorStatus() == 1) {
        // filter out the neutrinos and other invisibles
        if (!(abs(daughter->getPDG()) == 12 ||
              abs(daughter->getPDG()) == 14 ||
              abs(daughter->getPDG()) == 16 ||
              abs(daughter->getPDG()) == 1000022)) {
          //		  if(_nEvt<coutUpToEv || _nEvt==coutEv)
          // streamlog_out(DEBUG) <<"D vis "<<d<<" "<<daughter->getPDG()<<"
          // "<<daughter->getEnergy()<<endl;
          Evis += daughter->getEnergy();
          const double *mc_pvec = daughter->getMomentum();
          pxvis += mc_pvec[0];
          pyvis += mc_pvec[1];
          pzvis += mc_pvec[2];
        }
      }
    }
    // THIS WAS CHANGED - MAYBE WE ONLY HAVE TO LOOP OVER THE FINAL PARTICLES
    // WITHOUT DAUGHTERS, THUS PUTTING THESE BEFORE THE IF ABOVE? if
    // (daughter->getDaughters().size()) LoopDaughters(daughter, Evis, ptvis,
    // pvis);
  }
}

void MyEvaluateTauFinder::LoopDaughtersRelation(
    MCParticle *particle, LCRelationNavigator *relationNavigatorTau,
    LCRelationNavigator *relationNavigatorPFOMC, bool &relToTau) {
  for (unsigned int d = 0; d < particle->getDaughters().size(); d++) {
    MCParticle *daughter = particle->getDaughters()[d];
    // only particles visible in the detector
    if (daughter->hasLeftDetector() == false ||
        abs(daughter->getPDG()) == 13) {
      if (daughter->getGeneratorStatus() == 1) {
        // relation to the filled reconstructed particle
        EVENT::LCObjectVec relobjTO =
            relationNavigatorPFOMC->getRelatedFromObjects(daughter);
        for (unsigned int o = 0; o < relobjTO.size(); o++) {
          // relation to the reconstructed tau
          ReconstructedParticle *rec =
              static_cast<ReconstructedParticle *>(relobjTO[o]);
          EVENT::LCObjectVec relobj =
              relationNavigatorTau->getRelatedFromObjects(rec);
          if (relobj.size())
            relToTau = true;
        }
      }
    }
    if (relToTau)
      break;
    if (daughter->getDaughters().size())
      LoopDaughtersRelation(daughter, relationNavigatorTau,
                            relationNavigatorPFOMC, relToTau);
  }
}

int MyEvaluateTauFinder::FindDecayMode(MCParticle *tau) {
  /*
   *Decay modes for tau-:
   *0: -211, 16, 111 (pi-, pi0, nu-tau) (25.49%)
   *1: -211, 16 (pi- nu-tau) (10.82%)
   *2: -211, 16, 111, 111 (pi-, pi0x2, nu-tau) (9.26%)
   *3: -211, -211, 16, 211 (3-prong plus nu-tau) (8.99%)
   *4: -211, -211, 16, 111, 211 (3-prong plus pi0 and nu-tau) (2.74%)
   *5: -211, 16, 111, 111, 111 (3 pi0, pi-, nu-tau) (1.04%)
   *6: -12, 11, 16 (nu-tau, e, nu-e-bar) (17.82%)
   *7: -14, 13, 16 (nu-tau, mu, nu-mu-bar) (17.39%)

   *Decay modes for tau+:
   *0: 211, 111, -16 (pi+, pi0, nu-tau-bar) (25.49%)
   *1: 211, -16 (pi+, nu-tau-bar) (10.82%)
   *2: 211, 111, 111, -16 (pi+, pi0x2, nu-tau-bar) (9.26%)
   *3: 211, 211, -16, 211 (3-prong plus nu-tau-bar) (8.99%)
   *4: 211, 211, 111, -16, -211 (3-prong plus pi0 and nu-tau-bar) (2.74%)
   *5: 211, 111, 111, 111, -16 (3 pi0, pi+, nu-tau-bar) (1.04%)
   *6: 12, -11, -16 (nu-e, e+, nu-tau-bar) (17.82%)
   *7: 14, -13, -16 (nu-mu, mu+, nu-tau-bar) (17.39%)*/

  const int n_daughters = tau->getDaughters().size();

  std::vector<int> daughters_pdgs = {};

  /*
  for (int i = 0; i < n_daughters; i++)
  {
          daughters_pdgs.push_back(tau->getDaughters()[i]->getPDG());
  }*/

  for (auto const &daughter : tau->getDaughters())
    daughters_pdgs.push_back(daughter->getPDG());

  if (tau->getPDG() == 15) {
    std::sort(daughters_pdgs.begin(), daughters_pdgs.end());
    if (daughters_pdgs[0] == -211) {
      // hadronic decays
      if (daughters_pdgs[1] == -211) {
        // 3-prongs
        if (daughters_pdgs[3] == 211) {
          // 3-prong, no neutral pion
          return 3;
        } else if (daughters_pdgs[3] == 111) {
          // 3-prong with neutral pion
          return 4;
        }
      } else if (daughters_pdgs[1] == 16) {
        //
        if (n_daughters == 2) {
          return 1;
        } else if (daughters_pdgs[2] == 111) {
          if (n_daughters == 3) {
            return 0;
          } else if (n_daughters == 4) {
            return 2;
          } else if (n_daughters == 5) {
            return 5;
          }
        }
      }
    } else if (daughters_pdgs[0] == -12) {
      return 6;
    } else if (daughters_pdgs[0] == -14) {
      return 7;
    } else {
      return -1;
    }
  } else if (tau->getPDG() == -15) {
    std::sort(daughters_pdgs.begin(), daughters_pdgs.end(),
              std::greater<int>());
    if (daughters_pdgs[0] == +211) {
      // hadronic decays
      if (daughters_pdgs[1] == +211) {
        // 3-prongs
        if (daughters_pdgs[2] == -16) {
          // 3-prong, no neutral pion
          return 3;
        } else if (daughters_pdgs[2] == 111) {
          // 3-prong with neutral pion
          return 4;
        }
      } else if (daughters_pdgs[1] == 111) {
        //
        if (n_daughters == 2) {
          return 1;
        } else if (n_daughters == 3) {
          return 0;
        } else if (n_daughters == 4) {
          return 2;
        } else if (n_daughters == 5) {
          return 5;
        }
      }
    } else if (daughters_pdgs[0] == 12) {
      return 6;
    } else if (daughters_pdgs[0] == 14) {
      return 7;
    }
  }
  return -1;
}

void MyEvaluateTauFinder::check(LCEvent *) {
  // nothing to check here - could be used to fill checkplots in reconstruction
  // processor
}

void MyEvaluateTauFinder::end() {

  streamlog_out(DEBUG) << "MyEvaluateTauFinder::end()  " << name()
                       << " processed " << _nEvt << " events in " << _nRun
                       << " runs " << std::endl;

  //   evtuple->Write();
  //   tautuple->Write();
  //   mcmisstuple->Write();
  //   taumatchtuple->Write();
  //   tauexacttuple->Write();
  //   faketuple->Write();
  //   topofaketuple->Write();

  tree->Write();
  // Close File here
  rootfile->Write();
  delete rootfile;
  // rootfile->Close();
}