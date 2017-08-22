#include <memory>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <cstdint>

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"

#include "L1Trigger/L1TMuonEndCap/interface/TrackTools.h"


// From L1Trigger/L1TMuonEndCap/interface/MuonTriggerPrimitive.h
class TriggerPrimitive {
public:
  enum subsystem_type{kDT,kCSC,kRPC,kGEM,kNSubsystems};
};


// _____________________________________________________________________________
class RPCIntegration : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RPCIntegration(const edm::ParameterSet& iConfig);
  ~RPCIntegration();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  virtual void endJob() override;

  // Main functions
  void process();

  // Aux functions
  void getHandles(const edm::Event& iEvent);

  void bookHistograms();
  void writeHistograms();

  // Configurables
  const edm::InputTag emuHitTag_;
  const edm::InputTag emuTrackTag_;
  const edm::InputTag genPartTag_;
  const std::string outFileName_;
  int verbose_;

  // Member data
  edm::EDGetTokenT<l1t::EMTFHitCollection>   emuHitToken_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection> emuTrackToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPartToken_;

  l1t::EMTFHitCollection    emuHits_;
  l1t::EMTFTrackCollection  emuTracks_;
  reco::GenParticleCollection genParts_;

  std::map<TString, TH1F*> histograms_;
  std::map<TString, TH2F*> histogram2Ds_;
};


// _____________________________________________________________________________
//static std::random_device rd;
//static std::mt19937 genrd(rd());
static std::mt19937 genrd0(20230);
static std::mt19937 genrd1(20231);


RPCIntegration::RPCIntegration(const edm::ParameterSet& iConfig) :
    emuHitTag_    (iConfig.getParameter<edm::InputTag>("emuHitTag")),
    emuTrackTag_  (iConfig.getParameter<edm::InputTag>("emuTrackTag")),
    genPartTag_   (iConfig.getParameter<edm::InputTag>("genPartTag")),
    outFileName_  (iConfig.getParameter<std::string>  ("outFileName")),
    verbose_      (iConfig.getUntrackedParameter<int> ("verbosity"))
{
  usesResource("TFileService");

  emuHitToken_   = consumes<l1t::EMTFHitCollection>     (emuHitTag_);
  emuTrackToken_ = consumes<l1t::EMTFTrackCollection>   (emuTrackTag_);
  genPartToken_  = consumes<reco::GenParticleCollection>(genPartTag_);
}

RPCIntegration::~RPCIntegration() {}

void RPCIntegration::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  getHandles(iEvent);
  process();
}

// _____________________________________________________________________________
void RPCIntegration::getHandles(const edm::Event& iEvent) {
  // EMTF hits and tracks
  edm::Handle<decltype(emuHits_)>   emuHits_handle;
  edm::Handle<decltype(emuTracks_)> emuTracks_handle;

  if (!emuHitToken_.isUninitialized()) {
    iEvent.getByToken(emuHitToken_, emuHits_handle);
  }
  if (!emuHits_handle.isValid()) {
    edm::LogError("RPCIntegration") << "Cannot get the product: " << emuHitTag_;
    return;
  }

  if (!emuTrackToken_.isUninitialized()) {
    iEvent.getByToken(emuTrackToken_, emuTracks_handle);
  }
  if (!emuTracks_handle.isValid()) {
    edm::LogError("RPCIntegration") << "Cannot get the product: " << emuTrackTag_;
    return;
  }

  // Gen particles
  edm::Handle<decltype(genParts_)> genParts_handle;

  if (!iEvent.isRealData()) {
    if (!genPartToken_.isUninitialized()) {
      iEvent.getByToken(genPartToken_, genParts_handle);
    }
    if (!genParts_handle.isValid()) {
      edm::LogError("RPCIntegration") << "Cannot get the product: " << genPartTag_;
      return;
    }
  }

  // Object filters
  emuHits_.clear();
  for (const auto& hit : (*emuHits_handle)) {
    if (!(-1 <= hit.BX() && hit.BX() <= 1))  continue;  // only BX=[-1,+1]
    if (hit.Endcap() != 1)  continue;  // only positive endcap
    emuHits_.push_back(hit);
  }

  emuTracks_.clear();
  for (const auto& trk : (*emuTracks_handle)) {
    if (trk.BX() != 0)      continue;  // only BX=0
    if (trk.Endcap() != 1)  continue;  // only positive endcap
    emuTracks_.push_back(trk);
  }

  genParts_.clear();
  for (const auto& part : (*genParts_handle)) {
    //if (!(part.pt() >= 2.))     continue;  // only pT > 2
    if (!(1.24 <= part.eta() && part.eta() <= 2.4))  continue;  // only positive endcap
    genParts_.push_back(part);
  }

  // Reduced CSC efficiency
  {
    std::uniform_real_distribution<double> rng(0.,1.);

    l1t::EMTFHitCollection tmp_emuHits;
    for (const auto& hit : emuHits_) {
      bool keep = true;
      if (hit.Subsystem() == TriggerPrimitive::kCSC) {
        if (hit.Station() == 1) {
          //double eff = 0.8;
          //if (!(rng(genrd0) < eff))  keep = false;
        } else if (hit.Station() == 2) {
          //double eff = 0.8;
          //if (!(rng(genrd0) < eff))  keep = false;
        } else if (hit.Station() == 3) {
          //double eff = 0.8;
          //if (!(rng(genrd0) < eff))  keep = false;
        } else if (hit.Station() == 4) {
          //double eff = 0.8;
          //if (!(rng(genrd0) < eff))  keep = false;
        }
      }
      if (keep) {
        tmp_emuHits.push_back(hit);
      }
    }
    std::swap(emuHits_, tmp_emuHits);
  }
}

// _____________________________________________________________________________
void RPCIntegration::process() {
  TString hname;
  TH1F* h;

  auto get_mode_bin = [](const auto& trk) {
    int mode      = trk.Mode();
    assert(0 < mode && mode <= 15);
    if (mode == 15)                              return 3;
    if (mode == 11 || mode == 13 || mode == 14)  return 2;
    if (mode ==  7 || mode == 10 || mode == 12)  return 1;
    if (mode ==  3 || mode ==  5 || mode ==  6 || mode == 9)  return 0;
    return -1;
  };

  auto get_l1pt_bin = [](const auto& trk) {
    float pt      = trk.Pt();
    assert(pt >= 0.);
    if (pt >= 100.)  return 3;
    if (pt >=  22.)  return 2;
    if (pt >=  12.)  return 1;
    if (pt >=   0.)  return 0;
    return -1;
  };

  auto get_gen_eta_bin = [](const auto& part) {
    double absEta = std::abs(part.eta());
    if (absEta < 1.24)  return -1;
    if (absEta < 1.60)  return 0;
    if (absEta < 2.00)  return 1;
    if (absEta < 2.40)  return 2;
    return -1;
  };

  auto get_gen_eta_any_bin = [](const auto& part) {
    return 3;
  };

  // ___________________________________________________________________________
  bool keep_event = true;

  if (verbose_ > 0) {
    int ipart = 0;
    for (const auto& part : genParts_) {
      //std::cout << "part " << ipart++ << " " << part.pdgId() << " " << part.status() << " " << part.pt() << " " << part.eta() << " " << part.phi() << std::endl;
      std::cout << "part " << ipart++ << " " << part.pt() << " " << part.eta() << " " << part.phi() << std::endl;
    }

    int itrack = 0;
    for (const auto& trk : emuTracks_) {
      std::cout << "trk " << itrack++ << " " << trk.Pt() << " " << trk.GMT_eta() << " " << trk.GMT_phi() << " " << trk.Mode() << std::endl;
    }

    int ihit = 0;
    for (const auto& hit : emuHits_) {
      if (hit.Subsystem() == TriggerPrimitive::kCSC) {
        std::cout << "CSC hit " << ihit++ << " " << hit.Sector() << " " << hit.Station() << " " << hit.Ring() << " " << hit.Chamber() << " " << hit.BX() << std::endl;
      }
      if (hit.Subsystem() == TriggerPrimitive::kRPC) {
        std::cout << "RPC hit " << ihit++ << " " << hit.Sector() << " " << hit.Station() << " " << hit.Ring() << " " << hit.Subsector() << " " << hit.BX() << std::endl;
      }
      if (hit.Subsystem() == TriggerPrimitive::kGEM) {
        std::cout << "GEM hit " << ihit++ << " " << hit.Sector() << " " << hit.Station() << " " << hit.Ring() << " " << hit.Chamber() << " " << hit.BX() << std::endl;
      }
    }
  }

  if (genParts_.empty()) {
    keep_event = false;
  }

  if (keep_event) {
    const auto& part = genParts_.front();

    bool trigger = !emuTracks_.empty();

    //if (genParts_.size() != 1) {
    //  std::cout << "[WARNING] perche non uno? num of genParts: " << genParts_.size() << std::endl;
    //}

    //if (trigger && emuTracks_.size() != 1) {
    //  std::cout << "[WARNING] perche non uno? num of emuTracks: " << emuTracks_.size() << std::endl;
    //}

    int mode_bin    = trigger ? get_mode_bin(emuTracks_.front()) : -1;
    int l1pt_bin    = trigger ? get_l1pt_bin(emuTracks_.front()) : -1;
    int eta_bin     = get_gen_eta_bin(part);
    int eta_any_bin = get_gen_eta_any_bin(part);

    const int    charge = part.charge();
    const double pt     = part.pt();
    const double absEta = std::abs(part.eta());
    assert(charge == 1 || charge == -1);

    l1t::EMTFTrack ideal_track;
    ideal_track.set_mode(0);

    for (int sector = 1; sector <= 6; ++sector) {
      int mode = 0;
      for (const auto& hit : emuHits_) {
        assert(1 <= hit.PC_sector() && hit.PC_sector() <= 6);
        assert(1 <= hit.Station() && hit.Station() <= 4);
        if (hit.PC_sector() == sector) {
          mode |= (1<<(4-hit.Station()));
        }
      }
      if (ideal_track.Mode() < mode) {
        ideal_track.set_mode(mode);
      }
      //if (trigger && emuTracks_.front().Sector() == sector) {
      //  ideal_track.set_mode(mode);
      //}
    }
    assert(ideal_track.Mode() <= 15);

    int ideal_mode_bin = !emuHits_.empty() ? get_mode_bin(ideal_track) : -1;
    mode_bin = ideal_mode_bin;


    // _________________________________________________________________________
    // Fill histograms

    // Efficiency vs pT
    for (int i=0; i<4; ++i) {  // mode
      for (int j=0; j<4; ++j) {  // eta
        if (!(j == eta_bin || j == eta_any_bin)) continue;  // gen eta binning

        hname = Form("denom_vs_pt_mode%i_eta%i", i, j);
        h = histograms_.at(hname);
        h->Fill(pt);

        if (trigger) {
          if (!(i <= mode_bin)) continue;
          //if (!(2 <= l1pt_bin)) continue;  // trigger pT cut

          hname = Form("num_vs_pt_mode%i_eta%i", i, j);
          h = histograms_.at(hname);
          h->Fill(pt);
        }
      }
    }

    // Efficiency vs eta
    for (int i=0; i<4; ++i) {  // mode
      for (int j=0; j<4; ++j) {  // pT
        if (!(pt >= 20.)) continue;  // gen pT requirement

        hname = Form("denom_vs_eta_mode%i_l1pt%i", i, j);
        h = histograms_.at(hname);
        h->Fill(absEta);

        if (trigger) {
          if (!(i <= mode_bin)) continue;
          if (!(j <= l1pt_bin)) continue;

          hname = Form("num_vs_eta_mode%i_l1pt%i", i, j);
          h = histograms_.at(hname);
          h->Fill(absEta);
        }
      }
    }

    // Station efficiency vs eta
    for (int j=0; j<4; ++j) {  // station
      if (!(pt >= 20.)) continue;  // gen pT requirement

      hname = Form("denom_vs_eta_st%i", j);
      h = histograms_.at(hname);
      h->Fill(absEta);

      bool good_station = ideal_track.Mode() & (1<<(3-j));
      if (good_station) {
        hname = Form("num_vs_eta_st%i", j);
        h = histograms_.at(hname);
        h->Fill(absEta);
      }
    }

    // Deflection angles of GEM-CSC and RPC-CSC
    for (int j=0; j<4; ++j) {  // station
      if (trigger) {
        std::vector<l1t::EMTFHit> myhits0;
        const int sector = emuTracks_.front().Sector();
        std::copy_if(emuHits_.begin(), emuHits_.end(), std::back_inserter(myhits0), [sector](const auto& hit) { return (hit.PC_sector() == sector); });

        std::vector<l1t::EMTFHit> myhits1;
        std::vector<l1t::EMTFHit> myhits2;

        // For GEM-CSC
        if (j == 0) {  // station 1
          std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits1), [](const auto& hit) { return (hit.Station() == 1 && (hit.Ring() == 1 || hit.Ring() == 4) && hit.Subsystem() == TriggerPrimitive::kGEM); });
          std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits2), [](const auto& hit) { return (hit.Station() == 1 && (hit.Ring() == 1 || hit.Ring() == 4) && hit.Subsystem() == TriggerPrimitive::kCSC); });
        } else if (j == 1) {  // station 2
          std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits1), [](const auto& hit) { return (hit.Station() == 2 && hit.Ring() == 1 && hit.Subsystem() == TriggerPrimitive::kGEM); });
          std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits2), [](const auto& hit) { return (hit.Station() == 2 && hit.Ring() == 1 && hit.Subsystem() == TriggerPrimitive::kCSC); });
        } else if (j == 2) {  // station 3 (iRPC)
          std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits1), [](const auto& hit) { return (hit.Station() == 3 && hit.Ring() == 1 && hit.Subsystem() == TriggerPrimitive::kRPC); });
          std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits2), [](const auto& hit) { return (hit.Station() == 3 && hit.Ring() == 1 && hit.Subsystem() == TriggerPrimitive::kCSC); });
        } else if (j == 3) {  // station 4 (iRPC)
          std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits1), [](const auto& hit) { return (hit.Station() == 4 && hit.Ring() == 1 && hit.Subsystem() == TriggerPrimitive::kRPC); });
          std::copy_if(myhits0.begin(), myhits0.end(), std::back_inserter(myhits2), [](const auto& hit) { return (hit.Station() == 4 && hit.Ring() == 1 && hit.Subsystem() == TriggerPrimitive::kCSC); });
        } else {
          continue;
        }

        if (myhits1.size() > 0 && myhits2.size() > 0) {
          std::uniform_int_distribution<> index1(0, myhits1.size()-1);
          std::uniform_int_distribution<> index2(0, myhits2.size()-1);
          const l1t::EMTFHit& myhit1 = myhits1.at(index1(genrd1));
          const l1t::EMTFHit& myhit2 = myhits2.at(index2(genrd1));

          int ipt = -1;
          if ((1.0/2 - 0.01) < 1.0/pt && 1.0/pt <= (1.0/2)) {
            ipt = 0;
          } else if ((1.0/3 - 0.01) < 1.0/pt && 1.0/pt <= (1.0/3)) {
            ipt = 1;
          } else if ((1.0/5 - 0.01) < 1.0/pt && 1.0/pt <= (1.0/5)) {
            ipt = 2;
          } else if ((1.0/10 - 0.01) < 1.0/pt && 1.0/pt <= (1.0/10)) {
            ipt = 3;
          } else if ((1.0/20 - 0.01) < 1.0/pt && 1.0/pt <= (1.0/20)) {
            ipt = 4;
          } else if ((1.0/50 - 0.005) < 1.0/pt && 1.0/pt <= (1.0/50)) {
            ipt = 5;
          } else if ((1.0/100 - 0.002) < 1.0/pt && 1.0/pt <= (1.0/100)) {
            ipt = 6;
          } else if ((1.0/200 - 0.001) < 1.0/pt && 1.0/pt <= (1.0/200)) {
            ipt = 7;
          }

          // Find dphi
          double dphi = myhit1.Phi_sim() - myhit2.Phi_sim();
          //int dphi_int = myhit1.Phi_fp() - myhit2.Phi_fp();

          // Reverse dphi if positive muon
          if (charge > 0)  dphi = -dphi;

          if (ipt != -1) {
            hname = Form("deflection_gem_csc_st%i_pt%i", j, ipt);
            h = histograms_.at(hname);
            h->Fill(dphi);
          }
        }  // end if found two hits
      }  // end if trigger and good_station
    }  // end loop over station

  }  // end if keep_event

}

// _____________________________________________________________________________
void RPCIntegration::beginJob() {
  bookHistograms();
}

void RPCIntegration::endJob() {
  writeHistograms();
}

// _____________________________________________________________________________
void RPCIntegration::bookHistograms() {
  TString hname;
  TH1F* h;

  // Efficiency vs pT
  // Make [mode] x [eta] where
  //   mode 0,1,2,3 = MuOpen, DoubleMu, SingleMu, Mode15
  //   eta  0,1,2,3 = 1.2-1.6, 1.6-2.0, 2.0-2.4, inclusive

  const Double_t pt_bins[47] = {1., 1.5, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 35., 40., 45., 50., 55., 60., 70., 80., 100., 120., 140., 160., 180., 200., 250., 300., 350., 400., 450., 500., 550., 600., 700., 800., 1000.};

  for (int i=0; i<4; ++i) {  // mode
    for (int j=0; j<4; ++j) {  // eta
      for (int k=0; k<2; ++k) {
        if (k == 0)
          hname = Form("denom_vs_pt_mode%i_eta%i", i, j);
        else
          hname = Form("num_vs_pt_mode%i_eta%i", i, j);
        h = new TH1F(hname, "; gen p_{T} [GeV]; entries", 47-1, pt_bins);
        histograms_[hname] = h;
      }
    }
  }

  // Efficiency vs eta (requiring gen pT >= 20)
  // Make [mode] x [pT] where
  //   mode 0,1,2,3 = MuOpen, DoubleMu, SingleMu, Mode15
  //   pT   0,1,2,3 = >=0, >=12, >=22, >=100

  for (int i=0; i<4; ++i) {  // mode
    for (int j=0; j<4; ++j) {  // pT
      for (int k=0; k<2; ++k) {
        if (k == 0)
          hname = Form("denom_vs_eta_mode%i_l1pt%i", i, j);
        else
          hname = Form("num_vs_eta_mode%i_l1pt%i", i, j);
        h = new TH1F(hname, "; gen |#eta|; entries", 70, 1.1, 2.5);
        histograms_[hname] = h;
      }
    }
  }

  // Station efficiency vs eta (requiring gen pT >= 20)
  //   st 0,1,2,3 = st1, st2, st3, st4

  for (int j=0; j<4; ++j) {  // station
    for (int k=0; k<2; ++k) {
      if (k == 0)
        hname = Form("denom_vs_eta_st%i", j);
      else
        hname = Form("num_vs_eta_st%i", j);
      h = new TH1F(hname, "; gen |#eta|; entries", 70, 1.1, 2.5);
      histograms_[hname] = h;
    }
  }

  // Deflection angles of GEM-CSC and RPC-CSC
  //   st 0,1,2,3 = st1, st2, st3, st4
  //   pt 0..7 = 2, 3, 5, 10, 20, 50, 100, 200

  for (int j=0; j<4; ++j) {  // station
    for (int k=0; k<8; ++k) {  // pT
      hname = Form("deflection_gem_csc_st%i_pt%i", j, k);
      if (j == 0 || j == 1) {
        h = new TH1F(hname, "; GEM #phi - CSC #phi [deg]", 51, -1.54, 0.5);
      } else if (j == 2 || j == 3) {
        h = new TH1F(hname, "; iRPC #phi - CSC #phi [deg]", 51, -1.54, 0.5);
      }
      histograms_[hname] = h;
    }
  }

}

void RPCIntegration::writeHistograms() {
  // TFileService
  edm::Service<TFileService> fs;

  TH1F* h;
  TH2F* h2;

  TString hname;
  TH1F* denom;
  TH1F* num;
  TEfficiency* eff;

  for (const auto& kv : histograms_) {
    // Make TEfficiency
    if (kv.first.BeginsWith("denom_")) {
      hname = kv.first;
      denom = kv.second;
      hname.ReplaceAll("denom_", "num_");
      num = histograms_.at(hname);
      hname = kv.first;
      hname.ReplaceAll("denom_", "eff_");

      eff = fs->make<TEfficiency>(*num, *denom);
      eff->SetName(hname);
      eff->SetStatisticOption(TEfficiency::kFCP);
    }

    h = fs->make<TH1F>(*kv.second);
    if (h) {}
  }

  for (const auto& kv : histogram2Ds_) {
    h2 = fs->make<TH2F>(*kv.second);
    if (h2) {}
  }
}

// _____________________________________________________________________________
void RPCIntegration::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// _____________________________________________________________________________
// Define this as a plug-in
DEFINE_FWK_MODULE(RPCIntegration);
