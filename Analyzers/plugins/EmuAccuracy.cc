#include <memory>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"


// _____________________________________________________________________________
class EmuAccuracy : public edm::EDAnalyzer {
public:
  explicit EmuAccuracy(const edm::ParameterSet& iConfig);
  ~EmuAccuracy();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  virtual void endJob() override;

  // Member functions
  void getHandles(const edm::Event& iEvent);

  void findMatches();

  void printHits();
  void printTracks();

  void sitrep(const std::vector<int>& unp_matches, const std::vector<int>& emu_matches);

  void printSitrep();

  // Member data
  const edm::InputTag unpHitTag_;
  const edm::InputTag emuHitTag_;
  const edm::InputTag unpTrackTag_;
  const edm::InputTag emuTrackTag_;

  int verbose_;

  edm::EDGetTokenT<l1t::EMTFHitCollection>   unpHitToken_;
  edm::EDGetTokenT<l1t::EMTFHitCollection>   emuHitToken_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection> unpTrackToken_;
  edm::EDGetTokenT<l1t::EMTFTrackCollection> emuTrackToken_;

  l1t::EMTFHitCollection    unpHits_;
  l1t::EMTFHitCollection    emuHits_;
  l1t::EMTFTrackCollection  unpTracks_;
  l1t::EMTFTrackCollection  emuTracks_;

  std::map<int, int> sitrep_why_;
  std::map<int, int> sitrep_why_address_;

  unsigned nBadEvents_;
  unsigned nEvents_;
};

// _____________________________________________________________________________
EmuAccuracy::EmuAccuracy(const edm::ParameterSet& iConfig) :
    unpHitTag_   (iConfig.getParameter<edm::InputTag>("unpHitTag")),
    emuHitTag_   (iConfig.getParameter<edm::InputTag>("emuHitTag")),
    unpTrackTag_ (iConfig.getParameter<edm::InputTag>("unpTrackTag")),
    emuTrackTag_ (iConfig.getParameter<edm::InputTag>("emuTrackTag")),
    verbose_     (iConfig.getUntrackedParameter<int> ("verbosity"))
{
    unpHitToken_   = consumes<l1t::EMTFHitCollection>  (unpHitTag_);
    emuHitToken_   = consumes<l1t::EMTFHitCollection>  (emuHitTag_);
    unpTrackToken_ = consumes<l1t::EMTFTrackCollection>(unpTrackTag_);
    emuTrackToken_ = consumes<l1t::EMTFTrackCollection>(emuTrackTag_);

    nBadEvents_ = 0;
    nEvents_    = 0;
}

EmuAccuracy::~EmuAccuracy() {}

// _____________________________________________________________________________
struct TracksMatch {
  template<class T1, class T2>
  bool operator()(const T1& trk1, const T2& trk2) const {
    bool m = (
        std::make_tuple(
          trk1.BX(),
          trk1.Endcap(),
          trk1.Sector(),
          trk1.PtLUT().address,
          //trk1.GMT_pt(),
          trk1.Mode(),
          trk1.GMT_eta(),
          trk1.GMT_phi(),
          trk1.GMT_charge(),
          trk1.GMT_quality()
        ) == std::make_tuple(
          trk2.BX(),
          trk2.Endcap(),
          trk2.Sector(),
          trk2.PtLUT().address,
          //trk2.GMT_pt(),
          trk2.Mode(),
          trk2.GMT_eta(),
          trk2.GMT_phi(),
          trk2.GMT_charge(),
          trk2.GMT_quality()
        )
    );
    return m;
  }
};

struct TrackPrint {
  template<class T>
  void operator()(const T& trk) const {
    std::cout
        << "bx: "           << trk.BX()
        << " endcap: "      << trk.Endcap()
        << " sector: "      << trk.Sector()
        << " ptlut_addr: "  << trk.PtLUT().address
        << " gmt_pt: "      << trk.GMT_pt()
        << " mode: "        << trk.Mode()
        << " eta: "         << trk.GMT_eta()
        << " phi: "         << trk.GMT_phi()
        << " charge: "      << trk.GMT_charge()
        << " quality: "     << trk.GMT_quality()
        << std::endl;
    return;
  }
};

struct HitPrint {
  template<class T>
  void operator()(const T& h) const {
    int bx      = h.BX() + 3;
    int endcap  = (h.Endcap() == 1) ? 1 : 2;
    int sector  = h.PC_sector();
    int station = (h.PC_station() == 0 && h.Subsector() == 1) ? 1 : h.PC_station();
    int chamber = h.PC_chamber() + 1;
    int strip   = (h.Station() == 1 && h.Ring() == 4) ? h.Strip() + 128 : h.Strip();  // ME1/1a
    std::cout << bx << " " << endcap << " " << sector << " " << h.Subsector() << " " << station << " " << h.Valid() << " " << h.Quality() << " " << h.Pattern() << " " << h.Wire() << " " << chamber << " " << h.Bend() << " " << strip << std::endl;
    return;
  }
};

// See http://stackoverflow.com/a/21510185
namespace details {
  template <class T> struct _reversed {
    T& t; _reversed(T& _t): t(_t) {}
    decltype(t.rbegin()) begin() { return t.rbegin(); }
    decltype(t.rend()) end() { return t.rend(); }
  };
}
template <class T> details::_reversed<T> reversed(T& t) { return details::_reversed<T>(t); }


// _____________________________________________________________________________
void EmuAccuracy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  getHandles(iEvent);

  findMatches();
}

// _____________________________________________________________________________
void EmuAccuracy::getHandles(const edm::Event& iEvent) {

  edm::Handle<decltype(unpHits_)>    unpHits_handle;
  edm::Handle<decltype(emuHits_)>    emuHits_handle;
  edm::Handle<decltype(unpTracks_)>  unpTracks_handle;
  edm::Handle<decltype(emuTracks_)>  emuTracks_handle;

  if (iEvent.isRealData()) {
    if (!unpHitToken_.isUninitialized()) {
      iEvent.getByToken(unpHitToken_, unpHits_handle);
    }
    if (!unpHits_handle.isValid()) {
      edm::LogError("EmuAccuracy") << "Cannot get the product: " << unpHitTag_;
      return;
    }

    if (!unpTrackToken_.isUninitialized()) {
      iEvent.getByToken(unpTrackToken_, unpTracks_handle);
    }
    if (!unpTracks_handle.isValid()) {
      edm::LogError("EmuAccuracy") << "Cannot get the product: " << unpTrackTag_;
      return;
    }
  }

  if (!emuHitToken_.isUninitialized()) {
    iEvent.getByToken(emuHitToken_, emuHits_handle);
  }
  if (!emuHits_handle.isValid()) {
    edm::LogError("EmuAccuracy") << "Cannot get the product: " << emuHitTag_;
    return;
  }

  if (!emuTrackToken_.isUninitialized()) {
    iEvent.getByToken(emuTrackToken_, emuTracks_handle);
  }
  if (!emuTracks_handle.isValid()) {
    edm::LogError("EmuAccuracy") << "Cannot get the product: " << emuTrackTag_;
    return;
  }

  // Object filters

  // Exclude out-of-emu BX
  auto out_of_emu_bx = [](const auto& trk) {
    bool out = (trk.BX() < -1) || (+1 < trk.BX());
    //bool out = (trk.BX() != 0);
    return out;
  };

  // Exclude sector -6
  auto out_of_emu_sector = [](const auto& trk) {
    bool out = (trk.Endcap() == -1 && trk.Sector() == 6);
    return out;
  };

  unpHits_.clear();
  for (const auto& hit : (*unpHits_handle)) {
    unpHits_.push_back(hit);
  }

  emuHits_.clear();
  for (const auto& hit : (*emuHits_handle)) {
    emuHits_.push_back(hit);
  }

  unpTracks_.clear();
  for (const auto& trk : reversed(*unpTracks_handle)) {
    if (out_of_emu_bx(trk))      continue;
    if (out_of_emu_sector(trk))  continue;
    unpTracks_.push_back(trk);
  }

  emuTracks_.clear();
  for (const auto& trk : (*emuTracks_handle)) {
    if (out_of_emu_bx(trk))      continue;
    if (out_of_emu_sector(trk))  continue;
    emuTracks_.push_back(trk);
  }

}

// _____________________________________________________________________________
void EmuAccuracy::findMatches() {
  TracksMatch tracksMatch;

  // Check that unpacker tracks have matching emulator tracks
  std::vector<int> unp_matches(unpTracks_.size(), -99);
  bool unp_has_match = true;

  // Unpacker tracks
  if (!unpTracks_.empty()) {
    int i = 0;
    for (const auto& trk1 : unpTracks_) {
      // Compare emulator tracks
      int j = 0;
      for (const auto& trk2 : emuTracks_) {
        bool m = tracksMatch(trk1, trk2);
        if (m) {
          unp_matches.at(i) = j;
        }
        j++;
      }

      if (unp_matches.at(i) == -99) {
        unp_has_match = false;
      }
      i++;
    }
  }

  // Check that emulator tracks have matching unpacker tracks
  std::vector<int> emu_matches(emuTracks_.size(), -99);
  bool emu_has_match = true;

  // Emulator tracks
  if (!emuTracks_.empty()) {
    int i = 0;
    for (const auto& trk1 : emuTracks_) {
      // Compare unpacker tracks
      int j = 0;
      for (const auto& trk2 : unpTracks_) {
        bool m = tracksMatch(trk1, trk2);
        if (m) {
          emu_matches.at(i) = j;
        }
        j++;
      }

      if (emu_matches.at(i) == -99) {
        emu_has_match = false;
      }
      i++;
    }
  }

  bool no_match = (!unp_has_match || !emu_has_match);
  if (verbose_ || no_match) {  // if verbose, always print
    printHits();

    if (!unp_has_match) {
      std::cout << "ERROR: Unpacker tracks have no matching emulator tracks." << std::endl;
      for (unsigned i = 0; i < unp_matches.size(); ++i) {
        int j = unp_matches.at(i);
        std::cout << "  unp " << i << " --> emu " << j << std::endl;
      }
    }

    if (!emu_has_match) {
      std::cout << "ERROR: Emulator tracks have no matching unpacker tracks." << std::endl;
      for (unsigned i = 0; i < emu_matches.size(); ++i) {
        int j = emu_matches.at(i);
        std::cout << "  emu " << i << " --> unp " << j << std::endl;
      }
    }

    printTracks();
  }

  if (no_match) {
    nBadEvents_++;
  }
  nEvents_++;

  sitrep(unp_matches, emu_matches);
}

// _____________________________________________________________________________
void EmuAccuracy::printHits() {
  HitPrint hitPrint;

  // Emulator hits
  std::cout << "Num of emulator hits = " << emuHits_.size() << std::endl;
  std::cout << "bx e s ss st vf ql cp wg id bd hs" << std::endl;
  for (const auto& hit : emuHits_) {
    hitPrint(hit);
  }

  bool printUnpacker = false;
  if (printUnpacker) {
    // Unpacker hits
    std::cout << "Num of unpacker hits = " << unpHits_.size() << std::endl;
    std::cout << "bx e s ss st" << std::endl;
    for (const auto& hit : unpHits_) {
      std::cout << hit.BX() + 6 - 3 << " " << hit.Endcap() << " " << hit.Sector() << " " << hit.Subsector() << " " << hit.Station() << std::endl;
    }
  }
}

void EmuAccuracy::printTracks() {
  TrackPrint trackPrint;
  int i = 0;

  // Unpacker tracks
  std::cout << "Num of unpacker tracks = " << unpTracks_.size() << std::endl;
  i = 0;
  for (const auto& trk : unpTracks_) {
    std::cout << i++ << " ";
    trackPrint(trk);
  }

  // Emulator tracks
  std::cout << "Num of emulator tracks = " << emuTracks_.size() << std::endl;
  i = 0;
  for (const auto& trk : emuTracks_) {
    std::cout << i++ << " ";
    trackPrint(trk);
  }
}

// _____________________________________________________________________________
void EmuAccuracy::sitrep(const std::vector<int>& unp_matches, const std::vector<int>& emu_matches) {

  auto check_track_bx = [](const auto& trk1, const auto& trk2) {
    bool match = (
        std::make_tuple(trk1.BX(), trk1.Endcap(), trk1.Sector()) ==
        std::make_tuple(trk2.BX(), trk2.Endcap(), trk2.Sector())
    );
    return match;
  };

  auto check_track_mode = [](const auto& trk1, const auto& trk2) {
    bool match = (trk1.Mode() == trk2.Mode() && trk1.Mode_inv() == trk2.Mode_inv());
    return match;
  };

  auto check_track_address = [](const auto& trk1, const auto& trk2) {
    bool match = (trk1.PtLUT().address == trk2.PtLUT().address);
    return match;
  };

  auto check_track_pt = [](const auto& trk1, const auto& trk2) {
    bool match = (trk1.GMT_pt() == trk2.GMT_pt());
    return match;
  };

  auto check_track_eta = [](const auto& trk1, const auto& trk2) {
    bool match = (trk1.GMT_eta() == trk2.GMT_eta());
    return match;
  };

  auto check_track_phi = [](const auto& trk1, const auto& trk2) {
    bool match = (trk1.GMT_phi() == trk2.GMT_phi());
    return match;
  };

  auto check_track_charge = [](const auto& trk1, const auto& trk2) {
    bool match = (trk1.GMT_charge() == trk2.GMT_charge() && trk1.GMT_charge_valid() == trk2.GMT_charge_valid());
    return match;
  };

  auto check_track_quality = [](const auto& trk1, const auto& trk2) {
    bool match = (trk1.GMT_quality() == trk2.GMT_quality());
    return match;
  };

  auto check_me11_dupes = [](const auto& emuHits_) {
    typedef std::array<int, 5> reference_t;
    std::map<reference_t, int> counting;

    for (const auto& hit : emuHits_) {
      if (hit.Station() == 1 && (hit.Ring() == 1 || hit.Ring() == 4)) {  // ME1/1
        reference_t ref = {{hit.Endcap(), hit.Sector(), hit.Subsector(), hit.Station(), hit.CSC_ID()}};
        counting[ref]++;
      }
    }

    bool found = false;
    for (const auto& kv : counting)
      if (kv.second > 1)
        found = true;
    return found;
  };

  auto check_prim_match = [](const auto& emuHits_, const auto& trk) {
    std::array<int, 4> minima;
    std::array<int, 4> next_minima;

    minima.fill(999999);
    next_minima.fill(999999);

    const int bw_fph = 13;
    const int bpow = 7;
    int ph_pat = trk.Ph_num();

    for (const auto& hit : emuHits_) {
      if (
          (hit.Endcap() == trk.Endcap()) &&
          (hit.Sector() == trk.Sector()) &&
          (hit.Zone_code() & (1<<(trk.Zone()-1)))
      ) {
        int ph_seg = hit.Phi_fp();
        int ph_seg_red = ph_seg >> (bw_fph-bpow-1);
        int ph_diff = (ph_pat > ph_seg_red) ? (ph_pat - ph_seg_red) : (ph_seg_red - ph_pat);

        int istation = (hit.Station()-1);

        if (minima.at(istation) > ph_diff) {
          minima.at(istation) = ph_diff;
        } else if (next_minima.at(istation) > ph_diff) {
          next_minima.at(istation) = ph_diff;
        }
      }
    }

    bool found = false;
    for (int istation = 0; istation < 4; ++istation) {
      if (
        (minima.at(istation) != 999999) &&
        (minima.at(istation) == next_minima.at(istation))
      ) {
        found = true;
      }
    }
    return found;
  };


  bool empty           = false;
  bool mis_unp_ntracks = false;
  bool mis_emu_ntracks = false;
  bool mis_unp_2       = false;
  bool mis_emu_2       = false;
  bool mis_bx          = false;
  bool mis_mode        = false;
  bool mis_address     = false;
  bool mis_eta         = false;
  bool mis_phi         = false;
  bool mis_charge      = false;
  bool mis_quality     = false;
  bool mis_me11_dupes  = false;
  bool mis_prim_match  = false;

  empty = (unpTracks_.empty() && emuTracks_.empty());

  int cnt_mis_unp = (std::count(unp_matches.begin(), unp_matches.end(), -99));
  int cnt_mis_emu = (std::count(emu_matches.begin(), emu_matches.end(), -99));

  mis_unp_ntracks = cnt_mis_unp > cnt_mis_emu;
  mis_emu_ntracks = cnt_mis_unp < cnt_mis_emu;
  if (!mis_unp_ntracks && !mis_emu_ntracks) {
    mis_unp_2 = cnt_mis_unp > 1;
    mis_emu_2 = cnt_mis_emu > 1;

    if (cnt_mis_unp == 1)
      assert(cnt_mis_emu == 1);

    if (cnt_mis_unp == 1) {
      auto found_unp = std::find(unp_matches.begin(), unp_matches.end(), -99);
      auto found_emu = std::find(emu_matches.begin(), emu_matches.end(), -99);

      auto index_unp = std::distance(unp_matches.begin(), found_unp);
      auto index_emu = std::distance(emu_matches.begin(), found_emu);

      const auto& trk1 = unpTracks_.at(index_unp);
      const auto& trk2 = emuTracks_.at(index_emu);

      mis_bx = !check_track_bx(trk1, trk2);

      if (!mis_bx) {
        mis_mode    = !check_track_mode(trk1, trk2);
        mis_address = !check_track_address(trk1, trk2);
        mis_eta     = !check_track_eta(trk1, trk2);
        mis_phi     = !check_track_phi(trk1, trk2);
        mis_charge  = !check_track_charge(trk1, trk2);
        mis_quality = !check_track_quality(trk1, trk2);

        const auto& trk = emuTracks_.at(index_emu);
        mis_me11_dupes = check_me11_dupes(emuHits_);
        mis_prim_match = check_prim_match(emuHits_, trk);
      }
    }
  }

  int why = 0;
  if (empty) {
    why = 1;
  } else if (mis_unp_ntracks) {
    why = 2;
  } else if (mis_emu_ntracks) {
    why = 3;
  } else if (mis_unp_2 || mis_emu_2) {
    why = 4;
  } else if (mis_bx) {
    why = 5;
  } else if (mis_mode) {
    why = 6;
  } else if (mis_address) {
    why = 7;
  } else if (mis_eta) {
    why = 8;
  } else if (mis_phi) {
    why = 9;
  } else if (mis_charge) {
    why = 10;
  } else if (mis_quality) {
    why = 11;
  }

  int why_address = 0;
  if (!mis_bx && !mis_mode && mis_address) {
    if (mis_me11_dupes) {
      why_address = 1;
    } else if (mis_prim_match) {
      why_address = 2;
    }
  }

#if 0
  if (why == 0) {
    // Compare pT
    auto approx_equal = [] (float a, float b) { return std::abs(a-b) < 0.5+1e-3; };

    int i = 0;
    for (const auto& trk1 : unpTracks_) {
      if (unp_matches.at(i) >= 0) {
        const auto& trk2 = emuTracks_.at(unp_matches.at(i));
        if (!approx_equal(trk1.Pt(), trk2.Pt())) {
          // See DataFormats/L1TMuon/interface/EMTFTrack.h
          // See DataFormats/L1TMuon/interface/EMTF/SP.h
          std::cout << i << " Address  : " << trk1.Pt_LUT_addr()  << " vs " << trk2.Pt_LUT_addr()  << std::endl;
          std::cout << i << " Pt       : " << trk1.Pt()           << " vs " << trk2.Pt()           << std::endl;
          std::cout << i << " Pt_GMT   : " << trk1.Pt_GMT()       << " vs " << trk2.Pt_GMT()       << std::endl;
          std::cout << i << " Eta      : " << trk1.Eta()          << " vs " << trk2.Eta()          << std::endl;
          std::cout << i << " Eta_GMT  : " << trk1.Eta_GMT()      << " vs " << trk2.Eta_GMT()      << std::endl;
          std::cout << i << " Eta_LUT  : " << trk1.Eta_LUT()      << " vs " << trk2.Eta_LUT()      << std::endl;
          std::cout << i << " Phi_loc  : " << trk1.Phi_loc_int()  << " vs " << trk2.Phi_loc_int()  << std::endl;
          std::cout << i << " Phi_glob : " << trk1.Phi_glob_rad() << " vs " << trk2.Phi_glob_rad() << std::endl;
          std::cout << i << " Phi_GMT  : " << trk1.Phi_GMT()      << " vs " << trk2.Phi_GMT()      << std::endl;
          std::cout << i << " Mode_LUT : " << trk1.Mode_LUT()     << " vs " << trk2.Mode_LUT()     << std::endl;
          std::cout << i << " Quality  : " << trk1.Quality()      << " vs " << trk2.Quality()      << std::endl;
          std::cout << i << " Charge   : " << trk1.Charge()       << " vs " << trk2.Charge()       << std::endl;
          std::cout << i << " DPhi_12  : " << trk1.DPhi_12()      << " vs " << trk2.DPhi_12()      << std::endl;
          std::cout << i << " DPhi_23  : " << trk1.DPhi_23()      << " vs " << trk2.DPhi_23()      << std::endl;
          std::cout << i << " DTheta_12: " << trk1.DTheta_12()    << " vs " << trk2.DTheta_12()    << std::endl;
          std::cout << i << " DTheta_23: " << trk1.DTheta_23()    << " vs " << trk2.DTheta_23()    << std::endl;
          std::cout << i << " CLCT_1   : " << trk1.CLCT_1()       << " vs " << trk2.CLCT_1()       << std::endl;
          std::cout << i << " CLCT_2   : " << trk1.CLCT_2()       << " vs " << trk2.CLCT_2()       << std::endl;
          std::cout << i << " FR_1     : " << trk1.FR_1()         << " vs " << trk2.FR_1()         << std::endl;
          std::cout << i << " FR_2     : " << trk1.FR_2()         << " vs " << trk2.FR_2()         << std::endl;
        }
      }
      i += 1;
    }
  }
#endif

  sitrep_why_[why]++;
  sitrep_why_address_[why_address]++;

  if (why >= 2) {
    std::cout << ">> why        : " << why << std::endl;
    std::cout << ">> why_address: " << why_address << std::endl;
  }
}

void EmuAccuracy::printSitrep() {
  std::cout << "**************** SITREP ****************" << std::endl;
  for (const auto& kv : sitrep_why_) {
    std::cout << "why        : " << kv.first << " count: " << kv.second << std::endl;
  }
  for (const auto& kv : sitrep_why_address_) {
    std::cout << "why_address: " << kv.first << " count: " << kv.second << std::endl;
  }
}

// _____________________________________________________________________________
void EmuAccuracy::beginJob() {}

void EmuAccuracy::endJob() {
  std::cout << "Num of bad events: " << nBadEvents_ << "/" << nEvents_ << std::endl;
  printSitrep();
}

void EmuAccuracy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// _____________________________________________________________________________
// Define this as a plug-in
DEFINE_FWK_MODULE(EmuAccuracy);
