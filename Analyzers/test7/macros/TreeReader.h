#ifndef TREEREADER_H_
#define TREEREADER_H_

#include "TChain.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include <cstdint>
#include <memory>
#include <vector>

// _____________________________________________________________________________
// This is a simple wrapper around TChain. It sets the branch names, addresses
// etc. Its functions are essentially the same as the functions of TChain.

class TreeReader {
public:
  TreeReader(int verbose=1);
  ~TreeReader();

  void init(TString src);

  Long64_t loadTree(Long64_t entry) { return tchain_->LoadTree(entry); }

  Int_t getEntry(Long64_t entry) { return tchain_->GetEntry(entry); }

  Long64_t getEntries() const { return tchain_->GetEntries(); }

  //TChain* getChain() { return tchain_; }

  // Particle variables
  std::vector<float> *          vp_pt;         // particle pT [GeV]
  std::vector<float> *          vp_eta;        // particle eta
  std::vector<float> *          vp_phi;        // particle phi [rad]
  std::vector<float> *          vp_vx;         // particle vertex, x [cm]
  std::vector<float> *          vp_vy;         // particle vertex, y [cm]
  std::vector<float> *          vp_vz;         // particle vertex, z [cm]
  std::vector<int16_t> *        vp_q;          // particle charge
  std::vector<int16_t> *        vp_bx;         // particle BX
  std::vector<int32_t> *        vp_pdgid;      // particle PDG id

  // Hit variables
  std::vector<int16_t> *        vh_endcap;     // hit detector id: endcap
  std::vector<int16_t> *        vh_station;    // hit detector id: station
  std::vector<int16_t> *        vh_ring;       // hit detector id: ring
  std::vector<int16_t> *        vh_sector;     // hit detector id: sector
  std::vector<int16_t> *        vh_subsector;  // hit detector id: subsector
  std::vector<int16_t> *        vh_chamber;    // hit detector id: chamber
  std::vector<int16_t> *        vh_cscid;      // hit detector id: CSC id
  std::vector<int16_t> *        vh_bx;         // hit BX
  std::vector<int16_t> *        vh_type;       // hit subsystem type: DT=0,CSC=1,RPC=2,GEM=3,ME0=4
  std::vector<int16_t> *        vh_neighbor;   // hit from neighbor sector: false=0, true=1

  std::vector<int16_t> *        vh_strip;      // hit strip
  std::vector<int16_t> *        vh_wire;       // hit wire
  std::vector<int16_t> *        vh_roll;       // hit roll
  std::vector<int16_t> *        vh_quality;    // hit quality
  std::vector<int16_t> *        vh_pattern;    // hit pattern
  std::vector<int16_t> *        vh_bend;       // hit bend
  std::vector<int16_t> *        vh_time;       // hit time
  std::vector<int16_t> *        vh_fr;         // hit F/R: rear=0, front=1
  std::vector<int32_t> *        vh_emtf_phi;   // hit phi in integer EMTF unit
  std::vector<int32_t> *        vh_emtf_theta; // hit theta in integer EMTF unit

  std::vector<float> *          vh_sim_phi;    // hit position in global coordinate: phi [deg]
  std::vector<float> *          vh_sim_theta;  // hit position in global coordinate: theta [deg]
  std::vector<float> *          vh_sim_eta;    // hit position in global coordinate: eta
  std::vector<float> *          vh_sim_r;      // hit position in global coordinate: rho [cm]
  std::vector<float> *          vh_sim_z;      // hit position in global coordinate: z [cm]
  std::vector<int32_t> *        vh_sim_tp1;    // hit MC truth info [ignore this]
  std::vector<int32_t> *        vh_sim_tp2;    // hit MC truth info [ignore this]

protected:
  std::unique_ptr<TChain> tchain_;
  int treenumber_;
  const int verbose_;
};
#endif  // TREEREADER_H_

// _____________________________________________________________________________
// Implementation is included in the header file to simplify ROOT library generation
#define TREEREADER_CXX_
#ifdef  TREEREADER_CXX_

#include <cassert>
#include <stdexcept>

TreeReader::TreeReader(int verbose) :
  vp_pt               (0),
  vp_eta              (0),
  vp_phi              (0),
  vp_vx               (0),
  vp_vy               (0),
  vp_vz               (0),
  vp_q                (0),
  vp_bx               (0),
  vp_pdgid            (0),
  //
  vh_endcap           (0),
  vh_station          (0),
  vh_ring             (0),
  vh_sector           (0),
  vh_subsector        (0),
  vh_chamber          (0),
  vh_cscid            (0),
  vh_bx               (0),
  vh_type             (0),
  vh_neighbor         (0),
  //
  vh_strip            (0),
  vh_wire             (0),
  vh_roll             (0),
  vh_quality          (0),
  vh_pattern          (0),
  vh_bend             (0),
  vh_time             (0),
  vh_fr               (0),
  vh_emtf_phi         (0),
  vh_emtf_theta       (0),
  //
  vh_sim_phi          (0),
  vh_sim_theta        (0),
  vh_sim_eta          (0),
  vh_sim_r            (0),
  vh_sim_z            (0),
  vh_sim_tp1          (0),
  vh_sim_tp2          (0),
  //
  tchain_(),
  treenumber_(0),
  verbose_(verbose) {}

TreeReader::~TreeReader() {}

void TreeReader::init(TString src) {
  if (!src.EndsWith(".root") && !src.EndsWith(".txt")) {
    TString msg = "Input source must be either .root or .txt";
    throw std::invalid_argument(msg.Data());
  }

  //if (verbose_)  std::cout << "Opening " << src << std::endl;
  tchain_.reset(new TChain("ntupler/tree"));

  if (src.EndsWith(".root")) {
    if (tchain_->Add(src) ) {
      //if (verbose_)  std::cout << "Successfully read " << src << std::endl;
    } else {
      TString msg = "Failed to read " + src;
      throw std::invalid_argument(msg.Data());
    }
  } else if (src.EndsWith(".txt")) {
    TFileCollection fc("fileinfolist", "", src);
    if (tchain_->AddFileInfoList((TCollection*) fc.GetList()) ) {
      //if (verbose_)  std::cout << "Successfully read " << src << std::endl;
    } else {
      TString msg = "Failed to read " + src;
      throw std::invalid_argument(msg.Data());
    }
  }

  assert(tchain_ != 0);
  treenumber_ = tchain_->GetTreeNumber();

  tchain_->SetBranchAddress("vp_pt"         , &vp_pt        );
  tchain_->SetBranchAddress("vp_eta"        , &vp_eta       );
  tchain_->SetBranchAddress("vp_phi"        , &vp_phi       );
  tchain_->SetBranchAddress("vp_vx"         , &vp_vx        );
  tchain_->SetBranchAddress("vp_vy"         , &vp_vy        );
  tchain_->SetBranchAddress("vp_vz"         , &vp_vz        );
  tchain_->SetBranchAddress("vp_q"          , &vp_q         );
  tchain_->SetBranchAddress("vp_bx"         , &vp_bx        );
  tchain_->SetBranchAddress("vp_pdgid"      , &vp_pdgid     );
  //
  tchain_->SetBranchAddress("vh_endcap"     , &vh_endcap    );
  tchain_->SetBranchAddress("vh_station"    , &vh_station   );
  tchain_->SetBranchAddress("vh_ring"       , &vh_ring      );
  tchain_->SetBranchAddress("vh_sector"     , &vh_sector    );
  tchain_->SetBranchAddress("vh_subsector"  , &vh_subsector );
  tchain_->SetBranchAddress("vh_chamber"    , &vh_chamber   );
  tchain_->SetBranchAddress("vh_cscid"      , &vh_cscid     );
  tchain_->SetBranchAddress("vh_bx"         , &vh_bx        );
  tchain_->SetBranchAddress("vh_type"       , &vh_type      );
  tchain_->SetBranchAddress("vh_neighbor"   , &vh_neighbor  );
  //
  tchain_->SetBranchAddress("vh_strip"      , &vh_strip     );
  tchain_->SetBranchAddress("vh_wire"       , &vh_wire      );
  tchain_->SetBranchAddress("vh_roll"       , &vh_roll      );
  tchain_->SetBranchAddress("vh_quality"    , &vh_quality   );
  tchain_->SetBranchAddress("vh_pattern"    , &vh_pattern   );
  tchain_->SetBranchAddress("vh_bend"       , &vh_bend      );
  tchain_->SetBranchAddress("vh_time"       , &vh_time      );
  tchain_->SetBranchAddress("vh_fr"         , &vh_fr        );
  tchain_->SetBranchAddress("vh_emtf_phi"   , &vh_emtf_phi  );
  tchain_->SetBranchAddress("vh_emtf_theta" , &vh_emtf_theta);
  //
  tchain_->SetBranchAddress("vh_sim_phi"    , &vh_sim_phi   );
  tchain_->SetBranchAddress("vh_sim_theta"  , &vh_sim_theta );
  tchain_->SetBranchAddress("vh_sim_eta"    , &vh_sim_eta   );
  tchain_->SetBranchAddress("vh_sim_r"      , &vh_sim_r     );
  tchain_->SetBranchAddress("vh_sim_z"      , &vh_sim_z     );
  tchain_->SetBranchAddress("vh_sim_tp1"    , &vh_sim_tp1   );
  tchain_->SetBranchAddress("vh_sim_tp2"    , &vh_sim_tp2   );
}
#endif  // TREEREADER_CXX_
