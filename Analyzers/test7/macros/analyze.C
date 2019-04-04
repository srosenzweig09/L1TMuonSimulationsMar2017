#include "TreeReader.h"

#include <iostream>

TString filename = "ntuple_SingleMuon_Endcap_2GeV_add.5.root";

Long64_t maxEvents = 10;

// _____________________________________________________________________________
// To execute: root -l -b -q analyze.C+

void analyze() {
  // Open the file
  std::cout << "Opening file: " << filename << std::endl;

  TreeReader reader;
  reader.init(filename);

  // Get number of events
  std::cout << "Processing no of events: " << maxEvents << std::endl;
  Long64_t nentries = reader.getEntries();

  // Loop over events
  for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
    if (maxEvents != -1 && ievt == maxEvents)
      break;

    // Retrieve the event
    reader.getEntry(ievt);

    // Get the particle variables.
    // See TreeReader.h for more info.
    float   pt      = reader.vp_pt      ->front();
    float   eta     = reader.vp_eta     ->front();
    float   phi     = reader.vp_phi     ->front();
    int16_t q       = reader.vp_q       ->front();

    std::cout << "ievt: "     << ievt
              << " muon pt: " << pt
              << " q: "       << q
              << " phi: "     << phi
              << " eta: "     << eta
              << "\n";

    // Loop over hits
    unsigned nhits = reader.vh_endcap->size();

    for (unsigned ihit = 0; ihit < nhits; ++ihit) {
      // Get the hit variables.
      // See TreeReader.h for more info.
      int16_t hit_type      = reader.vh_type      ->at(ihit);
      int16_t hit_endcap    = reader.vh_endcap    ->at(ihit);
      int16_t hit_station   = reader.vh_station   ->at(ihit);
      int16_t hit_ring      = reader.vh_ring      ->at(ihit);
      int16_t hit_sector    = reader.vh_sector    ->at(ihit);
      int16_t hit_subsector = reader.vh_subsector ->at(ihit);
      int16_t hit_cscid     = reader.vh_cscid     ->at(ihit);
      int16_t hit_neighbor  = reader.vh_neighbor  ->at(ihit);
      int16_t hit_bx        = reader.vh_bx        ->at(ihit);
      int16_t hit_strip     = reader.vh_strip     ->at(ihit);
      int16_t hit_wire      = reader.vh_wire      ->at(ihit);
      int16_t hit_roll      = reader.vh_roll      ->at(ihit);
      float   hit_sim_phi   = reader.vh_sim_phi   ->at(ihit);
      float   hit_sim_theta = reader.vh_sim_theta ->at(ihit);

      std::cout << ".. ihit: "    << ihit
                << " type: "      << hit_type
                << " endcap: "    << hit_endcap
                << " station: "   << hit_station
                << " ring: "      << hit_ring
                << " sector: "    << hit_sector
                << " subsector: " << hit_subsector
                << " cscid: "     << hit_cscid
                << " neigh: "     << hit_neighbor
                << "\n  "
                << " bx: "        << hit_bx
                << " strip: "     << hit_strip
                << " wire: "      << hit_wire
                << " roll: "      << hit_roll
                << " phi: "       << hit_sim_phi
                << " theta: "     << hit_sim_theta
                << "\n";
    }  // end loop over hits

    std::cout << std::endl;

  }  // end loop over events

}  // end analyze()


#ifndef __CINT__
// Define main function if necessary
int main()
{
  analyze();
  return 0;
}
#endif
