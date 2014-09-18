#include "SusyntHlfv/EventFlags.h"

#include <sstream>      // std::ostringstream

using hlfv::EventFlags;

//-----------------------------------------
std::string EventFlags::str() const
{
  std::ostringstream oss;
  oss<<" grl: "             <<grl
     <<" larErr: "          <<larErr
     <<" tileErr: "         <<tileErr
     <<" ttcVeto: "         <<ttcVeto
     <<" goodVtx: "         <<goodVtx
     <<" tileTrip: "        <<tileTrip
     <<" lAr: "             <<lAr
     <<" badJet: "          <<badJet
     <<" deadRegions: "     <<deadRegions
     <<" badMuon: "         <<badMuon
     <<" cosmicMuon: "      <<cosmicMuon
     <<" hfor: "            <<hfor
     <<" ge2blep: "         <<ge2blep
     <<" eq2blep: "         <<eq2blep
     <<" mllMin: "          <<mllMin
     <<" l0eta: "           <<l0eta
     <<" l1eta: "           <<l1eta
     <<" tauVeto: "         <<tauVeto
     <<" dileptonTrigger: " <<dileptonTrigger
     <<" oppositeSign: "    <<oppositeSign
     <<" l0pt: "            <<l0pt;
  return oss.str();
}
//-----------------------------------------
