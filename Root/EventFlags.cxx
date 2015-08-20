#include "susynt-ss2lb/EventFlags.h"

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
     <<" eq2slep: "         <<eq2slep
     <<" mllMin: "          <<mllMin
     <<" l0eta: "           <<l0eta
     <<" l1eta: "           <<l1eta
     <<" tauVeto: "         <<tauVeto
     <<" forwardJetVeto: "  <<forwardJetVeto
     <<" bjetVeto: "        <<bjetVeto
     <<" dileptonTrigger: " <<dileptonTrigger
     <<" oppositeSign: "    <<oppositeSign
     <<" l0pt: "            <<l0pt;
  return oss.str();
}
//-----------------------------------------
std::string EventFlags::str_what_fails() const
{
    std::ostringstream oss;
    oss<<"failing: ";
    if(allTrue()){
        oss<<"--";
    } else {
        oss<<(!grl            ? " grl"             : "")
           <<(!larErr         ? " larErr"          : "")
           <<(!tileErr        ? " tileErr"         : "")
           <<(!ttcVeto        ? " ttcVeto"         : "")
           <<(!goodVtx        ? " goodVtx"         : "")
           <<(!tileTrip       ? " tileTrip"        : "")
           <<(!lAr            ? " lAr"             : "")
           <<(!badJet         ? " badJet"          : "")
           <<(!deadRegions    ? " deadRegions"     : "")
           <<(!badMuon        ? " badMuon"         : "")
           <<(!cosmicMuon     ? " cosmicMuon"      : "")
           <<(!hfor           ? " hfor"            : "")
           <<(!ge2blep        ? " ge2blep"         : "")
           <<(!eq2blep        ? " eq2blep"         : "")
           <<(!eq2slep        ? " eq2slep"         : "")
           <<(!mllMin         ? " mllMin"          : "")
           <<(!l0eta          ? " l0eta"           : "")
           <<(!l1eta          ? " l1eta"           : "")
           <<(!tauVeto        ? " tauVeto"         : "")
           <<(!forwardJetVeto ? " forwardJetVeto"  : "")
           <<(!bjetVeto       ? " bjetVeto"        : "")
           <<(!dileptonTrigger? " dileptonTrigger" : "")
           <<(!oppositeSign   ? " oppositeSign"    : "")
           <<(!l0pt           ? " l0pt"            : "");
    }
    return oss.str();
}
//-----------------------------------------
