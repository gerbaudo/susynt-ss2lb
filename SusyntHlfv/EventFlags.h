
// emacs -*- C++ -*-
#ifndef HLFV_EVENTFLAGS_H
#define HLFV_EVENTFLAGS_H

#include <string>

namespace hlfv
{
///  A holder to pass around the event-level flags; true = pass (including vetoes)
/**
   Imported from github.com/gerbaudo/SusyTest0/SusyTest0/
   davide.gerbaudo@gmail.com
   Feb 2014
*/
struct EventFlags {
  EventFlags() { reset(); }
  void reset() {
   grl = larErr = tileErr = ttcVeto = goodVtx = tileTrip = lAr = false;
   badJet = deadRegions = badMuon = cosmicMuon = hfor = ge2blep = eq2blep = mllMin = false;
   l0eta = l1eta = tauVeto = dileptonTrigger = oppositeSign = l0pt = false;
  }
  bool grl, larErr, tileErr, ttcVeto, goodVtx, tileTrip, lAr;
  bool badJet, deadRegions, badMuon, cosmicMuon, hfor, ge2blep, eq2blep, mllMin;
    bool l0eta, l1eta, tauVeto, dileptonTrigger, oppositeSign, l0pt;
  bool allTrue() {
    return  (grl && larErr && tileErr && ttcVeto && goodVtx && tileTrip && lAr &&
             badJet && deadRegions && badMuon && cosmicMuon && hfor && ge2blep && eq2blep && mllMin &&
             l0eta && l1eta && tauVeto && dileptonTrigger && oppositeSign && l0pt);
  }
  bool passAllEventCriteria() { return allTrue(); }
  bool failAny() { return !allTrue(); }
  std::string str() const;
};

} // namespace hlfv

#endif // end include guard
