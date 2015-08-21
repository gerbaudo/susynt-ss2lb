// emacs -*- C++ -*-
#ifndef SS3L_EVENTFLAGS_H
#define SS3L_EVENTFLAGS_H

#include <string>
#include <iostream>

namespace ss3l
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
        badJet = deadRegions = badMuon = cosmicMuon = hfor = ge2blep = eq2blep = eq2slep = mllMin = false;
        l0eta = l1eta = tauVeto = forwardJetVeto = bjetVeto = dileptonTrigger = oppositeSign = l0pt = false;
    }
    bool grl, larErr, tileErr, ttcVeto, goodVtx, tileTrip, lAr;
    bool badJet, deadRegions, badMuon, cosmicMuon, hfor, ge2blep, eq2blep, eq2slep, mllMin;
    bool l0eta, l1eta, tauVeto, forwardJetVeto, bjetVeto, dileptonTrigger, oppositeSign, l0pt;
    bool allTrue() const {
        return  (passAllEventCriteria() &&
                 ge2blep && eq2blep && eq2slep && mllMin &&
                 l0eta && l1eta && tauVeto && forwardJetVeto && bjetVeto && dileptonTrigger && oppositeSign && l0pt);
    }
    bool passAllEventCriteria() const {
        return (grl && larErr && tileErr && ttcVeto && goodVtx && tileTrip && lAr &&
                badJet && deadRegions && badMuon && cosmicMuon && hfor);
    }
    bool failAny() const { return !allTrue(); }
    std::string str() const;
    std::string str_what_fails() const;
};

} // namespace ss3l

#endif // end include guard
