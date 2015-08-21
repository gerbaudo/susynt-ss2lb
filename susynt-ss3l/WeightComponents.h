// emacs -*- C++ -*-
#ifndef HLVF_WEIGHTCOMPONENTS_H
#define HLVF_WEIGHTCOMPONENTS_H

#include <string>

namespace ss3l
{

/// A simple struct to pass around the event weight and all its components
/**
   Imported from github.com/gerbaudo/SusyTest0/SusyTest0/SusySelection.h
   The functions 'relative*' and 'replace*' are used to compute the
   relative factors that we need to fill the limit trees with
   systematic variations.

   davide.gerbaudo@gmail.com
   June 2014
*/
struct WeightComponents {
    WeightComponents() { reset(); }
    double product() const { return susynt * lepSf * btag * trigger * qflip * fake; } ///< product of all components
    void reset() { susynt = gen = pileup = norm = lepSf = btag = trigger = qflip = fake = 1.0; } ///< set everything to unity
    WeightComponents& replaceLepSf(double v) { lepSf = v; return *this; }
    WeightComponents& replaceBtag(double v) { btag = v; return *this; }
    WeightComponents& replaceTrig(double v) { trigger = v; return *this; }
    WeightComponents& replaceQflip(double v) { qflip = v; return *this; }
    WeightComponents& replaceFake(double v) { fake = v; return *this; }
    double relativeLepSf(double v) const { return (lepSf   ? v / lepSf   : 1.0); }
    double relativeBtag (double v) const { return (btag    ? v / btag    : 1.0); }
    double relativeTrig (double v) const { return (trigger ? v / trigger : 1.0); }
    double relativeQflip(double v) const { return (qflip   ? v / qflip   : 1.0); }
    double relativeFake (double v) const { return (fake    ? v / fake    : 1.0); }
    double relativeLepSf(const WeightComponents &w) const { return relativeLepSf(w.lepSf  ); }
    double relativeBtag (const WeightComponents &w) const { return relativeBtag (w.btag   ); }
    double relativeTrig (const WeightComponents &w) const { return relativeTrig (w.trigger); }
    double relativeQflip(const WeightComponents &w) const { return relativeQflip(w.qflip  ); }
    double relativeFake (const WeightComponents &w) const { return relativeFake (w.fake   ); }

    double susynt; ///< from SusyNtTools::getEventWeight: includes gen, pu, xs, lumi, sumw
    double gen, pileup, norm; ///< breakdown of the above; norm is xs*lumi/sumw
    double lepSf, btag, trigger, qflip, fake; ///< factors that we compute, not from upstream
    std::string str() const; ///< string representation of the weight factors
};

} // ss3l

#endif // HLVF_WEIGHTCOMPONENTS_H
