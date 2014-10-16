// emacs -*- C++ -*-
#ifndef HLVF_WEIGHTVARIATIONS_H
#define HLVF_WEIGHTVARIATIONS_H

#include <string>

namespace hlfv
{

/// A simple struct to pass around the systematic variations of the event weight
/**
   Imported from github.com/gerbaudo/SusyTest0/SusyTest0/HftFiller.h

   The variations are computed as the ratio to the nominal weight, and
   be stored in the nominal tree.

   davide.gerbaudo@gmail.com
   Oct 2014
*/

  struct WeightVariations {
    WeightVariations() { reset(); }
    WeightVariations& reset() {
      qflipUp = qflipDo = 1.0;
      elTrigUp = elTrigDo = 1.0;
      muTrigUp = muTrigDo = 1.0;
      elEffUp = elEffDo = 1.0;
      muEffUp = muEffDo = 1.0;
      bTagUp = bTagDo = 1.0;
      cTagUp = cTagDo = 1.0;
      lTagUp = lTagDo = 1.0;
      xsecUp = xsecDo = 1.0;
      mcgenUp = mcgenDo = 1.0;
      fakeElRealUp = fakeElRealDo = 1.0;
      fakeElFakeUp = fakeElFakeDo = 1.0;
      fakeMuRealUp = fakeMuRealDo = 1.0;
      fakeMuFakeUp = fakeMuFakeDo = 1.0;
      return *this;
    }
    double qflipUp, qflipDo;
    double elTrigUp, elTrigDo;
    double muTrigUp, muTrigDo;
    double elEffUp, elEffDo;
    double muEffUp, muEffDo;
    double bTagUp, bTagDo;
    double cTagUp, cTagDo;
    double lTagUp, lTagDo;
    double xsecUp, xsecDo;
    double mcgenUp, mcgenDo;
    double fakeElRealUp, fakeElRealDo;
    double fakeElFakeUp, fakeElFakeDo;
    double fakeMuRealUp, fakeMuRealDo;
    double fakeMuFakeUp, fakeMuFakeDo;
  };

} // hlfv

#endif // HLVF_WEIGHTVARIATIONS_H
