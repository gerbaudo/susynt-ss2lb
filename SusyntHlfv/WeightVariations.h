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
    /// set everything to 1
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
      pileupUp = pileupDo = 1.0;
      fakeElRealUp = fakeElRealDo = 1.0;
      fakeElFakeUp = fakeElFakeDo = 1.0;
      fakeMuRealUp = fakeMuRealDo = 1.0;
      fakeMuFakeUp = fakeMuFakeDo = 1.0;
      fakeElFracUp = fakeElFracDo = 1.0;
      fakeMuFracUp = fakeMuFracDo = 1.0;
      return *this;
    }
    /// HistFitter wants 'up' above 1, 'do' below 1
    WeightVariations& swapUpDoIfNecessary() {
      struct SwapFunc{ double tmp_; void operator()(double &a, double &b) {tmp_=a; a=b; b=tmp_;} } swap;
      if(qflipUp      < qflipDo     ) swap(qflipUp      , qflipDo     );
      if(elTrigUp     < elTrigDo    ) swap(elTrigUp     , elTrigDo    );
      if(muTrigUp     < muTrigDo    ) swap(muTrigUp     , muTrigDo    );
      if(elEffUp      < elEffDo     ) swap(elEffUp      , elEffDo     );
      if(muEffUp      < muEffDo     ) swap(muEffUp      , muEffDo     );
      if(bTagUp       < bTagDo      ) swap(bTagUp       , bTagDo      );
      if(cTagUp       < cTagDo      ) swap(cTagUp       , cTagDo      );
      if(lTagUp       < lTagDo      ) swap(lTagUp       , lTagDo      );
      if(xsecUp       < xsecDo      ) swap(xsecUp       , xsecDo      );
      if(mcgenUp      < mcgenDo     ) swap(mcgenUp      , mcgenDo     );
      if(pileupUp     < pileupDo    ) swap(pileupUp     , pileupDo    );
      if(fakeElRealUp < fakeElRealDo) swap(fakeElRealUp , fakeElRealDo);
      if(fakeElFakeUp < fakeElFakeDo) swap(fakeElFakeUp , fakeElFakeDo);
      if(fakeMuRealUp < fakeMuRealDo) swap(fakeMuRealUp , fakeMuRealDo);
      if(fakeMuFakeUp < fakeMuFakeDo) swap(fakeMuFakeUp , fakeMuFakeDo);
      if(fakeElFracUp < fakeElFracDo) swap(fakeElFracUp , fakeElFracDo);
      if(fakeMuFracUp < fakeMuFracDo) swap(fakeMuFracUp , fakeMuFracDo);
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
    double pileupUp, pileupDo;
    double fakeElRealUp, fakeElRealDo;
    double fakeElFakeUp, fakeElFakeDo;
    double fakeMuRealUp, fakeMuRealDo;
    double fakeMuFakeUp, fakeMuFakeDo;
    double fakeElFracUp, fakeElFracDo;
    double fakeMuFracUp, fakeMuFracDo;
  };

} // hlfv

#endif // HLVF_WEIGHTVARIATIONS_H
