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
      qflipUp_ = qflipDo_ = 1.0;
      elTrigUp_ = elTrigDo_ = 1.0;
      muTrigUp_ = muTrigDo_ = 1.0;
      elEffUp_ = elEffDo_ = 1.0;
      muEffUp_ = muEffDo_ = 1.0;
      bTagUp_ = bTagDo_ = 1.0;
      cTagUp_ = cTagDo_ = 1.0;
      lTagUp_ = lTagDo_ = 1.0;
      xsecUp_ = xsecDo_ = 1.0;
      mcgenUp_ = mcgenDo_ = 1.0;
      return *this;
    }
    float qflipUp_, qflipDo_;
    float elTrigUp_, elTrigDo_;
    float muTrigUp_, muTrigDo_;
    float elEffUp_, elEffDo_;
    float muEffUp_, muEffDo_;
    float bTagUp_, bTagDo_;
    float cTagUp_, cTagDo_;
    float lTagUp_, lTagDo_;
    float xsecUp_, xsecDo_;
    float mcgenUp_, mcgenDo_;
  };

} // hlfv

#endif // HLVF_WEIGHTVARIATIONS_H
