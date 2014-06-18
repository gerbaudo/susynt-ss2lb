#ifndef HLVF_PROGRESS_PRINTER_H
#define HLVF_PROGRESS_PRINTER_H

#include "Rtypes.h"

#include <iostream>

namespace hlfv
{
/*!
  A simple struct to print a timestamped line with the number of processed events.

  Details:
  Based on github.com/elaird/supy/steps/printer.py.

  davide.gerbaudo@gmail.com
  Oct 2013
*/
  struct ProgressPrinter {
    ProgressPrinter(int suppressionFactor=2, int suppressionOffset=300, bool quiet=false):
      m_suppressionFactor(suppressionFactor),
      m_suppressionOffset(suppressionOffset),
      m_eventCounter(0),
      m_intCounter(1),
      m_quiet(quiet) {};
    int m_suppressionFactor;
    int m_suppressionOffset;
    Long64_t m_eventCounter;
    Long64_t m_intCounter;
    bool m_quiet;
    void countAndPrint(std::ostream& oo);
  };
} // hlfv

#endif
