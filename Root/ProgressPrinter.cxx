#include "susynt-ss3l/ProgressPrinter.h"

#include <ctime>   // time_t, time
#include <iomanip> // setw 

using ss3l::ProgressPrinter;
//-----------------------------------------
void ProgressPrinter::countAndPrint(std::ostream& oo)
{
  m_eventCounter += 1;
  if(m_eventCounter!=m_intCounter) return;
  m_intCounter = m_suppressionFactor*m_intCounter;
  int colWidth(16), stampWidth(48);
  if(!m_quiet && (m_intCounter==m_suppressionFactor || m_intCounter>m_suppressionOffset)) {
    std::time_t t(std::time(NULL));
    oo<<""
      <<"Entry "<<std::setw(colWidth)<<m_eventCounter
      <<" "<<std::setw(stampWidth)<<std::asctime(std::localtime(&t))
      <<std::endl;
  }
}
