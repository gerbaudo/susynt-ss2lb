// emacs -*- C++ -*-
#ifndef HLVF_CUTFLOWCOUNTER_H
#define HLVF_CUTFLOWCOUNTER_H

#include <string>
#include <vector>

namespace hlfv
{

class WeightComponents;

/// A counter to keep track of the event counts in the cutflow.
/**
   Example usage:
   @code
   vector<string> cutnames;
   CutFlowCounter counter(cutnames);
   ...
   if(pass_cut) counter.pass(weightcomponents);
   else counter.setDebug(true).fail().setDebug(false);
   @endcode

   davide.gerbaudo@gmail.com
   June 2014
*/
class CutFlowCounter
{
public:
    CutFlowCounter();
    CutFlowCounter(const std::vector<std::string> &cutnames);
    CutFlowCounter& setDebug(bool v); ///< toggle debug
    CutFlowCounter& nextEvent(); ///< inform the counter that you are considering a new event
    CutFlowCounter& pass(); ///< increment the raw counts for the current cut stage
    CutFlowCounter& pass(const WeightComponents &w); ///< increment the raw and weighted counts for this cut stage
    ///
    /**
      You need to call this function only for a selection criterion
      that is non-sequential, i.e. for which you don't drop the event
      but you continue the cutflow. All it does it to increment the
      index of the internal cut. For example:
      @code
      if(hasSusyProp) counter.pass(); // just count, doesn't drop
      else counter.fail();
      @endcode
     */
    CutFlowCounter& fail();
    std::vector<std::string> cutNames() const { return m_cut_names; }
    std::vector<int> rawCounts() const { return m_raw_counts; }
    std::vector<double> weightedCounts() const { return m_weighted_counts; }
private:
    std::string defaultCutName(const size_t cutindex);
    void addCut();
    void addCut(const std::string &cutname);
    void addCutIfNecessary();
private:
    std::vector<std::string> m_cut_names; ///< names of the cuts; default is 'cut_N'
    std::vector<int> m_raw_counts; ///< counter for unweighted events (always filled)
    std::vector<double> m_weighted_counts; ///< counter for weighted events (filled if weight is provided)
    size_t m_iCut; ///< index of the cut that is going to be incremented
    size_t m_nErrors; ///< number of errors encountered (typically trying to increment beyond the number of known cuts)
    bool m_debug; ///< toggle to print debugging info
};

} // hlfv

#endif // HLVF_CUTFLOWCOUNTER_H
