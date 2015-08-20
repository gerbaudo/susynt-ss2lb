// emacs -*- C++ -*-
#ifndef HLVF_CUTFLOWCOUNTER_H
#define HLVF_CUTFLOWCOUNTER_H

#include <iostream>
#include <string>
#include <vector>


namespace hlfv
{

/// A counter to keep track of the event counts in the cutflow.
/**
   Example usage:
   @code
   CutFlowCounter counter;
   ...
   if(pass_cut) counter.pass(weightcomponents, "cut description");
   ...
   counter.printTableRaw();
   counter.printTableWeighted();
   @endcode

   davide.gerbaudo@gmail.com
   June 2014
*/
class CutFlowCounter
{
public:
    CutFlowCounter();
    CutFlowCounter& setDebug(bool v); ///< toggle debug
    /// increment the raw counts for the current cut stage
    CutFlowCounter& incrementRaw(const std::string &cutname);
    /// increment both raw and weighted counts for this cut stage
    CutFlowCounter& increment(const double &w, const std::string &cutname);
    std::vector<std::string> cutNames() const { return m_cut_names; }
    std::vector<int> rawCounts() const { return m_raw_counts; }
    std::vector<double> weightedCounts() const { return m_weighted_counts; }
    void printTableRaw     (std::ostream& oo) const;
    void printTableWeighted(std::ostream& oo) const;
private:
    bool hasCut(const std::string &cutname);
    size_t cutIndex(const std::string &cutname);
    void addCut(const std::string &cutname);
private:
    std::vector<std::string> m_cut_names; ///< names of the cuts; default is 'cut_N'
    std::vector<int> m_raw_counts; ///< counter for unweighted events (always filled)
    std::vector<double> m_weighted_counts; ///< counter for weighted events (filled if weight is provided)
    size_t m_index_last_cut; ///< index of the last cut that was incremented
    size_t m_nErrors; ///< number of errors encountered (typically trying to increment beyond the number of known cuts)
    bool m_debug; ///< toggle to print debugging info
};

} // hlfv

#endif // HLVF_CUTFLOWCOUNTER_H
