#include "susynt-ss3l/CutFlowCounter.h"

#include "susynt-ss3l/utils.h"

#include <algorithm>
#include <iomanip>
#include <stdio.h>

using ss3l::CutFlowCounter;

using std::cout;
using std::endl;
using std::string;
using std::vector;

//-----------------------------------------
CutFlowCounter::CutFlowCounter():
    m_index_last_cut(0),
    m_debug(false)
{
}
//-----------------------------------------
CutFlowCounter& CutFlowCounter::setDebug(bool v)
{
    m_debug = v;
    return *this;
}
//-----------------------------------------
CutFlowCounter& CutFlowCounter::incrementRaw(const std::string &cutname)
{
    size_t index=cutIndex(cutname);
    if(m_debug)
        cout<<"CutFlowCounter: increment"
            <<" raw["<<m_cut_names[index]<<"] ("<<m_raw_counts[index]<<"+1)"
            <<endl;
    m_raw_counts[index] += 1;
    m_index_last_cut = index;
    return *this;
}
//-----------------------------------------
CutFlowCounter& CutFlowCounter::increment(const double &w, const std::string &cutname)
{
    size_t index=cutIndex(cutname);
    if(m_debug)
        cout<<"CutFlowCounter: increment"
            <<" raw["<<m_cut_names[index]<<"] ("<<m_raw_counts[index]<<"+1)"
            <<endl;
    m_raw_counts[index] += 1;
    m_weighted_counts[index] += w;
    m_index_last_cut = index;
    return *this;
}
//-----------------------------------------
void CutFlowCounter::printTableRaw(std::ostream& oo) const
{
    int col0Width(3), col1Width(14), col2Width(24);
    oo<<std::setw(col1Width)<<"Cut "<<std::setw(col2Width)<<"raw counts"<<endl;
    if(m_raw_counts.size()==m_cut_names.size()){
        for(size_t i=0; i<m_raw_counts.size(); ++i)
            oo<<"Cut ["<<std::setw(col0Width)<<i<<"] "
              <<std::setw(col1Width)<<m_cut_names[i]
              <<std::setw(col2Width)<<m_raw_counts[i]<<endl;
    } else {
        oo<<"CutFlowCounter::printTableRaw: invalid vector size"
          <<" (names["<<m_cut_names.size()<<"]"
          <<", counts["<<m_raw_counts.size()<<"])"
          <<endl;
    }
}
//-----------------------------------------
void CutFlowCounter::printTableWeighted(std::ostream& oo) const
{
    int col0Width(3), col1Width(14), col2Width(24);
    oo<<std::setw(col1Width)<<"Cut "<<std::setw(col2Width)<<"weighted counts"<<endl;
    if(m_raw_counts.size()==m_cut_names.size()){
        for(size_t i=0; i<m_weighted_counts.size(); ++i)
            oo<<"Cut ["<<std::setw(col0Width)<<i<<"] "
              <<std::setw(col1Width)<<m_cut_names[i]
              <<std::setw(col2Width)<<m_weighted_counts[i]<<endl;
    } else {
        oo<<"CutFlowCounter::printTableWeighted: invalid vector size"
          <<" (names["<<m_cut_names.size()<<"]"
          <<", counts["<<m_weighted_counts.size()<<"])"
          <<endl;
    }
}
//-----------------------------------------
bool CutFlowCounter::hasCut(const std::string &cutname)
{
    bool has_cut=false;
    size_t best_guess = m_index_last_cut+1;
    if(best_guess<m_cut_names.size() &&
       m_cut_names[best_guess]==cutname){
        has_cut=true;
    } else {
        has_cut=ss3l::contains(m_cut_names, cutname);
    }
    return has_cut;
}
//-----------------------------------------
size_t CutFlowCounter::cutIndex(const std::string &cutname)
{
    size_t index=0;
    if(hasCut(cutname)){
        index = std::distance(m_cut_names.begin(), std::find(m_cut_names.begin(), m_cut_names.end(), cutname));
    } else {
        addCut(cutname);
        index = m_cut_names.size()-1;
        if(m_debug) cout<<"added cut["<<index<<"] '"<<cutname<<"'"<<endl;
    }
    return index;
}
//-----------------------------------------
void CutFlowCounter::addCut(const std::string &cutname)
{
    m_raw_counts.push_back(0);
    m_weighted_counts.push_back(0.0);
    m_cut_names.push_back(cutname);
}
//-----------------------------------------
