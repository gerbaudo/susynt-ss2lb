#include "SusyntHlfv/CutFlowCounter.h"

#include <iomanip>
#include <stdio.h>

using hlfv::CutFlowCounter;

using std::cout;
using std::endl;
using std::string;
using std::vector;

//-----------------------------------------
CutFlowCounter::CutFlowCounter():
    m_iCut(0),
    m_nErrors(0),
    m_debug(false)
{
}
//-----------------------------------------
CutFlowCounter::CutFlowCounter(const std::vector<std::string> &cutnames):
    m_iCut(cutnames.size()),
    m_nErrors(0),
    m_debug(false)
{
    m_cut_names = cutnames;
    m_raw_counts.resize(m_cut_names.size(), 0);
    m_weighted_counts.resize(m_cut_names.size(), 0.0);
}
//-----------------------------------------
CutFlowCounter& CutFlowCounter::setDebug(bool v)
{
    m_debug = v;
    return *this;
}
//-----------------------------------------
CutFlowCounter& CutFlowCounter::nextEvent()
{
    m_iCut = 0;
    if(m_debug)
        cout<<"CutFlowCounter::nextEvent: reset m_iCut"<<endl;
    return *this;
}
//-----------------------------------------
CutFlowCounter& CutFlowCounter::pass()
{
    addCutIfNecessary();
    if(m_debug)
        cout<<"CutFlowCounter: increment"
            <<" raw["<<m_cut_names[m_iCut]<<"] ("<<m_raw_counts[m_iCut]<<"+1)"
            <<endl;
    m_raw_counts[m_iCut] += 1;
    m_iCut++;
    return *this;
}
//-----------------------------------------
CutFlowCounter& CutFlowCounter::pass(const double &w)
{
    addCutIfNecessary();
    if(m_debug)
        cout<<"CutFlowCounter: increment"
            <<" raw["<<m_cut_names[m_iCut]<<"] ("<<m_raw_counts[m_iCut]<<"+1)"
            <<" and"
            <<" weighted["<<m_cut_names[m_iCut]<<"] ("<<m_weighted_counts[m_iCut]<<"+"<<w<<")"
            <<endl;
    m_raw_counts[m_iCut] += 1;
    m_weighted_counts[m_iCut] += w;
    m_iCut++;
    return *this;
}
//-----------------------------------------
CutFlowCounter& CutFlowCounter::fail()
{
    addCutIfNecessary();
    m_iCut++;
    return *this;
}
//-----------------------------------------
void CutFlowCounter::printTableRaw     (std::ostream& oo) const
{
    int col1Width(12), col2Width(24);
    oo<<std::setw(col1Width)<<"Cut "<<std::setw(col2Width)<<"raw counts"<<endl;
    if(m_raw_counts.size()==m_cut_names.size()){
        for(size_t i=0; i<m_raw_counts.size(); ++i)
            oo<<std::setw(col1Width)<<m_cut_names[i]<<std::setw(col2Width)<<m_raw_counts[i]<<endl;
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
    int col1Width(12), col2Width(24);
    oo<<std::setw(col1Width)<<"Cut "<<std::setw(col2Width)<<"weighted counts"<<endl;
    if(m_raw_counts.size()==m_cut_names.size()){
        for(size_t i=0; i<m_weighted_counts.size(); ++i)
            oo<<std::setw(col1Width)<<m_cut_names[i]<<std::setw(col2Width)<<m_weighted_counts[i]<<endl;
    } else {
        oo<<"CutFlowCounter::printTableWeighted: invalid vector size"
          <<" (names["<<m_cut_names.size()<<"]"
          <<", counts["<<m_weighted_counts.size()<<"])"
          <<endl;
    }
}
//-----------------------------------------
std::string CutFlowCounter::defaultCutName(const size_t cutindex)
{
    const int bufsize=512;
    unsigned long index = static_cast<unsigned long>(cutindex);
    char buf[bufsize];
    snprintf(buf, bufsize, "cut_%03lu", index);
    string cutname(buf);
    return cutname;
}
//-----------------------------------------
void CutFlowCounter::addCut()
{
    addCut(defaultCutName(m_iCut));
}
//-----------------------------------------
void CutFlowCounter::addCut(const std::string &cutname)
{
    m_raw_counts.push_back(0);
    m_weighted_counts.push_back(0.0);
    m_cut_names.push_back(cutname);
}
//-----------------------------------------
void CutFlowCounter::addCutIfNecessary()
{
    if(m_iCut >= m_raw_counts.size()) {
        if(m_debug)
            cout<<"CutFlowCounter: cut index "<<m_iCut<<", counters["<<m_raw_counts.size()<<"]"
                <<" : adding a counter for '"<<defaultCutName(m_iCut)<<"'"
                <<endl;
        addCut();
    }
}
