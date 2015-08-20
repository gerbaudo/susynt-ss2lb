// emacs -*- C++ -*-
#ifndef HLVF_FAKETUPLIZER_H
#define HLVF_FAKETUPLIZER_H

#include "susynt-ss2lb/Selector.h"
#include "susynt-ss2lb/TupleMaker.h"

namespace hlfv{

/// Select emu events and save the ntuple to be used to compute the fake matrix
/**
   Right now we're only writing out the ntuple used to compute the
   compositions. Later on I'll add also the ones for the scale factors
   and all the other steps.

   davide.gerbaudo@gmail.com
   July 2014
 */
class FakeTuplizer : public Selector
{
public:
    FakeTuplizer();
    virtual ~FakeTuplizer(){};
    virtual void    Begin(TTree *tree);      ///< called before looping on entries
    virtual void    Init(TTree *tree);       ///< called when the TChain is attached
    virtual void    Terminate();             ///< called after looping is finished
    virtual Bool_t  Process(Long64_t entry); ///< called at each event
protected:
    bool initTuples(); ///< initialize output ntuples
private:
    std::string tupleFilenameFromSamplename(const std::string &sampleName, const std::string &suffix) const;
protected:
    hlfv::TupleMaker m_tupleMakerEmu; ///< writer of the nutples for the emu inclusive selection
    std::string m_outTupleFileEmu; ///< name of the file where the nutple will be written
//--    // DG-2014-07-23 These will come later, with the other regions
//--    LeptonVector m_probes;            // Probe lepton vector
//--    LeptonVector m_tags;              // Tag Lepton vector
private:
    std::string m_outputdirectory; ///< where the ntuples will be written
    ClassDef(FakeTuplizer, 1);
};

} // hlfv

#endif // HLVF_FAKETUPLIZER_H
