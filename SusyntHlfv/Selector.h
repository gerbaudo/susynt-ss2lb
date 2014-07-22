// emacs -*- C++ -*-
#ifndef HLVF_SELECTOR_H
#define HLVF_SELECTOR_H

#include "SusyntHlfv/ProgressPrinter.h"
#include "SusyntHlfv/CutFlowCounter.h"
#include "SusyntHlfv/Systematic.h"
#include "SusyntHlfv/TupleMaker.h"

#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/EventlistHandler.h"


// fw decl
class DilTrigLogic;
class MCWeighter;
class JVFUncertaintyTool;
namespace Susy{ class Event; }

namespace hlfv{

class WeightComponents;
class EventFlags;
class DileptonVariables;

/// TSelector implementing the HLFV event selection
/**
   This TSelector reads in the SusyNtuple events and selects the ones
   passing the selection criteria we used for the Higgs lepton flavor
   violation (HLFV) study. It can also write out our analysis
   ntuples. For an example of its usage, see util/run_Selector.cxx.

   davide.gerbaudo@gmail.com
   June 2014
 */
class Selector : public SusyNtAna
{
public:
    Selector();
    virtual ~Selector(){};
    virtual void    Begin(TTree *tree);      ///< called before looping on entries
    virtual void    Init(TTree *tree);       ///< called when the TChain is attached
    virtual void    Terminate();             ///< called after looping is finished
    virtual Bool_t  Process(Long64_t entry); ///< called at each event
    Selector& setEventListFilename(const std::string filename);
    virtual void setDebug(int dbg); ///< overload SusyNtAna::setDebug
    /// toggle ouput ntuple option
    Selector& setWriteNtuple(bool val) { m_writeTuple = val; return *this; }
    Selector& setTupleFile(const std::string &name) { m_writeTuple = true; m_outTupleFile = name; return *this; }
    static std::vector<std::string> defaultCutNames(); ///< provide the list of default counter names (just labelling)
    /// provide the list of default counter names after emu/mue separation (just labelling)
    static std::vector<std::string> defaultCutNamesSplit();
protected:
    /// assign the weight components that depend only on event-level variables
    /**
       The output values are assigned to weightcomponents
       \todo using full lumi (LUMI_A_L); need to pass as parameter?
       \todo The input arguments should be const, but their getter
       methods are not declared const, so I cannot do that now.
     */
    static void assignStaticWeightComponents(/*const*/ Susy::SusyNtObject &ntobj,
                                             /*const*/ MCWeighter &weighter,
                                             hlfv::WeightComponents &weightComponents);
    /// assign the weight components that depend on the object-level variables
    /**
       The output values are assigned to vars and weightcomponents;
       need access to trigger and btag tools, so cannot be static.
       It's mainly weight components, but also need vars to keep track
       of the trigger matching.
    */
    bool assignNonStaticWeightComponents(const LeptonVector& leptons,
                                         const JetVector& jets,
                                         const hlfv::Systematic::Value sys,
                                         hlfv::DileptonVariables &vars,
                                         hlfv::WeightComponents &weightcomponents);
    /// compute the event-level flags
    /**
       Note that some of these quantities actually depend on the
       baseline objects, so you need to call it after selectObjects().
     */
    hlfv::EventFlags computeEventFlags();
    /// incremement the event-level counters
    void incrementEventCounters(const hlfv::EventFlags &f, const hlfv::WeightComponents &w);
    /// incremement the event counters after the emu/mue splitting
    void incrementEventSplitCounters(const hlfv::DileptonVariables &v, const hlfv::WeightComponents &w);
    /// dilepton trigger weight
    /**
       Need access to several internal variables, so cannot be static
     */
    double computeDileptonTriggerWeight(const LeptonVector &leptons, const hlfv::Systematic::Value sys);
    /// btag scale factor
    /**
       Need access to several internal variables, so cannot be static
     */
    double computeBtagWeight(const JetVector& jets, const Susy::Event* evt, const hlfv::Systematic::Value sys);
    /// select the jets we are interested in (central, high-pt)
    static JetVector filterJets(const JetVector &jets, JVFUncertaintyTool* jvfTool,
                                const hlfv::Systematic::Value sys,
                                AnalysisType anaType);
    /// lepton efficiency data/simulation scale factor
    static double computeLeptonEfficiencySf(const Susy::Lepton &lep, const hlfv::Systematic::Value sys);
    /// exactly two leptons
    static bool eventHasTwoLeptons(const LeptonVector &leptons);
    /// exactly one electron and one muon
    static bool eventIsEmu(const LeptonVector &leptons);
private:
    /// initialize weighter used for normalization
    bool initMcWeighter(TTree *tree);
    /// convention: we're using an event list if its filename was specified
    bool usingEventList() const { return m_eventListFilename.size()>0; }
    /// initialize event list
    /**
       To be called within Init(), after SusyNtAna::Begin() and after
       initMcWeighter().  Otherwise MCWeighter will get only a subset
       of the events and compute sumw incorrectly.
       Note to self: I am not sure that this feature would work on
       proof. We're not using proof, so who cares.
     */
    bool initEventList(TTree *tree);
protected:
    DilTrigLogic*       m_trigObj;      ///< trigger logic class
    MCWeighter*         m_mcWeighter;   ///< tool to determine the normalization
    hlfv::ProgressPrinter m_printer; ///< tool to print the progress
    hlfv::CutFlowCounter m_counter; ///< counters for cutflow
    hlfv::CutFlowCounter m_counterEmu; ///< counters for cutflow (after channel split, only emu)
    hlfv::CutFlowCounter m_counterMue; ///< counters for cutflow (after channel split, only mue)
    std::string m_eventListFilename; ///< name of the file with the eventlist (empty string means don't use this feature)
    bool m_useExistingList;        ///< to keep track of whether there is already an event list
    Susy::EventlistHandler m_eventList; ///< the actual event list
    hlfv::TupleMaker m_tupleMaker; ///< writer of our analysis nutples
    bool m_writeTuple; ///< whether we want to write the output ntuple
    std::string m_outTupleFile; ///< name of the file where the nutple will be written
    ClassDef(Selector, 1);
};

} // hlfv

#endif // HLVF_SELECTOR_H
