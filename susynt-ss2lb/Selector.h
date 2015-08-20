// emacs -*- C++ -*-
#ifndef HLVF_SELECTOR_H
#define HLVF_SELECTOR_H

#include "susynt-ss2lb/ProgressPrinter.h"
#include "susynt-ss2lb/CutFlowCounter.h"
#include "susynt-ss2lb/Systematic.h"
#include "susynt-ss2lb/TupleMaker.h"
#include "susynt-ss2lb/WeightVariations.h"

#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/EventlistHandler.h"

#include <string>

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
    Selector& useComputeSystematics(bool v) { m_computeSystematics = v; return *this; }
    virtual void setDebug(int dbg); ///< overload SusyNtAna::setDebug
    /// toggle ouput ntuple option
    Selector& setWriteNtuple(bool val) { m_writeTuple = val; return *this; }
    Selector& setTupleFile(const std::string &name) { m_writeTuple = true; m_outTupleFile = name; return *this; }
    /**
       @brief select events with two baseline leptons instead of two signal leptons

       This is to create the ntuples used for the fake estimate. All
       selection requirements (see `is_event_to_be_saved`) are the
       same as in the standard configuration, except for these two
       requirements that are dropped:

       \item the two leptons are baseline leptons rather than signal ones

       \item in simulated events, the leptons are not required to be true (prompt)

       The events in the resulting nutple are a superset of the events
       in the standard ntuple.
     */
    Selector& selectBaselineNonPromptLeptons() { m_saveBaselineNonPrompt = true; return *this; }
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
    /// compute both static and non-static weight systematic variations
   /**
      Mostly calls several times assignNonStaticWeightComponents
    */
    hlfv::WeightVariations computeSystematicWeightVariations(const Susy::Event &event,
                                                             const LeptonVector& leptons,
                                                             const JetVector& jets,
                                                             const hlfv::Systematic::Value sys,
                                                             const hlfv::WeightComponents &nominalWeightComponents);
    /// compute the event-level flags
    /**
       Note that some of these quantities actually depend on the
       baseline objects, so you need to call it after selectObjects().
     */
    hlfv::EventFlags computeEventFlags();
    /// incremement the event-level counters; return true if pass all event requirements
    bool incrementEventCounters(const hlfv::EventFlags &f, const hlfv::WeightComponents &w);
    /// incremement the object-level counters
    /**
       Increment the counters for selection criteria that depend on
       objects. Because here we keep several counters (inclusive, emu,
       mue), we also pass the counter as an input argument.
    */
    void incrementObjectCounters(const hlfv::DileptonVariables &v, const hlfv::WeightComponents &w,
                                 CutFlowCounter &counter);
    /// incremement the event counters after the emu/mue splitting
    /**
       Wraps incrementObjectCounters for the emu/mue specific cases.
     */
    void incrementObjectSplitCounters(const hlfv::DileptonVariables &v, const hlfv::WeightComponents &w);
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
    /// used to recompute the corrected iso when storing to ntuple (should refactor upstream SusyNtTools)
    float computeCorrectedEtCone(const Lepton *l);
    /// used to recompute the corrected iso when storing to ntuple (should refactor upstream SusyNtTools)
    float computeCorrectedPtCone(const Lepton *l);
    /// tuplemaker for a given systematic
    hlfv::TupleMaker& getTupleMaker(const Systematic::Value s);
    /// fill m_systematicsToProcess
    size_t initSystematicsList();
public:
    /// select the jets we are interested in (central, high-pt)
    static JetVector filterJets(const JetVector &jets, JVFUncertaintyTool* jvfTool,
                                const hlfv::Systematic::Value sys,
                                AnalysisType anaType);
    /// select forward jets
    static JetVector filterForwardJets(const JetVector &jets);
    /// select b jets
    static JetVector filterBtagJets(const JetVector &jets);
    /// lepton efficiency data/simulation scale factor
    static double computeLeptonEfficiencySf(const Susy::Lepton &lep, const hlfv::Systematic::Value sys);
    /// same as computeLeptonEfficiencySf, but for two leptons
    static double computeDileptonEfficiencySf(const Susy::Lepton &l0, const Susy::Lepton &l1, const hlfv::Systematic::Value sys);
    /// exactly one electron and one muon
    static bool eventIsEmu(const LeptonVector &leptons);
    /// two opposite-sign leptons
    static bool eventIsOppositeSign(const LeptonVector &leptons);
    /// two same-sign leptons
    static bool eventIsSameSign(const LeptonVector &leptons);

protected:
    /// initialize the 2L trig logic
    bool initDilTrigLogic();
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
    /// initialize output trees
    bool initTupleWriters();
    /// close output trees
    bool closeTupleWriters();
    /// remove from m_baseLeptons and m_signalLeptons the muons with |eta|>2.4
    /**
       These muons cannot be triggered, so in a dilepton event they
       will not be considered because they will fail the trigger
       match.  However, in a trilepton event they can cause the
       nlep==2 requirement to fail. This inconsistency is due to the
       fact that the third lepton veto (nlep==2) is applied before the
       trigger match. In fact, I think that this criterion should be
       applied within SusyNt::selectObjects; here it's a bit of a
       hack, but that's what we have for now. DG 2015-03-13
     */
    void removeForwardMuons();
    DilTrigLogic*       m_trigObj;      ///< trigger logic class
    MCWeighter*         m_mcWeighter;   ///< tool to determine the normalization
    hlfv::ProgressPrinter m_printer; ///< tool to print the progress
    hlfv::CutFlowCounter m_counter; ///< counters for cutflow
    hlfv::CutFlowCounter m_counterEmu; ///< counters for cutflow (after channel split, only emu)
    hlfv::CutFlowCounter m_counterMue; ///< counters for cutflow (after channel split, only mue)
    std::string m_eventListFilename; ///< name of the file with the eventlist (empty string means don't use this feature)
    bool m_useExistingList;        ///< to keep track of whether there is already an event list
    Susy::EventlistHandler m_eventList; ///< the actual event list
    bool m_computeSystematics; ///< whether the syst (weights for now) should be filled
    hlfv::TupleMaker m_tupleMaker; ///< writer of our analysis nutples
    std::vector<hlfv::TupleMaker*> m_systTupleMakers; ///< same as m_tupleMaker, but one for each object syst
    std::vector<hlfv::Systematic::Value> m_systematicsToProcess; ///< object systematics requiring a call to selectObjects()
    bool m_writeTuple; ///< whether we want to write the output ntuple
    std::string m_outTupleFile; ///< name of the file where the nutple will be written
    bool m_saveBaselineNonPrompt; ///< consider baseline and nonprompt in addition to signal leptons
    ClassDef(Selector, 2);
};

} // hlfv

#endif // HLVF_SELECTOR_H
