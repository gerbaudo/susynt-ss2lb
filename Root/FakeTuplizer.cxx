#include "susynt-ss2lb/FakeTuplizer.h"
#include "susynt-ss2lb/WeightComponents.h"
#include "susynt-ss2lb/EventFlags.h"
#include "susynt-ss2lb/DileptonVariables.h"


using namespace std;
using hlfv::FakeTuplizer;
using hlfv::Selector;
using hlfv::WeightComponents;
using hlfv::EventFlags;
using hlfv::Systematic;
using hlfv::DileptonVariables;

//-----------------------------------------
FakeTuplizer::FakeTuplizer() :
    m_tupleMakerEmu("", "")
{
  setAnaType(Ana_2LepWH);
  setSelectTaus(true);
}
//-----------------------------------------
void FakeTuplizer::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  initDilTrigLogic();
  initTuples();
}
//-----------------------------------------
void FakeTuplizer::Init(TTree* tree)
{
    Selector::Init(tree);
}
//-----------------------------------------
Bool_t FakeTuplizer::Process(Long64_t entry)
{
    m_printer.countAndPrint(cout);
    GetEntry(entry);
    m_chainEntry++; // SusyNtAna counter
    clearObjects();
    WeightComponents weightComponents;
    assignStaticWeightComponents(nt, *m_mcWeighter, weightComponents);
    m_counter.increment(weightComponents.product(), "input");
    bool removeLepsFromIso(false);
    selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium); // always select with nominal? (to compute event flags)
    EventFlags eventFlags = computeEventFlags();
    incrementEventCounters(eventFlags, weightComponents);
    if(eventFlags.passAllEventCriteria()) {
        const Systematic::Value sys = Systematic::CENTRAL; // syst loop will go here
        const JetVector&   bj = m_baseJets; // why are we using basejets and not m_signalJets2Lep?
        const LeptonVector& l = m_signalLeptons;
        if(l.size()==2 && eventIsEmu(l)) { // several vars cannot be computed if we don't have 2 lep
            const JetVector cljets(Selector::filterJets(m_signalJets2Lep, m_jvfTool, sys, m_anaType));
            DileptonVariables vars = computeDileptonVariables(l, m_met, cljets, m_signalJets2Lep, m_signalTaus);
            assignNonStaticWeightComponents(l, bj, sys, vars, weightComponents);
            incrementObjectSplitCounters(vars, weightComponents);
            if(usingEventList() && !m_useExistingList) m_eventList.addEvent(entry);
            if(m_writeTuple) {
                double weight(weightComponents.product());
                unsigned int run(nt.evt()->run), event(nt.evt()->event);
                const Lepton &l0 = *m_signalLeptons[0];
                const Lepton &l1 = *m_signalLeptons[1];
                m_tupleMaker
                    .setTriggerBits(nt.evt()->trigFlags)
                    .fill(weight, run, event, l0, l1, *m_met);
            }
        }
    }
  // m_debugThisEvent = susy::isEventInList(nt.evt()->event);
    return kTRUE;
}
//-----------------------------------------
void FakeTuplizer::Terminate()
{
    if(m_writeTuple) m_tupleMaker.close();
    SusyNtAna::Terminate();
    m_counter.printTableRaw     (cout);
    m_counter.printTableWeighted(cout);
    cout<<"--- emu ---"<<endl;
    m_counterEmu.printTableRaw     (cout);
    m_counterEmu.printTableWeighted(cout);
    cout<<"--- mue ---"<<endl;
    m_counterMue.printTableRaw     (cout);
    m_counterMue.printTableWeighted(cout);
    if(m_mcWeighter) delete m_mcWeighter;
}
//-----------------------------------------
bool FakeTuplizer::initTuples()
{
    return true;
}
//-----------------------------------------
std::string FakeTuplizer::tupleFilenameFromSamplename(const std::string &sampleName,
                                                      const std::string &suffix) const
{
    string tupleFname = m_outputdirectory+/*you are here*/suffix+".root";
    #warning todo
    return tupleFname;
}
//-----------------------------------------
