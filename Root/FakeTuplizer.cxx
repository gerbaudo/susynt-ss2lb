#include "susynt-ss3l/FakeTuplizer.h"
#include "susynt-ss3l/WeightComponents.h"
#include "susynt-ss3l/EventFlags.h"
#include "susynt-ss3l/DileptonVariables.h"

using namespace std;
using ss3l::FakeTuplizer;
using ss3l::Selector;
using ss3l::WeightComponents;
using ss3l::EventFlags;
using ss3l::Systematic;
using ss3l::DileptonVariables;

//-----------------------------------------
FakeTuplizer::FakeTuplizer() :
    m_tupleMakerEmu("", "")
{
  nttools().setAnaType(Susy::AnalysisType::Ana_2LepWH);
  setSelectTaus(true);
}
//-----------------------------------------
void FakeTuplizer::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
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
    selectObjects(Susy::NtSys::NOM, removeLepsFromIso, TauID_medium); // always select with nominal? (to compute event flags)
    EventFlags eventFlags = computeEventFlags();
    incrementEventCounters(eventFlags, weightComponents);
    if(eventFlags.passAllEventCriteria()) {
        const Systematic::Value sys = Systematic::CENTRAL; // syst loop will go here
        const JetVector&   bj = m_baseJets; // why are we using basejets and not m_signalJets2Lep?
        const LeptonVector& l = m_signalLeptons;
        if(l.size()==2 && eventIsEmu(l)) { // several vars cannot be computed if we don't have 2 lep
            const JetVector cljets; // TODO (Selector::filterJets(m_signalJets2Lep, m_jvfTool, sys, nttools().getAnaType()));
            DileptonVariables vars = computeDileptonVariables(l, m_met, cljets, m_signalJets2Lep, m_signalTaus);
            assignNonStaticWeightComponents(l, bj, sys, vars, weightComponents);
            incrementObjectSplitCounters(vars, weightComponents);
            if(usingEventList() && !m_useExistingList) m_eventList.addEvent(entry);
            if(m_writeTuple) {
                double weight(weightComponents.product());
                unsigned int run(nt.evt()->run), event(nt.evt()->eventNumber);
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
