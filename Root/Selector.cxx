#include "SusyntHlfv/Selector.h"

#include "SusyntHlfv/EventFlags.h"
#include "SusyntHlfv/criteria.h"
#include "SusyntHlfv/kinematic.h"
#include "SusyntHlfv/utils.h"

#include "SusyNtuple/MCWeighter.h"

#include "TSystem.h"
#include "TVector2.h"

#include <cassert>
#include <cmath> // isnan
#include <cfloat> // FLT_MAX, FLT_MIN
#include <iomanip> // setw, setprecision
#include <sstream>      // std::ostringstream

using namespace std;
using hlfv::Selector;

//-----------------------------------------
Selector::Selector() :
  m_trigObj(NULL),
  m_mcWeighter(NULL),
{
  setAnaType(Ana_2Lep);
  setSelectTaus(true);
}
//-----------------------------------------
void Selector::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  string period = "Moriond";
  bool useReweightUtils = false;
  m_trigObj = new DilTrigLogic(period, useReweightUtils);
//  if(m_useMCTrig) m_trigObj->useMCTrigger();
}
//-----------------------------------------
void Selector::Init(TTree* tree)
{
    SusyNtAna::Init(tree);
    initMcWeighter(tree);
}
//-----------------------------------------
Bool_t Selector::Process(Long64_t entry)
{
    m_printer.countAndPrint(cout);
    GetEntry(entry);
    clearObjects();
    cacheStaticWeightComponents();
//  increment(n_readin, m_weightComponents);
    bool removeLepsFromIso(false), allowQflip(true);
    selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium);
    swh::EventFlags eventFlags(computeEventFlags());
    incrementCounters(eventFlags, m_weightComponents);
  // if(eventFlags.failAny()) return kTRUE;
  // m_debugThisEvent = susy::isEventInList(nt.evt()->event);

  // const JetVector&   bj = m_baseJets;
  // const LeptonVector& l = m_signalLeptons;
  // if(l.size()>1) assignNonStaticWeightComponents(computeNonStaticWeightComponents(l, bj, susy::wh::WH_CENTRAL));
  // else return false;
  // VarFlag_t varsFlags = computeSsFlags(m_signalLeptons, m_signalTaus, m_signalJets2Lep, m_met, susy::wh::WH_CENTRAL, allowQflip);
  // const SsPassFlags &ssf = varsFlags.second;
  // incrementSsCounters(ssf, m_weightComponents);
  // if(ssf.lepPt) {
  //     if(m_writeTuple) {
  //         double weight(m_weightComponents.product());
  //         unsigned int run(nt.evt()->run), event(nt.evt()->event);
  //         LeptonVector anyLep(getAnyElOrMu(nt));
  //         LeptonVector lowPtLep(subtract_vector(anyLep, m_baseLeptons));
  //         const Lepton *l0 = m_signalLeptons[0];
  //         const Lepton *l1 = m_signalLeptons[1];
  //         const JetVector clJets(Selector::filterClJets(m_signalJets2Lep, m_jvfTool, NtSys_NOM, m_anaType));
  //         m_tupleMaker.fill(weight, run, event, *l0, *l1, *m_met, lowPtLep, clJets);
  //     }
  // }
    return kTRUE;
}
//-----------------------------------------
void Selector::Terminate()
{
    SusyNtAna::Terminate();
    // dumpEventCounters();
    if(m_mcWeighter) delete m_mcWeighter;
}
//-----------------------------------------
bool Selector::initMcWeighter(TTree *tree)
{
    bool success=false;
    if(tree){
        string xsecDir = gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc12_8TeV/");
        m_mcWeighter = new MCWeighter(tree, xsecDir);
        bool isPmssmSample(contains(sampleName(), "Herwigpp_UEEE3_CTEQ6L1_DGnoSL_TB10"));
        if(isPmssmSample) m_mcWeighter->setLabelBinCounter("Initial").clearAndRebuildSumwMap(m_tree);
        if(m_dbg) cout<<"Selector: MCWeighter has been initialized"<<endl;
    } else {
        cout<<"Selector::initMcWeighter: error, invalid input tree, cannot initialize Mcweighter"<<endl;
    }
    return success;
}
//-----------------------------------------
