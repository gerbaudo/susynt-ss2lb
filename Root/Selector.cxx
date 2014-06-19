#include "SusyntHlfv/Selector.h"
#include "SusyntHlfv/WeightComponents.h"
#include "SusyntHlfv/EventFlags.h"

// #include "SusyntHlfv/EventFlags.h"
// #include "SusyntHlfv/criteria.h"
// #include "SusyntHlfv/kinematic.h"
// #include "SusyntHlfv/utils.h"

#include "SusyNtuple/MCWeighter.h"
#include "SusyNtuple/DilTrigLogic.h"
#include "SusyNtuple/string_utils.h"

#include "TSystem.h"
#include "TVector2.h"

#include <cassert>
#include <cmath> // isnan
#include <cfloat> // FLT_MAX, FLT_MIN
#include <iomanip> // setw, setprecision
#include <sstream>      // std::ostringstream

using namespace std;
using hlfv::Selector;
using hlfv::WeightComponents;
using hlfv::EventFlags;

//-----------------------------------------
Selector::Selector() :
  m_trigObj(NULL),
  m_mcWeighter(NULL),
  m_useExistingList(false)
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
    // note: the event list must be initialized after sumw has been computed by mcweighter
    if(usingEventList()) initEventList(tree);
}
//-----------------------------------------
Bool_t Selector::Process(Long64_t entry)
{
    m_counter.nextEvent();
    m_printer.countAndPrint(cout);
    GetEntry(entry);
    m_chainEntry++; // SusyNtAna counter
    clearObjects();
    WeightComponents weightComponents;
    assignStaticWeightComponents(nt, *m_mcWeighter, weightComponents);
    m_counter.pass(weightComponents.product());
    bool removeLepsFromIso(false);
    selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium);
    hlfv::EventFlags eventFlags = computeEventFlags();
    incrementEventCounters(eventFlags, weightComponents);
    bool isSelectedEvent = !eventFlags.failAny();
    if(isSelectedEvent) {
        if(usingEventList() && !m_useExistingList)
            m_eventList.addEvent(entry);
    }
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
    m_counter.printTableRaw     (cout);
    m_counter.printTableWeighted(cout);
    if(m_mcWeighter) delete m_mcWeighter;
}
//-----------------------------------------
bool Selector::initMcWeighter(TTree *tree)
{
    bool success=false;
    if(tree){
        string xsecDir = gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc12_8TeV/");
        m_mcWeighter = new MCWeighter(tree, xsecDir);
        bool isPmssmSample(susy::utils::contains(sampleName(), "Herwigpp_UEEE3_CTEQ6L1_DGnoSL_TB10"));
        if(isPmssmSample) m_mcWeighter->setLabelBinCounter("Initial").clearAndRebuildSumwMap(m_tree);
        if(m_dbg) cout<<"Selector: MCWeighter has been initialized"<<endl;
    } else {
        cout<<"Selector::initMcWeighter: error, invalid input tree, cannot initialize Mcweighter"<<endl;
    }
    return success;
}
//-----------------------------------------
bool Selector::initEventList(TTree *tree)
{
   bool success = false;
   // check if the event list is there; if so, fetch it and loop only on those events
   m_useExistingList = m_eventList.cacheDoesExists();
   if(m_useExistingList) {
       tree->SetEventList(m_eventList.fetchEventList());
       if(m_dbg) cout<<"using existing event list from "<<m_eventList.cacheFilename()<<endl;
       success = true;
   }
   return success;
}
//-----------------------------------------
void Selector::assignStaticWeightComponents(/*const*/ Susy::SusyNtObject &ntobj,
                                            /*const*/ MCWeighter &weighter,
                                            WeightComponents &weightComponents)
{
    if(ntobj.evt()->isMC) {
        weightComponents.gen = ntobj.evt()->w;
        weightComponents.pileup = ntobj.evt()->wPileup;
        // for now just nom since we're handling the syst variation when filling trees
        const MCWeighter::WeightSys wSys = MCWeighter::Sys_NOM;
        weightComponents.susynt = weighter.getMCWeight(ntobj.evt(), LUMI_A_L, wSys);
        // getMCWeight provides gen * pu * xsec * lumi / sumw, so norm is xsec * lumi / sumw = susynt/(gen*pu)
        float genpu(weightComponents.gen*weightComponents.pileup);
        weightComponents.norm = (genpu != 0.0 ? weightComponents.susynt/genpu : 1.0);
    }
}
//-----------------------------------------
bool Selector::passEventCriteria()
{
    bool pass=true;

    return pass;
}
//-----------------------------------------
Selector& Selector::setEventListFilename(const std::string filename)
{
    m_eventListFilename = filename;
    if(usingEventList()) {
        m_eventList.setListName("hlfv_event_list");
        m_eventList.setCacheFilename(filename);
    }
    return *this;
}
//-----------------------------------------
void Selector::setDebug(int dbg)
{
    SusyNtAna::setDebug(dbg);
    m_eventList.setVerbose(dbg>0);
}
//-----------------------------------------
hlfv::EventFlags Selector::computeEventFlags()
{
    EventFlags f;
    if(m_dbg) cout<<"Selector::computeEventFlags"<<endl;
    int flag = nt.evt()->cutFlags[NtSys_NOM];
    const LeptonVector &bleps = m_baseLeptons;
    const JetVector     &jets = m_baseJets;
    const JetVector    &pjets = m_preJets;
    const Susy::Met      *met = m_met;
    uint run = nt.evt()->run;
    bool mc = nt.evt()->isMC;
    float mllMin(20);
    bool has2lep(bleps.size()>1 && bleps[0] && bleps[1]);
    float mll(has2lep ? (*bleps[0] + *bleps[1]).M() : 0.0);
    const int killHfor(4); // inheriting hardcoded magic values from HforToolD3PD.cxx
    if(passGRL        (flag           ))  f.grl         = true;
    if(passLarErr     (flag           ))  f.larErr      = true;
    if(passTileErr    (flag           ))  f.tileErr     = true;
    if(passTTCVeto    (flag           ))  f.ttcVeto     = true;
    if(passGoodVtx    (flag           ))  f.goodVtx     = true;
    if(passTileTripCut(flag           ))  f.tileTrip    = true;
    if(passLAr        (flag           ))  f.lAr         = true;
    if(!hasBadJet     (jets           ))  f.badJet      = true;
    if(passDeadRegions(pjets,met,run,mc)) f.deadRegions = true;
    if(!hasBadMuon    (m_preMuons     ))  f.badMuon     = true;
    if(!hasCosmicMuon (m_baseMuons    ))  f.cosmicMuon  = true;
    if(nt.evt()->hfor != killHfor      )  f.hfor        = true;
    if(bleps.size() >= 2               )  f.ge2blep     = true;
    if(bleps.size() == 2               )  f.eq2blep     = true;
    if(mll>mllMin                      )  f.mllMin      = true;
    return f;
}
//-----------------------------------------
void Selector::incrementEventCounters(const hlfv::EventFlags &f, const hlfv::WeightComponents &w)
{
    double weight = w.product();
    if(f.grl        ) m_counter.pass(weight); else return;
    if(f.larErr     ) m_counter.pass(weight); else return;
    if(f.tileErr    ) m_counter.pass(weight); else return;
    if(f.ttcVeto    ) m_counter.pass(weight); else return;
    if(f.goodVtx    ) m_counter.pass(weight); else return;
    if(f.tileTrip   ) m_counter.pass(weight); else return;
    if(f.lAr        ) m_counter.pass(weight); else return;
    if(f.badJet     ) m_counter.pass(weight); else return;
    if(f.deadRegions) m_counter.pass(weight); else return;
    if(f.badMuon    ) m_counter.pass(weight); else return;
    if(f.cosmicMuon ) m_counter.pass(weight); else return;
    if(f.hfor       ) m_counter.pass(weight); else return;
    if(f.ge2blep    ) m_counter.pass(weight); else return;
    if(f.eq2blep    ) m_counter.pass(weight); else return;
    if(f.mllMin     ) m_counter.pass(weight); else return;
}
//-----------------------------------------
