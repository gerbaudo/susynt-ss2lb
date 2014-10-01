#include "SusyntHlfv/Selector.h"
#include "SusyntHlfv/WeightComponents.h"
#include "SusyntHlfv/EventFlags.h"
#include "SusyntHlfv/DileptonVariables.h"

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
using hlfv::Systematic;
using hlfv::DileptonVariables;

//-----------------------------------------
Selector::Selector() :
  m_trigObj(NULL),
  m_mcWeighter(NULL),
  m_counter(Selector::defaultCutNames()),
  m_counterEmu(Selector::defaultCutNamesSplit()),
  m_counterMue(Selector::defaultCutNamesSplit()),
  m_useExistingList(false),
  m_tupleMaker("",""),
  m_writeTuple(false),
  m_outTupleFile("")
{
  setAnaType(Ana_2Lep);
  setSelectTaus(true);
}
//-----------------------------------------
void Selector::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  initDilTrigLogic();
  if(m_writeTuple) {
      if(susy::utils::endswith(m_outTupleFile, ".root") &&
         m_tupleMaker.init(m_outTupleFile, "hlfv_tuple"))
          cout<<"initialized ntuple file "<<m_outTupleFile<<endl;
      else {
          cout<<"cannot initialize ntuple file '"<<m_outTupleFile<<"'"<<endl;
          m_writeTuple = false;
      }
  }
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
    m_counterEmu.nextEvent();
    m_counterMue.nextEvent();
    m_printer.countAndPrint(cout);
    GetEntry(entry);
    m_chainEntry++; // SusyNtAna counter
    clearObjects();
    WeightComponents weightComponents;
    assignStaticWeightComponents(nt, *m_mcWeighter, weightComponents);
    m_counter.pass(weightComponents.product());
    bool removeLepsFromIso(false);
    selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium); // always select with nominal? (to compute event flags)
    EventFlags eventFlags = computeEventFlags();
    incrementEventCounters(eventFlags, weightComponents);
    if(eventFlags.passAllEventCriteria()) {
        const Systematic::Value sys = Systematic::CENTRAL; // syst loop will go here
        const JetVector&   bj = m_baseJets; // why are we using basejets and not m_signalJets2Lep?
        const LeptonVector& l = m_signalLeptons;
        if(eventHasTwoLeptons(l)) { // several vars cannot be computed if we don't have 2 lep
            const JetVector jets(Selector::filterJets(m_signalJets2Lep, m_jvfTool, sys, m_anaType));
            DileptonVariables vars = computeDileptonVariables(l, m_met, jets);
            assignNonStaticWeightComponents(l, bj, sys, vars, weightComponents);
            incrementObjectCounters(vars, weightComponents, m_counter);
            incrementObjectSplitCounters(vars, weightComponents);
            // m_tupleMaker.fill(weight, run, event, *l0, *l1, *m_met, jets); // todo (just re-use the one from wh)
            bool is_event_to_be_saved = eventIsEmu(l);
            if(is_event_to_be_saved){
                if(usingEventList() && !m_useExistingList) m_eventList.addEvent(entry);
                if(m_writeTuple) {
                    double weight(weightComponents.product());
                    unsigned int run(nt.evt()->run), event(nt.evt()->event);
                    const Lepton &l0 = *m_signalLeptons[0];
                    const Lepton &l1 = *m_signalLeptons[1];
                    m_tupleMaker.fill(weight, run, event, l0, l1, *m_met);
                }
            } // is_event_to_be_saved
        } // eventHasTwoLeptons
    } // passAllEventCriteria
    // m_debugThisEvent = susy::isEventInList(nt.evt()->event);
    return kTRUE;
}
//-----------------------------------------
void Selector::Terminate()
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
bool Selector::initDilTrigLogic()
{
  string period = "Moriond";
  bool useReweightUtils = false;
  m_trigObj = new DilTrigLogic(period, useReweightUtils);
//  if(m_useMCTrig) m_trigObj->useMCTrigger);
  return (m_trigObj!=NULL);
}
//-----------------------------------------
bool Selector::initMcWeighter(TTree *tree)
{
    bool success=false;
    if(tree){
        string xsecDir = gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc12_8TeV/");
        m_mcWeighter = new MCWeighter(tree, xsecDir);
        bool isPmssmSample(susy::utils::contains(sampleName(), "Herwigpp_UEEE3_CTEQ6L1_DGnoSL_TB10"));
        m_mcWeighter->parseAdditionalXsecFile("${ROOTCOREBIN}/data/SusyntHlfv/LFV.txt", m_dbg);
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
                                            hlfv::WeightComponents &weightComponents)
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

bool Selector::assignNonStaticWeightComponents(const LeptonVector& leptons,
                                               const JetVector& jets,
                                               const hlfv::Systematic::Value sys,
                                               hlfv::DileptonVariables &vars,
                                               hlfv::WeightComponents &weightcomponents)
{
    bool success=false;
    WeightComponents &wc = weightcomponents;
    if(leptons.size()>1) {
        vars.hasFiredTrig = m_trigObj->passDilEvtTrig  (leptons, m_met->Et, nt.evt());
        vars.hasTrigMatch = m_trigObj->passDilTrigMatch(leptons, m_met->Et, nt.evt());
        const Lepton &l0 = *(leptons[0]);
        const Lepton &l1 = *(leptons[1]);
        if(nt.evt()->isMC) {
            wc.lepSf   = (computeLeptonEfficiencySf(l0, sys)*
                          computeLeptonEfficiencySf(l1, sys));
            wc.trigger = computeDileptonTriggerWeight(leptons, sys);
            wc.btag    = computeBtagWeight(jets, nt.evt(), sys);
        }
        success = true;
    }
    return success;
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
    int flag = nt.evt()->cutFlags[NtSys_NOM];
    const LeptonVector &bleps = m_baseLeptons;
    const JetVector     &jets = m_baseJets;
    const JetVector    &pjets = m_preJets;
    const TauVector     &taus = m_signalTaus;
    const Susy::Met      *met = m_met;
    uint run = nt.evt()->run;
    bool mc = nt.evt()->isMC;
    float mllMin(20);
    bool has2lep(bleps.size()>1 && bleps[0] && bleps[1]);
    float mll(has2lep ? (*bleps[0] + *bleps[1]).M() : 0.0);
    const int killHfor(4); // inheriting hardcoded magic values from HforToolD3PD.cxx
    bool pass_hfor(nt.evt()->hfor != killHfor);
    if(pass_hfor)                         f.hfor        = true;
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
    if(bleps.size() >= 2               )  f.ge2blep     = true;
    if(bleps.size() == 2               )  f.eq2blep     = true;

    // const LeptonVector& leptons, DilTrigLogic *dtl, float met, Event* evt
    // dtl->passDilEvtTrig(leptons, met, evt);
    // dtl->passDilTrigMatch(leptons, met, evt);

    if(mll>mllMin                      )  f.mllMin      = true;
    if(taus.size()==0                  )  f.tauVeto     = true;
    return f;
}
//-----------------------------------------
void Selector::incrementEventCounters(const hlfv::EventFlags &f, const hlfv::WeightComponents &w)
{
    double weight = w.product();
    if(f.hfor       ) m_counter.pass(weight); else return;
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
    if(f.ge2blep    ) m_counter.pass(weight); else return;
    if(f.eq2blep    ) m_counter.pass(weight); else return;
    if(f.mllMin     ) m_counter.pass(weight); else return;
    if(f.tauVeto    ) m_counter.pass(weight); else return;
}
//-----------------------------------------
void Selector::incrementObjectCounters(const hlfv::DileptonVariables &v, const hlfv::WeightComponents &w,
                                       CutFlowCounter &counter)
{
    double weight = w.product();
    if(true                 ) counter.pass(weight); else return; // flavor
    if(v.hasFiredTrig       ) counter.pass(weight); else return;
    if(v.hasTrigMatch       ) counter.pass(weight); else return;
    if(v.isOs()             ) counter.pass(weight); else return;
    if(v.pt0 > 45.0         ) counter.pass(weight); else return;
    if(bool passpt1=true    ) counter.pass(weight); else return;
    if(v.numCentralLightJets==0) counter.pass(weight); else return;
    if(bool passDphiLl=true ) counter.pass(weight); else return;
    if(bool passDphiL1Met=true) counter.pass(weight); else return;
}
//-----------------------------------------
void Selector::incrementObjectSplitCounters(const hlfv::DileptonVariables &v, const hlfv::WeightComponents &w)
{
    if(v.isEmu() || v.isMue()) {
        CutFlowCounter &counter = (v.isEmu() ? m_counterEmu : m_counterMue);
        incrementObjectCounters(v, w, counter);
    }
}
//-----------------------------------------
double Selector::computeDileptonTriggerWeight(const LeptonVector &leptons, const hlfv::Systematic::Value sys)
{
    double trigW = 1.0;
    if(leptons.size()==2){

        trigW = m_trigObj->getTriggerWeight(leptons, nt.evt()->isMC, m_met->Et,
                                            m_signalJets2Lep.size(),
                                            nt.evt()->nVtx, hlfv::sys2ntsys(sys));
        bool twIsInvalid(isnan(trigW) || trigW<0.0);
        if(twIsInvalid){
            if(m_dbg) cout<<"SusySelection::getTriggerWeight: invalid weight "<<trigW<<", using 0.0"<<endl;
            trigW = (twIsInvalid ? 0.0 : trigW);
        }
        assert(!twIsInvalid); // is this still necessary? DG-2014-06-24
    }
    return trigW;
}
//-----------------------------------------
double Selector::computeBtagWeight(const JetVector& jets, const Susy::Event* evt, const hlfv::Systematic::Value sys)
{
    JetVector taggableJets = SusyNtTools::getBTagSFJets2Lep(jets);
    return SusyNtTools::bTagSF(evt, taggableJets, evt->mcChannel, hlfv::sys2ntbsys(sys));
}
//-----------------------------------------
double Selector::computeLeptonEfficiencySf(const Susy::Lepton &lep, const hlfv::Systematic::Value sys)
{
    float effFactor = 1.0;
    float sf(lep.effSF), delta(0.0);
    if     (lep.isEle() && sys==hlfv::Systematic::ESFUP   ) delta = (+lep.errEffSF);
    else if(lep.isEle() && sys==hlfv::Systematic::ESFDOWN ) delta = (-lep.errEffSF);
    else if(lep.isMu()  && sys==hlfv::Systematic::MEFFUP  ) delta = (+lep.errEffSF);
    else if(lep.isMu()  && sys==hlfv::Systematic::MEFFDOWN) delta = (-lep.errEffSF);
    effFactor = (sf + delta);
    return effFactor;
}
//-----------------------------------------
bool Selector::eventHasTwoLeptons(const LeptonVector &leptons)
{
    return leptons.size()==2;
}
//-----------------------------------------
bool Selector::eventIsEmu(const LeptonVector &leptons)
{
    bool isEmu = false;
    if(leptons.size()==2) {
        const Lepton &l0 = *(leptons[0]);
        const Lepton &l1 = *(leptons[1]);
        isEmu = ((l0.isEle() && l1.isMu()) ||
                 (l0.isMu() && l1.isEle()) );
    }
    return isEmu;
}
//-----------------------------------------
std::vector<std::string> Selector::defaultCutNames()
{ // note to self: these labels must match the calls of
  // m_counter.pass() (i.e. the first one and the subsequent ones in
  // incrementEventCounters
    vector<string> labels;
    labels.push_back("input"      );
    labels.push_back("hfor"       );
    labels.push_back("grl"        );
    labels.push_back("larErr"     );
    labels.push_back("tileErr"    );
    labels.push_back("ttcVeto"    );
    labels.push_back("goodVtx"    );
    labels.push_back("tileTrip"   );
    labels.push_back("lAr"        );
    labels.push_back("badJet"     );
    labels.push_back("deadRegions");
    labels.push_back("badMuon"    );
    labels.push_back("cosmicMuon" );
    labels.push_back("ge2blep"    );
    labels.push_back("eq2blep"    );
    labels.push_back("mllMin"     );
    labels.push_back("tauVeto"    );
    labels.push_back("flavor");
    labels.push_back("trigger-bit");
    labels.push_back("trigger-match");
    labels.push_back("opp-sign");
    labels.push_back("pt0");
    labels.push_back("pt1");
    labels.push_back("jet-veto");
    labels.push_back("dphi-l0l1");
    labels.push_back("dphi-l1met");
    return labels;
}
//-----------------------------------------
std::vector<std::string> Selector::defaultCutNamesSplit()
{
    vector<string> labels;
    labels.push_back("flavor");
    labels.push_back("trigger-bit");
    labels.push_back("trigger-match");
    labels.push_back("opp-sign");
    labels.push_back("pt0");
    labels.push_back("pt1");
    labels.push_back("jet-veto");
    labels.push_back("dphi-l0l1");
    labels.push_back("dphi-l1met");
    return labels;
}
//-----------------------------------------
JetVector Selector::filterJets(const JetVector &jets, JVFUncertaintyTool* jvfTool,
                               const hlfv::Systematic::Value sys,
                               AnalysisType anaType)
{
    JetVector outjets;
    for(size_t i=0; i<jets.size(); ++i){
        if(SusyNtTools::isCentralLightJet(jets[i], jvfTool, hlfv::sys2ntsys(sys), anaType))
            outjets.push_back(jets[i]);
    }
    return outjets;
}
//-----------------------------------------
