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
    if(incrementEventCounters(eventFlags, weightComponents)){
        const Systematic::Value sys = Systematic::CENTRAL; // syst loop will go here
        const JetVector&   bj = m_baseJets; // why are we using basejets and not m_signalJets2Lep?
        const JetVector&  jets= m_signalJets; // shouldn't we use m_signalJets2Lep?
        const LeptonVector& l = m_signalLeptons;
        if(eventHasTwoLeptons(l)) { // several vars cannot be computed if we don't have 2 lep
            const JetVector cljets(Selector::filterJets(jets, m_jvfTool, sys, m_anaType));
            DileptonVariables vars = computeDileptonVariables(l, m_met, cljets, jets, m_signalTaus);
            assignNonStaticWeightComponents(l, bj, sys, vars, weightComponents);
            incrementObjectCounters(vars, weightComponents, m_counter);
            incrementObjectSplitCounters(vars, weightComponents);
            bool is_event_to_be_saved = (vars.numTaus==0 &&
                                         vars.numBtagJets==0 &&
                                         vars.numForwardJets==0 &&
                                         eventFlags.mllMin &&
                                         vars.hasFiredTrig &&
                                         vars.hasTrigMatch &&
                                         eventIsEmu(l));
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
    const Susy::Met      *met = m_met;
    uint run = nt.evt()->run, event(nt.evt()->event);
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
    if(mll>mllMin                      )  f.mllMin      = true;
    if(m_signalLeptons.size()==2)         f.eq2slep     = true;
    if(m_signalTaus.size()==0          )  f.tauVeto     = true;
    return f;
}
//-----------------------------------------
bool Selector::incrementEventCounters(const hlfv::EventFlags &f, const hlfv::WeightComponents &w)
{
    double weight = w.product();
    if(f.hfor       ) m_counter.increment(weight, "hfor"       ); else return false;
    if(f.grl        ) m_counter.increment(weight, "grl"        ); else return false;
    if(f.larErr     ) m_counter.increment(weight, "larErr"     ); else return false;
    if(f.tileErr    ) m_counter.increment(weight, "tileErr"    ); else return false;
    if(f.ttcVeto    ) m_counter.increment(weight, "ttcVeto"    ); else return false;
    if(f.goodVtx    ) m_counter.increment(weight, "goodVtx"    ); else return false;
    if(f.tileTrip   ) m_counter.increment(weight, "tileTrip"   ); else return false;
    if(f.lAr        ) m_counter.increment(weight, "lAr"        ); else return false;
    if(f.badJet     ) m_counter.increment(weight, "badJet"     ); else return false;
    if(f.deadRegions) m_counter.increment(weight, "deadRegions"); else return false;
    if(f.badMuon    ) m_counter.increment(weight, "badMuon"    ); else return false;
    if(f.cosmicMuon ) m_counter.increment(weight, "cosmicMuon" ); else return false;
    if(f.ge2blep    ) m_counter.increment(weight, "ge2blep"    ); else return false;
    if(f.eq2blep    ) m_counter.increment(weight, "eq2blep"    ); else return false;
    if(f.mllMin     ) m_counter.increment(weight, "mllMin"     ); else return false; // todo: this should go in DileptonVariables
    if(f.eq2slep    ) m_counter.increment(weight, "eq2slep"); else return false;
    if(f.tauVeto    ) m_counter.increment(weight, "tauVeto"); else return false; // todo: this should go in DileptonVariables
    return true;
}
//-----------------------------------------
void Selector::incrementObjectCounters(const hlfv::DileptonVariables &v, const hlfv::WeightComponents &w,
                                       CutFlowCounter &counter)
{
    double weight = w.product();
    if(abs(v.eta0)<2.4      )    counter.increment(weight, "l0_eta<2.4"   ); else return;
    if(abs(v.eta1)<2.4      )    counter.increment(weight, "l1_eta<2.4"   ); else return;
    if(v.hasFiredTrig       )    counter.increment(weight, "trigger-bit"  ); else return;
    if(v.hasTrigMatch       )    counter.increment(weight, "trigger-match"); else return;
    if(v.isOs()             )    counter.increment(weight, "opp-sign"     ); else return;
    if(v.isOf()             )    counter.increment(weight, "opp-flav"     ); else return;
    if(v.pt0 > 45.0         )    counter.increment(weight, "pt0"          ); else return;
    if(bool passpt1=true    )    counter.increment(weight, "pt1"          ); else return;
    if(v.numCentralLightJets==0) counter.increment(weight, "jet-veto"     ); else return;
    if(bool passDphiLl=true )    counter.increment(weight, "dphi-l0l1"    ); else return;
    if(bool passDphiL1Met=true)  counter.increment(weight, "dphi-l1met"   ); else return;
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
JetVector Selector::filterForwardJets(const JetVector &jets)
{
    JetVector outjets;
    for(size_t i=0; i<jets.size(); ++i){
        if(SusyNtTools::isForwardJet(jets[i]))
            outjets.push_back(jets[i]);
    }
    return outjets;
}
//-----------------------------------------
JetVector Selector::filterBtagJets(const JetVector &jets)
{
    JetVector outjets;
    for(size_t i=0; i<jets.size(); ++i){
        if(SusyNtTools::isCentralBJet(jets[i]))
            outjets.push_back(jets[i]);
    }
    return outjets;
}
//-----------------------------------------
