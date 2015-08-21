#include "susynt-ss3l/Selector.h"
#include "susynt-ss3l/WeightComponents.h"
#include "susynt-ss3l/EventFlags.h"
#include "susynt-ss3l/DileptonVariables.h"
#include "susynt-ss3l/LeptonTruthType.h"
#include "susynt-ss3l/NtUtils.h"
#include "susynt-ss3l/utils.h"

// #include "susynt-ss3l/EventFlags.h"
// #include "susynt-ss3l/criteria.h"
// #include "susynt-ss3l/kinematic.h"
// #include "susynt-ss3l/utils.h"

#include "SusyNtuple/MCWeighter.h"
#include "SusyNtuple/string_utils.h"

#include "TString.h"
#include "TSystem.h"
#include "TVector2.h"

#include <algorithm> // remove_if
#include <cassert>
#include <cmath> // isnan
#include <cfloat> // FLT_MAX, FLT_MIN
#include <iomanip> // setw, setprecision
#include <iterator> // distance
#include <sstream>      // std::ostringstream

using namespace std;
using ss3l::Selector;
using ss3l::WeightComponents;
using ss3l::WeightVariations;
using ss3l::EventFlags;
using ss3l::Systematic;
using ss3l::DileptonVariables;

//-----------------------------------------
Selector::Selector() :
  m_mcWeighter(NULL),
  m_useExistingList(false),
  m_computeSystematics(false),
  m_tupleMaker("",""),
  m_writeTuple(false),
  m_outTupleFile(""),
  m_saveBaselineNonPrompt(false)
{
  nttools().setAnaType(Susy::AnalysisType::Ana_2LepWH);
  setSelectTaus(true);
}
//-----------------------------------------
void Selector::Begin(TTree* /*tree*/)
{
  SusyNtAna::Begin(0);
  initSystematicsList();
  if(m_writeTuple) initTupleWriters();
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
    selectObjects(Susy::NtSys::NOM, removeLepsFromIso, TauID_medium); // always select with nominal? (to compute event flags)
    removeForwardMuons(); // do this before computing the event flags (can affect the 2l test)
    EventFlags eventFlags = computeEventFlags();
    if(m_saveBaselineNonPrompt) eventFlags.eq2slep = eventFlags.eq2blep;
    if(incrementEventCounters(eventFlags, weightComponents)){
        for(size_t iSys=0; iSys<m_systematicsToProcess.size(); ++iSys){
            const Systematic::Value sys = m_systematicsToProcess[iSys];
            bool isNominal = sys==Systematic::CENTRAL;
            selectObjects(ss3l::sys2ntsys(sys), removeLepsFromIso, TauID_medium);
            const JetVector&   bj = m_baseJets; // these are just used to compute the btag weight
            const JetVector&  jets= m_signalJets2Lep;
            const LeptonVector& l = m_saveBaselineNonPrompt ? m_baseLeptons : m_signalLeptons;
            if(l.size()==2) { // several vars cannot be computed if we don't have 2 lep
                const JetVector cljets; // TODO (Selector::filterJets(jets, m_jvfTool, sys, m_anaType));
                DileptonVariables vars = computeDileptonVariables(l, m_met, cljets, jets, m_signalTaus);
                assignNonStaticWeightComponents(l, bj, sys, vars, weightComponents);
                if(isNominal){
                    incrementObjectCounters(vars, weightComponents, m_counter);
                    incrementObjectSplitCounters(vars, weightComponents);
                }
                bool is_data(!nt.evt()->isMC);
                bool two_mc_prompt = m_saveBaselineNonPrompt ? true : vars.hasTwoPromptLeptons;
                bool is_e_mu(eventIsEmu(l));
                bool is_same_sign(eventIsSameSign(l));
                bool is_event_to_be_saved = (vars.numTaus==0 &&
                                             (is_data || two_mc_prompt) &&
                                             eventFlags.mllMin &&
                                             vars.hasFiredTrig &&
                                             vars.hasTrigMatch &&
                                             abs(vars.eta0)<2.4 &&
                                             abs(vars.eta1)<2.4 &&
                                             (is_e_mu || is_same_sign));
                if(is_event_to_be_saved){
                    if(usingEventList() && isNominal && !m_useExistingList) m_eventList.addEvent(entry);
                    if(m_writeTuple) {
                        TupleMaker &tupleMaker = getTupleMaker(sys);
                        double weight(weightComponents.product());
                        unsigned int run(nt.evt()->run), event(nt.evt()->eventNumber), nVtx(nt.evt()->nVtx);
                        bool isMc = nt.evt()->isMC;
                        const Lepton &l0 = *l[0];
                        const Lepton &l1 = *l[1];
                        LeptonTruthType::Value l0Source = (isMc ? ss3l::getLeptonSource(l0) : LeptonTruthType::Unknown);
                        LeptonTruthType::Value l1Source = (isMc ? ss3l::getLeptonSource(l1) : LeptonTruthType::Unknown);
                        bool l0IsTight(nttools().isSignalLepton(&l0, m_baseElectrons, m_baseMuons, nVtx, isMc));
                        bool l1IsTight(nttools().isSignalLepton(&l1, m_baseElectrons, m_baseMuons, nVtx, isMc));
                        bool computeWeightVariations = (m_computeSystematics && sys==Systematic::CENTRAL);
                        WeightVariations wv = (computeWeightVariations ?
                                               computeSystematicWeightVariations(*nt.evt(), l, bj, sys, weightComponents) :
                                               WeightVariations());
                        tupleMaker
                            .setTriggerBits(nt.evt()->trigFlags)
                            .setWeightVariations(wv)
                            .setNumFjets(vars.numForwardJets)
                            .setNumBjets(vars.numBtagJets)
                            .setL0IsTight(l0IsTight).setL0Source(l0Source)
                            .setL1IsTight(l1IsTight).setL1Source(l1Source)
                            .setL0EtConeCorr(computeCorrectedEtCone(&l0))
                            .setL0PtConeCorr(computeCorrectedPtCone(&l0))
                            .setL1EtConeCorr(computeCorrectedEtCone(&l1))
                            .setL1PtConeCorr(computeCorrectedPtCone(&l1))
                            .fill(weight, run, event, l0, l1, *m_met, cljets);
                    } // m_writeTuple
                } // is_event_to_be_saved
            } // l.size()==2
        } // for(iSys)
    } // passAllEventCriteria
    // m_debugThisEvent = susy::isEventInList(nt.evt()->event);
    return kTRUE;
}
//-----------------------------------------
void Selector::Terminate()
{
    if(m_writeTuple) closeTupleWriters();
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
bool Selector::initMcWeighter(TTree *tree)
{
    bool success=false;
    // TODO drop
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
size_t Selector::initSystematicsList()
{
    m_systematicsToProcess.push_back(Systematic::CENTRAL);
    if(m_computeSystematics){
        m_systematicsToProcess.push_back(Systematic::EERUP      );
        m_systematicsToProcess.push_back(Systematic::EERDOWN    );
        m_systematicsToProcess.push_back(Systematic::EESZUP     );
        m_systematicsToProcess.push_back(Systematic::EESZDOWN   );
        m_systematicsToProcess.push_back(Systematic::EESLOWUP   );
        m_systematicsToProcess.push_back(Systematic::EESLOWDOWN );
        m_systematicsToProcess.push_back(Systematic::EESMATUP   );
        m_systematicsToProcess.push_back(Systematic::EESMATDOWN );
        m_systematicsToProcess.push_back(Systematic::EESPSUP    );
        m_systematicsToProcess.push_back(Systematic::EESPSDOWN  );
        m_systematicsToProcess.push_back(Systematic::MIDUP      );
        m_systematicsToProcess.push_back(Systematic::MIDDOWN    );
        m_systematicsToProcess.push_back(Systematic::JESUP      );
        m_systematicsToProcess.push_back(Systematic::JESDOWN    );
        m_systematicsToProcess.push_back(Systematic::JER        );
        m_systematicsToProcess.push_back(Systematic::RESOST     );
        m_systematicsToProcess.push_back(Systematic::MESUP      );
        m_systematicsToProcess.push_back(Systematic::MESDOWN    );
        m_systematicsToProcess.push_back(Systematic::SCALESTUP  );
        m_systematicsToProcess.push_back(Systematic::SCALESTDOWN);
    }
    return m_systematicsToProcess.size();
}
//-----------------------------------------
bool is_forward_muon(const Lepton *l)
{
    const float max_mu_eta_trigger = 2.4;
    return (l && l->isMu() && abs(l->Eta()) > max_mu_eta_trigger);
}
void Selector::removeForwardMuons()
{
    std::remove_if(m_baseLeptons.begin(),   m_baseLeptons.end(),   is_forward_muon);
    std::remove_if(m_signalLeptons.begin(), m_signalLeptons.end(), is_forward_muon);
}
//-----------------------------------------
void Selector::assignStaticWeightComponents(/*const*/ Susy::SusyNtObject &ntobj,
                                            /*const*/ MCWeighter &weighter,
                                            ss3l::WeightComponents &weightComponents)
{
    if(ntobj.evt()->isMC) {
        weightComponents.gen = ntobj.evt()->w;
        weightComponents.pileup = ntobj.evt()->wPileup;
        // TODO
        // for now just nom since we're handling the syst variation when filling trees
        // const MCWeighter::WeightSys wSys = MCWeighter::Sys_NOM;
        // weightComponents.susynt = 1.0; // TODO weighter.getMCWeight(ntobj.evt(), LUMI_A_L, wSys);
        // // getMCWeight provides gen * pu * xsec * lumi / sumw, so norm is xsec * lumi / sumw = susynt/(gen*pu)
        // float genpu(weightComponents.gen*weightComponents.pileup);
        // weightComponents.norm = (genpu != 0.0 ? weightComponents.susynt/genpu : 1.0);
    }
}
//-----------------------------------------
bool Selector::assignNonStaticWeightComponents(const LeptonVector& leptons,
                                               const JetVector& jets,
                                               const ss3l::Systematic::Value sys,
                                               ss3l::DileptonVariables &vars,
                                               ss3l::WeightComponents &weightcomponents)
{
    bool success=false;
    WeightComponents &wc = weightcomponents;
    if(leptons.size()>1) {
        const Lepton &l0 = *(leptons[0]);
        const Lepton &l1 = *(leptons[1]);
        if(nt.evt()->isMC) {
            // wc.lepSf   = computeDileptonEfficiencySf(l0, l1, sys);
            // wc.trigger = computeDileptonTriggerWeight(leptons, sys);
            // wc.btag    = computeBtagWeight(jets, nt.evt(), sys);
        }
        success = true;
    } else {
        cout<<"Selector::assignNonStaticWeightComponents: warning"<<endl
            <<" without at least two leptons several of these variables are not well defined."
            <<endl;
    }
    return success;
}
//-----------------------------------------
Selector& Selector::setEventListFilename(const std::string filename)
{
    m_eventListFilename = filename;
    if(usingEventList()) {
        m_eventList.setListName("ss3l_event_list");
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
ss3l::EventFlags Selector::computeEventFlags()
{
    EventFlags f;
    int flag = nt.evt()->cutFlags[Susy::NtSys::NOM];
    const LeptonVector &bleps = m_baseLeptons;
    const JetVector     &jets = m_baseJets;
    const JetVector    &pjets = m_preJets;
    const Susy::Met      *met = m_met;
    uint run = nt.evt()->run, event(nt.evt()->eventNumber);
    bool mc = nt.evt()->isMC;
    float mllMin(20);
    bool has2lep(bleps.size()>1 && bleps[0] && bleps[1]);
    float mll(has2lep ? (*bleps[0] + *bleps[1]).M() : 0.0);
    const int killHfor(4); // inheriting hardcoded magic values from HforToolD3PD.cxx
    bool pass_hfor(nt.evt()->hfor != killHfor);
    SusyNtTools &nt = nttools();
    if(pass_hfor                           )  f.hfor        = true;
    if(nt.passGRL        (flag            ))  f.grl         = true;
    if(nt.passLarErr     (flag            ))  f.larErr      = true;
    if(nt.passTileErr    (flag            ))  f.tileErr     = true;
    if(nt.passTTCVeto    (flag            ))  f.ttcVeto     = true;
    if(nt.passGoodVtx    (flag            ))  f.goodVtx     = true;
    if(nt.passTileTripCut(flag            ))  f.tileTrip    = true;
    if(nt.passLAr        (flag            ))  f.lAr         = true;
    if(!nt.hasBadJet     (jets            ))  f.badJet      = true;
    if(nt.passDeadRegions(pjets,met,run,mc))  f.deadRegions = true;
    if(!nt.hasBadMuon    (m_preMuons      ))  f.badMuon     = true;
    if(!nt.hasCosmicMuon (m_baseMuons     ))  f.cosmicMuon  = true;
    if(bleps.size() >= 2                   )  f.ge2blep     = true;
    if(bleps.size() == 2                   )  f.eq2blep     = true;
    if(mll>mllMin                          )  f.mllMin      = true;
    if(m_signalLeptons.size()==2           )  f.eq2slep     = true;
    if(m_signalTaus.size()==0              )  f.tauVeto     = true;
    return f;
}
//-----------------------------------------
bool Selector::incrementEventCounters(const ss3l::EventFlags &f, const ss3l::WeightComponents &w)
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
    if(f.eq2slep    ) m_counter.increment(weight, "eq2slep"    ); else return false;
    if(f.tauVeto    ) m_counter.increment(weight, "tauVeto"    ); else return false; // todo: this should go in DileptonVariables
    return true;
}
//-----------------------------------------
void Selector::incrementObjectCounters(const ss3l::DileptonVariables &v, const ss3l::WeightComponents &w,
                                       CutFlowCounter &counter)
{
    double weight = w.product();
    bool is_data = !nt.evt()->isMC;
    bool two_mc_prompt = is_data || v.hasTwoPromptLeptons;
    if(abs(v.eta0)<2.4      )    counter.increment(weight, "l0_eta<2.4"   ); else return;
    if(abs(v.eta1)<2.4      )    counter.increment(weight, "l1_eta<2.4"   ); else return;
    if(v.hasFiredTrig       )    counter.increment(weight, "trigger-bit"  ); else return;
    if(v.hasTrigMatch       )    counter.increment(weight, "trigger-match"); else return;
    if(two_mc_prompt        )    counter.increment(weight, "two-mc-propt" ); else return;
    if(v.isOs()             )    counter.increment(weight, "opp-sign"     ); else return;
    if(v.isOf()             )    counter.increment(weight, "opp-flav"     ); else return;
    if(v.pt0 > 45.0         )    counter.increment(weight, "pt0"          ); else return;
    if(bool passpt1=true    )    counter.increment(weight, "pt1"          ); else return;
    if(v.numCentralLightJets==0) counter.increment(weight, "jet-veto"     ); else return;
    if(bool passDphiLl=true )    counter.increment(weight, "dphi-l0l1"    ); else return;
    if(bool passDphiL1Met=true)  counter.increment(weight, "dphi-l1met"   ); else return;
}
//-----------------------------------------
void Selector::incrementObjectSplitCounters(const ss3l::DileptonVariables &v, const ss3l::WeightComponents &w)
{
    if(v.isEmu() || v.isMue()) {
        CutFlowCounter &counter = (v.isEmu() ? m_counterEmu : m_counterMue);
        incrementObjectCounters(v, w, counter);
    }
}
//-----------------------------------------
double Selector::computeDileptonTriggerWeight(const LeptonVector &leptons, const ss3l::Systematic::Value sys)
{
    double trigW = 1.0;
// TODO
    return trigW;
}
//-----------------------------------------
double Selector::computeBtagWeight(const JetVector& jets, const Susy::Event* evt, const ss3l::Systematic::Value sys)
{
    // TODO
    return 1.0;
}
//-----------------------------------------
double Selector::computeLeptonEfficiencySf(const Susy::Lepton &lep, const ss3l::Systematic::Value sys)
{
    float effFactor = 1.0;
    // TODO
    return effFactor;
}
//-----------------------------------------
double Selector::computeDileptonEfficiencySf(const Susy::Lepton &l0, const Susy::Lepton &l1, const ss3l::Systematic::Value sys)
{
    return (computeLeptonEfficiencySf(l0, sys)*
            computeLeptonEfficiencySf(l1, sys));
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
bool Selector::eventIsOppositeSign(const LeptonVector &leptons)
{
    return (leptons.size()==2 && (leptons[0]->q*leptons[1]->q <0));
}
//-----------------------------------------
bool Selector::eventIsSameSign(const LeptonVector &leptons)
{
    return (leptons.size()==2 && (leptons[0]->q*leptons[1]->q >0));
}
//-----------------------------------------
JetVector Selector::filterJets(const JetVector &jets, JVFUncertaintyTool* jvfTool,
                               const ss3l::Systematic::Value sys,
                               AnalysisType anaType)
{
    JetVector outjets;
    for(size_t i=0; i<jets.size(); ++i){
        if(nttools().m_jetSelector.isCentralLightJet(jets[i]))
            outjets.push_back(jets[i]);
    }
    return outjets;
}
//-----------------------------------------
JetVector Selector::filterForwardJets(const JetVector &jets)
{
    JetVector outjets;
    for(size_t i=0; i<jets.size(); ++i){
        if(nttools().m_jetSelector.isForwardJet(jets[i]))
            outjets.push_back(jets[i]);
    }
    return outjets;
}
//-----------------------------------------
JetVector Selector::filterBtagJets(const JetVector &jets)
{
    JetVector outjets;
    for(size_t i=0; i<jets.size(); ++i){
        if(nttools().m_jetSelector.isCentralBJet(jets[i]))
            outjets.push_back(jets[i]);
    }
    return outjets;
}
//----------------------------------------------------------
float Selector::computeCorrectedEtCone(const Lepton *l)
{
    float correctedEtCone = 0.0;
    if(l){
        uint nVtx(nt.evt()->nVtx);
        bool isMC(nt.evt()->isMC);
        if(l->isEle()) {
            if(const Electron* e = static_cast<const Electron*>(l))
                correctedEtCone = nttools().elEtTopoConeCorr(e, m_baseElectrons, m_baseMuons, nVtx, isMC);
        } else if(l->isMu()) {
            if(const Muon* m = static_cast<const Muon*>(l))
                correctedEtCone = nttools().muEtConeCorr(m, m_baseElectrons, m_baseMuons, nVtx, isMC);
        }
    }
    return correctedEtCone;
}
//----------------------------------------------------------
float Selector::computeCorrectedPtCone(const Lepton *l)
{
    float correctedPtCone = 0.0;
    if(l){
        uint nVtx(nt.evt()->nVtx);
        bool isMC(nt.evt()->isMC);
        if(l->isEle()) {
            if(const Electron* e = static_cast<const Electron*>(l))
                correctedPtCone = nttools().elPtConeCorr(e, m_baseElectrons, m_baseMuons, nVtx, isMC);
        } else if(l->isMu()) {
            if(const Muon* m = static_cast<const Muon*>(l))
                correctedPtCone = nttools().muPtConeCorr(m, m_baseElectrons, m_baseMuons, nVtx, isMC);
        }
    }
    return correctedPtCone;
}
//----------------------------------------------------------
WeightVariations Selector::computeSystematicWeightVariations(const Susy::Event &event,
                                                             const LeptonVector& leptons,
                                                             const JetVector& jets,
                                                             const ss3l::Systematic::Value sys,
                                                             const ss3l::WeightComponents &nominalWeightComponents)
{
    WeightVariations wv;
    ss3l::DileptonVariables dummyVars; // just used to retrieve the output from assignNonStaticWeightComponents (lepSf, trig, btag)
    ss3l::WeightComponents dummyWeightComps; // same as above
    // just shorter names
    const LeptonVector &l = leptons;
    const JetVector    &j = jets;
    const Lepton      &l0 = *l[0];
    const Lepton      &l1 = *l[1];
    ss3l::WeightComponents &dwc = dummyWeightComps;
    const ss3l::WeightComponents &nwc = nominalWeightComponents;
    dwc.trigger = computeDileptonTriggerWeight(l, Systematic::ETRIGREWUP  ); wv.elTrigUp = nwc.relativeTrig (dwc);
    dwc.trigger = computeDileptonTriggerWeight(l, Systematic::ETRIGREWDOWN); wv.elTrigDo = nwc.relativeTrig (dwc);
    dwc.trigger = computeDileptonTriggerWeight(l, Systematic::MTRIGREWUP  ); wv.muTrigUp = nwc.relativeTrig (dwc);
    dwc.trigger = computeDileptonTriggerWeight(l, Systematic::MTRIGREWDOWN); wv.muTrigDo = nwc.relativeTrig (dwc);
    dwc.btag    = computeBtagWeight   (j, &event, Systematic::BJETUP      ); wv.bTagUp   = nwc.relativeBtag (dwc);
    dwc.btag    = computeBtagWeight   (j, &event, Systematic::BJETDOWN    ); wv.bTagDo   = nwc.relativeBtag (dwc);
    dwc.btag    = computeBtagWeight   (j, &event, Systematic::CJETUP      ); wv.cTagUp   = nwc.relativeBtag (dwc);
    dwc.btag    = computeBtagWeight   (j, &event, Systematic::CJETDOWN    ); wv.cTagDo   = nwc.relativeBtag (dwc);
    dwc.btag    = computeBtagWeight   (j, &event, Systematic::BMISTAGUP   ); wv.lTagUp   = nwc.relativeBtag (dwc);
    dwc.btag    = computeBtagWeight   (j, &event, Systematic::BMISTAGDOWN ); wv.lTagDo   = nwc.relativeBtag (dwc);
    dwc.lepSf   = computeDileptonEfficiencySf(l0, l1, Systematic::ESFUP   ); wv.elEffUp  = nwc.relativeLepSf(dwc);
    dwc.lepSf   = computeDileptonEfficiencySf(l0, l1, Systematic::ESFDOWN ); wv.elEffDo  = nwc.relativeLepSf(dwc);
    dwc.lepSf   = computeDileptonEfficiencySf(l0, l1, Systematic::MEFFUP  ); wv.muEffUp  = nwc.relativeLepSf(dwc);
    dwc.lepSf   = computeDileptonEfficiencySf(l0, l1, Systematic::MEFFDOWN); wv.muEffDo  = nwc.relativeLepSf(dwc);
    wv.pileupUp = event.wPileup ? (event.wPileup_up / event.wPileup) : 1.0;
    wv.pileupDo = event.wPileup ? (event.wPileup_dn / event.wPileup) : 1.0;
    wv.swapUpDoIfNecessary();
    return wv;
}
//----------------------------------------------------------
ss3l::TupleMaker& Selector::getTupleMaker(const Systematic::Value s)
{
    if(s==Systematic::CENTRAL) return m_tupleMaker;
    else{
        vector<Systematic::Value>::iterator it = std::find(m_systematicsToProcess.begin(), m_systematicsToProcess.end(), s);
        if(it==m_systematicsToProcess.end()){
            cout<<"We do not have a tuple maker for "<<syst2str(s)<<endl
                <<"Stopping here"<<endl;
            assert(false);
        }
        size_t nominalOffset = 1;
        size_t systIndex = std::distance(m_systematicsToProcess.begin(), it)-nominalOffset;
        return *(m_systTupleMakers[systIndex]);
    }
}
//----------------------------------------------------------
bool Selector::initTupleWriters()
{
    bool success = false;
    if(!Susy::utils::endswith(m_outTupleFile, ".root")){
        cout<<"You must provide a root output file for the ntuple"<<endl;
    } else {
        string filenameNom = m_outTupleFile;
        string treename = "ss3l_tuple";
        success = m_tupleMaker.init(filenameNom, treename);
        if(success) cout<<"initialized ntuple file "<<m_tupleMaker.filename()<<endl;
        size_t nominalOffset=1;
        for(size_t iSys=nominalOffset; iSys<m_systematicsToProcess.size(); ++iSys) {
            const Systematic::Value sys = m_systematicsToProcess[iSys];
            string filename(filenameNom);
            // ss3l::replace(filename, string(".root"), string("_"+syst2str(sys)+".root")); // why doesn't it work?? todo
            filename = TString(filenameNom.c_str()).ReplaceAll(TString(".root"),
                                                               TString("_"+syst2str(sys)+".root")).Data();
            TupleMaker *tm = new TupleMaker("", "");
            m_systTupleMakers.push_back(tm);
            success &= tm->init(filename, treename);
            if(success) cout<<"initialized ntuple file "<<tm->filename()<<endl;
        }
    }
    if(!success){
        cout<<"failed to initialize some of the ntuples; disabling m_writeTuple"<<endl;
        m_writeTuple = false;
    }
    return success;
}
//----------------------------------------------------------
bool Selector::closeTupleWriters()
{
    cout<<m_tupleMaker.summary()<<endl;
    m_tupleMaker.close();
    for(size_t iTm=0; iTm<m_systTupleMakers.size(); ++iTm){
        cout<<m_systTupleMakers[iTm]->summary()<<endl;
        m_systTupleMakers[iTm]->close();
        delete m_systTupleMakers[iTm];
    }
    m_systTupleMakers.clear();
    return true; // find a criterion to decide whether we've closed everything properly
}
//----------------------------------------------------------
