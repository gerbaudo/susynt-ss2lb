#include "SusyntHlfv/MatrixPrediction.h"

#include "SusyntHlfv/WeightComponents.h"
#include "SusyntHlfv/EventFlags.h"
#include "SusyntHlfv/DileptonVariables.h"
#include "SusyntHlfv/WeightVariations.h"
#include "SusyntHlfv/utils.h"

#include "DileptonMatrixMethod/Systematic.h"
#include "SusyNtuple/SusyDefs.h"

#include <iomanip>
#include <sstream>      // std::ostringstream

using namespace std;
using namespace Susy;
namespace sf = susy::fake;
using hlfv::MatrixPrediction;
using hlfv::Selector;
using hlfv::WeightComponents;
using hlfv::WeightVariations;
using hlfv::EventFlags;
using hlfv::DileptonVariables;
using sf::Parametrization;


//----------------------------------------------------------
MatrixPrediction::MatrixPrediction() :
   Selector(),
   m_dbg(false),
   m_matrix(0),
   m_use2dparametrization(false),
   m_allconfigured(false),
   m_computeSystematics(false)
{
}
//----------------------------------------------------------
void MatrixPrediction::Begin(TTree* /*tree*/)
{
  if(m_dbg) cout << "MatrixPrediction::Begin" << endl;
  Selector::Begin(0);
  m_allconfigured = initMatrixTool();
}
//----------------------------------------------------------
bool leptonIsFromPv(const Lepton &l)
{
    bool from_pv = false;
    if(!l.isMu() && !l.isEle()){
        cout<<"leptonIsFromPv, invalid lepton type"<<endl;
    } else {
        float d0Sig = fabs(l.d0Sig(true));
        float z0SinTheta = fabs(l.z0SinTheta(true));
        float d0Sig_max = l.isEle() ? ELECTRON_D0SIG_CUT : MUON_D0SIG_CUT;
        float z0SinTheta_max = l.isEle() ? ELECTRON_Z0_SINTHETA_CUT : MUON_Z0_SINTHETA_CUT;
        from_pv = ((d0Sig < d0Sig_max) && (z0SinTheta < z0SinTheta_max));
    }
    return from_pv;
}
//----------------------------------------------------------
Bool_t MatrixPrediction::Process(Long64_t entry)
{
#warning todo rename smm hmp


    m_printer.countAndPrint(cout);
    GetEntry(entry);
    m_chainEntry++; // SusyNtAna counter
    clearObjects();
    WeightComponents weightComponents;
    assignStaticWeightComponents(nt, *m_mcWeighter, weightComponents);
    m_counter.increment(weightComponents.product(), "input");
    bool removeLepsFromIso(false);
    selectObjects(NtSys_NOM, removeLepsFromIso, TauID_medium); // always select with nominal? (to compute event flags)
    removeForwardMuons();
    EventFlags eventFlags = computeEventFlags();
    eventFlags.eq2slep = true; // for fakes we only need baseline leptons
    // incrementEventCounters(eventFlags, weightComponents);
    // cout<<eventFlags.str()<<endl;
    // if(eventFlags.passAllEventCriteria()) {
    if(incrementEventCounters(eventFlags, weightComponents)) {
        // cout<<"passAllEventCriteria"<<endl;
        const JetVector& jets = m_signalJets2Lep;
        const JetVector&   bj = m_baseJets; // these are just used to compute the btag weight
        const LeptonVector& l = m_baseLeptons;
        const Met*          m = m_met;
        if(l.size()==2) { // several vars cannot be computed if we don't have 2 lep
            const JetVector cljets(Selector::filterJets(jets, m_jvfTool, Systematic::CENTRAL, m_anaType));
            DileptonVariables vars = computeDileptonVariables(l, m_met, cljets, m_signalJets2Lep, m_signalTaus);
            double gev=1.0;
            unsigned int run(nt.evt()->run), event(nt.evt()->event);
            const Lepton &l0 = *l[0];
            const Lepton &l1 = *l[1];
            assignNonStaticWeightComponents(l, bj, Systematic::CENTRAL, vars, weightComponents);
            incrementObjectCounters(vars, weightComponents, m_counter);
            incrementObjectSplitCounters(vars, weightComponents);
            bool is_e_mu(eventIsEmu(l)), is_same_sign(eventIsSameSign(l));
            bool is_event_to_be_saved = (vars.numTaus==0 &&
                                         eventFlags.mllMin &&
                                         vars.hasFiredTrig &&
                                         abs(vars.eta0)<2.4 &&
                                         abs(vars.eta1)<2.4 &&
                                         vars.hasTrigMatch &&
                                         (is_e_mu || is_same_sign));
            if(is_event_to_be_saved){
                uint nVtx = nt.evt()->nVtx;
                bool isMC = nt.evt()->isMC;
                float metRel = getMetRel(m, l, jets);
                bool l0IsSig(SusyNtTools::isSignalLepton(&l0, m_baseElectrons, m_baseMuons, nVtx, isMC));
                bool l1IsSig(SusyNtTools::isSignalLepton(&l1, m_baseElectrons, m_baseMuons, nVtx, isMC));
                string regionName="emuInc";
//            regionName = "CR_SSInc1j";
                sf::Systematic::Value sys = sf::Systematic::SYS_NOM;
                size_t iRegion = m_matrix->getIndexRegion(regionName);
                sf::Lepton fl0(l0IsSig, l0.isEle(), l0.Pt()*gev, l0.Eta());
                sf::Lepton fl1(l1IsSig, l1.isEle(), l1.Pt()*gev, l1.Eta());
                double weight = m_matrix->getTotalFake(fl0, fl1, iRegion, sys);
                WeightVariations wv = computeSystematicWeights(fl0, fl1, iRegion);
                m_tupleMaker
                    .setTriggerBits(nt.evt()->trigFlags)
                    .setWeightVariations(wv)
                    .setNumFjets(vars.numForwardJets)
                    .setNumBjets(vars.numBtagJets)
                    .setL0IsTight(l0IsSig)//.setL0Source(l0Source) // not available in data
                    .setL1IsTight(l1IsSig)//.setL1Source(l1Source)
                    .setL0EtConeCorr(computeCorrectedEtCone(&l0))
                    .setL0PtConeCorr(computeCorrectedPtCone(&l0))
                    .setL1EtConeCorr(computeCorrectedEtCone(&l1))
                    .setL1PtConeCorr(computeCorrectedPtCone(&l1))
                    .fill(weight, run, event, l0, l1, *m_met, cljets);
                // const JetVector clJets(SusySelection::filterClJets(m_signalJets2Lep));
                // m_tupleMaker.fill(weight, run, event, *l0, *l1, *m, lowPtLep, m_signalJets2Lep);
            } // is_event_to_be_saved
        } // l.size()==2
    } // if(emu)
    return true;
}
//----------------------------------------------------------
void MatrixPrediction::Terminate()
{
    Selector::Terminate();
    if(m_dbg) cout << "MatrixPrediction::Terminate" << endl;
    delete m_matrix;
}
//----------------------------------------------------------
float MatrixPrediction::getFakeWeight(const LeptonVector &baseLeps,
                                      std::string &regionName,
                                      float metRel,
                                      susy::fake::Systematic::Value)
{
/*
  if(baseLeps.size() != 2) return 0.0;
  uint nVtx = nt.evt()->nVtx;
  bool isMC = nt.evt()->isMC;
  float gev2mev(1000.);
  //m_matrix->setDileptonType(baseLeps[0]->isEle(), baseLeps[1]->isEle());
  const Susy::Lepton *l0=baseLeps[0], *l1=baseLeps[1];
  bool l0IsSig(SusyNtTools::isSignalLepton(l0, m_baseElectrons, m_baseMuons, nVtx, isMC));
  bool l1IsSig(SusyNtTools::isSignalLepton(l1, m_baseElectrons, m_baseMuons, nVtx, isMC));
  return m_matrix->getTotalFake(l0IsSig, l0->isEle(), l0->Pt(), l0->Eta(),
                                l1IsSig, l1->isEle(), l1->Pt(), l1->Eta(),
                                region, metRel, sys);
*/
    return 0.0;
}
//----------------------------------------------------------
float MatrixPrediction::getRFWeight(const LeptonVector &baseLeps,
                                    std::string &regionName,
                                    float metRel,
                                    susy::fake::Systematic::Value)
{
/*
  if(baseLeps.size() != 2) return 0.0;
  uint nVtx = nt.evt()->nVtx;
  bool isMC = nt.evt()->isMC;
  return m_matrix->getRF( isSignalLepton(baseLeps[0],m_baseElectrons, m_baseMuons,nVtx,isMC),
                          baseLeps[0]->isEle(),
                          baseLeps[0]->Pt(),
                          baseLeps[0]->Eta(),
                          isSignalLepton(baseLeps[1],m_baseElectrons, m_baseMuons,nVtx,isMC),
                          baseLeps[1]->isEle(),
                          baseLeps[1]->Pt(),
                          baseLeps[1]->Eta(),
                          regionName,
                          metRel,
                          sys);
*/
    return 0.0;
}
//----------------------------------------------------------
MatrixPrediction& MatrixPrediction::setMatrixFilename(const std::string filename)
{
    if(!hlfv::fileExists(filename))
    if(true)
    cout<<"MatrixPrediction::setMatrixFilename: invalid file '"<<filename<<"'"<<endl
        <<"\t"<<"something will go wrong"<<endl;
  m_matrixFilename = filename;
  return *this;
}
//----------------------------------------------------------
bool MatrixPrediction::initMatrixTool()
{
    m_matrix = new sf::DileptonMatrixMethod();
    sf::Parametrization::Value p = (m_use2dparametrization ? sf::Parametrization::PT_ETA : sf::Parametrization::PT);
    std::vector<std::string> regions;
    regions.push_back("emuInc");
    // regions.push_back("CR_SSInc1j");
    return m_matrix->configure(m_matrixFilename, regions, p, p, p, p);
}
//----------------------------------------------------------
std::string MatrixPrediction::dilepDetails(const Susy::Event &event,
                                           const DiLepEvtType &ll,
                                           const LeptonVector &ls)
{
  bool ee(ll==ET_ee), mm(ll==ET_mm);
  const Lepton *l0(ls.size()>0 ? ls[0] : NULL), *l1(ls.size()>1 ? ls[1] : NULL);
  float l0pt(l0 ? l0->Pt() : 0.0), l0eta(l0 ? l1->Eta() : 0.0);
  float l1pt(l1 ? l1->Pt() : 0.0), l1eta(l1 ? l1->Eta() : 0.0);
  std::ostringstream oss;
  oss<<"run "<<event.run
     <<" evt "<<event.event
     <<" "<<(ee?"ee":(mm?"mm":"em"))
     <<" l0: pt="<<l0pt<<" eta="<<l0eta
     <<" l1: pt="<<l1pt<<" eta="<<l1eta;
  return oss.str();
}
//----------------------------------------------------------
std::string MatrixPrediction::eventDetails(bool passSrSs, const Susy::Event &event,
                                           const DiLepEvtType &ll,
                                           const LeptonVector &ls)
{
  std::ostringstream oss;
  oss<<"MatrixPrediction passSrSs("<<(passSrSs?"true":"false")<<")"
     <<" "<<dilepDetails(event, ll, ls)
//     <<" weight="<<m_weightComponents.fake
      ;
  return oss.str();
}
//----------------------------------------------------------
hlfv::WeightVariations MatrixPrediction::computeSystematicWeights(const sf::Lepton &l0, const sf::Lepton &l1, size_t regionIndex)
{
    using sf::Systematic;
    WeightVariations wv;
    size_t &ri = regionIndex;
    if(m_computeSystematics){
        double nominal = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_NOM);
        if(nominal!=0){
            double in = 1.0/nominal;
            wv.fakeElRealUp = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_EL_RE_UP   ) * in;
            wv.fakeElRealDo = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_EL_RE_DOWN ) * in;
            wv.fakeElFakeUp = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_EL_FR_UP   ) * in;
            wv.fakeElFakeDo = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_EL_FR_DOWN ) * in;
            wv.fakeMuRealUp = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_MU_RE_UP   ) * in;
            wv.fakeMuRealDo = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_MU_RE_DOWN ) * in;
            wv.fakeMuFakeUp = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_MU_FR_UP   ) * in;
            wv.fakeMuFakeDo = m_matrix->getTotalFake(l0, l1, ri, Systematic::SYS_MU_FR_DOWN ) * in;
        }
    }
    return wv;
}
//----------------------------------------------------------
