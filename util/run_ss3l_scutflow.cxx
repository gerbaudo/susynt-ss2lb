/**
   run a basic ss3l cutflow based on Superflow

   davide.gerbaudo@gmail.com
   Aug 2015
*/



#include <cstdlib>
#include <string>

#include <cstdlib>  // todo cleanup
#include <cmath>
#include <iostream>
#include <string>

#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"

#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"

using namespace std;
using namespace sflow; // TODO drop me

void usage(const char *exeName) {
  cout<<"Usage:"<<endl
      <<exeName<<" options"<<endl
      <<"\t"<<"-n [--num-event]   n_events (default -1, all)"<<endl
      <<"\t"<<"-k [--num-skip]    n_skip_events (default 0)" <<endl
      <<"\t"<<"-i [--input]       (file, list, or dir)"      <<endl
      <<"\t"<<"-s [--sample]      samplename"                <<endl
      <<"\t"<<"-d [--debug]       debug (>0 print stuff)"    <<endl
      <<"\t"<<"-h [--help]        print help"                <<endl
      <<endl;
}


bool passTrigger(const TriggerTools &t, const TBits & b)
{
    return (t.passTrigger(b, "HLT_2e12_loose_L12EM10VH"  ) ||
            t.passTrigger(b, "HLT_2e12_lhloose_L12EM10VH") ||
            t.passTrigger(b, "HLT_e17_lhloose_mu14"      ) ||
            t.passTrigger(b, "HLT_e17_loose_mu14"        ) ||
            t.passTrigger(b, "HLT_mu18_mu8noL1"          ) );
}

bool filterEvent(int run, long event)
{
    const int &r = run;
    const long &e = event;
    // return true;
    // davide not huan
    return (
        (r==270806 && e==5600684) ||
        (r==270806 && e==7314718) ||
        (r==270806 && e==8178077) ||
        (r==270806 && e==8468143) ||
        (r==270806 && e==8537765) ||
        (r==270806 && e==8282939) ||
        (r==270806 && e==8903766)
        );
    // huan not davide
    return (
        (r==270806 && e==13508931) ||
        (r==270806 && e==14466502) ||
        (r==270806 && e==15098331) ||
        (r==270806 && e==21461126) ||
        (r==270806 && e==30045453)
        );
}

int main(int argc, char** argv)
{
    int n_events = -1;
    int n_skip_events = 0;
    int dbg = 0;
    const double lumiPeriodC = 55.4;
    string sample;
    string input;

    int optind(1);
    while ((optind < argc)) {
        std::string sw = argv[optind];
        if(sw[0]!='-') { if(dbg) cout<<"skip "<<sw<<endl; optind++; continue; }
        else if(sw=="-n"||sw=="--num-event"  ) { n_events = atoi(argv[++optind]); }
        else if(sw=="-k"||sw=="--num-skip"   ) { n_skip_events = atoi(argv[++optind]); }
        else if(sw=="-d"||sw=="--debug"      ) { dbg = atoi(argv[++optind]); }
        else if(sw=="-i"||sw=="--input"      ) { input = argv[++optind]; }
        else if(sw=="-s"||sw=="--sample"     ) { sample = argv[++optind]; }
        else if(sw=="-h"||sw=="--help"       ) { usage(argv[0]); return 0; }
        else cout<<"Unknown switch "<<sw<<endl;
        optind++;
    } // end while(optind<argc)
    cout<<"flags:"                           <<endl
        <<"  sample          "<<sample       <<endl
        <<"  n_events        "<<n_events     <<endl
        <<"  n_skip_events   "<<n_skip_events<<endl
        <<"  dbg             "<<dbg          <<endl
        <<"  input           "<<input        <<endl
        <<endl;

    if(!input.length()){
        cout<<"--input is a required option"<<endl;
        return 1;
    }

    TChain* chain = new TChain("susyNt");
    bool verbose(dbg>0);
    bool printEvent = true;
    ChainHelper::addInput(chain, input, verbose);
    Long64_t nEntries = chain->GetEntries();
    n_events = (n_events<0 ? nEntries : n_events);
    if(dbg) chain->ls();

    SuperflowRunMode run_mode = SuperflowRunMode::nominal;
    // SusyNtSys nt_sys = NtSys::NOM;
    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setAnaName("ss3l cutflow");
    cutflow->setAnaType(AnalysisType::Ana_SS3L);
    cutflow->nttools().triggerTool().m_dbg = true;
    cutflow->nttools().initTriggerTool(ChainHelper::firstFile(input, verbose));
    cutflow->setSelectTaus(false);
    cutflow->setSampleName(sample);
    cutflow->setRunMode(run_mode);
    cutflow->setCountWeights(true); // print the weighted cutflows
    cutflow->setLumi(lumiPeriodC);
    cutflow->setChain(chain);
    cout<<"Number of chain entries: "<<chain->GetEntries()<<endl;

    // begin cutflow
    // *cutflow << CutName("read in") << [](Superlink*) -> bool { return true; };

    *cutflow << CutName("read in") << [&](Superlink* sl) -> bool {
        return filterEvent(sl->nt->evt()->run, sl->nt->evt()->eventNumber);
        if(printEvent)
            cout<<"read in"
                <<" run "<<sl->nt->evt()->run
                <<" event "<<sl->nt->evt()->eventNumber
                <<endl;
        return filterEvent(sl->nt->evt()->run, sl->nt->evt()->eventNumber);
        return true;
    };


    *cutflow << CutName("atleast-one-baseline") << [&](Superlink* sl) -> bool {
        ElectronVector &preEls = *sl->preElectrons;
        MuonVector & preMus = *sl->preMuons;
        size_t pre10lep = 0;
        for(auto e : preEls) if(e->Pt()>10.0 && fabs(e->Eta())<2.47 && e->looseLH) pre10lep++;
        for(auto m : preMus) if(m->Pt()>10.0 && fabs(m->Eta())<2.5 && m->medium ) pre10lep++;
        return pre10lep>0;
        // return (sl->preElectrons->size() + sl->preMuons->size())>0;
    };

    int cutflags = 0;

    *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
        // cout<<"run "<<sl->nt->evt()->run
        //     <<" event "<<sl->nt->evt()->eventNumber
        //     <<endl;
        cutflags = sl->nt->evt()->cutFlags[sl->nt_sys];
        if(printEvent)
            cout<<"---------------------------------------"
                <<endl
                <<"pass GRL : "<<bool(cutflags & ECut_GRL)
                <<" run "<<sl->nt->evt()->run
                <<" event "<<sl->nt->evt()->eventNumber
                <<endl;
        return (cutflags & ECut_GRL);
    };

    *cutflow << CutName("pass good vertex") << [&](Superlink*) -> bool {
        if(printEvent) cout<<"pass good vertex : "<<bool(cutflags & ECut_GoodVtx)<<endl;
        return (cutflags & ECut_GoodVtx);
    };

    *cutflow << CutName("Trigger") << [&] (Superlink *sl) -> bool {
        if(printEvent) cout<<"pass Trigger : "<<bool(passTrigger(cutflow->nttools().triggerTool(), sl->nt->evt()->trigBits))<<endl;
        return passTrigger(cutflow->nttools().triggerTool(), sl->nt->evt()->trigBits);
    };

    *cutflow << CutName("LAr error") << [&](Superlink*) -> bool {
        if(printEvent) cout<<"pass LAr error: "<<bool(cutflags & ECut_LarErr)<<endl;
        return (cutflags & ECut_LarErr);
    };

    *cutflow << CutName("Tile error") << [&](Superlink*) -> bool {
        if(printEvent) cout<<"pass Tile error "<<bool(cutflags & ECut_TileErr)<<endl;
        return (cutflags & ECut_TileErr);
    };

    *cutflow << CutName("TTC veto") << [&](Superlink*) -> bool {
        if(printEvent) cout<<"pass TTC veto: "<<bool(cutflags & ECut_TTC)<<endl;
        return (cutflags & ECut_TTC);
    };

    // *cutflow << CutName("at least 1 signal lepton") << [](Superlink* sl) -> bool {
    //     return sl->leptons->size() >= 1;
    // };


    // base leptons (after overlap removal)
    LeptonVector   baseLeptons;
    ElectronVector baseElectrons;
    MuonVector     baseMuons;
    JetVector      baseJets;
    // signal leptons
    LeptonVector   signalLeptons;
    ElectronVector signalElectrons;
    MuonVector     signalMuons;
    JetVector      signalJets;

    *cutflow << [&](Superlink* sl, var_void*) { baseLeptons     = *sl->baseLeptons; };
    *cutflow << [&](Superlink* sl, var_void*) { baseElectrons   = *sl->baseElectrons; };
    *cutflow << [&](Superlink* sl, var_void*) { baseMuons       = *sl->baseMuons; };
    *cutflow << [&](Superlink* sl, var_void*) { baseJets        = *sl->jets; };
    *cutflow << [&](Superlink* sl, var_void*) { signalLeptons   = *sl->leptons; };
    *cutflow << [&](Superlink* sl, var_void*) { signalElectrons = *sl->electrons; };
    *cutflow << [&](Superlink* sl, var_void*) { signalMuons     = *sl->muons; };

    *cutflow << CutName("bad muon veto") << [&](Superlink*) -> bool {
        if(printEvent) cout<<"pass bad muon veto: "<<bool(cutflags & ECut_BadMuon)<<endl;
        return (cutflags & ECut_BadMuon);
    };

    *cutflow << CutName(">=1 jet after OR") << [&](Superlink *sl) -> bool {
        baseJets = *sl->baseJets;
        if(printEvent) cout<<"pass >=1 jet after OR : "<<bool(baseJets.size()>0)<<endl;
        return baseJets.size()>0;
    };

    *cutflow << CutName("Bad jet") << [&](Superlink*) -> bool {
        if(printEvent) cout<<"pass Bad jet: "<<bool(cutflags & ECut_BadJet)<<endl;
        return (cutflags & ECut_BadJet);
    };

    *cutflow << CutName(">=1 signal jet") << [&](Superlink *sl) -> bool {
        // for(auto jet : baseJets) { // not necessary, use sl->jets; just leave as an example
        //     if(cutflow->nttools().getJetSelector().isSignalJet(jet)) {
        //         signalJets.push_back(jet);
        //     }
        // }
        signalJets = *sl->jets;
        if(printEvent) cout<<"pass >=1 signal jet : "<<bool(signalJets.size()>0)<<endl;
        return signalJets.size()>0;
    };

    *cutflow << CutName("pass cosmic veto") << [&](Superlink *sl) -> bool {
        if(printEvent) cout<<"pass cosmic veto : "<<bool(cutflags & ECut_Cosmic)<<endl;
        return (cutflags & ECut_Cosmic);
    };

    *cutflow << CutName(">=2 baseline leptons") << [&](Superlink *sl) -> bool {
        baseLeptons   = *sl->baseLeptons;
        baseMuons     = *sl->baseMuons;
        baseElectrons = *sl->baseElectrons;
        if(printEvent) cout<<"pass >=2 baseline leptons : "<<bool(baseLeptons.size()>1)<<endl;
        if(printEvent && (baseLeptons.size()<2)) {
            cout<<"failed pass11: "
                <<" run "<<sl->nt->evt()->run
                <<" event "<<sl->nt->evt()->eventNumber
                <<" now print leptons"
                <<endl;
            cutflow->dumpBaselineObjects();
            cutflow->dumpSignalObjects();
        }
        return baseLeptons.size()>1;
    };

    *cutflow << CutName(">=2 signal leptons print") << [&](Superlink *sl) -> bool {
        cout<<"after pass11: "
            <<" run "<<sl->nt->evt()->run
            <<" event "<<sl->nt->evt()->eventNumber
            <<" now print leptons"
            <<endl;
        cutflow->dumpBaselineObjects();
        cutflow->dumpSignalObjects();
        return true;
    };

    *cutflow << CutName(">=2 signal leptons") << [&](Superlink *sl) -> bool {
        signalLeptons   = *sl->leptons;
        signalMuons     = *sl->muons;
        signalElectrons = *sl->electrons;
        return signalLeptons.size()>1;
    };

    const Lepton *l0 = nullptr;
    const Lepton *l1 = nullptr;
    *cutflow << CutName("same-sign 20-20") << [&](Superlink *sl) -> bool {
        l0 = nullptr;
        l1 = nullptr;
        size_t nplus = 0;
        size_t nminus = 0;
        for(auto l : signalLeptons) {
            if(l->Pt()>20.0) {
                if(l->q>0) nplus++;
                if(l->q<0) nminus++;
            }
        }
        bool samesign = (nplus>1 || nminus>1);
        if(samesign) {
            for(auto l : signalLeptons) {
                if(l->Pt()>20.0) {
                    bool is_right_sign = ((nplus>1 && l->q>0) || (nminus>1 && l->q<0));
                    if     (!l0 && is_right_sign) l0 = l;
                    else if(!l1 && is_right_sign) l1 = l;
                    else if(l0 && l1) break;
                }
            }
        }
        // // bool twolep = signalLeptons.size()==2;
        // // l0 = twolep ? signalLeptons[0] : nullptr;
        // // l1 = twolep ? signalLeptons[1] : nullptr;
        // // bool samesign = twolep && (l0->q * l1->q > 0);
        // // bool pt2020 = twolep && (l0->Pt()>20.0 && l1->Pt()>20.0);
        // // return samesign && pt2020;
        return (l0 && l1);
    };
    string channel_label;
    *cutflow << CutName("same-sign channel") << [&](Superlink *sl) -> bool {
        channel_label = (string(l0->isEle() ? "e" : (l0->isMu() ? "m" : "?"))+
                         string(l1->isEle() ? "e" : (l1->isMu() ? "m" : "?")));
        cout<<"ll : "<<channel_label<<endl;
        return true;
    };

    size_t n_bjets(0), n_cljets(0), n_jets(0), n_jets50(0);
    *cutflow << CutName(">=1 b-jets>20") << [&](Superlink *sl) -> bool {
        n_bjets = cutflow->nttools().jetSelector().count_CB_jets(signalJets);
        cout<<">=1 b : "<<channel_label<<" "<<(n_bjets > 0 ? "pass":"fail")<<endl;
        return n_bjets > 0;
    };

    *cutflow << CutName("==4 jets>50") << [&](Superlink *sl) -> bool {
        n_cljets = cutflow->nttools().jetSelector().count_CL_jets(signalJets);
        n_jets = n_bjets + n_cljets;
        for(auto jet : signalJets) {
            if(jet->Pt()>50.0)
                n_jets50++;
        }
        cout<<"n_jets "<<n_jets<<" n_cljets "<<n_cljets<<" n_bjets "<<n_bjets<<" n_jets50 "<<n_jets50<<endl;
        cout<<"== 4j : "<<channel_label<<" "<<(n_jets50==4 ? "pass":"fail")<<endl;
        return n_jets50==4;
    };

    *cutflow << CutName("MET>150") << [&](Superlink *sl) -> bool {
        double met = sl->met->lv().Pt();
        cout<<"met>150 : "<<channel_label<<" "<<(met>150.0 ? "pass":"fail")<<endl;
        return (met > 150.0);
    };

    *cutflow << [&](Superlink* sl, var_void*) { baseLeptons.clear();     };
    *cutflow << [&](Superlink* sl, var_void*) { baseElectrons.clear();   };
    *cutflow << [&](Superlink* sl, var_void*) { baseMuons.clear();       };
    *cutflow << [&](Superlink* sl, var_void*) { baseJets.clear();        };
    *cutflow << [&](Superlink* sl, var_void*) { signalLeptons.clear();   };
    *cutflow << [&](Superlink* sl, var_void*) { signalElectrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { signalMuons.clear();     };

    chain->Process(cutflow, sample.c_str(), n_events, n_skip_events);

    delete chain;
    return 0;
}

