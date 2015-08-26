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
    ChainHelper::addInput(chain, input, verbose);
    Long64_t nEntries = chain->GetEntries();
    n_events = (n_events<0 ? nEntries : n_events);
    if(dbg) chain->ls();

    SuperflowRunMode run_mode = SuperflowRunMode::nominal;
    SusyNtSys nt_sys = NtSys::NOM;
    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setAnaName("ss3l cutflow");
    cutflow->setAnaType(AnalysisType::Ana_SS3L);
    cutflow->setSelectTaus(false);
    cutflow->setSampleName(sample);
    cutflow->setRunMode(run_mode);
    cutflow->setCountWeights(true); // print the weighted cutflows
    cutflow->setLumi(lumiPeriodC);
    cutflow->setChain(chain);
    cout<<"Number of chain entries: "<<chain->GetEntries()<<endl;

    // begin cutflow
    *cutflow << CutName("read in") << [](Superlink*) -> bool { return true; };

    int cutflags = 0;

    *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
        cutflags = sl->nt->evt()->cutFlags[sl->nt_sys];
        return (cutflags & ECut_GRL);
    };

    *cutflow << CutName("pass good vertex") << [&](Superlink*) -> bool {
        return (cutflags & ECut_GoodVtx);
    };

    *cutflow << CutName("LAr error") << [&](Superlink*) -> bool {
        return (cutflags & ECut_LarErr);
    };

    *cutflow << CutName("Tile error") << [&](Superlink*) -> bool {
        return (cutflags & ECut_TileErr);
    };

    *cutflow << CutName("TTC veto") << [&](Superlink*) -> bool {
        return (cutflags & ECut_TTC);
    };

    *cutflow << CutName("bad muon veto") << [&](Superlink*) -> bool {
        return (cutflags & ECut_BadMuon);
    };

    *cutflow << CutName("jet cleaning") << [&](Superlink*) -> bool {
        return (cutflags & ECut_BadJet);
    };

    *cutflow << CutName("pass cosmic veto") << [&](Superlink*) -> bool {
        return (cutflags & ECut_Cosmic);
    };

    chain->Process(cutflow, sample.c_str(), n_events, n_skip_events);

    delete chain;
    return 0;
}

/*
//Pre-processing Def
#undef __DEBUG__
#define __DEBUG__  cout << "DEBUG " << __FILE__ << ":" << __FUNCTION__ << " line: " << __LINE__ << endl;

#define __INTLUMI__ 55.4  //Period C

// Functions
void usage(std::string progName="makeSS_miniNt");
void exec_options(int argc, char* argv[], TChain* chain, int& num_events, string& sample,
                  SuperflowRunMode& run_mode, bool& fake_mode, SusyNtSys& nt_sysnt);


//========================================//
// Fixe fake_mode - ain't working w/ -b
//========================================//

///////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    int n_events  = -1;
    int n_skip_events  = 0;
    string sample;

    SuperflowRunMode run_mode = SuperflowRunMode::nominal;
    bool fake_mode = false;
    SusyNtSys nt_sys = NtSys::NOM;
    TChain* chain = new TChain("susyNt");
    chain->SetDirectory(0);

    ////////////////////////////////////////////////////////////////////////
    // Read command-line options to excecutable
    ////////////////////////////////////////////////////////////////////////
    exec_options(argc, argv, chain, n_events, sample, run_mode, fake_mode, nt_sys);

    ////////////////////////////////////////////////////////////////////////
    // Initialize and configure the analysis
    ////////////////////////////////////////////////////////////////////////
    Superflow* cutflow = new Superflow(); // initialize the cutflow
    cutflow->setAnaName("SS_miniNt");
    cutflow->setAnaType(AnalysisType::Ana_SS3L);
    cutflow->setSelectTaus(false);
    cutflow->setSampleName(sample);
    cutflow->setRunMode(run_mode);
    cutflow->setCountWeights(true); // print the weighted cutflows
    cutflow->setLumi(__INTLUMI__);  // set the MC normalized to lumi
    cutflow->setChain(chain);
    if (run_mode == SuperflowRunMode::single_event_syst) cutflow->setSingleEventSyst(nt_sys);
    cout << "Analysis    Total Entries: " << chain->GetEntries() << endl;


    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //
    //  Superflow methods [BEGIN]
    //
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    *cutflow << CutName("read in") << [](Superlink*) -> bool { return true; };

    ////////////////////////////////////////////////////////////
    //  Cleaning Cuts
    ////////////////////////////////////////////////////////////
    int cutflags = 0;

    *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
        cutflags = sl->nt->evt()->cutFlags[sl->nt_sys];
        return (cutflags & ECut_GRL);
    };

    *cutflow << CutName("LAr error") << [&](Superlink*) -> bool {
        return (cutflags & ECut_LarErr);
    };

    *cutflow << CutName("Tile error") << [&](Superlink*) -> bool {
        return (cutflags & ECut_TileErr);
    };

    *cutflow << CutName("TTC veto") << [&](Superlink*) -> bool {
        return (cutflags & ECut_TTC);
    };

    *cutflow << CutName("bad muon veto") << [&](Superlink*) -> bool {
        return (cutflags & ECut_BadMuon);
    };

    *cutflow << CutName("jet cleaning") << [&](Superlink*) -> bool {
        return (cutflags & ECut_BadJet);
    };

    *cutflow << CutName("pass good vertex") << [&](Superlink*) -> bool {
        return (cutflags & ECut_GoodVtx);
    };

    *cutflow << CutName("pass cosmic veto") << [&](Superlink*) -> bool {
        return (cutflags & ECut_Cosmic);
    };


    ////////////////////////////////////////////////////////////
    //  Analysis Cuts
    ////////////////////////////////////////////////////////////

    // *cutflow << CutName("tau veto") << [](Superlink* sl) -> bool {
    //        return sl->taus->size() == 0;
    //    };
    *cutflow << CutName("at least 2 base lepton") << [](Superlink* sl) -> bool {
        return sl->baseLeptons->size() >= 1;
    };

    *cutflow << CutName("at least 2 signal lepton") << [](Superlink* sl) -> bool {
        return sl->leptons->size() >= 1;
    };
    // *cutflow << CutName("at least 2 signal lepton pt>20") << [](Superlink* sl) -> bool {
    //     return sl->tools->hasNLeptonsPtThreshold(*sl->leptons, 2, 20.);
    // };

    ////////////////////////////////////////////////////////////
    //  Output Ntuple Setup
    //      > Ntuple variables
    ////////////////////////////////////////////////////////////

    //
    // Event Vars
    //
    *cutflow << NewVar("Event run number"); {
        *cutflow << HFTname("runNumber");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->run; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Event number"); {
        *cutflow << HFTname("eventNumber");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->eventNumber; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("is Monte Carlo"); {
        *cutflow << HFTname("isMC");
        *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("mcChannel (dsid)"); {
        *cutflow << HFTname("dsid");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->isMC ? sl->nt->evt()->mcChannel : 0;};
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("event weight"); {
        *cutflow << HFTname("eventweight");
        *cutflow << [](Superlink* sl, var_double*) -> double { return sl->weights->product(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Number of good vertices"); {
        *cutflow << HFTname("nVtx");
        *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->nVtx; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("Average dnumber of interactions"); {
        *cutflow << HFTname("avgMu");
        *cutflow << [](Superlink* sl, var_float*) -> float { return sl->nt->evt()->avgMu; };
        *cutflow << SaveVar();
    }

    //
    // Triggers
    //


  // ee HLT_2e12_loose_L12EM10VH HLT_2e12_lhloose_L12EM10VH
  // mm HLT_mu18_mu8noL1
  // em HLT_e17_loose_mu14 HLT_e17_lhloose_mu14
  // Met HLT_xe100



    Trigger* trigtool = new Trigger(chain, true);
    TBits triggerBits;
    *cutflow << NewVar("pass trigger 2e12_loose"); {
        *cutflow << HFTname("trig_2e12_loose");
        *cutflow << [&](Superlink* sl, var_bool*) -> bool { return trigtool->passTrigger(triggerBits, "HLT_2e12_loose_L12EM10VH");};
        *cutflow << SaveVar();
    }



    //
    // Leptons
    //
    LeptonVector   baseLeptons;
    ElectronVector baseElectrons;
    MuonVector     baseMuons;
    LeptonVector   signalLeptons;
    ElectronVector signalElectrons;
    MuonVector     signalMuons;

    //baseline or signal depending on running mode - default is signal
    LeptonVector   storedLeptons;
    ElectronVector storedElectrons;
    MuonVector     storedMuons;

    LeptonVector   SSLeptons;//highest leading pT SS pair of stored leptons

    *cutflow << [&](Superlink* sl, var_void*) { baseLeptons     = *sl->baseLeptons; };
    *cutflow << [&](Superlink* sl, var_void*) { baseElectrons   = *sl->baseElectrons; };
    *cutflow << [&](Superlink* sl, var_void*) { baseMuons       = *sl->baseMuons; };
    *cutflow << [&](Superlink* sl, var_void*) { signalLeptons   = *sl->leptons; };
    *cutflow << [&](Superlink* sl, var_void*) { signalElectrons = *sl->electrons; };
    *cutflow << [&](Superlink* sl, var_void*) { signalMuons     = *sl->muons; };

    *cutflow << [&](Superlink* sl, var_void*) { storedLeptons   = fake_mode ? baseLeptons : signalLeptons; };
    *cutflow << [&](Superlink* sl, var_void*) { storedElectrons = fake_mode ? baseElectrons : signalElectrons; };
    *cutflow << [&](Superlink* sl, var_void*) { storedMuons     = fake_mode ? baseMuons : signalMuons; };

    *cutflow << [&](Superlink* sl, var_void*) { SSLeptons = sl->tools->getSSLeptonPair(storedLeptons); };

    *cutflow << NewVar("number of stored leptons"); {
        *cutflow << HFTname("nLep");
        *cutflow << [&](Superlink* sl, var_int*) -> int { return storedLeptons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of stored Electrons"); {
        *cutflow << HFTname("nEle");
        *cutflow << [&](Superlink* sl, var_int*) -> int {return storedElectrons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of stored Muons"); {
        *cutflow << HFTname("nMuo");
        *cutflow <<[&](Superlink* sl, var_int*) -> int {return storedMuons.size(); };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("has 2 leptons pT>20"); {
        *cutflow << HFTname("nLepPt20");
        *cutflow <<[&](Superlink* sl, var_bool*) -> int {return sl->tools->hasNLeptonsPtThreshold(*sl->leptons, 2, 20.); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("has SS "); {
        *cutflow << HFTname("hasSS");
        *cutflow <<[&](Superlink* sl, var_bool*) -> int {return SSLeptons.size()>0; };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("SS lepton flavor (0: ee, 1:em/me, 2: mm, 3 unknown)"); {
        *cutflow << HFTname("ss_flav");
        *cutflow <<[&](Superlink* sl, var_int*) -> int {
            int type=3;
            if(SSLeptons.size()>1){
                DiLepEvtType dil_type = getDiLepEvtType(SSLeptons);
                if(dil_type == ET_ee ) type=0;
                if(dil_type == ET_me || dil_type == ET_em) type=1;
                if(dil_type == ET_mm ) type=2;
            }
            return type;
        };
        *cutflow << SaveVar();
    }


    *cutflow << NewVar("lepton flavor (0: e, 1: m)"); {
        *cutflow << HFTname("l_flav");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> lep_flav;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                lep_flav.push_back(storedLeptons.at(i)->isEle() ? 0 : 1);
            }
            return lep_flav;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton pt"); {
        *cutflow << HFTname("l_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> lep_pt;
            for(uint i = 0; i< storedLeptons.size(); i++) {
                lep_pt.push_back(storedLeptons.at(i)->Pt());
            }
            return lep_pt;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton eta"); {
        *cutflow << HFTname("l_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> lep_eta;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                lep_eta.push_back(storedLeptons.at(i)->Eta());
            }
            return lep_eta;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton phi"); {
        *cutflow << HFTname("l_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> lep_phi;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                lep_phi.push_back(storedLeptons.at(i)->Phi());
            }
            return lep_phi;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton d0sig (BSCorr)"); {
        *cutflow << HFTname("l_d0sig");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> d0sigBSCorr;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                d0sigBSCorr.push_back(storedLeptons.at(i)->d0sigBSCorr);
            }
            return d0sigBSCorr;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton z0"); {
        *cutflow << HFTname("l_z0");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                out.push_back(storedLeptons.at(i)->z0);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton z0sinTheta"); {
        *cutflow << HFTname("l_z0sinTheta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> z0;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                z0.push_back(storedLeptons.at(i)->z0SinTheta());
            }
            return z0;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton charge"); {
        *cutflow << HFTname("l_q");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                out.push_back(storedLeptons.at(i)->q);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton is iso gradient loose"); {
        *cutflow << HFTname("l_isGradientLoose");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                if(storedLeptons.at(i)->isEle()) out.push_back(((Electron*) storedLeptons.at(i))->isoGradientLoose);
                else                           out.push_back(((Muon*) storedLeptons.at(i))->isoGradientLoose);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton ptvarcone20"); {
        *cutflow << HFTname("l_ptvarcone20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                out.push_back(storedLeptons.at(i)->ptvarcone20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton ptvarcone30"); {
        *cutflow << HFTname("l_ptvarcone30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                out.push_back(storedLeptons.at(i)->ptvarcone30);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton etconetopo20"); {
        *cutflow << HFTname("l_etconetopo20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                out.push_back(storedLeptons.at(i)->etconetopo20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton etconetopo30"); {
        *cutflow << HFTname("l_etconetopo30");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                out.push_back(storedLeptons.at(i)->etconetopo30);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("lepton is signal (SUSYTools 'signal' flag)"); {
        *cutflow << HFTname("l_STsignal");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                out.push_back(storedLeptons.at(i)->isSignal);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton is loose"); {
        *cutflow << HFTname("l_isLoose");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                if(storedLeptons.at(i)->isEle()) out.push_back(((Electron*) storedLeptons.at(i))->looseLH);
                else                       out.push_back(((Muon*) storedLeptons.at(i))->loose);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton is medium"); {
        *cutflow << HFTname("l_isMedium");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                if(storedLeptons.at(i)->isEle()) out.push_back(((Electron*) storedLeptons.at(i))->mediumLH);
                else                       out.push_back(((Muon*) storedLeptons.at(i))->medium);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton is tight"); {
        *cutflow << HFTname("l_isTight");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                if(storedLeptons.at(i)->isEle()) out.push_back(((Electron*) storedLeptons.at(i))->tightLH);
                else                       out.push_back(((Muon*) storedLeptons.at(i))->tight);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton is loose no d0"); {
        *cutflow << HFTname("l_isLoose_noD0");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                if(storedLeptons.at(i)->isEle()) out.push_back(((Electron*) storedLeptons.at(i))->looseLH_nod0);
                else                       out.push_back(false);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton is medium no d0"); {
        *cutflow << HFTname("l_isMedium_noD0");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                if(storedLeptons.at(i)->isEle()) out.push_back(((Electron*) storedLeptons.at(i))->mediumLH_nod0);
                else                       out.push_back(false);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton is tight no d0"); {
        *cutflow << HFTname("l_isTight_noD0");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(uint i = 0; i < storedLeptons.size(); i++) {
                if(storedLeptons.at(i)->isEle()) out.push_back(((Electron*) storedLeptons.at(i))->tightLH_nod0);
                else                       out.push_back(false);
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    //
    // Jets
    //
    JetVector sjets; //signal jets
    JetVector bjets;

    *cutflow << [&](Superlink* sl, var_void*) {
        JetVector susyJets = *sl->jets;
        for(int i = 0; i < susyJets.size(); i++) {
            Jet* j = susyJets.at(i);
            sjets.push_back(j);
            if(sl->tools->m_jetSelector.isCentralBJet(sjets.at(i))) { bjets.push_back(j); }
        }
    };


    //Signal Jets
    *cutflow << NewVar("number of sjets"); {
        *cutflow << HFTname("nSJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return sjets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet pt"); {
        *cutflow << HFTname("sj_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Pt());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet eta"); {
        *cutflow << HFTname("sj_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Eta());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet phi"); {
        *cutflow << HFTname("sj_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->Phi());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet jvt"); {
        *cutflow << HFTname("sj_jvt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->jvt);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet mv2c20"); {
        *cutflow << HFTname("sj_mv2c20");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < sjets.size(); i++) {
                out.push_back(sjets.at(i)->mv2c20);
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("sjet isBJet"); {
        *cutflow << HFTname("sj_isBJet");
        *cutflow << [&](Superlink* sl, var_bool_array*) -> vector<bool> {
            vector<bool> out;
            for(int i = 0; i < sjets.size(); i++) {
                if(sl->tools->m_jetSelector.isCentralBJet(sjets.at(i))) { out.push_back(true); }
                else { out.push_back(false); }
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    //b-jets
    *cutflow << NewVar("number of bjets"); {
        *cutflow << HFTname("nBJets");
        *cutflow << [&](Superlink* sl, var_int*) -> int {
            return bjets.size();
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet pt"); {
        *cutflow << HFTname("bj_pt");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Pt());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet eta"); {
        *cutflow << HFTname("bj_eta");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Eta());
            }
            return out;
        };
        *cutflow << SaveVar();
    }
    *cutflow << NewVar("bjet phi"); {
        *cutflow << HFTname("bj_phi");
        *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
            vector<double> out;
            for(int i = 0; i < bjets.size(); i++) {
                out.push_back(bjets.at(i)->Phi());
            }
            return out;
        };
        *cutflow << SaveVar();
    }

    //
    // Met
    //
    Met met;

    *cutflow << [&](Superlink* sl, var_void*) { met = *sl->met; };
    *cutflow << NewVar("transverse missing energy (Etmiss)"); {
        *cutflow << HFTname("met");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.lv().Pt(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("phi coord. of Etmiss"); {
        *cutflow << HFTname("metPhi");
        *cutflow << [&](Superlink* sl, var_float*) -> double { return met.lv().Phi(); };
        *cutflow << SaveVar();
    }

    *cutflow << NewVar("sumet"); {
        *cutflow << HFTname("sumet");
        *cutflow << [&](Superlink* sl, var_float*) -> double {
            return met.sumet;
        };
        *cutflow << SaveVar();
    }


    //
    // Cleanup container
    //

    *cutflow << [&](Superlink* sl, var_void*) { baseLeptons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { baseElectrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { baseMuons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { signalLeptons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { signalElectrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { signalMuons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { storedLeptons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { storedElectrons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { storedMuons.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { SSLeptons.clear(); };
    // *cutflow << [&](Superlink* sl, var_void*) { jets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { sjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { bjets.clear(); };
    *cutflow << [&](Superlink* sl, var_void*) { met.clear(); };



    ////////////////////////////////////////////////////////////////////////
    // Initialize the cutflow and start the event loop.
    ////////////////////////////////////////////////////////////////////////

    chain->Process(cutflow, sample.c_str(), n_events, n_skip_events);
    //delete cutflow;
    //delete chain;

    cout << "Done." << endl;
    exit(0);

}

///////////////////////////////////////////////////////////////////////
// Usage
///////////////////////////////////////////////////////////////////////
void usage(std::string progName)
{
    printf("=================================================================\n");
    printf("%s [options]\n",progName.c_str());
    printf("=================================================================\n");
    printf("Options:\n");
    printf("-h        Print this help\n");
    printf("-n        Number of events to be processed (default: -1)\n");
    printf("-f        Input file as *.root, list of *.root in a *.txt,\n");
    printf("          or a DIR/ containing *.root (default: none)\n");
    printf("-c        Nominal weighting\n");
    printf("-e        Nomianl and single systematic \n");
    printf("-a        All systemactics\n");
    printf("-b        Do ana with baseline for fake estimate (default:false)\n");
    printf("-i        Input file list (.root or .txt or /)\n");
    printf("-s        systematic name \n");
    printf("=================================================================\n");
}

///////////////////////////////////////////////////////////////////////
// exec_options
///////////////////////////////////////////////////////////////////////
void exec_options(int argc, char* argv[], TChain* chain, int& num_events, string& sample,
                  SuperflowRunMode& run_mode, bool& fake_mode, SusyNtSys& nt_sys)
{
    bool nominal = false;
    bool nominal_and_weight_syst = false;
    bool all_syst = false;
    bool single_event_syst = false;

    string systematic = "undefined";
    string input;

    // Read inputs to program
    opterr = 0;
    int c;
    while ((c = getopt (argc, argv, "n:c:e:a:b:i:s:f:h")) != -1)
        switch (c){
        case 'n':
            num_events = atoi(optarg);
            break;
        case 'c':
            nominal = true;
            break;
        case 'e':
            single_event_syst = true;
            break;
        case 'a':
            all_syst = true;
            break;
        case 'b':
            fake_mode = true;
            break;
        case 'i':
            input = string(optarg);
            break;
        case 's':
            systematic = optarg;
            break;
        case 'f':
            input = optarg;
            break;
        case 'h':
            usage();
            break;
        default:
            break;
        }
    //--------------//

    // Catch problems or cast
    for (int index = optind; index < argc; index++)
        printf ("make_miniNt \t Non-option argument %s\n", argv[index]);
    if (input.size()==0) {
        printf("make_miniNt\t An input file must be provided with option -f (a list, a DIR or single file)\n");
        exit(EXIT_FAILURE);
    }


    bool inputIsFile = Susy::utils::endswith(input, ".root");
    bool inputIsList = Susy::utils::endswith(input, ".txt");
    bool inputIsDir  = Susy::utils::endswith(input, "/");
    bool validInput(inputIsFile || inputIsList || inputIsDir);
    if (!validInput) {
        cout << "Analysis    invalid input '" << input << "'" << endl;
        exit(EXIT_FAILURE);
    }
    if (inputIsFile) {
        ChainHelper::addFile(chain, input);
        cout << "Analysis    file: " << input << endl;
        cout << "Analysis    file: " << input << endl;
        cout << "Analysis    file: " << input << endl;
        sample = input;
    }
    if (inputIsList) {
        ChainHelper::addFileList(chain, input);
        cout << "Analysis    list: " << input << endl;
        cout << "Analysis    list: " << input << endl;
        cout << "Analysis    list: " << input << endl;
        ifstream infile(input.c_str());
        if (infile.good()) {
            string sLine;
            getline(infile, sLine);
            sample = sLine;
        }
        else {
            sample = input;
        }
        infile.close();
    }
    if (inputIsDir) {
        ChainHelper::addFileDir(chain, input);
        cout << "Analysis    dir: " << input << endl;
        cout << "Analysis    dir: " << input << endl;
        cout << "Analysis    dir: " << input << endl;
        sample = input;
    }
    Long64_t tot_num_events = chain->GetEntries();
    num_events = (num_events < 0 ? tot_num_events : num_events);
    // if (debug) chain->ls();

    if (nominal) {
        run_mode = SuperflowRunMode::nominal;
        cout << "Analysis    run mode: SuperflowRunMode::nominal" << endl;
    }
    if (nominal_and_weight_syst) {
        run_mode = SuperflowRunMode::nominal_and_weight_syst;
        cout << "Analysis    run mode: SuperflowRunMode::nominal_and_weight_syst" << endl;
    }
    if (single_event_syst) {
        run_mode = SuperflowRunMode::single_event_syst;
        cout << "Analysis    run mode: SuperflowRunMode::single_event_syst" << endl;
    }

    if (all_syst) {
        run_mode = SuperflowRunMode::all_syst;
        cout << "Analysis    run mode: SuperflowRunMode::all_syst" << endl;
    }
    if (fake_mode) {
        cout << "Analysis    fake mode TRUE " << endl;
    }
    else cout << "Analysis    fake mode FALSE " << endl;

    map <string, SusyNtSys> event_syst_map;
      // event_syst_map["EESZUP"] = NtSys::EES_Z_UP;
      // event_syst_map["EESZDOWN"] = NtSys::EES_Z_DN;
      // event_syst_map["EESMATUP"] = NtSys::EES_MAT_UP;
      // event_syst_map["EESMATDOWN"] = NtSys::EES_MAT_DN;
      // event_syst_map["EESPSUP"] = NtSys::EES_PS_UP;
      // event_syst_map["EESPSDOWN"] = NtSys::EES_PS_DN;
      // event_syst_map["EESLOWUP"] = NtSys::EES_LOW_UP;
      // event_syst_map["EESLOWDOWN"] = NtSys::EES_LOW_DN;
      // event_syst_map["EERUP"] = NtSys::EER_UP;
      // event_syst_map["EERDOWN"] = NtSys::EER_DN;
      // event_syst_map["MSUP"] = NtSys::MS_UP;
      // event_syst_map["MSDOWN"] = NtSys::MS_DN;
      // event_syst_map["IDUP"] = NtSys::ID_UP;
      // event_syst_map["IDDOWN"] = NtSys::ID_DN;
      // event_syst_map["JESUP"] = NtSys::JES_UP;
      // event_syst_map["JESDOWN"] = NtSys::JES_DN;
      // event_syst_map["JER"] = NtSys::JER;
      // event_syst_map["SCALESTUP"] = NtSys::SCALEST_UP;
      // event_syst_map["SCALESTDOWN"] = NtSys::SCALEST_DN;
      // event_syst_map["RESOST"] = NtSys::RESOST;
      // event_syst_map["TRIGSFELUP"] = NtSys::TRIGSF_EL_UP;
      // event_syst_map["TRIGSFELDN"] = NtSys::TRIGSF_EL_DN;
      // event_syst_map["TRIGSFMUUP"] = NtSys::TRIGSF_MU_UP;
      // event_syst_map["TRIGSFMUDN"] = NtSys::TRIGSF_MU_DN;
      // event_syst_map["TESUP"] = NtSys::TES_UP;
      // event_syst_map["TESDOWN"] = NtSys::TES_DN;
      // event_syst_map["JVFUP"] = NtSys::JVF_UP;
      // event_syst_map["JVFDOWN"] = NtSys::JVF_DN;
      if (single_event_syst) {
        if (event_syst_map.count(systematic) == 1) {
            nt_sys = event_syst_map[systematic];
        }
        else {
            cout << "Analysis" << "    ERROR (fatal): Event systematic option /s "
                 << systematic << " -> not found." << endl;
            exit(EXIT_FAILURE);
        }
    }
}

*/
