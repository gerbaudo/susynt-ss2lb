/**

   Main executable for hlfv::Selector

   davide.gerbaudo@gmail.com
   June 2014
*/

#include "SusyntHlfv/Selector.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"

#include "TChain.h"
#include "Cintex/Cintex.h"


#include <cstdlib>
#include <string>
#include <getopt.h>

using namespace std;

//----------------------------------------------------------
void print_usage(const char *exeName) {
    cout<<"Usage: "<<endl
        <<exeName<<" options"<<endl
        <<"\t -h [--help]                   : print this message \n"
        <<"\t -v [--verbose]                : toggle verbose \n"
        <<"\t -d [--debug]                  : toggle debug \n"
        <<"\t -i [--input]   <file.txt>     : input root file (or filelist or dir) \n"
        <<"\t -s [--sample]  <samplename>   : sample name \n"
        <<"\t -e [--event-list] <file.root> : file where the eventlist is cached \n"
        <<"Example: \n"
        <<exeName<<"\n"
        <<" -i filelist/Sherpa_CT10_lllnu_WZ_MassiveCB.txt \n"
        <<" -s Sherpa_CT10_lllnu_WZ_MassiveCB \n"
        <<" -e out/selection/eventlist_files \n"
        <<endl;
}
//----------------------------------------------------------
int main(int argc, char** argv) {
    ROOT::Cintex::Cintex::Enable();
    int nEvt = -1;
    int nSkip = 0;
    bool verbose=false;
    bool debug=false;
    string sample;
    string input;
    string output;
    string eventlist;

    int opt=0;
    static struct option long_options[] = {
        {"help",       no_argument,       0, 'h'},
        {"verbose",    no_argument,       0, 'v'},
        {"debug",      no_argument,       0, 'd'},
        {"input",      required_argument, 0, 'i'},
        {"sample",     required_argument, 0, 's'},
        {"event-list", required_argument, 0, 'e'},
        {0, 0, 0, 0}
    };
    while ((opt=getopt_long (argc, argv, "hvd:i:e:s:", long_options, NULL)) != -1) {
        switch (opt) {
        case 'h' : print_usage(argv[0]); exit(0);
        case 'v' : verbose=true; break;
        case 'd' : debug=true; break;
        case 'i' : input = optarg; break;
        case 's' : sample = optarg; break;
        case 'e' : eventlist = optarg; break;
        default: cout<<"unknown option "<<argv[optind]<<endl; exit(1);
        }
    }
    if (optind < argc) {
        cout<<"The following arguments are not expected; they will be discarded:"<<endl;
        while (optind < argc) cout<<argv[optind++]<<endl;
    }
    if(verbose)
        cout<<"verbose   '"<<verbose  <<"'"<<endl
            <<"debug     '"<<debug    <<"'"<<endl
            <<"input     '"<<input    <<"'"<<endl
            <<"sample    '"<<sample   <<"'"<<endl
            <<"eventlist '"<<eventlist<<"'"<<endl;

    TChain* chain = new TChain("susyNt");
    bool inputIsFile = susy::utils::endswith(input, ".root");
    bool inputIsList = susy::utils::endswith(input, ".txt");
    bool inputIsDir  = susy::utils::endswith(input, "/");
    bool validInput(inputIsFile||inputIsList||inputIsDir);
    if(!validInput) {
        cout<<"invalid input '"<<input<<"'"<<endl;
        print_usage(argv[0]); return 1;
    }
    if(!sample.size()) {
        cout<<" sample '"<<sample<<"'"<<endl;
        print_usage(argv[0]); return 1;
    }
    if(inputIsFile) ChainHelper::addFile    (chain, input);
    if(inputIsList) ChainHelper::addFileList(chain, input);
    if(inputIsDir ) ChainHelper::addFileDir (chain, input);
    Long64_t nEntries = chain->GetEntries();
    nEvt = (nEvt<0 ? nEntries : nEvt);
    if(debug) chain->ls();

    hlfv::Selector selector;
    if(debug) selector.setDebug(1);
    selector.setSampleName(sample);
    cout<<"eventlist : "<<eventlist<<endl;
    selector.setEventListFilename(eventlist); // setting the event list to '' is disabling it
    chain->Process(&selector, sample.c_str(), nEvt, nSkip);

    delete chain;
    return 0;
}
