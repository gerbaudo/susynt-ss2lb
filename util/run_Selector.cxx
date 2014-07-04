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
        <<"\t -n [--num-events] <int>       : number of events to be processed \n"
        <<"\t -s [--sample]  <samplename>   : sample name \n"
        <<"\t -t [--tuple-out] fname.root (out ntuple file)\n"
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
    int num_events = -1;
    int nSkip = 0;
    bool verbose=false;
    bool debug=false;
    string sample;
    string input;
    string output;
    string eventlist;
    bool write_tuple=false;
    string tuple_out;

    int opt=0;
    static struct option long_options[] = {
        {"help",       no_argument,       0, 'h'},
        {"verbose",    no_argument,       0, 'v'},
        {"debug",      no_argument,       0, 'd'},
        {"input",      required_argument, 0, 'i'},
        {"num-events", required_argument, 0, 'n'},
        {"sample",     required_argument, 0, 's'},
        {"tuple-out",  required_argument, 0, 't'},
        {"event-list", required_argument, 0, 'e'},
        {0, 0, 0, 0}
    };
    while ((opt=getopt_long (argc, argv, "hvd:i:n:e:s:t:", long_options, NULL)) != -1) {
        switch (opt) {
        case 'h' : print_usage(argv[0]); exit(0);
        case 'v' : verbose=true; break;
        case 'd' : debug=true; break;
        case 'i' : input = optarg; break;
        case 'n' : num_events = atoi(optarg); break;
        case 's' : sample = optarg; break;
        case 't' : write_tuple = true; tuple_out = optarg; break;
        case 'e' : eventlist = optarg; break;
        default: cout<<"unknown option "<<argv[optind]<<endl; exit(1);
        }
    }
    if (optind < argc) {
        cout<<"The following arguments are not expected; they will be discarded:"<<endl;
        while (optind < argc) cout<<argv[optind++]<<endl;
    }
    if(verbose)
        cout<<"verbose   '"<<verbose   <<"'"<<endl
            <<"debug     '"<<debug     <<"'"<<endl
            <<"input     '"<<input     <<"'"<<endl
            <<"num-events'"<<num_events<<"'"<<endl
            <<"sample    '"<<sample    <<"'"<<endl
            <<"tuple-out '"<<tuple_out <<"'"<<endl
            <<"eventlist '"<<eventlist <<"'"<<endl;

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
    Long64_t tot_num_events = chain->GetEntries();
    num_events = (num_events<0 ? tot_num_events : num_events);
    if(debug) chain->ls();

    hlfv::Selector selector;
    if(debug) selector.setDebug(1);
    selector.setSampleName(sample);
    cout<<"eventlist : "<<eventlist<<endl;
    selector.setEventListFilename(eventlist); // setting the event list to '' is disabling it
    if(write_tuple) selector.setTupleFile(tuple_out);
    chain->Process(&selector, sample.c_str(), num_events, nSkip);

    delete chain;
    return 0;
}
