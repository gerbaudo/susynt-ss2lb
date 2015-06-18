
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"

#include "SusyntHlfv/MatrixPrediction.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyntHlfv/utils.h"

using namespace std;

/**

   run the fake matrix estimate with the MatrixPrediction

*/

void usage(const char *exeName, const char *defaultMatrixFile) {
  cout<<"Usage:"<<endl
      <<exeName<<" options"<<endl
      <<"\t"<<"-m [--matrix-file] input matrix file"     <<endl
      <<"\t"<<" (default "<<defaultMatrixFile<<")"       <<endl
      <<"\t"<<"-n [--num-event]   nEvt (default -1, all)"<<endl
      <<"\t"<<"-k [--num-skip]    nSkip (default 0)"     <<endl
      <<"\t"<<"-i [--input]       (file, list, or dir)"  <<endl
      <<"\t"<<"-o [--output]      output file (required)"<<endl
      <<"\t"<<"-s [--sample]      samplename"            <<endl
      <<"\t"<<"-S [--systematics] : store also syst variations"<<endl
      <<"\t"<<"-e [--etapt]       : use eta-pt (default pt only)"<<endl
      <<"\t"<<"-d [--debug]       : debug (>0 print stuff)"<<endl
      <<"\t"<<"-h [--help]        : print help"            <<endl
      <<endl;
}

int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();
  int nEvt = -1;
  int nSkip = 0;
  int dbg = 0;
  string sample;
  string input;
  string output;
  bool etapt = false;
  bool systematics = false;
  string matrixFile(hlfv::getRootCoreDir()+"/data/DileptonMatrixMethod/FakeMatrix_Jul_26.root");
  int optind(1);
  while ((optind < argc)) {
    std::string sw = argv[optind];
    if(sw[0]!='-') { if(dbg) cout<<"skip "<<sw<<endl; optind++; continue; }
    if     (sw=="-m"||sw=="--matrix-file") { matrixFile = argv[++optind]; }
    else if(sw=="-n"||sw=="--num-event"  ) { nEvt = atoi(argv[++optind]); }
    else if(sw=="-k"||sw=="--num-skip"   ) { nSkip = atoi(argv[++optind]); }
    else if(sw=="-e"||sw=="--etapt"      ) { etapt = true;}
    else if(sw=="-d"||sw=="--debug"      ) { dbg = atoi(argv[++optind]); }
    else if(sw=="-i"||sw=="--input"      ) { input = argv[++optind]; }
    else if(sw=="-o"||sw=="--output"     ) { output = argv[++optind]; }
    else if(sw=="-s"||sw=="--sample"     ) { sample = argv[++optind]; }
    else if(sw=="-S"||sw=="--systematics") { systematics = true; }
    else if(sw=="-h"||sw=="--help"       ) { usage(argv[0], matrixFile.c_str()); return 0; }
    else cout<<"Unknown switch "<<sw<<endl;
    optind++;
  } // end while(optind<argc)
  cout<<"flags:"                <<endl
      <<"  sample  "<<sample    <<endl
      <<"  nEvt    "<<nEvt      <<endl
      <<"  nSkip   "<<nSkip     <<endl
      <<"  etapt   "<<etapt     <<endl
      <<"  dbg     "<<dbg       <<endl
      <<"  input   "<<input     <<endl
      <<"  output  "<<output    <<endl
      <<"  systematics "<<systematics<<endl
      <<endl;

  if(!output.length()){
      cout<<"--output is a required option"<<endl;
      return 1;
  }
  TChain* chain = new TChain("susyNt");
  bool verbose(dbg>0);
  ChainHelper::addInput(chain, input, verbose);
  Long64_t nEntries = chain->GetEntries();
  nEvt = (nEvt<0 ? nEntries : nEvt);
  if(dbg) chain->ls();

  hlfv::MatrixPrediction fakePred;
  fakePred.setMatrixFilename(matrixFile);
  fakePred.setDebug(dbg);
  fakePred.useComputeSystematics(systematics);
  fakePred.setSampleName(sample);
  if(etapt) fakePred.use2dParametrization();
  fakePred.setTupleFile(output);
  fakePred.buildSumwMap(chain);
  chain->Process(&fakePred, sample.c_str(), nEvt, nSkip);

  delete chain;
  return 0;
}
