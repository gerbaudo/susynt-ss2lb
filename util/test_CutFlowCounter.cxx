/**

   Unit test and example usage for hlfv::CutFlowCounter

   davide.gerbaudo@gmail.com
   June 2014
*/


#include "SusyntHlfv/CutFlowCounter.h"

#include <iostream>
using namespace std;

// end should go in SusyNtuple/vec_utils
// from http://stackoverflow.com/questions/4268886/initialize-a-vector-array-of-strings
template<typename T, size_t N>
T * end(T (&ra)[N]) {
    return ra + N;
}
//----------------------------------------------------------
int main(int argc, char** argv) {
    
    const string tmp_cutnames[] = {"FirstCut", "SecondCut", "ThirdCut"}; // just needed to initialize a vec w/out push_back
    vector<string> cutnames(tmp_cutnames, end(tmp_cutnames));
    hlfv::CutFlowCounter counter(cutnames);
    counter.setDebug(true);

    for(size_t i=0; i<10; ++i){
        counter.nextEvent();
        counter.pass(0.1);
        if(i%2) counter.pass();
        else    counter.fail();
        if(i%3) {
            counter.pass(0.1); // test weighted counter
            counter.pass(); // test automatic extension of internal vectors
        }
    }
    cout<<"Summary of counts: "<<endl;
    size_t nCuts = counter.rawCounts().size();
    string header = "Cut name\traw\tweighted";
    cout<<header<<endl;
    for(size_t iCut=0; iCut<nCuts; ++iCut)
        cout<<counter.cutNames()[iCut]<<":\t"<<counter.rawCounts()[iCut]<<"\t"<<counter.weightedCounts()[iCut]<<endl;
    return 0;
}
