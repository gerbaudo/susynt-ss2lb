/**

   Unit test and example usage for hlfv::CutFlowCounter

   davide.gerbaudo@gmail.com
   June 2014
*/


#include "susynt-ss2lb/CutFlowCounter.h"

#include <iostream>
using namespace std;

//----------------------------------------------------------
int main(int argc, char** argv) {
    
    hlfv::CutFlowCounter counter;
    counter.setDebug(true);

    for(size_t i=0; i<10; ++i){
        cout<<i<<")"<<endl;
        counter.increment(0.1, "cut1");
        if(i%2==0) counter.incrementRaw("cut2"); // test raw counter only
        if(i%3==0) counter.increment(0.1, "cut3"); // test raw+weighted counter
    }
    cout<<"Summary of counts: "<<endl;
    size_t nCuts = counter.rawCounts().size();
    string header = "Cut name\traw\tweighted";
    cout<<header<<endl;
    for(size_t iCut=0; iCut<nCuts; ++iCut)
        cout<<counter.cutNames()[iCut]<<":\t"<<counter.rawCounts()[iCut]<<"\t"<<counter.weightedCounts()[iCut]<<endl;
    return 0;
}
