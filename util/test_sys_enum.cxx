/*
  Unit test for Systematic.h (validate names and self-consistentcy)

  davide.gerbaudo@gmail.com
  Feb 2014
*/

#include <iostream>
#include <string>

#include "susynt-ss3l/Systematic.h"
#include "SusyNtuple/SusyDefs.h"

using std::cout;
using std::endl;
using std::string;
using ss3l::Systematic;
using Susy::NtSys::SusyNtSys;

int main(int argc, char** argv)
{
    size_t failing=0;
    cout<<"\tSusyNtSys\tswh::Systematic"<<endl;
    for(int i=0; i<Susy::NtSys::SYS_UNKNOWN; ++i){
        SusyNtSys sNt = static_cast<SusyNtSys>(i);
        Systematic::Value sWh = ss3l::ntsys2sys(sNt);
        SusyNtSys sNt2 = ss3l::sys2ntsys(sWh);
        if(sNt!=sNt2) failing++;
        cout<<"["<<i<<"]"
            <<" : "<<ss3l::syst2str(sNt)
            <<"\t -> "<<ss3l::syst2str(sWh)
            <<" ; "<<sNt<<" : "<<sNt2
            <<(sNt==sNt2 ? "" : "\t <<<<<<<")
            <<endl;
    }
    if(failing) cout<<"\nFailing "<<failing<<" cases."<<endl;
    else        cout<<"All fine."<<endl;
    return 0;
}
