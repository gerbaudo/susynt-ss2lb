/**
  Unit test for NtUtils.h

  davide.gerbaudo@gmail.com
  Dec 2014
*/

#include "SusyntHlfv/NtUtils.h"
#include "SusyNtuple/SusyNt.h"

#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using Susy::Electron;
using Susy::Lepton;
using Susy::Muon;

using namespace hlfv;

struct validate_enum_value{
    bool operator()(int value, int expected, const char* msg){
        bool pass= value==expected;
        if(!pass)
            cout<<msg<<": got "<<value<<", expected "<<expected<<endl;
        return pass;
    }
} validate;

int main(int argc, char** argv)
{
    bool pass_all=true;
    { // pdg functions
        bool pass_pdg=true;
        Muon muplus;      muplus.q  = +1.0;
        Muon muminus;     muminus.q = -1.0;
        Electron elplus;  elplus.q  = +1.0;
        Electron elminus; elminus.q = -1.0;
        pass_pdg &= validate(pdgIdFromLep(elplus),  LeptonPdg::Ael, "positron");
        pass_pdg &= validate(pdgIdFromLep(elminus), LeptonPdg::Pel, "electron");
        pass_pdg &= validate(pdgIdFromLep(muplus),  LeptonPdg::Amu, "antimuon");
        pass_pdg &= validate(pdgIdFromLep(muminus), LeptonPdg::Pmu, "muon");
        if(!pass_pdg)
            cout<<"failed some of the pdg tests"<<endl;
        pass_all &= pass_pdg;
    }
    cout<<endl
        <<(pass_all ? "passed all" : "failed some")
        <<endl;
    return 0;
}
