// emacs -*- C++ -*-
#ifndef SUSYNTHLFV_UTILS_H
#define SUSYNTHLFV_UTILS_H
/*
  SusyNt utility functions

  davide.gerbaudo@gmail.com
  Dec 2014
 */

namespace Susy
{
class Lepton;
}

namespace hlfv
{
/// pdg id for Susy::Lepton
struct LeptonPdg{
    /// particles have positive codes, see doi:10.1146/annurev.ns.25.120175.003011
    enum Value {
        Pel=+11, /// electron
        Ael=-11, /// positron
        Pmu=+13, /// muon
        Amu=-13, /// antimuon
        Unknown=0
    };
};

/// pdg id from Lepton attributes (isEle/isMu and charge)
int pdgIdFromLep(const Susy::Lepton &l);
} // hlvf
#endif
