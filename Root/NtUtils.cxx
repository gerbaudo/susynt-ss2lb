#include "susynt-ss2lb/NtUtils.h"

#include "SusyNtuple/SusyNt.h"

using Susy::Lepton;

namespace hlfv
{
//-----------------------------------------
int pdgIdFromLep(const Lepton &l)
{
    if     (l.isEle()) return (l.q < 0 ? LeptonPdg::Pel : LeptonPdg::Ael);
    else if(l.isMu() ) return (l.q < 0 ? LeptonPdg::Pmu : LeptonPdg::Amu);
    else                return LeptonPdg::Unknown;
}
//-----------------------------------------

} // hlfv
