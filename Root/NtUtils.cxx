#include "susynt-ss3l/NtUtils.h"

#include "SusyNtuple/SusyNt.h"

using Susy::Lepton;

namespace ss3l
{
//-----------------------------------------
int pdgIdFromLep(const Lepton &l)
{
    if     (l.isEle()) return (l.q < 0 ? LeptonPdg::Pel : LeptonPdg::Ael);
    else if(l.isMu() ) return (l.q < 0 ? LeptonPdg::Pmu : LeptonPdg::Amu);
    else                return LeptonPdg::Unknown;
}
//-----------------------------------------

} // ss3l
