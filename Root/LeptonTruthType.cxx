#include "SusyntHlfv/LeptonTruthType.h"

#include "SusyNtuple/SusyNt.h"

using Susy::Lepton;
using hlfv::LeptonTruthType;

namespace hlfv
{
//------------------------------------------------------------------------------
std::string source2string(const LeptonTruthType::Value &v)
{
    std::string s;
    switch(v) {
    case LeptonTruthType::Prompt      : s = "Prompt"     ; break;
    case LeptonTruthType::Conversion  : s = "Conversion" ; break;
    case LeptonTruthType::HeavyFlavor : s = "HeavyFlavor"; break;
    case LeptonTruthType::LightFlavor : s = "LightFlavor"; break;
    case LeptonTruthType::Unknown     : s = "Unknown"    ; break;
    // no default, so that the compiler will warn us of un-handled cases
    }
    return s;
}
//------------------------------------------------------------------------------
LeptonTruthType::Value getLeptonSource(const Susy::Lepton &l)
{
    // note to self (DG 2015-01-29)
    // Matt used to have a special treatment for some of the samples,
    // see MeasureFakeRate2::isRealLepton(). However, that implementation
    // is not clear, it might be outdated, and it has several question marks.
    // For now try to keep it simple.
    return int2source(l.truthType);
}
//------------------------------------------------------------------------------

} // hlfv
