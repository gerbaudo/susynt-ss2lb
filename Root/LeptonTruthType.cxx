#include "SusyntHlfv/LeptonTruthType.h"

using hlfv::LeptonTruthType;
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
