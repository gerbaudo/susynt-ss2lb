// emacs -*- C++ -*-
#ifndef SS3L_LEPTONTRUTHTYPE_H
#define SS3L_LEPTONTRUTHTYPE_H

#include <string>

namespace Susy {
class Lepton;
}

namespace ss3l {

/// Possible lepton sources from simulation truth
/**
   These are the values from

   @RecoTruthMatch::type
   from
   svn+ssh://svn.cern.ch/reps/atlasphys/Physics/SUSY/Analyses/WeakProduction/LeptonTruthTools/tags/LeptonTruthTools-00-01-06

   Note to self 1.
   Enclose enum in a struct to avoid collisions; see
   explanation at
   <A HREF="http://stackoverflow.com/questions/7090130/enum-in-a-namespace">stackoveflow</A>.

   Note to self 2.
   reasons for duplicating this enum:
   (1) the original enum is defined inside the reconstruction class
   (2) we don't need all the vlaues from the original enum, only the
       ones returned by RecoTruthMatch::fakeType(), see
       SusyNtMaker::fillElectronVars(), SusyNtMaker::fillMuonVars()
   (3) we want to define a few additional converter/repr

   davide.gerbaudo@gmail.com
   Oct 2014
 */
struct LeptonTruthType{
    enum Value {
        Prompt,
        Conversion,
        HeavyFlavor,
        LightFlavor,
        Unknown
    };
    inline static LeptonTruthType::Value first() { return LeptonTruthType::Prompt; }
    inline static LeptonTruthType::Value last() { return LeptonTruthType::Unknown; }
};

inline bool isValid(const LeptonTruthType::Value &v) { return v>=LeptonTruthType::first() && v<=LeptonTruthType::last(); }
inline bool isValid(const int &v) { return isValid(static_cast<LeptonTruthType::Value>(v)); }
/// safe conversion; out-of-range values are turned into unknown
LeptonTruthType::Value int2source(const int &v);
/// string representation
std::string source2string(const LeptonTruthType::Value &v);

/// helper functions
LeptonTruthType::Value getLeptonSource(const Susy::Lepton &l);

} // ss3l
#endif
