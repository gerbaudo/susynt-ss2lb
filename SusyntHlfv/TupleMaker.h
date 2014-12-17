// emacs -*- C++ -*-
#ifndef HLVF_TUPLEMAKER_H
#define HLVF_TUPLEMAKER_H

#include "SusyntHlfv/TupleMakerObjects.h"
#include "SusyntHlfv/WeightVariations.h"

#include <string>
#include <vector>

class TTree;
class TFile;
namespace Susy
{
class Lepton;
class Jet;
class Met;
}
// LeptonVector is defined in SusyDefs.h, but that's a huge include just for one def...refactor
typedef std::vector<Susy::Lepton*> LeptonVector;
typedef std::vector<Susy::Jet*>    JetVector;

namespace hlfv
{
/// A class to save the information from SusyNt to a simpler ntuple.
/**
  Details:
  This class is meant to create small ntuples for faster turnaround.
  The nutples store the information relative to the following objects:
  - two leading leptons
  - met
  - jets
  - other leptons
  - event variables
  The information is converted from SusyNt classes to smaller and
  simpler objects (see TupleMakerObjects.h)

  davide.gerbaudo@gmail.com
  November 2013
*/
class TupleMaker {
public:
    /// c'tor.
    /**
       If delayInit is false, the output tree is created right away;
       otherwise will need to call TupleMaker::init()
    */
    TupleMaker(const std::string &outFilename, const std::string &treename, bool delayInit=true);
    /// d'tor
    /**
       If necessary, takes care of closing the output file and
       cleaning things up.
     */
    ~TupleMaker();
    /// initialize output file (and tree)
    bool init(const std::string &outFilename, const std::string &treename);
    /// close output file
    bool close();
    /// fill variables for event parameters, dilepton, and met
    bool fill(const double weight, const unsigned int run, const unsigned int event,
              const Susy::Lepton &l0, const Susy::Lepton &l1, const Susy::Met &met);
    /// fill as above, + jets
    bool fill(const double weight, const unsigned int run, const unsigned int event,
              const Susy::Lepton &l0, const Susy::Lepton &l1, const Susy::Met &met,
              const JetVector &jets);
    /// fill as above, + other leptons and jets
    bool fill(const double weight, const unsigned int run, const unsigned int event,
              const Susy::Lepton &l0, const Susy::Lepton &l1, const Susy::Met &met,
              const LeptonVector &otherLeptons, const JetVector &jets);
    //! methods to assign the pieces of info that are not accessible from Lepton (mostly fake-related)
    TupleMaker& setNumTaus(size_t n) { eventPars_.numTaus = n; return *this; }
    TupleMaker& setNumFjets(size_t n) { eventPars_.numFjets = n; return *this; }
    TupleMaker& setNumBjets(size_t n) { eventPars_.numBjets = n; return *this; }
    TupleMaker& setL0IsTight(bool v) { l0_.setIsTight(v); return *this; }
    TupleMaker& setL1IsTight(bool v) { l1_.setIsTight(v); return *this; }
    TupleMaker& setL0Source(int v) { l0_.setSource(v); return *this; }
    TupleMaker& setL1Source(int v) { l1_.setSource(v); return *this; }
    TupleMaker& setL0EtConeCorr(float v) { l0_.setEtConeCorr(v); return *this; }
    TupleMaker& setL1EtConeCorr(float v) { l1_.setEtConeCorr(v); return *this; }
    TupleMaker& setL0PtConeCorr(float v) { l0_.setPtConeCorr(v); return *this; }
    TupleMaker& setL1PtConeCorr(float v) { l1_.setPtConeCorr(v); return *this; }
    TupleMaker& setQflipWeight(const double &v) { eventPars_.setQflipWeight(v); return *this; }
    TupleMaker& setWeightVariations(const hlfv::WeightVariations &w) { weightVariations_ = w; return *this; }

private: // rule of three
    TupleMaker(const TupleMaker&);
    TupleMaker& operator=(const TupleMaker&);
private:
    bool initFile(const std::string &outFilename); ///< init output file
    bool initTree(const std::string &treename); ///< init output tree
    bool initTreeBranches(); ///< create branches
private:
    TFile *file_;
    TTree *tree_;
    hlfv::FourMom l0_, l1_, met_;
    std::vector<FourMom> jets_, lowptLepts_;
    EventParameters eventPars_;
    hlfv::WeightVariations weightVariations_;
}; // end TupleMaker

} // namespace hlfv

#endif // end include guard
