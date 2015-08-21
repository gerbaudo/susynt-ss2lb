// emacs -*- C++ -*-
#ifndef SusyAna_MatrixPrediction_h
#define SusyAna_MatrixPrediction_h


#include "susynt-ss3l/Selector.h"
#include "DileptonMatrixMethod/DileptonMatrixMethod.h"

#include <string>

namespace susy {
namespace fake {
class Lepton;
}
}

namespace ss3l
{

const int nFakeSys = 21;
static std::string FAKESYSNames[nFakeSys] =
{
  "NONE",
  "ALL_UP",
  "ALL_DN",
  "EL_RE_UP",
  "EL_RE_DOWN",
  "EL_FR_UP",
  "EL_FR_DOWN",
  "MU_RE_UP",
  "MU_RE_DOWN",
  "MU_FR_UP",
  "MU_FR_DOWN",
  "EL_SSOS_UP",
  "EL_SSOS_DOWN",
  "EL_WJETS",
  "EL_TTBAR",
  "EL_METREL",
  "MU_SSOS_UP",
  "MU_SSOS_DOWN",
  "MU_WJETS",
  "MU_TTBAR",
  "MU_METREL"
};

/**
    @brief Looper to run over data and make fake lepton prediction

    It applies the LL/LT/TT event weight from the matrix method.

    \todo cleanup

    davide.gerbaudo@gmail.com
    Sep 2014
*/
 class MatrixPrediction : public Selector
{

 public:
  MatrixPrediction();
  virtual void    Begin(TTree *tree); //!< called before looping on entries
  virtual void    Terminate(); //!< called after looping is finished
  virtual Bool_t  Process(Long64_t entry);
    // Get the fake event weight given a signal region
    float getFakeWeight(const LeptonVector &baseLeps,
                        std::string &regionName,
                        float metRel,
                        susy::fake::Systematic::Value);
    float getRFWeight(const LeptonVector &baseLeps,
                      std::string &regionName,
                      float metRel,
                      susy::fake::Systematic::Value);
    MatrixPrediction& setMatrixFilename(const std::string filename); ///< to be called before Begin
    MatrixPrediction& use2dParametrization() { m_use2dparametrization = true; return *this; }
    static std::string dilepDetails(const Susy::Event &event, const DiLepEvtType &ll,
                                    const LeptonVector &ls);
    std::string eventDetails(bool passSrSs, const Susy::Event &event, const DiLepEvtType &ll,
                             const LeptonVector &ls);
    bool m_dbg;
    ClassDef(MatrixPrediction, 3);

  protected:
    std::vector<uint> m_matrixSys;      ///< systematics to process
    susy::fake::DileptonMatrixMethod* m_matrix;
    std::string m_matrixFilename;
    bool m_use2dparametrization;
    bool m_allconfigured;
    bool initMatrixTool();
    /// compute all the weight variations for the fake systematics
    ss3l::WeightVariations computeSystematicWeights(const susy::fake::Lepton &l0, const susy::fake::Lepton &l1, size_t regionIndex);
};
} // ss3l
#endif
