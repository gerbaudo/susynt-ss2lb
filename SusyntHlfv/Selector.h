// emacs -*- C++ -*-
#ifndef HLVF_SELECTOR_H
#define HLVF_SELECTOR_H



#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyDefs.h"

#include "SusyntHlfv/ProgressPrinter.h"

#include <utility> // std::pair


// fw decl
class DilTrigLogic;
class MCWeighter;

namespace hlfv{

class Selector : public SusyNtAna
{
 public:
    Selector();
    virtual ~Selector(){};
    virtual void    Begin(TTree *tree);      ///< called before looping on entries
    virtual void    Init(TTree *tree);       ///< called when the TChain is attached
    virtual void    Terminate();             ///< called after looping is finished
    virtual Bool_t  Process(Long64_t entry); ///< called at each event
    bool initMcWeighter(TTree *tree);
 protected:
    DilTrigLogic*       m_trigObj;      ///< trigger logic class
    MCWeighter*         m_mcWeighter;   ///< tool to determine the normalization
    hlvf::ProgressPrinter m_printer;

    ClassDef(Selector, 1);
};

} // hlvf

#endif // HLVF_SELECTOR_H
