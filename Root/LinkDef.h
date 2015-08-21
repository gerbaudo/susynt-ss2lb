#include "susynt-ss3l/Selector.h"
#include "susynt-ss3l/FakeTuplizer.h"
#include "susynt-ss3l/TupleMakerObjects.h"
#include "susynt-ss3l/WeightVariations.h"
#include "susynt-ss3l/LeptonTruthType.h"
/* #include "susynt-ss3l/kinematic.h" */

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#pragma link C++ namespace ss3l;
#pragma link C++ class ss3l::Selector+;
#pragma link C++ class ss3l::FakeTuplizer+;
#pragma link C++ struct ss3l::FourMom+;
#pragma link C++ struct ss3l::LeptonTruthType+;
#pragma link C++ struct ss3l::EventParameters+;
#pragma link C++ struct ss3l::WeightVariations+;
#pragma link C++ enum ss3l::Systematic;
#pragma link C++ class vector<ss3l::FourMom>+;

#pragma link C++ function ss3l::source2string(const ss3l::LeptonTruthType::Value &);
#pragma link C++ function ss3l::int2source(const ss3l::LeptonTruthType::Value &);

/* #pragma link C++ namespace susy::wh::kin; */
/* #pragma link C++ struct susy::wh::kin::DilepVars+; */
#endif
