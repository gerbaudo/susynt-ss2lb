#include "susynt-ss2lb/Selector.h"
#include "susynt-ss2lb/FakeTuplizer.h"
#include "susynt-ss2lb/TupleMakerObjects.h"
#include "susynt-ss2lb/WeightVariations.h"
#include "susynt-ss2lb/MatrixPrediction.h"
#include "susynt-ss2lb/LeptonTruthType.h"
/* #include "susynt-ss2lb/kinematic.h" */

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#pragma link C++ namespace hlfv;
#pragma link C++ class hlfv::Selector+;
#pragma link C++ class hlfv::FakeTuplizer+;
#pragma link C++ class hlfv::MatrixPrediction+;
#pragma link C++ struct hlfv::FourMom+;
#pragma link C++ struct hlfv::LeptonTruthType+;
#pragma link C++ struct hlfv::EventParameters+;
#pragma link C++ struct hlfv::WeightVariations+;
#pragma link C++ enum hlfv::Systematic;
#pragma link C++ class vector<hlfv::FourMom>+;

#pragma link C++ function hlfv::source2string(const hlfv::LeptonTruthType::Value &);
#pragma link C++ function hlfv::int2source(const hlfv::LeptonTruthType::Value &);

/* #pragma link C++ namespace susy::wh::kin; */
/* #pragma link C++ struct susy::wh::kin::DilepVars+; */
#endif
