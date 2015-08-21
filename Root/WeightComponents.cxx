#include "susynt-ss3l/WeightComponents.h"

#include <sstream>      // std::ostringstream

using ss3l::WeightComponents;
using namespace std;

//-----------------------------------------
std::string WeightComponents::str() const
{
  ostringstream oss;
  oss<<" susynt: "<<susynt
     <<" lepSf: "<<lepSf
     <<" btag: "<<btag
     <<" trigger: "<<trigger
     <<" qflip: "<<qflip
     <<" fake: "<<fake;
  return oss.str();
}
//-----------------------------------------
