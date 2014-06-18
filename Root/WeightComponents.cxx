#include "SusyntHlfv/WeightComponents.h"

#include <sstream>      // std::ostringstream

using hlfv::WeightComponents;
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
