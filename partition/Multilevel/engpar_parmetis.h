#ifndef __ENGPAR_PARMETIS__
#define __ENGPAR_PARMETIS__


namespace agi {
  class Migration;
}
namespace engpar {
  class SplitInput;
  agi::Migration* EnGPar_ParMETIS(SplitInput*,int target_parts, bool isLocal);

}
#endif
