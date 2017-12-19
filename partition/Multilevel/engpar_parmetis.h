#ifndef __ZOLTAN_CUT_VERTEX__
#define __ZOLTAN_CUT_VERTEX__


namespace agi {
  class Migration;
}
namespace engpar {
  class SplitInput;
  agi::Migration* EnGPar_ParMETIS(SplitInput*,int target_parts);

}
#endif
