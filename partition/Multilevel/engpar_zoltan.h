#ifndef __ENGPAR_ZOLTAN__
#define __ENGPAR_ZOLTAN__

namespace agi {
  class Migration;
}
namespace engpar {
  class SplitInput;
  agi::Migration* EnGPar_Zoltan(SplitInput*, int target_parts, bool isLocal);
}
#endif
