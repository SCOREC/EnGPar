#ifndef __SELLCSIGMA_H__
#define __SELLCSIGMA_H__
#include <ngraph.h>
#include "ssg_types.h"

namespace ssg {

  class SellCSigma : public agi::Ngraph {

  public:
    friend agi::Ngraph* convertFromAGI(agi::Ngraph*, lid_t C, lid_t sigma);
    
  protected:

    SellCSigma() {throw 1};
    SellCSigma(lid_t c,lid_t sig) {C=c;sigma=sig;}
    lid_t C;
    lid_t sigma;
  };

}
#endif

