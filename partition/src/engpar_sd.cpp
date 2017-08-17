#include "engpar_sd.h"
#include <cassert>

namespace {
  unsigned int order = 2;
  const double ws[] = {-3./2,2.,-1./2};
}

namespace engpar {
  SDSlope::SDSlope() {
    max_len = order+1;
    cur_len = 0;
    pos = 0;
    vals = new double[max_len];
  }
  SDSlope::~SDSlope() {
    delete [] vals;
  }

  bool SDSlope::isFull() {
    return cur_len == max_len;
  }

  double SDSlope::get(unsigned int i) {
    unsigned int p = pos - i;
    p = (p + cur_len) % cur_len;
    return vals[p];
  }

  void SDSlope::push(double v) {
    if (!isFull())
      cur_len++;
    vals[pos++] = v;
    pos %= max_len;
  }

  double SDSlope::slope() {
    assert(isFull());
    double s = 0;
    for (unsigned int i = 0; i < max_len; i++)
      s += get(i) * ws[i];
    return s;
  }
  double SDSlope::average() {
    assert(isFull());
    double a = 0;
    for (unsigned int i = 0; i < max_len; i++)
      a += get(i);
    return a / max_len;

  }

}
