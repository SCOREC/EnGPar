#ifndef ENGPAR_TYPES_H__
#define ENGPAR_TYPES_H__
namespace engpar {
  typedef agi::wgt_t wgt_t;
  typedef agi::part_t part_t;

  enum SPLIT_METHOD {
    GLOBAL_PARMETIS,
    LOCAL_PARMETIS
  };

}

#endif
