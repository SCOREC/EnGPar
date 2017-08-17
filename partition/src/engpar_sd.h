#ifndef ENGPAR_SD_H__
#define ENGPAR_SD_H__

namespace engpar {
  class SDSlope {
  public:
    SDSlope();
    ~SDSlope();

    bool isFull();
    double get(unsigned int i);
    void push(double v);

    double slope();
    double average();
    
  private:
    unsigned int max_len;
    unsigned int cur_len;
    unsigned int pos;
    double* vals;

  };
  
}

#endif
