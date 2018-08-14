
namespace agi {

  /*This class allows creating an array on each process to store
    some kind of data on either vertices or edges in the graph

    Note: Currently there is no error checking for passing in something
          that doesnt correspond to this tag or if looking at garbage
   */
  class Info {
  public:
    Info(lid_t s,int b) {data = new char[b*s];}

    virtual ~Info() {delete [] data;}

    int getInt(lid_t i) {return reinterpret_cast<int*>(data)[i];}
    void setInt(lid_t i, int val) {reinterpret_cast<int*>(data)[i] = val;}
    double getDouble(lid_t i) {return reinterpret_cast<double*>(data)[i];}
    void setDouble(lid_t i, double val) {reinterpret_cast<double*>(data)[i] = val;}
    long getLong(lid_t i) {return reinterpret_cast<long*>(data)[i];}
    void setLong(lid_t i, long val) {reinterpret_cast<long*>(data)[i] = val;}
  protected:
    char* data;
  };

}
