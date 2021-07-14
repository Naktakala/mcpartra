#ifndef MCPARTRA_GRIDTALLYBLOCK_H
#define MCPARTRA_GRIDTALLYBLOCK_H

#include <vector>

namespace mcpartra
{
//######################################################### Tally struct
/**Tally block.*/
class GridTallyBlock
{
private:
  bool is_empty = true;
public:
  std::vector<double> tally_local;
  std::vector<double> tally_global;
  std::vector<double> tally_sqr_local;
  std::vector<double> tally_sqr_global;

  std::vector<double> tally_sigma;
  std::vector<double> tally_relative_sigma;

  void Resize(size_t tally_size);
  GridTallyBlock& operator=(const GridTallyBlock& that);
  void ZeroOut();
  GridTallyBlock& operator+=(const GridTallyBlock& that);
  GridTallyBlock& operator*=(double value);
  bool empty() const {return is_empty;}

//  void Resize(size_t tally_size)
//  {
//    tally_local         .resize(tally_size, 0.0);
//    tally_global        .resize(tally_size, 0.0);
//    tally_sqr_local     .resize(tally_size, 0.0);
//    tally_sqr_global    .resize(tally_size, 0.0);
//    tally_sigma         .resize(tally_size, 0.0);
//    tally_relative_sigma.resize(tally_size, 0.0);
//
//    is_empty = false;
//  }
//
//  GridTallyBlock& operator=(const GridTallyBlock& that)
//  {
//    tally_local          = that.tally_local;
//    tally_global         = that.tally_global;
//    tally_sqr_local      = that.tally_sqr_local;
//    tally_sqr_global     = that.tally_sqr_global;
//
//    tally_sigma          = that.tally_sigma;
//    tally_relative_sigma = that.tally_relative_sigma;
//
//    return *this;
//  }
//
//  void ZeroOut()
//  {
//    size_t tally_size = tally_local.size();
//
//    tally_local          .assign(tally_size, 0.0);
//    tally_global         .assign(tally_size, 0.0);
//    tally_sqr_local      .assign(tally_size, 0.0);
//    tally_sqr_global     .assign(tally_size, 0.0);
//
//    tally_sigma          .assign(tally_size, 0.0);
//    tally_relative_sigma .assign(tally_size, 0.0);
//  }
//
//  GridTallyBlock& operator+=(const GridTallyBlock& that)
//  {
//    size_t tally_size = tally_local.size();
//
//    for (size_t i=0; i<tally_size; ++i)
//    {
//      tally_local         [i] += that.tally_local         [i];
//      tally_global        [i] += that.tally_global        [i];
//      tally_sqr_local     [i] += that.tally_sqr_local     [i];
//      tally_sqr_global    [i] += that.tally_sqr_global    [i];
//      tally_sigma         [i] += that.tally_sigma         [i];
//      tally_relative_sigma[i] += that.tally_relative_sigma[i];
//    }
//
//    return *this;
//  }
//
//  GridTallyBlock& operator*=(const double value)
//  {
//    size_t tally_size = tally_local.size();
//
//    for (size_t i=0; i<tally_size; ++i)
//    {
//      tally_local         [i] *= value;
//      tally_global        [i] *= value;
//      tally_sqr_local     [i] *= value;
//      tally_sqr_global    [i] *= value;
//      tally_sigma         [i] *= value;
//      tally_relative_sigma[i] *= value;
//    }
//
//    return *this;
//  }

//  bool empty() const {return is_empty;}
};
}


#endif //MCPARTRA_GRIDTALLYBLOCK_H