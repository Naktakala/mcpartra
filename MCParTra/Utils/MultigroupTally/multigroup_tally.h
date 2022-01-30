#ifndef MCPARTRA_GRIDTALLYBLOCK_H
#define MCPARTRA_GRIDTALLYBLOCK_H

#include <vector>
#include <cstdint>

namespace mcpartra
{
//######################################################### Tally struct
/**Tally block.*/
class MultigroupTally
{
private:
  bool is_empty = true;
public:
  std::vector<uint64_t> counter_local;
  std::vector<uint64_t> counter_global;
  std::vector<double> tally_local;
  std::vector<double> tally_global;
  std::vector<double> tally_sqr_local;
  std::vector<double> tally_sqr_global;

  std::vector<double> tally_sigma;
  std::vector<double> tally_relative_sigma;

  void Resize(size_t tally_size);
  MultigroupTally& operator=(const MultigroupTally& that);
  void ZeroOut();
  MultigroupTally& operator+=(const MultigroupTally& that);
  MultigroupTally& operator*=(double value);
  bool empty() const {return is_empty;}
};
}


#endif //MCPARTRA_GRIDTALLYBLOCK_H