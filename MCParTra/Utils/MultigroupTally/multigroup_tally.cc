#include "multigroup_tally.h"

//###################################################################
/**Resizes all the underlying data structures.*/
void mcpartra::MultigroupTally::Resize(size_t tally_size)
{
  tally_local         .resize(tally_size, 0.0);
  tally_global        .resize(tally_size, 0.0);
  tally_sqr_local     .resize(tally_size, 0.0);
  tally_sqr_global    .resize(tally_size, 0.0);
  tally_sigma         .resize(tally_size, 0.0);
  tally_relative_sigma.resize(tally_size, 0.0);

  is_empty = false;
}

//###################################################################
/**Assignment operator.*/
mcpartra::MultigroupTally& mcpartra::MultigroupTally::
  operator=(const MultigroupTally& that)
{
  tally_local          = that.tally_local;
  tally_global         = that.tally_global;
  tally_sqr_local      = that.tally_sqr_local;
  tally_sqr_global     = that.tally_sqr_global;

  tally_sigma          = that.tally_sigma;
  tally_relative_sigma = that.tally_relative_sigma;

  return *this;
}

//###################################################################
/**Zeroes all the data structures leaving the dimension the same.*/
void mcpartra::MultigroupTally::ZeroOut()
{
  size_t tally_size = tally_local.size();

  tally_local          .assign(tally_size, 0.0);
  tally_global         .assign(tally_size, 0.0);
  tally_sqr_local      .assign(tally_size, 0.0);
  tally_sqr_global     .assign(tally_size, 0.0);

  tally_sigma          .assign(tally_size, 0.0);
  tally_relative_sigma .assign(tally_size, 0.0);
}

//###################################################################
/**Addition operator.*/
mcpartra::MultigroupTally& mcpartra::MultigroupTally::
  operator+=(const MultigroupTally& that)
{
  size_t tally_size = tally_local.size();

  for (size_t i=0; i<tally_size; ++i)
  {
    tally_local         [i] += that.tally_local         [i];
    tally_global        [i] += that.tally_global        [i];
    tally_sqr_local     [i] += that.tally_sqr_local     [i];
    tally_sqr_global    [i] += that.tally_sqr_global    [i];
    tally_sigma         [i] += that.tally_sigma         [i];
    tally_relative_sigma[i] += that.tally_relative_sigma[i];
  }

  return *this;
}

//###################################################################
/**Multiplication by scalar operator.*/
mcpartra::MultigroupTally& mcpartra::MultigroupTally::
  operator*=(const double value)
{
  size_t tally_size = tally_local.size();

  for (size_t i=0; i<tally_size; ++i)
  {
    tally_local         [i] *= value;
    tally_global        [i] *= value;
    tally_sqr_local     [i] *= value;
    tally_sqr_global    [i] *= value;
    tally_sigma         [i] *= value;
    tally_relative_sigma[i] *= value;
  }

  return *this;
}