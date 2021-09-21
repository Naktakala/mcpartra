#ifndef MCPARTRA_MG_SCATTERING_CDFS_H
#define MCPARTRA_MG_SCATTERING_CDFS_H

#include <vector>

namespace chi_physics
{
  class TransportCrossSections;
}

namespace mcpartra
{

class MultigroupScatteringCDFs
{
private:
  typedef std::vector<std::pair<double,double>> Tvecdbl_vecdbl;
  std::vector<std::vector<double>>         cdf_gprime_g;
  std::vector<std::vector<Tvecdbl_vecdbl>> scat_angles_gprime_g;

public:
  MultigroupScatteringCDFs(const chi_physics::TransportCrossSections& mgxs,
                           unsigned int L_parameter);

  unsigned int Sample_gprime(unsigned int gp, double rn) const;
  double SampleMu_gprime_g(unsigned int gprime,
                           unsigned int g,
                           double random_number, bool isotropic) const;
};

}//namespace mcpartra

#endif //MCPARTRA_MG_SCATTERING_CDFS_H
