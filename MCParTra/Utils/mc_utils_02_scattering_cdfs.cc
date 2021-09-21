#include "mcpartra.h"

#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"
#include "ChiMath/GolubFischer/GolubFischer.h"

#include <cassert>

//###################################################################
/**Uses a transfer matrix to determine group-to-group scattering
 * probabilities between G amount of groups. The result is a
 * vector of vectors. The first index is the initial group (gprime)
 * and the second index the final group (g).*/
MatDbl mcpartra::ComputeGroupToGroupScatteringCDFs(
  const size_t G,
  const chi_math::SparseMatrix &isotropic_transfer_matrix)
{
  //============================================= Make dense matrix from
  //                                              sparse matrix
  MatDbl S_0(G, VecDbl(G, 0.0)); //Dense transfer matrix
  for (size_t g=0; g < G; g++)
  {
    size_t num_transfer = isotropic_transfer_matrix.rowI_indices[g].size();
    for (size_t j=0; j<num_transfer; j++)
    {
      size_t gp = isotropic_transfer_matrix.rowI_indices[g][j];
      S_0[g][gp] = isotropic_transfer_matrix.rowI_values[g][j];
    }//for j
  }//for g

  //============================================= Build CDFs for scattering
  //                                              from gprime to g
  int G_int = static_cast<int>(G);

  //====================================== Compute column sum of S[0]
  VecDbl S0_column_sum(G,0.0);
  for (int gp=0; gp<G_int; ++gp)
    for (int g=0; g<G_int; ++g)
      S0_column_sum[gp] += S_0[g][gp];

  //====================================== Initialize cdf with a running
  //                                       total along the columns
  // Note: this will result in a cdf
  // that is not yet normalized
  MatDbl cdf_g_gprime(G, VecDbl(G,0.0));
  for (int gp=0; gp<G_int; ++gp)
    for (int g=0; g<G_int; ++g)
      cdf_g_gprime[g][gp] = cdf_g_gprime[std::max(g-1,0)][gp] + S_0[g][gp];

  //====================================== Normalize cdf
  for (int gp=0; gp<G_int; ++gp)
    for (int g=0; g<G_int; ++g)
      cdf_g_gprime[g][gp] /= S0_column_sum[gp];

  return chi_math::Transpose(cdf_g_gprime);
}

//###################################################################
/**Compute Discrete Scattering Angles that support the required number
 * of moments (num_moments_to_support or L). The required number of moments
 * to support, L, can be less than or equal to the number of transfer-matrices
 * supplied (this rule will automatically be enforced). The function returns
 * a vector of a vector containing a CDF. The first index is gprime, the
 * second index is g. The CDF contains pairs where the first value is the
 * cosine and the second value is the associated cumulative probability.*/
std::vector<mcpartra::VecCosineCDFs> mcpartra::ComputeDiscreteScatteringAngles(
  size_t G,
  const std::vector<chi_math::SparseMatrix>& transfer_matrices,
  size_t num_moments_to_support)
{
  assert(not transfer_matrices.empty());
  assert(G == transfer_matrices[0].NumRows());

  std::vector<VecCosineCDFs> cosine_cdf_gprime_g(G, VecCosineCDFs(G));

  size_t num_moments = std::min(transfer_matrices.size(), num_moments_to_support);
  if (num_moments == 0) return cosine_cdf_gprime_g;

  //============================================= Make dense matrices from
  //                                              sparse matrices
  std::vector<MatDbl> S(num_moments);            //Dense transfer matrices
  for (size_t m=0; m<num_moments; ++m)
  {
    MatDbl S_m(G, VecDbl(G,0.0));
    for (size_t g=0; g < G; g++)
    {
      size_t num_transfer = transfer_matrices[m].rowI_indices[g].size();
      for (size_t j=0; j<num_transfer; j++)
      {
        size_t gp = transfer_matrices[m].rowI_indices[g][j];
        S_m[g][gp] = transfer_matrices[m].rowI_values[g][j];
      }//for j
    }//for g

    S[m] = std::move(S_m);
  }//for m

  //============================================= Build cosine-CDFs for
  //                                              scattering from gprime to g
  typedef std::pair<double,double> MuW;  //Cosine weight pairs
  typedef std::vector<MuW>  ScatterMuWs; //Scattering cosines and weights


  //====================================== Make reindexed transfer matrices
  // For a given gprime to g, the
  // Golub-Fischer algorithm uses just
  // the moments
  std::vector<MatDbl> S_gprime_g_m(G, MatDbl(G, VecDbl(num_moments, 0.0)));
  for (size_t gp=0; gp<G; ++gp)
    for (size_t g=0; g<G; ++g)
      for (size_t m=0; m<num_moments; ++m)
        S_gprime_g_m[gp][g][m] = S[m][g][gp];

  //====================================== Make cosine-cdfs for each
  //                                       gprime-g pair
  for (size_t gp=0; gp<G; ++gp)
    for (size_t g=0; g<G; ++g)
    {
      GolubFischer gb;
      if (S[0][g][gp] < 1.0e-16 or num_moments_to_support == 0) continue;
      ScatterMuWs scatter_mu_ws = gb.GetDiscreteScatAngles(S_gprime_g_m[gp][g]);

      double total_weight = 0.0;
      for (auto mu_w : scatter_mu_ws)
        total_weight += mu_w.second;

      CosineCDF cosine_cdf;
      cosine_cdf.reserve(scatter_mu_ws.size());
      double running_weight_total = 0.0;
      for (auto mu_w : scatter_mu_ws)
      {
        double mu = mu_w.first;
        double w  = mu_w.second;

        running_weight_total += w;

        double cp = running_weight_total / total_weight; //Cumulative probability
        cosine_cdf.emplace_back(mu, cp);
      }//for mu_w

      cosine_cdf_gprime_g[gp][g] = std::move(cosine_cdf);
    }//for g

  return cosine_cdf_gprime_g;
}