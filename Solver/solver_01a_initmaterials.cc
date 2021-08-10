#include"solver_montecarlon.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "ChiMath/GolubFischer/GolubFischer.h"
#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

#include "ChiLog/chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include <assert.h>

//###################################################################
/**Initialize materials.*/
void mcpartra::Solver::InitMaterials()
{
  chi_log.Log() << "MCParTra: Initializing Materials.";

  typedef chi_physics::TransportCrossSections TrXS;

  //=================================== Check materials exist
  if (chi_physics_handler.material_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::Solver::Initialize : No materials found.";
    exit(EXIT_FAILURE);
  }
  size_t num_mat = chi_physics_handler.material_stack.size();

  //=================================== Check that cells all map to a
  //                                    valid material
  size_t invalid_cell_mats = 0;
  for (const auto& cell : grid->local_cells)
    if (cell.material_id < 0 or cell.material_id >= static_cast<int>(num_mat))
      invalid_cell_mats++;

  if (invalid_cell_mats > 0)
    throw std::invalid_argument(
      "MCParTra: " + std::to_string(invalid_cell_mats) + " cells found with "
      "material-IDs pointing to invalid materials.");

  //=================================== Initialize Materials and make property
  //                                    mappings
  matid_xs_map.resize(num_mat,-1);
  matid_q_map.resize(num_mat,-1);
  for (size_t m=0; m<num_mat; m++)
  {
    auto cur_mat = chi_physics_handler.material_stack[m];
    bool material_xs_mapped = false;

    //======================= Only first xs will be used
    size_t num_props = cur_mat->properties.size();
    for (size_t p=0; p<num_props; p++)
    {
      auto& property = cur_mat->properties[p];
      auto property_type = cur_mat->properties[p]->Type();

      if (property_type == chi_physics::PropertyType::TRANSPORT_XSECTIONS and
        (not material_xs_mapped))
      {
        auto transp_xs = std::static_pointer_cast<TrXS>(property);

        transp_xs->ComputeDiscreteScattering(options.scattering_order);

        if (transp_xs->num_groups > num_groups)
          num_groups = transp_xs->num_groups;

        matid_xs_map[m] = static_cast<int>(p);
        material_xs_mapped = true;
      }

      if (property_type == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
        matid_q_map[m] = static_cast<int>(p);
    }//for prop

    using namespace std;

    if (matid_xs_map[m]<0)
      throw invalid_argument(
        "The mapping of material " + to_string(m) + " to a cross section "
        "property failed. This indicates that the given material might not "
        "have a transport cross section property.");
  }//for mat

  chi_log.Log() << "MCParTra: Number of groups = " << num_groups;
  MPI_Barrier(MPI_COMM_WORLD);
}


//###################################################################
/**Compute Discrete Scattering Angles.*/
void mcpartra::ComputeDiscreteScatteringAngles(
  size_t G,
  const std::vector<chi_math::SparseMatrix>& transfer_matrices)
{
  assert(not transfer_matrices.empty());
  assert(G == transfer_matrices[0].NumRows());

  chi_log.Log(LOG_0) << "Creating Discrete scattering angles.";

  size_t num_moments = transfer_matrices.size();

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

  //============================================= Build CDFs for scattering
  //                                              from gprime to g
  MatDbl cdf_gprime_g(G,VecDbl(G,0.0));
  {
    int G_int = static_cast<int>(G);

    //====================================== Compute column sum of S[0]
    VecDbl S0_column_sum(G,0.0);
    for (int gp=0; gp<G_int; ++gp)
      for (int g=0; g<G_int; ++g)
        S0_column_sum[gp] += S[0][g][gp];

    //====================================== Initialize cdf with a running
    //                                       total along the columns
    // Note: this will result in a cdf
    // that is not yet normalized
    MatDbl cdf_g_gprime(G, VecDbl(G,0.0));
    for (int gp=0; gp<G_int; ++gp)
      for (int g=0; g<G_int; ++g)
        cdf_g_gprime[g][gp] = cdf_g_gprime[std::max(g-1,0)][gp] + S[0][g][gp];

    //====================================== Normalize cdf
    for (int gp=0; gp<G_int; ++gp)
      for (int g=0; g<G_int; ++g)
        cdf_g_gprime[g][gp] /= S0_column_sum[gp];

    //====================================== Make the transpose for easier
    //                                       sampling
    cdf_gprime_g = chi_math::Transpose(cdf_g_gprime);
  }

  //============================================= Build cosine-CDFs for
  //                                              scattering from gprime to g
  typedef std::pair<double,double> MuCP; //Cosine and cumulative probability pairs
  typedef std::pair<double,double> MuW;  //Cosine weight pairs
  typedef std::vector<MuW>  ScatterMuWs; //Scattering cosines and weights
  typedef std::vector<MuCP> CosineCDF;   //CDF of cosines
  typedef std::vector<CosineCDF> VecCosineCDFs; //Vector of CDFs

  std::vector<std::vector<CosineCDF>> cosine_cdf_gprime_g(G, VecCosineCDFs(G));
  {
    //====================================== Make reindexed transfer matrices
    // For a given gprime to g, the
    // Golub-Fischer algorithm uses just
    // the moments
    std::vector<MatDbl> S_gprime_g_m(G, MatDbl(G, VecDbl(num_moments, 0.0)));
    for (size_t gp=0; gp<G; ++gp)
      for (size_t g=0; g<G; ++g)
        for (size_t m=0; m<num_moments; ++m)
          S_gprime_g_m[gp][g][m] = S[m][g][gp];

    //====================================== Make cosine-cdfs for each gprime-g
    //                                       pair
    for (size_t gp=0; gp<G; ++gp)
      for (size_t g=0; g<G; ++g)
      {
        GolubFischer gb;
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
  }

  chi_log.Log(LOG_0) << "Done creating Discrete scattering angles.";
}