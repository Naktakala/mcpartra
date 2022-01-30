#include "sdsolver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiMath/RandomNumberGeneration/random_number_generator.h"

#include "mcpartra.h"

void mcpartra::SourceDrivenSolver::TestCode()
{
  chi_log.Log() << "\nExecuting test code.\n\n";

  typedef chi_mesh::Vector3 Vec3;
  typedef chi_mesh::Matrix3x3 Mat3x3;
  const Vec3 i_hat(1.0,0.0,0.0);

  chi_math::RandomNumberGenerator rng;

  //================================================== Define source
  const double S = 1.0;
  const auto   x_s = Vec3(-1.0,0.0,0.0);

  for (int k=0; k<1; ++k)
  {
    //=========================================== Define Scale
    const double scale = 0.1;

    //=========================================== Define Rotation matrix
    const auto vec_R = SampleRandomDirection(rng);
//    const auto vec_R = Vec3(0.0,0.0,1.0);
    Mat3x3 R = Mat3x3::MakeRotationMatrixFromVector(vec_R);

    std::cout << R.PrintS() << std::endl;

    //=========================================== Define Translation
    const auto t = Vec3(15.0, 0.0, 0.0);

    std::cout << "Translation: " << t.PrintStr() << std::endl;

    //=========================================== Transform mesh and determine
    //                                            biasing cosine
    double mu_b = 1.0;
    const auto& cell = grid->local_cells[0];
    const double cell_volume = fv->MapFeView(cell.local_id)->volume*scale*scale*scale;

    for (uint64_t vid : cell.vertex_ids)
    {
      const Vec3 v = grid->vertices[vid];
      grid->vertices[vid] = scale * (R * v) + t;

      std::cout << grid->vertices[vid].PrintStr() << std::endl;

      const Vec3 vx_s_hat = (grid->vertices[vid] - x_s).Normalized();

      mu_b = std::min(mu_b, vx_s_hat.Dot(i_hat));
    }
    grid->local_cells[0].RecomputeCentroidsAndNormals(*grid);
    chi_log.Log() << "Biasing cosine: " << mu_b;

    //=========================================== Shoot src particles
    const int N_p = 1000000;
    int num_score = 0;
    int num_lost = 0;
    double total_tlw     = 0.0;
    double total_tlw_sqr = 0.0;

    double avg_D = 0.0;
    for (int p=0; p<N_p; ++p)
    {
      const double dmu = (1.0 - mu_b);
      const double mu = mu_b + dmu*rng.Rand();
      const double w  = 0.5*dmu;

      const double varphi = 2.0*M_PI*rng.Rand();
      const double theta  = acos(mu);

      const auto omega = Vec3(cos(theta),
                              sin(theta)*sin(varphi),
                              -sin(theta)*cos(varphi));

      chi_mesh::RayTracerOutputInformation ray_info =
        default_raytracer.TraceIncidentRay(cell, x_s, omega);

      if (not ray_info.particle_lost)
      {
        ++num_score;
        auto omega_t = omega;
        auto track_pos_i = ray_info.pos_f;

        avg_D += (track_pos_i - x_s).Norm();

        chi_mesh::RayTracerOutputInformation track_info =
          default_raytracer.TraceRay(cell, track_pos_i, omega_t);

        if (not track_info.particle_lost)
        {
          const auto track_pos_f = track_info.pos_f;
          const double tl  = (track_pos_f - track_pos_i).Norm();
          const double tlw = tl*w;

          total_tlw     += tlw;
          total_tlw_sqr += tlw*tlw;
        }
        else
          ++num_lost;
      }

    }//for p
    avg_D /= num_score;
    const double avg_flux = total_tlw/(N_p*cell_volume);

    chi_log.Log() << "Num score: " << num_score << " (" << double(num_score)/N_p
                  << ") " << 1/(0.25*M_PI*(sqrt(2.0)*sqrt(2.0)));
    chi_log.Log() << "Num lost: " << num_lost;


    chi_log.Log() << "Average D: " << avg_D;
    chi_log.Log() << "Average flux: " << avg_flux;

    const double D = (t - x_s).Norm();
    const double psi_analytic = S/(4.0*M_PI*D*D);

    chi_log.Log() << "Analytical flux: " << psi_analytic << " " << psi_analytic/avg_flux;
  }


  chi_log.Log() << "\nDone executing test code.\n\n";
}