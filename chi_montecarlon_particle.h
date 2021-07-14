#ifndef MCPARTRA_PARTICLE_H
#define MCPARTRA_PARTICLE_H

#include <ChiMesh/chi_mesh.h>
#include <chi_mpi.h>

/**Structure for storing particle information.*/
struct mcpartra::Particle final
{
  //Pos
  chi_mesh::Vector3 pos;

  //Dir
  chi_mesh::Vector3 dir = chi_mesh::Vector3(0.0,0.0,1.0);

  double w = 1.0; //Weight
  int egrp = 0; //Energy group

  uint64_t cur_cell_global_id = 0;
  uint64_t pre_cell_global_id = 0;

  uint64_t cur_cell_local_id = 0;
  uint64_t pre_cell_local_id = 0;

  int ray_trace_method = 0; //STANDARD
  int tally_method = 0;     //STANDARD
  int tally_mask = (1 << 0) | //DEFAULT_FVTALLY
                   (1 << 1);  //DEFAULT_PWLTALLY

  double cur_cell_importance = 1.0;
  double pre_cell_importance = 1.0;

  bool alive = true;
  bool banked = false;

  double moment_values[64] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


  Particle() {}

  static Particle MakeDeadParticle()
  {
    Particle new_particle;
    new_particle.alive = false;
    return  new_particle;
  }

  //=================================== Copy operator
  Particle& operator=(const Particle& that)
  {
    this->pos.x = that.pos.x;  this->dir.x = that.dir.x;
    this->pos.y = that.pos.y;  this->dir.y = that.dir.y;
    this->pos.z = that.pos.z;  this->dir.z = that.dir.z;

    this->w = that.w;
    this->egrp = that.egrp;
    this->cur_cell_global_id = that.cur_cell_global_id;
    this->pre_cell_global_id = that.pre_cell_global_id;

    this->cur_cell_local_id = that.cur_cell_local_id;
    this->pre_cell_local_id = that.pre_cell_local_id;

    this->ray_trace_method = that.ray_trace_method;
    this->tally_method = that.tally_method;
    this->tally_mask       = that.tally_mask;

    this->cur_cell_importance = that.cur_cell_importance;
    this->pre_cell_importance = that.pre_cell_importance;

    this->alive = that.alive;
    this->banked = that.banked;

    for (int i=0; i<64; ++i)
      this->moment_values[i] = that.moment_values[i];

    return *this;
  }

  Particle(const Particle& that):
    pos(that.pos),
    dir(that.dir),
    w(that.w),
    egrp(that.egrp),
    cur_cell_global_id(that.cur_cell_global_id),
    pre_cell_global_id(that.pre_cell_global_id),
    cur_cell_local_id(that.cur_cell_local_id),
    pre_cell_local_id(that.pre_cell_local_id),
    ray_trace_method(that.ray_trace_method),
    tally_method(that.tally_method),
    tally_mask(that.tally_mask),
    cur_cell_importance(that.cur_cell_importance),
    pre_cell_importance(that.pre_cell_importance),
    alive(that.alive),
    banked(that.banked)
  {
    for (int i=0; i<64; ++i)
      this->moment_values[i] = that.moment_values[i];
  }

  static void BuildMPIDatatype(MPI_Datatype& prtcl_data_type)
  {
    int block_lengths[] = {3,  //pos
                           3,  //dir
                           1,  //w
                           1,  //egrp
                           1,  //cur_cell_global_id
                           1,  //pre_cell_global_id
                           1,  //cur_cell_local_id
                           1,  //pre_cell_local_id
                           1,  //ray_trace_method
                           1,  //tally_method
                           1,  //tally_mask
                           1,  //cur_cell_importance
                           1,  //pre_cell_importance
                           1,  //alive
                           1,  //banked
                           64};//moments

    MPI_Aint block_disp[] = {
      offsetof(Particle, pos),
      offsetof(Particle, dir),
      offsetof(Particle, w),
      offsetof(Particle, egrp),
      offsetof(Particle, cur_cell_global_id),
      offsetof(Particle, pre_cell_global_id),
      offsetof(Particle, cur_cell_local_id),
      offsetof(Particle, pre_cell_local_id),
      offsetof(Particle, ray_trace_method),
      offsetof(Particle, tally_method),
      offsetof(Particle, tally_mask),
      offsetof(Particle, cur_cell_importance),
      offsetof(Particle, pre_cell_importance),
      offsetof(Particle, alive),
      offsetof(Particle, banked),
      offsetof(Particle, moment_values)};

    MPI_Datatype block_types[] = {
      MPI_DOUBLE,                //pos
      MPI_DOUBLE,                //dir
      MPI_DOUBLE,                //w
      MPI_INT,                   //egrp
      MPI_UNSIGNED_LONG_LONG,    //cur_cell_global_id
      MPI_UNSIGNED_LONG_LONG,    //pre_cell_global_id
      MPI_UNSIGNED_LONG_LONG,    //cur_cell_local_id
      MPI_UNSIGNED_LONG_LONG,    //pre_cell_local_id
      MPI_INT,                   //ray_trace_method
      MPI_INT,                   //tally_method
      MPI_INT,                   //tally_mask
      MPI_DOUBLE,                //cur_cell_importance
      MPI_DOUBLE,                //pre_cell_importance
      MPI_BYTE,                  //alive
      MPI_BYTE,                  //banked
      MPI_DOUBLE                 //moment_values
    };

    MPI_Type_create_struct(
      16,             // Number of blocks
      block_lengths,  // Block lengths
      block_disp,     // Block displacements
      block_types,    // Block types
      &prtcl_data_type);

    MPI_Type_commit(&prtcl_data_type);
  }


};


#endif //MCPARTRA_PARTICLE_H
