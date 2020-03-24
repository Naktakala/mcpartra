#ifndef _chi_montecarlon_particle_h
#define _chi_montecarlon_particle_h

#include <ChiMesh/chi_mesh.h>
#include <chi_mpi.h>

/**Structure for storing particle information.*/
struct chi_montecarlon::Particle
{
  //Pos
  chi_mesh::Vector3 pos;

  //Dir
  chi_mesh::Vector3 dir = chi_mesh::Vector3(0.0,0.0,1.0);

  double w = 1.0; //Weight
  int egrp = 0; //Energy group

  int cur_cell_global_id = -1;
  int pre_cell_global_id = -1;

  int cur_cell_local_id = -1;
  int pre_cell_local_id = -1;

  bool alive = true;
  bool banked = false;


  Particle() {}

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

    this->alive = that.alive;
    this->banked = that.banked;

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
    alive(that.alive),
    banked(that.banked)
  {}

  static void BuildMPIDatatype(MPI_Datatype& prtcl_data_type)
  {
    int block_lengths[] = {3, //pos
                           3, //dir
                           1, //w
                           1, //egrp
                           1, //cur_cell_global_id
                           1, //pre_cell_global_id
                           1, //cur_cell_local_id
                           1, //pre_cell_local_id
                           1, //alive
                           1};//banked

    MPI_Aint block_disp[] = {
      offsetof(Particle, pos),
      offsetof(Particle, dir),
      offsetof(Particle, w),
      offsetof(Particle, egrp),
      offsetof(Particle, cur_cell_global_id),
      offsetof(Particle, pre_cell_global_id),
      offsetof(Particle, cur_cell_local_id),
      offsetof(Particle, pre_cell_local_id),
      offsetof(Particle, alive),
      offsetof(Particle, banked)};

    MPI_Datatype block_types[] = {
      MPI_DOUBLE,
      MPI_DOUBLE,
      MPI_DOUBLE,
      MPI_INT,
      MPI_INT,
      MPI_INT,
      MPI_INT,
      MPI_INT,
      MPI_BYTE,
      MPI_BYTE
    };

    MPI_Type_create_struct(
      10,             // Number of blocks
      block_lengths,  // Block lengths
      block_disp,     // Block displacements
      block_types,    // Block types
      &prtcl_data_type);

    MPI_Type_commit(&prtcl_data_type);
  }


};


#endif
