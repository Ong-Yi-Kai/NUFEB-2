/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef GRID_CLASS

GridStyle(nufeb/gsmm,GridVecGSMM)

#else

#ifndef LMP_GRID_VEC_GSMM_H
#define LMP_GRID_VEC_GSMM_H

#include "grid_vec_chemostat.h"

namespace LAMMPS_NS {

class GridVecGSMM : public GridVecChemostat {
 public:
   GridVecGSMM(class LAMMPS *);
   void init();
   void grow(int);
   int pack_comm(int, int *, double *);
   void unpack_comm(int, int *, double *);
   int pack_exchange(int, int *, double *);
   void unpack_exchange(int, int *, double *);
   void set(int, char **);
   void set_grid(int, double, double);

 private:

};

}

#endif
#endif
