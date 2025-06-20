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

#include "grid_vec_gsmm.h"

#include "grid.h"
#include "force.h"
#include "error.h"
#include "memory.h"
#include "grid_masks.h"
#include "comm.h"
#include "atom.h"
#include "group.h"

#include <iostream>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridVecGSMM::GridVecGSMM(LAMMPS *lmp) : GridVecChemostat(lmp) {
  grid->gsmm_flag = 1; 
}

/* ---------------------------------------------------------------------- */

void GridVecGSMM::init()
{
  // Call parent initialization first
  GridVecChemostat::init();
  
  // Any additional GSMM-specific initialization can go here
}

/* ---------------------------------------------------------------------- */

void GridVecGSMM::grow(int n)
{
  GridVecChemostat::grow(n);
}

/* ---------------------------------------------------------------------- */

int GridVecGSMM::pack_comm(int n, int *cells, double *buf)
{
  return GridVecChemostat::pack_comm(n, cells, buf);
}

/* ---------------------------------------------------------------------- */

void GridVecGSMM::unpack_comm(int n, int *cells, double *buf)
{
  GridVecChemostat::unpack_comm(n, cells, buf);
}

/* ---------------------------------------------------------------------- */

int GridVecGSMM::pack_exchange(int n, int *cells, double *buf)
{
  return GridVecChemostat::pack_exchange(n, cells, buf);
}

/* ---------------------------------------------------------------------- */

void GridVecGSMM::unpack_exchange(int n, int *cells, double *buf)
{
  GridVecChemostat::unpack_exchange(n, cells, buf);
}

/* ---------------------------------------------------------------------- */

void GridVecGSMM::set(int narg, char **arg)
{
  GridVecChemostat::set(narg, arg);
}


/* ---------------------------------------------------------------------- */

void GridVecGSMM::set_grid(int isub, double domain, double bulk)
{
  GridVecChemostat::set_grid(isub, domain, bulk);
}
