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


#ifndef LMP_FIX_REACTOR_H
#define LMP_FIX_REACTOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixReactor : public Fix {
 public:
  int compute_flag;

  FixReactor(class LAMMPS *, int, char **);
  virtual ~FixReactor() {}
  int modify_param(int, char **);
  virtual int setmask();
  virtual void post_integrate();
  virtual void compute();
};

}

#endif

/* ERROR/WARNING messages:
*/
