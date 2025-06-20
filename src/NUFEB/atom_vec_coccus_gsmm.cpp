#include "atom_vec_coccus_gsmm.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "modify.h"

#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecGSMM::AtomVecGSMM(LAMMPS *lmp) : AtomVecCoccus(lmp){
  atom->gsmm_flag = 1;
}

void AtomVecGSMM::create_atom_post(int ilocal)
{
  AtomVecCoccus::create_atom_post(ilocal);
}