#ifdef ATOM_CLASS

AtomStyle(gsmm,AtomVecGSMM)

#else

#ifndef LMP_ATOM_VEC_COCCUSGSMM_H
#define LMP_ATOM_VEC_COCCUSGSMM_H

#include "atom_vec_coccus.h"

namespace LAMMPS_NS {

class AtomVecGSMM : public AtomVecCoccus {
 public:
  AtomVecGSMM(class LAMMPS *);
  void create_atom_post(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid radius in Atoms section of data file

Radius must be >= 0.0.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

*/
