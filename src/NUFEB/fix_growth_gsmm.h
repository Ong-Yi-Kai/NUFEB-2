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
#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#include <sbml/packages/fbc/extension/FbcExtension.h>
#include <sbml/packages/fbc/extension/FbcModelPlugin.h>
#include <glpk.h>
#include <string>
#include <vector>
#include <map>
#include "fix_growth.h"

#ifdef FIX_CLASS

FixStyle(nufeb/growth/gsmm,FixGrowthGSMM)

#else

#ifndef LMP_FIX_GROWTH_GSMM_H
#define LMP_FIX_GROWTH_GSMM_H


namespace LAMMPS_NS {

class FixGrowthGSMM: public FixGrowth {
    public:
        FixGrowthGSMM(class LAMMPS *, int, char **);
        FixGrowthGSMM(const std::string& fpath, 
                    const std::vector<std::string>& ex_mets);
        ~FixGrowthGSMM();

        virtual void update_atoms();
        virtual void update_cells();

    private:
        Model* model;
        glp_prob* lp;
        std::map<std::string, int> met_indices;
        std::map<std::string, int> rxn_indices;
        std::vector<std::string> exchange_mets;
        std::vector<int> ex_rxn_indices;

        void setExchangeConc(const std::string& ex_met_id, const double& conc);
        void solve();
        double getBiomass();
        double getExchangeFlux(const std::string& ex_met_id);


   
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
