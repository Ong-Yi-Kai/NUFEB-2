/* ----------------------------------------------------------------------
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
#include <sbml/packages/fbc/extension/FbcReactionPlugin.h>
#include <glpk.h>

#include "fix_growth_gsmm.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"
#include "comm.h"



using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

bool endsWith(const std::string& mainStr, const std::string& suffix) {
    if (mainStr.length() < suffix.length()) {
        return false;
    }
    return mainStr.substr(mainStr.length() - suffix.length()) == suffix;
}

bool startsWith(const std::string& mainStr, const std::string& prefix) {
    if (mainStr.length() < prefix.length()) {
        return false;
    }
    return mainStr.substr(0, prefix.length()) == prefix;
}

void splitByCommas(std::string str, std::vector<std::string>& result) {
    size_t pos = 0;
    while ((pos = str.find(',')) != std::string::npos) {
        std::string token = str.substr(0, pos);
        if (!token.empty()) result.push_back(token);
        str.erase(0, pos + 1);
    }
    if (!str.empty()) result.push_back(str); 
}

FixGrowthGSMM::FixGrowthGSMM(LAMMPS *lmp, int narg, char **arg) :
    FixGrowth(lmp, narg, arg){

    if (narg < 7) error->all(FLERR, "Illegal fix nufeb/growth/gsmm command");
    if (grid->gsmm_flag == 0) error->all(FLERR, "fix nufeb/growth/gsmm requires grid_style nufeb/gsmm");
    if (atom->gsmm_flag == 0) error->all(FLERR, "fix nufeb/growth/gsmm requires atom_style gsmm");
    
    // parse arguments
    std::string fpath;
    std::vector<std::string> ex_mets;

    for (int i = 3; i < narg; i+=2){
        if (strcmp(arg[i], "fpath") == 0) fpath = arg[i+1];
        else if (strcmp(arg[i], "exchange_mets") == 0) splitByCommas(arg[i+1], ex_mets);
        else error->all(FLERR, "Unknown argument in fix nufeb/growth/gsmm command, need fpath, or exchange_mets");
    }

    // initialize SBML document
    SBMLReader reader;
    SBMLDocument * doc = reader.readSBML(fpath);
    doc->enablePackage(FbcExtension::getXmlnsL3V1V2(), "fbc", true);
    doc->setPackageRequired("fbc", false);
    if (doc->getNumErrors() > 0){
        std::string error_msg;
        for (unsigned int i = 0; i < doc->getNumErrors(); ++i) {
        error_msg += doc->getError(i)->getMessage() + "\n";
        }
        std::cerr << "Error reading SBML file: " << fpath << "\n" << error_msg;
        return;
    }

    model = doc->getModel();
    FbcModelPlugin * fbcModel = static_cast<FbcModelPlugin*>(model->getPlugin("fbc"));
    
    int nMets = model->getNumSpecies();
    int nReactions = model->getNumReactions();
    std::cout << "Number of metabolites:\t" << nMets << "\n";
    std::cout << "Number of reactions:\t" << nReactions << "\n";

    lp = glp_create_prob();

    // populate exhange metabolites and their indices
    for (const auto& ex_met_id : ex_mets){
        bool found = false;
        for (unsigned int i = 0; i < model->getNumSpecies(); ++i){
            if (model->getSpecies(i)->getId() == ex_met_id){
                exchange_mets.push_back(ex_met_id);
                found = true;
                break;
            }
        }
        if (!found) std::cout << "Warning: " << ex_met_id << " not found in the model." << std::endl;
    }
    ex_rxn_indices.resize(exchange_mets.size(), 0);
    
    // Create LP rows that correspond to metabolites
    glp_add_rows(lp, model->getNumSpecies());
    for (unsigned int i = 0; i < model->getNumSpecies(); ++i) {
        Species* species = model->getSpecies(i);
        glp_set_row_name(lp, i+1, species->getId().c_str());
        glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
        met_indices[species->getId()] = i+1; 
    }

    // stoichiometric matrix variables
    std::vector<int> ia, ja;
    std::vector<double> ar;
    ia.push_back(0);
    ja.push_back(0);
    ar.push_back(0.0);
    int rowIndex;

    // create reaction columns and propagate stochiometric matrix
    Reaction * rxn;
    FbcReactionPlugin * fbcRxn;
    Parameter *lbparam, *ubparam;
    double lowerBound, upperBound;
    glp_add_cols(lp, model->getNumReactions());
    for (unsigned int i = 1; i <= model->getNumReactions(); ++i){
        rxn = model->getReaction(i-1);
        if (rxn==nullptr) {
            std::cout << "Warning: Reaction at index " << i << " is null. Skipping." << std::endl;
            continue;
        }
        fbcRxn = static_cast<FbcReactionPlugin*>(rxn->getPlugin("fbc"));
        glp_set_col_name(lp, i, rxn->getId().c_str());
        rxn_indices[rxn->getId()] = i; 

        lbparam = model->getParameter(fbcRxn->getLowerFluxBound());
        lowerBound = lbparam ? lbparam->getValue() : 0.0;
        ubparam = model->getParameter(fbcRxn->getUpperFluxBound());
        upperBound = ubparam ? ubparam->getValue() : 0.0;

        if (lowerBound == upperBound) glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
        else glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);

        // populate the stoichiometric matrix with reactants
        for (unsigned int j=0; j<rxn->getNumReactants(); ++j) {
            SpeciesReference* reactant = rxn->getReactant(j);
            std::string metId = reactant->getSpecies();
            if (met_indices.find(metId) != met_indices.end()) {
                rowIndex = met_indices[metId];
                ia.push_back(rowIndex);
                ja.push_back(i);
                ar.push_back(-reactant->getStoichiometry());
            }
        }
        
        // populate the stoichiometric matrix with products
        for (unsigned int j=0; j<rxn->getNumProducts(); ++j) {
            SpeciesReference* product = rxn->getProduct(j);
            std::string metId = product->getSpecies();
            if (met_indices.find(metId) != met_indices.end()) {
                rowIndex = met_indices[metId];
                ia.push_back(rowIndex);
                ja.push_back(i);
                ar.push_back(product->getStoichiometry());
            }
        }

        // Find reactions that contain exchange metabolites
        if (startsWith(rxn->getId(), "R_EX")){
            std::string reactant = rxn->getReactant(0)->getSpecies();
            for (long unsigned int j = 0; j < exchange_mets.size(); ++j) {
                if (exchange_mets[j] == reactant){
                    if (ex_rxn_indices[j] != 0 ) std::cout << "Warning: Multiple exchange reactions found for metabolite " 
                                  << reactant << ". Using the first one.\n";
                    else ex_rxn_indices[j] = i ;
                }
            }
        }
    }

    glp_load_matrix(lp, ia.size() - 1, &ia[0], &ja[0], &ar[0]);

    // Set the objective function
    if (fbcModel->getNumObjectives() == 0) throw std::runtime_error("No objectives set in the model.");
    Objective* obj = fbcModel->getObjective(0);
    const ListOfFluxObjectives* fluxObjectives = obj->getListOfFluxObjectives();
    const FluxObjective* fluxObj;
    std::string reactionId;
    double coeff;
    for (unsigned int i = 0; i < fluxObjectives->size(); ++i){
        fluxObj = fluxObjectives->get(i);
        reactionId = fluxObj->getReaction();
        coeff = fluxObj->getCoefficient();
        if (rxn_indices.find(reactionId) != rxn_indices.end()) {
            int colIndex = rxn_indices[reactionId];
            glp_set_obj_coef(lp, colIndex, coeff);
        } else {
            std::cerr << "Warning: Reaction " << reactionId << " not found in the model.\n";
        }
    }
    glp_set_obj_dir(lp, GLP_MAX);
    
    // solve the LP problem
    glp_simplex(lp, NULL);

    // print objective value
    std::cout << "Objective value: " << glp_get_obj_val(lp) << std::endl;
    
    // find exchange reactions and print their fluxes
    double flux;
    for (long unsigned int i = 0; i < exchange_mets.size(); ++i) {
        if (ex_rxn_indices[i] != 0) {
            flux = glp_get_col_prim(lp, ex_rxn_indices[i]);
            std::cout << "Exchange reaction for " << exchange_mets[i] << ": " << flux << std::endl;
        } else {
            std::cout << "No exchange reaction found for " << exchange_mets[i] << std::endl;
        }
    }

    delete doc;
    doc = nullptr;
}

void FixGrowthGSMM::setExchangeConc(const std::string& ex_met_id, const double& conc) {
    auto it = find(exchange_mets.begin(), exchange_mets.end(), ex_met_id);
    if (it == exchange_mets.end()) {
        std::cerr << "Error: Exchange metabolite " << ex_met_id << " not found in the model." << std::endl;
        return;
    }
    int rxn_index = ex_rxn_indices[std::distance(exchange_mets.begin(), it)];
    if (rxn_index == 0) std::cerr << "Error: No exchange reaction found for metabolite " << ex_met_id << "." << std::endl;
    else glp_set_col_bnds(lp, rxn_index, GLP_FX, -conc, 1000);
        
}

void FixGrowthGSMM::solve(){
    glp_simplex(lp, NULL);
}

double FixGrowthGSMM::getBiomass() {
    return glp_get_obj_val(lp);
}

double FixGrowthGSMM::getExchangeFlux(const std::string& ex_met_id) {
    auto it = find(exchange_mets.begin(), exchange_mets.end(), ex_met_id);
    if (it == exchange_mets.end()) {
        std::cerr << "Error: Exchange metabolite " << ex_met_id << " not found in the model." << std::endl;
        return 0.0;
    }
    int rxn_index = ex_rxn_indices[std::distance(exchange_mets.begin(), it)];
    if (rxn_index == 0) {
        std::cerr << "Error: No exchange reaction found for metabolite " << ex_met_id << "." << std::endl;
        return 0.0;
    }
    else return glp_get_col_prim(lp, rxn_index);
}

FixGrowthGSMM::~FixGrowthGSMM() {
    if (lp) {
        glp_delete_prob(lp);
        lp = nullptr;
    }
}

void FixGrowthGSMM::update_atoms(){
    std::cout << "FixGrowthGSMM::update_atoms() NOT IMPLEMENTED!" << std::endl;
}

void FixGrowthGSMM::update_cells(){
    std::cout << "FixGrowthGSMM::update_cells() NOT IMPLEMENTED!" << std::endl;
}