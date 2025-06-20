#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#include <sbml/packages/fbc/extension/FbcExtension.h>
#include <sbml/packages/fbc/extension/FbcModelPlugin.h>
#include <sbml/packages/fbc/extension/FbcReactionPlugin.h>
#include <sbml/packages/fbc/sbml/Objective.h>
#include <glpk.h>

#include <iostream>

#include "MetModel.h"

MetModel::MetModel(const std::string &fname, 
        const std::vector<std::string> &eComp, 
        const std::vector<std::string> &cComp) {
    SBMLDocument* document;
    SBMLReader reader;
    document = reader.readSBML(fname);
    if (document->getNumErrors() > 0) {
        document->printErrors(std::cerr);
        return;
    }

    // Enable FBC package
    document->enablePackage(FbcExtension::getXmlnsL3V1V2(), "fbc", true);
    document->setPackageRequired("fbc", false);
    if (document->getNumErrors() > 0) {
        document->printErrors(std::cerr);
        return;
    }

    // Get FBC plugin
    Model* model = document->getModel();
    FbcModelPlugin* fbcModel = static_cast<FbcModelPlugin*>(model->getPlugin("fbc"));
    if (fbcModel != nullptr) 
        std::cout << "FBC Objectives: " << fbcModel->getNumObjectives() << std::endl;

    // printModelInfo(model);
    
    parseMetabolites(model, eComp, cComp);
    std::cout << "Number of metabolites: " << metabolites.size() << std::endl;
    
    parseReactions(model);
    std::cout << "Number of reactions: " << reactions.size() << std::endl;

    parseObjective(fbcModel, 0);
    std::cout << "Num reactions in objective: " << objectiveRxnIds.size() << std::endl;


    // create glpk instance
    ia.push_back(0); // GLPK uses 1-based indexing, so we start with 0
    ja.push_back(0); 
    ar.push_back(0.0); 

    lp = glp_create_prob();
    glp_set_prob_name(lp, model->getName().c_str());
    createLP();
    printLPInfo();
    glp_simplex(lp, NULL);

    // print simplex solved flux values
    std::cout << "Optimal value: " << glp_get_obj_val(lp) << std::endl;
    std::cout << "Simplex solved flux values:" << std::endl;
    for (int i = 0; i < reactions.size(); ++i) {
        std::shared_ptr<Reactn> rxn = reactions[i];
        int rxnIdx = i + 1; // GLPK indices are 1-based
        double fluxValue = glp_get_col_prim(lp, rxnIdx);
        std::cout << "Reaction: " << rxn->id 
                  << ", Flux: " << fluxValue << std::endl;
    }
    
    
    delete document;

}


void MetModel::createLP(){
    // create metabolite rows
    glp_add_rows(lp, metabolites.size());
    std::cout << "Creating " << metabolites.size() 
              << " metabolite rows in GLPK." << std::endl;
    for (int i=0; i<metabolites.size(); ++i){
        glp_set_row_name(lp, i+1, metabolites[i]->id.c_str());
        glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
    }

    // create reaction columns
    glp_add_cols(lp, reactions.size());
    double ub, lb;
    for (int i=0; i<reactions.size(); ++i){
        glp_set_col_name(lp, i+1, reactions[i]->id.c_str());
        lb = reactions[i]->lowerBound;
        ub = reactions[i]->upperBound;
        if (lb == ub){
            glp_set_col_bnds(lp, i+1, GLP_FX, lb, ub);
        } else{
            glp_set_col_bnds(lp, i+1, GLP_DB, lb, ub);
        }
    }

    // fill in the stoichiometric matrix
    int cnt = 0;
    std::shared_ptr<Reactn> rxn;
    for (int i=0; i<reactions.size(); ++i){
        rxn = reactions[i];

        // reactants of this reactions
        for (unsigned int j = 0; j < rxn->reactants.size(); ++j) {
            std::shared_ptr<Metabolite> met = rxn->reactants[j];
            double coeff = rxn->reactantCoeffs[j];            
            if (met == nullptr) continue;
            int metIdx = getMetIdx(met->id);
            if (metIdx < 0) {
                std::cerr << "Error: Metabolite " << met->id << " not found in the model." << std::endl;
                continue;
            }
            ia.push_back(metIdx + 1); // row index (metabolite)
            ja.push_back(i + 1); // column index (reaction)
            ar.push_back(-coeff); // stoichiometric coefficient for reactants
            cnt++;
        }

        // products of this reaction
        for (unsigned int j = 0; j < rxn->products.size(); ++j) {
            std::shared_ptr<Metabolite> met = rxn->products[j];
            double coeff = rxn->productCoeffs[j];
            if (met == nullptr) continue;
            int metIdx = getMetIdx(met->id);
            if (metIdx < 0) {
                std::cerr << "Error: Metabolite " << met->id << " not found in the model." << std::endl;
                continue;
            }
            ia.push_back(metIdx + 1); // row index (metabolite)
            ja.push_back(i + 1); // column index (reaction)
            ar.push_back(coeff); // stoichiometric coefficient for reactants
            cnt++;
        }
    }
    // load the stoichiometric matrix into GLPK
    glp_load_matrix(lp, cnt, ia.data(), ja.data(), ar.data());

    // set the objective function
    glp_set_obj_dir(lp, isObjectiveMaximize ? GLP_MAX : GLP_MIN);
    for (int i=0; i<objectiveRxnIds.size(); ++i) {
        int rxnIdx = objectiveRxnIds[i];
        double coeff = objectiveCoeffs[i];
        std::cout << coeff << " * " 
                  << reactions[rxnIdx]->id << std::endl;
        if (rxnIdx < 0 || rxnIdx >= reactions.size()) {
            std::cerr << "Error: Invalid reaction index " << rxnIdx << " for objective." << std::endl;
            continue;
        }
        glp_set_obj_coef(lp, rxnIdx + 1, coeff);
    }

}

void MetModel::parseObjective(FbcModelPlugin* fbcModel, int n) {
    Objective* obj = fbcModel->getObjective(n);
    if (obj == nullptr) {
        std::cerr << "Number of objectives in model: " << fbcModel->getListOfObjectives()->getNumObjectives()
     << ". Got objective number " << n << std::endl;
     return;
    }

    std::cout << "Objective direction: " << obj->getType() << std::endl;
    if (strcmp(obj->getType().c_str(), "maximize") == 0) {
        isObjectiveMaximize = true;
    } else if (strcmp(obj->getType().c_str(), "minimize") == 0) {
        isObjectiveMaximize = false;
    } else {
        std::cerr << "Unknown objective type: " << obj->getType() << std::endl;
        return;
    }

    // Get the objective reactions and coefficients
    const ListOfFluxObjectives* fluxObjectives = obj->getListOfFluxObjectives();
    for (unsigned int i = 0; i < fluxObjectives->size(); i++) {
        const FluxObjective* fluxObj = fluxObjectives->get(i);
        std::string reactionId = fluxObj->getReaction();
        double coefficient = fluxObj->getCoefficient();
        
        int rxnIdx = getRxnIdx(reactionId);
        if (rxnIdx < 0) {
            std::cerr << "Error: Reaction " << reactionId 
                      << " not found in the model." << std::endl;
            continue;
        }
        objectiveRxnIds.push_back(rxnIdx); // store reaction id
        objectiveCoeffs.push_back(coefficient); // store coefficient
        std::cout << "Objective reaction: " << reactionId 
                  << " with coefficient: " << coefficient << std::endl;
    }
    
}

int MetModel::getMetIdx(const std::string &metId){
    for (int i = 0; i < metabolites.size(); ++i) {
        if (metabolites[i]->id == metId) {
            return i; // return index of metabolite
        }
    }
    return -1;
}

int MetModel::getRxnIdx(const std::string &rxnId) {
    for (int i = 0; i < reactions.size(); ++i) {
        if (reactions[i]->id == rxnId) {
            return i; // return index of reaction
        }
    }
    return -1;
}


void MetModel::parseMetabolites(Model *model, 
                            const std::vector<std::string> &eComp, 
                            const std::vector<std::string> &cComp) {

    metabolites.reserve(model->getNumSpecies());
    for (unsigned int i = 0; i < model->getNumSpecies(); ++i) {
        Species* species = model->getSpecies(i);
        if (species == nullptr) continue;

        std::string id = species->getId();
        std::string name = species->getName();
        Comp compartment = (species->getCompartment() == eComp[0]) ? E : C;
    
        std::shared_ptr<Metabolite> met = std::make_shared<Metabolite>();
        met->id = id;
        met->name = name;
        met->compartment = compartment;

        metabolites.push_back(met); 
        metId2Met[id] = met; // store in map for quick access
    }
}

void MetModel::parseReactions(Model* model) {

    reactions.reserve(model->getNumReactions());
    for (int i = 0; i < model->getNumReactions(); ++i) {

        std::shared_ptr<Reactn> rxn = std::make_shared<Reactn>();

        Reaction* reaction = model->getReaction(i);
        if (reaction == nullptr) continue;
        FbcReactionPlugin* fbcReaction = static_cast<FbcReactionPlugin*>(reaction->getPlugin("fbc"));
        if (fbcReaction == nullptr) continue;

        rxn->id = reaction->getId();
        rxn->name = reaction->getName();
        rxn->reversible = reaction->getReversible();

        double lb = 0;
        Parameter* lbparam = model->getParameter(fbcReaction->getLowerFluxBound());
        if (lbparam != nullptr) lb = lbparam->getValue();
        rxn->lowerBound = lb;

        double ub = 0;
        Parameter* ubparam = model->getParameter(fbcReaction->getUpperFluxBound());
        if (ubparam != nullptr) ub = ubparam->getValue();
        rxn->upperBound = ub;

        SpeciesReference* ref;
        std::string metId;

        int numReactants = reaction->getNumReactants();
        rxn->reactants = std::vector<std::shared_ptr<Metabolite>>(numReactants, nullptr);
        for (unsigned int j = 0; j < numReactants; ++j) {
            ref = reaction->getReactant(j);
            if (ref == nullptr) continue;
            metId = ref->getSpecies();
            if (metId2Met.find(metId) != metId2Met.end()) {
                rxn->reactants[j] = metId2Met[metId];
                rxn->reactantCoeffs.push_back(ref->getStoichiometry());
            } else {
                std::cerr << "Warning: Metabolite " << metId << " not found in metId2met." << std::endl;
            }
        }

        int numProducts = reaction->getNumProducts();
        rxn->products = std::vector<std::shared_ptr<Metabolite>>(numProducts, nullptr);
        for (int j = 0; j < numProducts; ++j) {
            ref = reaction->getProduct(j);
            if (ref == nullptr) continue;
            metId = ref->getSpecies();
            if (metId2Met.find(metId) != metId2Met.end()) {
                rxn->products[j] = metId2Met[metId];
                rxn->productCoeffs.push_back(ref->getStoichiometry());
            } else {
                std::cerr << "Warning: Metabolite " << metId << " not found in metId2met." << std::endl;
            }
        }
        reactions.push_back(rxn);

        if (isExRxn(rxn)){
            exchangeRxnIds.push_back(reactions.size() - 1); // store index of exchange reaction
            // std::cout << "Exchange reaction found (" << rxn->id << "): "
            //     << rxn->reactants[0]-> id << "->" << rxn->products[0]->id << std::endl;
        }
    }
}   

bool MetModel::isExRxn(std::shared_ptr<Reactn> rxn) {
    if (rxn->reactants.size() == 1 && rxn->products.size() == 1){
        if (rxn->reactants[0]->compartment != rxn->products[0]->compartment) {
            return true; 
        } 
    }
    return false;
}

void MetModel::printModelInfo(Model *model) const{
    std::cout << "functionDefinitions: " << model->getNumFunctionDefinitions() << std::endl;
    std::cout << "    unitDefinitions: " << model->getNumUnitDefinitions    () << std::endl;
    std::cout << "   compartmentTypes: " << model->getNumCompartmentTypes   () << std::endl;
    std::cout << "        specieTypes: " << model->getNumSpeciesTypes       () << std::endl;
    std::cout << "       compartments: " << model->getNumCompartments       () << std::endl;
    std::cout << "            species: " << model->getNumSpecies            () << std::endl;
    std::cout << "         parameters: " << model->getNumParameters         () << std::endl;
    std::cout << " initialAssignments: " << model->getNumInitialAssignments () << std::endl;
    std::cout << "              rules: " << model->getNumRules              () << std::endl;
    std::cout << "        constraints: " << model->getNumConstraints        () << std::endl;
    std::cout << "          reactions: " << model->getNumReactions          () << std::endl;
    std::cout << "             events: " << model->getNumEvents             () << std::endl;
}

void MetModel::printLPInfo() {
    std::cout << "==============================================" << std::endl;
    for (int i = 1; i <= reactions.size(); i++) {
        std::cout << "Reaction " << reactions[i-1]->id 
                << " (col " << i << "):" << std::endl;
        
        // Check bounds
        double lb = glp_get_col_lb(lp, i);
        double ub = glp_get_col_ub(lp, i);
        std::cout << "  Bounds: [" << lb << ", " << ub << "]" << std::endl;
        
        // Check coefficients
        int len = glp_get_mat_col(lp, i, NULL, NULL);
        std::vector<int> ind(len + 1);
        std::vector<double> val(len + 1);
        glp_get_mat_col(lp, i, ind.data(), val.data());
        
        std::cout << "  Coefficients:" << std::endl;
        for (int j = 1; j <= len; j++) {
            std::cout << "    Row " << ind[j] << ": " << val[j] << std::endl;
        }
    }
    std::cout << "==============================================" << std::endl;
}

MetModel::~MetModel() {
    if (lp != nullptr) {
        glp_delete_prob(lp); // delete GLPK problem instance
        lp = nullptr;
    }
}

void MetModel::computeFluxDistribution() {}
void MetModel::getFluxDistribution(std::vector<double> &fluxes) const {}