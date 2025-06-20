#pragma once

#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#include <sbml/packages/fbc/extension/FbcExtension.h>
#include <sbml/packages/fbc/extension/FbcModelPlugin.h>
#include <glpk.h>

#include <vector>
#include <string>
#include <map>
#include <memory>

// declare enums for compartment types
enum Comp{C, E};

struct Metabolite{
    std::string id;
    std::string name;
    Comp compartment;
};

struct Reactn{
    std::string id;
    std::string name;
    std::vector<std::shared_ptr<Metabolite>> reactants;
    std::vector<double> reactantCoeffs; 
    std::vector<std::shared_ptr<Metabolite>> products;
    std::vector<double> productCoeffs; 
    int lowerBound;
    int upperBound;
    bool reversible;
};

class MetModel{
    public:
        std::string *name;
        std::vector<std::shared_ptr<Metabolite>> metabolites;
        std::vector<std::shared_ptr<Reactn>> reactions;
        std::vector<int> exchangeRxnIds; // indices of exchange reactions 
        std::vector<int> objectiveRxnIds; // indices of objective reactions
        std::vector<double> objectiveCoeffs; // coefficients for objective reactions
        bool isObjectiveMaximize;

        glp_prob *lp;
        std::vector<int> ia, ja;
        std::vector<double> ar;

        MetModel(const std::string &, 
            const std::vector<std::string> &eComp, // name of extracellular compartment
            const std::vector<std::string> &cComp); // name of cytosolic compartment
        ~MetModel();

        void computeFluxDistribution();
        void getFluxDistribution(std::vector<double> &fluxes) const;

    private:
        std::map<std::string, std::shared_ptr<Metabolite>> metId2Met;
        std::vector<std::vector<int>> stoichMat;
        std::vector<int> objVec;
        std::vector<double> flux;

        void printModelInfo(Model *) const;
        void parseMetabolites(Model *model, 
                        const std::vector<std::string> &eComp, 
                        const std::vector<std::string> &cComp);
        void parseReactions(Model *model);
        bool isExRxn(std::shared_ptr<Reactn> rxn);
        void parseObjective(FbcModelPlugin* fbcModel, int n);
        int getMetIdx(const std::string &metId);
        int getRxnIdx(const std::string &rxnId);
        void createLP();
        void printLPInfo();
};
