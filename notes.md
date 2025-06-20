# Architecture 

This is akin to a model-view-controller pattern. There are four main classes `grid`, `atom_vec`, `fix` and `compute`. `grid` and `atom_vec` should be seen as the model, their attributes encode the state of the model. `fix` should be seen as a controller, it uses and modifies attributes in `grid` and `atom_vec`. `compute` reads attributes from `grid` and `atom_vec` and outputs the data. (I still don't exactly know how computes work)  

For the purpose of integrating GSMM into NUFEB we will create a new `grid`, `atom_vec` and `fix` class. `grid_vec_gsmm` will inherit from `grid_vec_chemostat` to include attributes that track the new nutrients. `fix_growth_gsmm` will inherit from `fix_growth` and contain the metabolic model corresponding to a species. This design choice requires each new species to have it's own fix with the associated metabolic model. Lastly, we will have `atom_vec_coccus_gsmm` to inherit from `atom_vec_coccus` and would additionally include information for gene knockout. 

Echange reactions in the GSMM .xml have reaction id `R_EX_<met_id>_` with stoichiometry `<met_id>_e <->` where the metabolite to be supplies is a reactant.  

TODO (Jun 18):
1. Incorporate FBA computes into `fix_growth_gsmm`
2. Specify nutrients to be exchanged
