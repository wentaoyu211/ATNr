Beta version of the ATNr package that propose solutions to estimate populations dynamics in food webs with Allometric Trophic Networks

# To do

There might be some opportunities to optimise execution time of the armadillo models, especially regarding bracketing of expressions based on matrix algebra. Not important to investigate, but could be nice to have a look at some point. Need a good knowledge of the approach though

Some parameters now use the same values for all nodes it refers to. Maybe we should change that to vectors to make them node specific and having something more general: 
* In Schneider: D (turnover rate of the nutrients) is a single value common to all nutrient pools. 
* In all models, interspecific competition is a scalar. 
 
# Added models

Two models were added to the source code (scr), namely Unscaled_nuts_plant and Unscaled_nuts_faci. The former only modeling plants competition for nutrients (without animals on top), with the possibility to add facilitation schemes to the plants. The latter is the UnScaled_nuts_plant model with animals feeding interactions on top of plants. In the case that the facilitation interaction matrix is 0 in Unscaled_nuts_faci, this model is equal to Unscaled_nuts. 
