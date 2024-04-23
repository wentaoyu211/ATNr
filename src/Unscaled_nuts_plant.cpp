#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
using namespace arma;

// generate documentation with:
// Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
// roxygen2::roxygenize(roclets="rd")  # this updates the documentation based on roxygen comments

//' @name Unscaled_nuts_plant
//' @title Store parameters and functions associated to the unscaled version of ATN including nutrient dynamics
//' @description Type the name of the class to see its methods
//' @field nb_b Number of basal species (only plants)
//' @field nb_n Number of nutrient pool
//' @field X vector of metabolic rates (length = number of species)
//' @field K matrix of plant nutrient efficiencies (dim = number of nutrients * number of plants)
//' @field V matrix of plant relative nutrient content (dim = number of nutrients * number of plants)
//' @field S Vector of maximum nutrient concentration (length = number of plants)
//' @field r Vector of maximum growth rate of plant species (length = number of plant species)
//' @field BM Vector of body masses (length = number of species)
//' @field dB Vector of local derivatives (length = number of species)
//' @field faci Matrix of plant-plant facilitation schemes
//' @field ext Extinction threshold for species
//' @field theta the exponent of neighbour tree biomass to affect the realized growth
//' @field ODE Calculate the derivatives for the scaled version of the ATN model \itemize{
//' \item Parameter: bioms -  Local species biomasses
//' \item Parameter: t - Integration time point
//' \item Returns a vector of growth rate for each species at time t
//' }

class Unscaled_nuts_plant{ // adding facilitation process to G:realized growth rate
public:
  int nb_b; // number of basal species
  int nb_n; // number of nutrients
  int n_tot; // bool prefs MORE PRECISE?
  // global nutrient turn over rate (rate of replenishment)
  // used in calculating change in nutrient concentration
  double D;

  double ext;
  
  double theta; // (first try the empirical exponent 0.15?)

  vec X; // metabolic rates
  vec total_X; // population metabolism
  vec r; // growth rates of plants
  vec S; // maximal nutrient level
  //vec out_fluxes; // out fluxes for all species
  // body masses
  vec BM ;

  
  // biomasses
  vec bioms;
  // r*G for plants, as there is no need to compute it at each ODE call
  
  // NumericVector test;
  // vec p;
  // Coltor of derivatives
  vec dB;
  

  // plants nutrient uptake efficiency (K(i,j): plant i on nutrient j)
  mat K; //!!!!!!!! change to same row*col in inits in R
  // relative content in the plant species' biomass
  mat V; //rows = nuts, cols = plants
  
  // adding  a facilitation matrix with dimension nb_b*nb_b
  mat faci;

  // internal variables for optimisation
  vec G; // species specific growth factor Gi
  vec G_faci; // realized growth factor with facilitation
  
  //vec low;
  vec uptake;
  // vec bioms_non_nut;
  vec theta_bioms;


  // iterators
  int res;
  vec::iterator res_end;

  // slicers
  uvec plants;
  uvec nut;
  //uvec non_nut;
  uvec extinct;
  
  Unscaled_nuts_plant(int nb, int nn): // should probably add facilitation matrix to the argument???
     nb_b(nb), nb_n(nn) {
      n_tot = nb_b + nb_n;

	  
      // initialise vectors to 0
      X.zeros(nb_b);
      r.zeros(nb_b);
      S.zeros(nb_n);
      G.zeros(nb_b);
	  G_faci.zeros(nb_b);
      dB.zeros(n_tot);
      uptake.zeros(nb_b);
	  theta_bioms.zeros(nb_b);
      BM.ones(nb_b);


      // initialise matrices
      K.zeros(nb_n, nb_b);
      V.zeros(nb_n, nb_b);
	  faci.zeros(nb_b, nb_b);

      // scalars
      D = 0.0;
      theta = 0.0;

      // iterator
     // res_end = G.end();

      // slicers TODO check optional argument N, =100 by default
      // but inconsistant with this: https://stackoverflow.com/questions/25790634/how-to-create-a-vector-from-1n-in-c-armadillo
      plants = linspace<uvec>(nb_n, nb_n + nb_b-1, nb_b);
      nut = linspace<uvec>(0, nb_n-1, nb_n);
      //non_nut = linspace<uvec>(nb_n, n_tot-1, nb_s);
    }
// the following 5 lines are for error detection
//template<typename T>
//T mod(T a, int n)
//{
//    return a - floor(a/n)*n;
//} 
// ----

  void initialisations(){

  }
   
   
  void print(){
    Rcpp::Rcout << "nb_n:"  << std::endl << nb_n << std::endl;
    Rcpp::Rcout << "nb_b:"  << std::endl << nb_b << std::endl; 
    Rcpp::Rcout << "bioms: " << bioms << std::endl; 
    Rcpp::Rcout << "G: " << G << std::endl; 
    // Rcout << " prey" << prey << std::endl;
  }
  
  
  
  // NumericVector ODE(double t, NumericVector bioms, NumericVector p){  // for sundials
  vec ODE(vec bioms, double t){ // for odeintr
    // here, the matricies of attack rates, handling times, feeding rates ...
    // are of diemnsions: nb species * nb consumers
    // so non square matrices
    // this is because I removed all the elements relative to a plant on a plant
    // it implies less claculations for the feeding rates
    // but more subsetting afterwards. 
    // not sure of the best option though 
    // (would it not be faster to use square matrices everywhere?)
    
    // Rcpp::Rcout << "aaa " << bioms.n_elem << std::endl;

    extinct = find(bioms < ext);
    uvec not_nut_extinct = find(extinct >= nb_n);
    extinct = extinct(not_nut_extinct);
    bioms.elem(extinct).zeros();
	
    // Rcpp::Rcout << bioms.t()  << std::endl;
    // format the output 
	//bioms_non_nut = bioms.elem(non_nut);
    //pow_bioms.each_col() = bioms_non_nut;
    //pow_bioms = pow(pow_bioms.each_row(), q.t());
	 theta_bioms = pow(bioms(plants), theta);

    // calculate values for feeding rates
    // F contains first the upper part of the feeding rates
    // F = wb_mat % pow_bioms;

    // wbh_mat*bioms: gives for each consumer i the sum over prey j of
    // wij*hij*bij*Bj
    //low = sum(wbh_mat%pow_bioms,0).t() + c%bioms(animals) + 1;
    //low = low%BM(animals-nb_n);

    //F.each_row() /=low.t();

    // out_fluxes: sum of out flux for each resource species, col vector
    //out_fluxes = F*bioms(animals);

    // realised met. rate
    total_X = X%bioms(plants);

    // species specific growth factor for each plant
    // iterators are pointers, so not good here as I have to access to 
    // the ith element of different vectors
    for (res = 0; res != nb_b; ++res){
      
	  G(res) = min(bioms(nut) / (K.col(res) + bioms(nut)));	  
	  
    }
    // G could be calculated using matrix inversion, 
    // not sure there is much to win here though
    // Cholesky algs are in general in O(n^3)
    // KandBioms = k.each_col() + bioms(nut);
    // G = min(inv(KandBioms.each_col() / bioms(nut)), DIMENSION)
     
	// the total facilitation effect is the sum of all the facilitator's 
	// facilitation capacity (reflected by facilitation matrix) * its own biomass to an exponent theta
	// if elements in faci matrix are all 0, it's reduced to nutrient model
  	G_faci = G % (1+ faci * theta_bioms);
	 //Rcpp::Rcout << bioms.t()  << std::endl;
	 //Rcpp::Rcout << t  << std::endl;
    //if (mod(t, 100) <2){
      
    //}

    // plant uptake (use updated realised G_faci:
    uptake = r%bioms(plants)%G_faci;
	
    // derivatives for non nutrients
    dB(plants) = uptake - total_X(plants-nb_n);
 
     
    // note: e may be preintegrated as constant over time

    dB.elem(extinct).zeros();

    // derivate for nutrients
    dB(nut) = D * (S - bioms(nut)) - V*uptake;
    
    return(dB);
  }
  
};



RCPP_MODULE(Unscaled_nuts_plantModule){
  using namespace Rcpp;
  class_<Unscaled_nuts_plant>("Unscaled_nuts_plant")
    .constructor<int, int>("constructor") //constructor
    .method("print", &Unscaled_nuts_plant::print)
    .method("ODE", &Unscaled_nuts_plant::ODE)
    .method("initialisations", &Unscaled_nuts_plant::initialisations)
    .field("nb_b", &Unscaled_nuts_plant::nb_b)
    .field("nb_n", &Unscaled_nuts_plant::nb_n)
    .field("BM", &Unscaled_nuts_plant::BM)
	.field("faci", &Unscaled_nuts_plant::faci)
	.field("theta", &Unscaled_nuts_plant::theta)
    .field("K", &Unscaled_nuts_plant::K)
    .field("D", &Unscaled_nuts_plant::D)
    .field("S", &Unscaled_nuts_plant::S)
    .field("r", &Unscaled_nuts_plant::r)
    .field("X", &Unscaled_nuts_plant::X)
    .field("V", &Unscaled_nuts_plant::V)
    .field("dB", &Unscaled_nuts_plant::dB)
    .field("ext", &Unscaled_nuts_plant::ext)
    ;  
}

