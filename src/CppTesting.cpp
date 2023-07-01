#include <stdio.h>
#include <math.h>
#include <string.h>
#include <RcppArmadillo.h>
#include <random>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp17")]]


// IntegerVector mywhich(const NumericVector v){
//   IntegerVector out = Rcpp::match(v, 2.0);
//   return(out)
// }

//pair list to store frequency table for factors
typedef std::pair<double, int>  ptype;


LogicalVector table_cpp_(const Rcpp::NumericVector & v, const NumericVector full){
  std::vector<double> data = as<std::vector<double>>(v);
  unsigned int nTot = v.size();
  // Create a map
  Rcpp::CharacterVector tempLevs = full.attr("names");//need to made sure new vectors has same levels in same order
  std::vector<double> levels(tempLevs.size());
  for(int i = 0; i < tempLevs.size(); i++){
    levels[i] = std::strtod(tempLevs[i], NULL);
  }

  std::map<double, int> Elt;
  Elt.clear();

  // Initialize with all zero - this might be slow
  for (unsigned int i = 0; i != levels.size(); ++i) {
    Elt[levels[i]] = 0;
  }
  //count frequencies
  for (int i = 0; i != v.size(); ++i) {
    Elt[data[i]] += 1;
  }
  unsigned int n_obs = Elt.size();

  std::vector<ptype> sorted_Elt(Elt.begin(), Elt.end());
  Rcpp::NumericVector result_vals(n_obs);

  unsigned int count = 0;
  double temp;
  //Use iterators to access objects in map
  for(std::vector<ptype>::iterator it = sorted_Elt.begin(); it != sorted_Elt.end(); ++it){
    temp = it->second;
    result_vals(count) = temp/nTot;
    count++;
  }
  // Rcout << "Names: " << tempLevs << "\n";
  // Rcout << "Small: " << result_vals << "\n";
  // Rcout << "Full: " << full << "\n";
  //Rcout << "Curr Facts: " << result_vals << "\n" << "Orig: " << full << "\n";

  //result_vals = abs(result_vals - full);
  result_vals = result_vals - full;
  int maxv = which_max(result_vals);
  double oversampled = strtod(tempLevs[maxv], NULL);
  
  LogicalVector isbad (v.size());
  for (int i = 0; i < v.size(); i++) {
    if(data[i] == oversampled){
      isbad[i] = 1;
    }
  }
  //which samples contribute to oversampled values?
  return (isbad);
}


