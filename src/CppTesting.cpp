#include <stdio.h>
#include <math.h>
#include <string.h>
#include <RcppArmadillo.h>
#include <random>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp17")]]


// [[Rcpp::export]]
std::vector<int> calc_dist(double min_dist, arma::mat ll_curr) {
  //Now calculate distances if applicable 
  double xdist, ydist;
  double totdist;
  //int numOver;
  int ll_len = ll_curr.n_rows;
  std::vector<int> dist_all(ll_len,0);
  Rcout << "Length: " << ll_len << "\n";
  if(ll_len > 1){
    for(int i = 0; i < ll_len; i++){
      for(int j = i+1; j < ll_len; j++){
        //Rcout << "i: " << i << " j: " << j << "\n";
        xdist = ll_curr(i,0) - ll_curr(j,0);
        ydist = ll_curr(i,1) - ll_curr(j,1);
        totdist = sqrt(pow(xdist,2) + pow(ydist,2));
        if(totdist < min_dist){
          dist_all[i]++;
        }
      }
    }
  }
  return(dist_all);
}

