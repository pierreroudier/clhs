/*
 * Kiri Daust, August 2020
 * 
 * This program contains the main function (and helper functions)
 * used for the clhs optimisation if use.cpp = T
 * 
 * The method is based on the R code by Pierre Roudier
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <RcppArmadillo.h>
#include <random>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp17")]]

//difference in vectors (i.e. setdiff in R)
std::vector<int> vector_diff( const std::vector<int>& model, const std::vector<int>& pattern ){
  std::set<int> s_model( model.begin(), model.end() );
  std::set<int> s_pattern( pattern.begin(), pattern.end() );
  std::vector<int> result;
  
  std::set_difference( s_model.begin(), s_model.end(), s_pattern.begin(), s_pattern.end(),
                       std::back_inserter( result ) );
  
  return result;
}

//The following functions are used to quickly calculate the covariance matrix
//Note that only the lower triangular part is returned
struct asset_info {
  double sum, sum2, stdev;
};

// cor = ( n * sXY - sX * sY ) / ( sqrt(n * sX2 - sX^2) * sqrt(n * sY2 - sY^2) )
inline asset_info compute_asset_info(const NumericMatrix& mat, 
                                     const int icol, const int rstart, const int rend) {
  double sum, sum2;
  sum = sum2 = 0;
  
  for (int r = rstart; r < rend; r++) {
    double d = mat(r, icol);
    sum += d;
    sum2 += pow(d,2);
  }
  
  asset_info res;
  res.sum = sum;
  res.sum2 = sum2;
  res.stdev = sqrt((rend-rstart) * sum2 - pow(sum, 2));
  return res;
}

inline NumericMatrix c_cor_helper(const NumericMatrix& mat, const int rstart, const int rend) {
  int nc = mat.ncol();
  int nperiod = rend - rstart;
  double tempCor;
  NumericMatrix rmat(nc, nc);
  
  vector<asset_info> info(nc);
  for (int c = 0; c < nc; c++)
    info[c] = compute_asset_info(mat, c, rstart, rend);
  
  for (int c1 = 0; c1 < nc; c1++) {
    for (int c2 = 0; c2 < c1; c2++) {
      double sXY = 0;
      
      for (int r = rstart; r < rend; r++)
        sXY += mat(r, c1) * mat(r, c2);
      
      tempCor = (nperiod * sXY - info[c1].sum * info[c2].sum) / (info[c1].stdev * info[c2].stdev);
      if(isnan(tempCor)){
        tempCor = 0;
      }
      rmat(c1, c2) = tempCor;
    }
  }
  
  return rmat;
}

//c_cor is called to calculate the correlation matrix of mat
NumericMatrix c_cor(NumericMatrix mat) {
  return c_cor_helper(mat, 0, mat.nrow());
}

//pair list to store frequency table for factors
typedef std::pair<double, int>  ptype;

//this function creates a frequency table of sampled factor data
//it takes a vector v of sampled data, and a frequency table for the full factor
//it returns the absolute values of the original table - the sampled frequencies
NumericVector table_cpp(const Rcpp::NumericVector & v, const NumericVector full){
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
  
  result_vals = abs(result_vals - full);
  return (result_vals);
}

//bincount - basically just a histogram function
//takes vector of data and breaks, returns a vector of counts
IntegerVector hist(NumericVector x, NumericVector breaks){ //based on C_bincount from graphics package
    int n = x.length();
    int nb = breaks.length();
    int nb1 = nb-1;
    int i,lo,hi,newVal;
    
    IntegerVector counts(nb1);
    //Rcout << "In Hist function \n";
    
    for(i = 0 ; i < n ; i++){
      lo = 0;
      hi = nb1;
      if(breaks[lo] <= x[i] &&
         (x[i] < breaks[hi] || (x[i] == breaks[hi]))) {
        while(hi-lo >= 2) {
          newVal = (hi+lo)/2;
          if(x[i] > breaks[newVal]){
            lo = newVal;
          }else{
            hi = newVal;
          }
        }
        counts[lo]++;
      }
    }
    return(counts);
  }

//structure to store objective function result
struct objResult {
  double objRes;
  std::vector<double> obj_cont_res;
};

//objective function - calculates objective value for current sample
//returns objResult struct
objResult obj_fn(arma::mat x, NumericMatrix strata, arma::mat include, bool factors, 
                 arma::uvec i_fact, NumericMatrix cor_full, Rcpp::List fact_full, 
                 double wCont, double wFact, double wCorr, arma::mat etaMat){
  int nsamps = x.n_rows;
  arma::mat x_all = join_vert(x,include);//join with include - does nothing if no include
  NumericMatrix fact_all;
  
  //if there are factors, remove them from the continuous data
  if(factors){
    arma::mat tempMat = x_all.cols(i_fact);
    fact_all = wrap(tempMat);
    x_all.shed_cols(i_fact);
  }
  
  int num_vars = x_all.n_cols;
  int num_obs = strata.nrow();
  NumericVector hist_cnt;
  IntegerMatrix hist_out(num_obs - 1, num_vars);
  NumericVector data;
  NumericVector data_orig;
  NumericVector strata_curr;
  NumericVector hist_temp;
  NumericVector obj_cont;
  NumericMatrix t2;
  std::vector<double> obj_cont2;
  
  //send each column to bincount function
  for(int i = 0; i < num_vars; i++){
    //Rcout << "Hist idex is " << i << "\n";
    data = wrap(x_all.col(i));
    strata_curr = strata(_,i);
    hist_out(_,i) = hist(data,strata_curr);
  }
  
  //Rcout << "Histout: " << hist_out << "\n";
  //convert to arma mat because subtraction is faster
  arma::mat hist2 = as<arma::mat>(hist_out);
  t2 = wrap(arma::abs(hist2 - etaMat));//subtract eta - either input matrix, or all 1

  obj_cont = rowSums(t2);
  //Rcout << "Full ObjCont: " << obj_cont << "\n";
  
  //send factor data to get tabulated
  NumericVector factRes(num_vars);
  if(factors){
    int num_vars2 = fact_full.size();
    NumericVector temp(x_all.n_rows);
    double total;
    for(int i = 0; i < num_vars2; i++){
      temp = fact_full[i];
      total = sum(temp);
      temp = temp/total;
      factRes[i] = sum(table_cpp(fact_all(_,i),temp));
    }
  }
  
  //Rcout << "FactRes: " << factRes << "\n";
  //correlation matrix for current sample
  NumericMatrix cor_new = c_cor(wrap(x_all));

  double obj_cor = sum(abs(cor_full - cor_new));
  //combined objective values - since corr_mat is lower tri, have to multiply by 2
  double objFinal = sum(obj_cont)*wCont + obj_cor*2*wCorr + sum(factRes)*wFact;
  //Rcout << "FinalObj" << objFinal << "\n";
  obj_cont = obj_cont[Rcpp::Range(0,nsamps-1)];
  obj_cont2 = as<std::vector<double>>(obj_cont);
  struct objResult out = {objFinal, obj_cont2};
  return(out);
}


//' This is the internal Cpp function used to run the metropolis hasting algorithm if use.cpp = T. 
//' In general, it shouldn't be used as a stand alone function, because some preprocessing is done in R
//' 
//' @param xA matrix of data - must be numeric (factors are converted to numeric in R)
//' @param cost cost vector (0 if no cost)
//' @param strata matrix of continuous strata
//' @param include matrix of included data
//' @param idx integer vector of rows from which sampling is allowed
//' @param factors boolean factor flag
//' @param i_fact indices of factors in xA
//' @param nsample number of samples
//' @param cost_mode bool cost flag
//' @param iter number of iterations
//' @param wCont continuous weight
//' @param wFact factor weights
//' @param wCorr correlation weights
//' @param etaMat eta matrix - either all 1, or user input
//' @param temperature initial temperature
//' @param tdecrease temperature decrease every length_cycle iterations
//' @param length_cycle number of iterations between temperature decrease
//' @return list with sampled data, indices, objective values, cost value, and final continuous objectives for each strata

// [[Rcpp::export]]
List CppLHS(arma::mat xA, NumericVector cost, NumericMatrix strata, 
            arma::mat include, std::vector<int> idx, bool factors, arma::uvec i_fact, 
            int nsample, bool cost_mode, int iter, double wCont,
            double wFact, double wCorr, arma::mat etaMat,
            double temperature, double tdecrease, int length_cycle){
  
  //initialise objects
  double prev_obj;
  double obj;
  double delta_obj;
  double metropolis = 1.0;
  double prev_opCost = 0;
  double opCost = 0;
  double metropolis_cost;
  double delta_cost;
  IntegerVector prev_sampled;
  IntegerVector prev_unsampled;
  NumericVector prev_contObj;
  arma::mat x_curr;
  
  std::vector<double> delta_cont;
  std::vector<double> delta_cont_prev;
  IntegerVector idx_removed;
  IntegerVector idx_added;
  int spl_removed;
  std::vector<int> i_sampled;
  std::vector<int> i_sampled_prev;
  std::vector<int> i_unsampled;
  std::vector<int> i_unsampled_prev;
  //std::vector<int> idx(ndata);
  //std::iota(idx.begin(),idx.end(),0);//populate idx with seq(0:ndata)
  std::vector<double>::iterator it_worse;
  int i_worse;
  IntegerVector idx2;
  NumericVector temp;
  NumericVector obj_values(iter);
  
  //initial sample
  idx2 = wrap(idx);
  i_sampled = as<std::vector<int>>(Rcpp::sample(idx2,nsample,false));
  i_unsampled = vector_diff(idx,i_sampled);
  arma::uvec arm_isamp = arma::conv_to<arma::uvec>::from(i_sampled);//index needs to be aram::uvec to extract rows
  x_curr = xA.rows(arm_isamp); // is this efficient?
  //Rcout << "Finish initial sample \n";
  arma::mat x_cont = xA;
  if(factors){
    x_cont.shed_cols(i_fact);// remove factors for correlation
  }
  NumericMatrix cor_full = c_cor(wrap(x_cont));//full correlation matrix
  
  //setup table list for factors
  Rcpp::List factTab(i_fact.size()); 
  if(factors){
    arma::mat xFact = xA.cols(i_fact);
    NumericMatrix xFact2 = wrap(xFact);
    for(unsigned int i = 0; i < i_fact.size(); i++){
      factTab[i] = Rcpp::table(xFact2(_,i));//we use Rcpp here because it gives us name attributes
    }
  }
  
  //initial objective value
  struct objResult res = obj_fn(x_curr,strata,include,factors,i_fact,cor_full,factTab,wCont,wFact,wCorr,etaMat);
  obj = res.objRes;
  delta_cont = res.obj_cont_res;
  //Rcout << "Finished function; obj = "<< obj << "\n";
  
  
  NumericVector cost_values(iter, NA_REAL);//store cost
  
  if(cost_mode){
    idx2 = wrap(i_sampled);
    temp = cost[idx2];
    opCost = sum(temp);
  }
  
  //start metropolis hasting iterations
  for(int i = 0; i < iter; i++){
    if(i % 50 == 0)
      Rcpp::checkUserInterrupt();
    prev_obj = obj;
    i_sampled_prev = i_sampled;
    i_unsampled_prev = i_unsampled;
    delta_cont_prev = delta_cont;
    //Rcout << "Deltacont: " << wrap(delta_cont) << "\n";
    
    if(cost_mode){
      prev_opCost = opCost;
    }
    
    if(Rcpp::runif(1,0,1)[0] < 0.5){
      //Rcout << "In Simple Swap \n";
      idx_removed = Rcpp::sample(i_sampled.size()-1, 1, false);
      spl_removed = i_sampled[idx_removed[0]];
      i_sampled.erase(i_sampled.begin() + idx_removed[0]);
      idx_added = Rcpp::sample(i_unsampled.size()-1, 1, false);
      i_sampled.push_back(i_unsampled[idx_added[0]]);
      i_unsampled.erase(i_unsampled.begin()+idx_added[0]);
      i_unsampled.push_back(spl_removed);
    }else{
      //Rcout << "In remove worst \n";
      it_worse = std::max_element(delta_cont.begin(),delta_cont.end()); //returns max element
      i_worse = std::distance(delta_cont.begin(), it_worse); //find location of worst
      //NumericVector temp42 = wrap(delta_cont);
      //NumericVector temp43 = wrap(i_sampled);
      // Rcout << "i_sampled: " << temp43 << "\n";
      // Rcout << "Delta vec: " << temp42 << "\n";
      // Rcout << "Delta Worst: " << i_worse << "\n";
      spl_removed = i_sampled[i_worse]; //this is the problem
      //Rcout << "spl_removed " << spl_removed << "\n";
      idx_added = Rcpp::sample(i_unsampled.size()-1, 1, false);
      i_sampled.erase(i_sampled.begin()+i_worse);
      i_sampled.push_back(i_unsampled[idx_added[0]]);
      i_unsampled.erase(i_unsampled.begin()+idx_added[0]);
      i_unsampled.push_back(spl_removed);
      // temp43 = wrap(i_sampled);
      // Rcout << "i_samplednew: " << temp43 << "\n";
    }
    //Rcout << "sending to function \n";
    arm_isamp = arma::conv_to<arma::uvec>::from(i_sampled);
    x_curr = xA.rows(arm_isamp); // is this efficient?
    //calculate objective value
    res = obj_fn(x_curr,strata,include,factors,i_fact,cor_full,factTab,wCont,wFact,wCorr,etaMat); //test this
    obj = res.objRes;
    delta_cont = res.obj_cont_res;
    NumericVector dcTemp = wrap(delta_cont);
    //update variables
    delta_obj = obj - prev_obj;
    metropolis = exp(-1*delta_obj/temperature);
    
    if(cost_mode){//update cost variables
      idx2 = wrap(i_sampled);
      temp = cost[idx2];
      opCost = sum(temp);
      delta_cost = opCost - prev_opCost;
      metropolis_cost = exp(-1*delta_cost/temperature);
    }else{
      metropolis_cost = R_PosInf;
    }
    
    //Revert Change
    if((delta_obj > 0 && runif(1,0,1)[0] >= metropolis) || runif(1,0,1)[0] >= metropolis_cost){
      //Rcout << "In revert change \n";
      i_sampled = i_sampled_prev;
      i_unsampled = i_unsampled_prev;
      obj = prev_obj;
      delta_cont = delta_cont_prev;
      if(cost_mode){
        opCost = prev_opCost;
      }
    }
    obj_values[i] = obj;
    if(cost_mode){
      cost_values[i] = opCost;
    }
    if(i % length_cycle == 0){
      temperature = temperature*tdecrease;
    }
  }
  arm_isamp = arma::conv_to<arma::uvec>::from(i_sampled);
  x_curr = xA.rows(arm_isamp);
  Rcout << "vroom vroom \n";
  //create list to return
  return List::create(_["sampled_data"] = x_curr,
                      _["index_samples"] = i_sampled,
                      _["obj"] = obj_values,
                      _["cost"] = cost_values,
                      _["final_obj_continuous"] = delta_cont);
}