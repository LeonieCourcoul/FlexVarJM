#include <RcppArmadillo.h>
#include <stdexcept>
#include <math.h>
#include <stdlib.h> /* srand, rand */
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]

double log_llh_ind(bool variability_hetero, arma::rowvec Otime_i, arma::vec Wtime_i,
             arma::mat Os_i, arma::mat Ws_i, arma::vec omega, arma::mat b_om, List nb_points_int,
             double alpha_sigma, bool left_trunc, arma::mat Os_0_i, arma::mat Ws_0_i, bool competing_risk,
             double alpha_sigma_CR, arma::vec beta, arma::rowvec Xtime_i, arma::vec Utime_i, arma::mat b_al,
             arma::mat Xs_i, arma::mat Us_i, arma::mat Xs_0_i, arma::mat Us_0_i, double alpha_current,arma::vec sharedtype, arma::vec sharedtype_CR,
             double alpha_current_CR, arma::vec beta_slope, arma::rowvec Xslope_i, arma::mat b_al_slope, arma::vec Uslope_i, arma::mat Xs_slope_i, 
             arma::mat Us_slope_i,arma::mat Xs_slope_0_i, arma::mat Us_slope_0_i, double alpha_slope, double alpha_slope_CR, List HB, 
             arma::vec wk, double Time_i, arma::vec st_i, arma::vec st_0_i, double shape, arma::vec gamma, arma::vec B_i, arma::mat Bs_i, arma::mat Bs_0_i,
             arma::rowvec Z_i, arma::vec alpha, double P_i, double P_0_i, 
             double shape_CR, arma::vec gamma_CR, arma::vec B_i_CR, arma::mat Bs_i_CR, arma::mat Bs_0_i_CR,
             arma::rowvec Z_i_CR, arma::vec alpha_CR, arma::mat X_base_i, arma::mat O_base_i, arma::mat W_base_i, arma::mat U_i, arma::vec y_i, int n_row_X, double sigma_epsilon,List event){
  
  //do something
  const std::string& hazard_baseline = HB[0];
  const std::string& hazard_baseline_CR = HB[1];
  int S = nb_points_int[0];
  int nb_pointsGK= nb_points_int[1];
  int event1_i = event[0];
  int event2_i = event[1];
  bool dep_current = sharedtype[0];
  bool dep_slope = sharedtype[1];
  bool dep_variability = sharedtype[2];
  bool dep_current_CR = sharedtype_CR[0];
  bool dep_slope_CR = sharedtype_CR[1];
  bool dep_variability_CR = sharedtype_CR[2];
  arma::vec h(S,fill::ones);
  arma::vec h_CR(S,fill::ones);
  //double h = 1;
  double etaBaseline = 0;
  arma::mat survLong(S,nb_pointsGK, fill::zeros);
  arma::mat survLong_0(S,nb_pointsGK, fill::zeros);
  arma::mat Sigma_CV;
  arma::mat Sigma_current_GK;
  arma::mat Sigma_current_GK_0;
  arma::mat CV;
  arma::mat current_GK;
  arma::mat current_GK_0;
  arma::mat slope;
  arma::mat slope_GK;
  arma::mat slope_GK_0;
  double etaBaseline_CR = 0;
  double etaBaseline_0 =0;
  double etaBaseline_0_CR =0;
  arma::mat survLong_CR(S,nb_pointsGK, fill::zeros);
  arma::mat survLong_0_CR(S,nb_pointsGK, fill::zeros);
  if(dep_variability || dep_variability_CR){
    Sigma_CV = exp(arma::dot(omega, Otime_i)+b_om*Wtime_i);
    Sigma_current_GK = exp(arma::repmat(omega.t()*Os_i.t(),S,1) + b_om*Ws_i.t());
    }
  if(dep_variability){
    h= h%exp(alpha_sigma*Sigma_CV);
    survLong = survLong + alpha_sigma*Sigma_current_GK;
    if(left_trunc){
      Sigma_current_GK_0 = exp(arma::repmat(omega.t()*Os_0_i.t(),S,1) + b_om*Ws_0_i.t());
      survLong_0 = survLong_0 + alpha_sigma*Sigma_current_GK_0;
    }
  }
  if(competing_risk){
    if(dep_variability_CR){
      h_CR= h_CR%exp(alpha_sigma_CR*Sigma_CV);
      survLong_CR = survLong_CR + alpha_sigma_CR*Sigma_current_GK;
      if(left_trunc){
        survLong_0_CR = survLong_0_CR + alpha_sigma_CR*Sigma_current_GK_0;
      }
    }
  }
  if((dep_current) ||
     (competing_risk && dep_current_CR)){
    CV =arma::dot(beta, Xtime_i)+b_al*Utime_i;
    current_GK = arma::repmat(beta.t()*Xs_i.t(),S,1) + b_al*Us_i.t();
    if(left_trunc){
      current_GK_0 = arma::repmat(beta.t()*Xs_0_i.t(),S,1) + b_al*Us_0_i.t();
    }
    if(dep_current){
      h = h%exp(alpha_current*CV);
      survLong = survLong + alpha_current*current_GK;
      if(left_trunc){
        survLong_0 = survLong_0 + alpha_current*current_GK_0;
      }
    }
    if(competing_risk && dep_current_CR){
      h_CR = h_CR%exp(alpha_current_CR*CV);
      survLong_CR = survLong_CR + alpha_current_CR*current_GK;
      if(left_trunc){
        survLong_0_CR = survLong_0_CR + alpha_current_CR*current_GK_0;
      }
    }
  }
  if((dep_slope) ||
     (competing_risk && dep_slope_CR)){
    slope = arma::dot(beta_slope, Xslope_i)+b_al_slope*Uslope_i;
    slope_GK = arma::repmat(beta_slope.t()*Xs_slope_i.t(),S,1) + b_al_slope*Us_slope_i.t();
    if(left_trunc){
      slope_GK_0 = arma::repmat(beta_slope.t()*Xs_slope_0_i.t(),S,1) + b_al_slope*Us_slope_0_i.t();
    }
    if(dep_slope){
      //Rcout << "The value of v : \n" << 21 << "\n";
      h = h%exp(alpha_slope*slope);
     // Rcout << "The value of v : \n" << 21 << "\n";
      survLong = survLong + alpha_slope*slope_GK;
      if(left_trunc){
        survLong_0 = survLong_0 + alpha_slope*slope_GK_0;
      }
    }
    if(competing_risk && dep_slope_CR){
      h_CR = h_CR%exp(alpha_slope_CR*slope);
      survLong_CR = survLong_CR + alpha_slope_CR*slope_GK;
      if(left_trunc){
        survLong_0_CR = survLong_0_CR + alpha_slope_CR*slope_GK_0;
      }
    }
  }
 // Rcout << "The value of v : \n" << 2 << "\n";
  //////h0
  double h_0;
  arma::vec h_0_GK;
  arma::vec h_0_GK_0;
  double h_0_CR;
  arma::vec h_0_GK_CR;
  arma::vec h_0_GK_0_CR;
  
  if(hazard_baseline == "Exponential"){
    h_0 = 1;
    h_0_GK = wk;
    if(left_trunc){
      h_0_GK_0 = wk;
    }
  }
  if(hazard_baseline == "Weibull"){
    h_0 = shape*(pow(Time_i,(shape-1)));
    h_0_GK = shape*(pow(st_i,shape-1))%wk;
    if(left_trunc){
      h_0_GK_0 = shape*(pow(st_0_i,shape-1))%wk;
    }
  }
  
  if(hazard_baseline == "Splines"){
    h_0 = exp(arma::dot(gamma,B_i));
    h_0_GK = wk%exp(Bs_i*gamma);
    if(left_trunc){
      h_0_GK_0 = wk%exp(Bs_0_i*gamma);
    }
  }
  double predsurv;
  if(Z_i.is_empty()){
    predsurv = 0;
  }
  else{
    predsurv = arma::dot(alpha, Z_i);
  }
  h = h_0*exp(predsurv)*h;
  etaBaseline = etaBaseline + predsurv;
  if(left_trunc){
    etaBaseline_0 = etaBaseline_0 + predsurv;
  }
  //###GK integration
  survLong = exp(survLong);
  survLong = survLong*h_0_GK;
  arma::vec Surv;
  Surv = (-exp(etaBaseline)*P_i*survLong);
  arma::vec Surv_0;
  arma::vec Surv_CR;
  arma::vec Surv_0_CR;
  if(left_trunc){
    survLong_0 = exp(survLong_0);
    survLong_0 = survLong_0*h_0_GK_0;
    Surv_0 = exp((-exp(etaBaseline_0)*P_0_i*survLong_0));
  }
  if(competing_risk){
    if(hazard_baseline_CR == "Exponential"){
      h_0_CR = 1;
      h_0_GK_CR = wk;
      if(left_trunc){
        h_0_GK_0_CR = wk;
      }
    }
    if(hazard_baseline_CR == "Weibull"){
      h_0_CR = shape_CR*(pow(Time_i,(shape_CR-1)));
      h_0_GK_CR = shape_CR*(pow(st_i,shape_CR-1))%wk;
      if(left_trunc){
        h_0_GK_0_CR = shape_CR*(pow(st_0_i,shape_CR-1))%wk;
      }
    }
    
    if(hazard_baseline_CR == "Splines"){
      h_0_CR = exp(arma::dot(gamma_CR,B_i_CR));
      h_0_GK_CR = wk%exp(Bs_i_CR*gamma_CR);
      if(left_trunc){
        h_0_GK_0_CR = wk%exp(Bs_0_i_CR*gamma_CR);
      }
    }
    double predsurv_CR;
    if(Z_i_CR.is_empty()){
      predsurv_CR = 0;
    }
    else{
      predsurv_CR = arma::dot(alpha_CR, Z_i_CR);
    }
    h_CR = h_0_CR*exp(predsurv_CR)*h_CR;
    etaBaseline_CR = etaBaseline_CR + predsurv_CR;
    if(left_trunc){
      etaBaseline_0_CR = etaBaseline_0_CR + predsurv_CR;
    }
    
    //###GK integration
    survLong_CR = exp(survLong_CR);
    survLong_CR = survLong_CR*h_0_GK_CR;
    Surv_CR = (-exp(etaBaseline_CR)*P_i*survLong_CR);
    if(left_trunc){
      survLong_0_CR = exp(survLong_0_CR);
      survLong_0_CR = survLong_0_CR*h_0_GK_0_CR;
      Surv_0_CR = exp((-exp(etaBaseline_0_CR)*P_0_i*survLong_0_CR));
    }
  }
  // Longitudinal part
  arma::vec f_Y_b_sigma(S,fill::zeros);
  arma::vec sigma_long;
  if(n_row_X == 0){
    if(variability_hetero){
      sigma_long = exp(dot(omega,O_base_i) + b_om*W_base_i );
    }
    else{
      sigma_long = sigma_epsilon;
    }
    CV  = dot(beta,X_base_i) + b_al*U_i ;
    f_Y_b_sigma = log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i-CV)/sigma_long,2);
  }
  else{
    for(int k=0; k<n_row_X; k++){
      if(variability_hetero){
        
        sigma_long = exp(dot(omega,O_base_i.row(k)) + b_om*W_base_i.row(k).t());
        
      }
      else{
        sigma_long = sigma_epsilon;
      }
      CV = dot(beta,X_base_i.row(k)) + b_al*U_i.row(k).t(); 
      f_Y_b_sigma = f_Y_b_sigma + log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i(k)-CV)/sigma_long,2);
     // Rcout << "The value of v : \n" << f_Y_b_sigma << "\n";
    }
  }
  arma::vec log_dens_int;
  double Clogexp;
  double log_dens;
  if(competing_risk){
    log_dens_int = f_Y_b_sigma + log(pow(h,event1_i))+log(pow(h_CR,event2_i))+Surv+Surv_CR;
    //Rcout << "The value of v : \n" << 2 << "\n";
    Clogexp = max(log_dens_int) - 500;
   // Rcout << "The value of v : \n" << 3 << "\n";
    log_dens_int = log_dens_int - Clogexp;
    //Rcout << "The value of v : \n" << 4 << "\n";
    log_dens = Clogexp + log(sum(exp(log_dens_int))) - log(S);
    //Rcout << "The value of v : \n" << 5 << "\n";
  }
  else{
    log_dens_int = f_Y_b_sigma + log(pow(h,event1_i))+Surv;
    Clogexp = max(log_dens_int) - 500;
    log_dens_int = log_dens_int - Clogexp;
    log_dens = Clogexp  +log(sum(exp(log_dens_int))) - log(S);
    
  }
  double den;
 //Rcout << "The value of v : \n" << 1 << "\n";
  if(left_trunc){
    if(competing_risk){
      den = log(arma::dot(Surv_0,Surv_0_CR))-log(S);
      //Rcout << "The value of v : \n" << den << "\n";
      //Rcout << "The value of v : \n" << log(sum(exp(log(Surv_0)+log(Surv_0_CR))))-log(S) << "\n";
    }
    else{
      den = log(sum(Surv_0))-log(S);
    }
    log_dens = log_dens - den;
  }
  
  return log_dens; // Remplacez cette valeur par le résultat que vous souhaitez retourner
}