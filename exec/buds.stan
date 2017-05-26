/*
  STAN implementation of Bayesian Unidimensional Scaling for
  modeling pseudotemporal ordering of samples.  
  Author: Lan Huong
*/
  
data {
  int<lower=1> N;                                     // number of pairwise dissimilarities
  int<lower=1> Npairs;                                // number of pairwise dissimilarities
  int<lower=1> i_idx[Npairs];                  
  int<lower=1> j_idx[Npairs];                     
  row_vector<lower=0>[Npairs] dvec;                   // vector of original pairwise distances
  row_vector<lower=0>[Npairs] rel_sd;                 // vector of weights for the distances

  // Hyperparameters
  real<lower = 0> gamma_tau;                          // Cauchy scale param for shift term
  real<lower = 0> gamma_bias;                         // Cauchy scale param for shift term
  real<lower = 0> gamma_rho;                          // Cauchy scale param for scale term
  real<lower = 0> gamma_epsilon;                      // Cauchy scale param for error term
  real<lower = 0> min_sigma;                          // Cauchy scale param for error term
}

parameters {
  real<lower = 0> tau_shape1;                         // shape1 hyperparameter for tau                                         
  real<lower = 0> tau_shape2;                         // shape2 hyperparameter for tau                               
  real<lower = 0> bias;                               // global shift for 1D latent distances                                         
  real<lower=0> rho_sq;                               // global scale for 1D latent distance      
  real<lower = min_sigma> meansd;                     // mean standard deviation parameter                                
  real<lower = 0, upper = 1> tau[N];                  // latent 1D coordinate
}

model {
  real sig_sq;
  real alpha;
  real beta;
  row_vector[Npairs] delta;                           // latent distances for 1D-coordinates
  
  tau_shape1 ~ cauchy(1, gamma_tau);    
  tau_shape2 ~ cauchy(1, gamma_tau);      
  bias ~ cauchy(0, gamma_bias);                           
  rho_sq ~ cauchy(1, gamma_rho);                      
  meansd ~ cauchy(0, gamma_epsilon);  
  
  for(n in 1:N) {
    tau[n] ~ beta(tau_shape1, tau_shape2);
  }
  
  for(kdx in 1:Npairs){
    delta[kdx] = bias + rho_sq * fabs(tau[i_idx[kdx]] - tau[j_idx[kdx]]);
    sig_sq = pow(fmax(meansd * rel_sd[kdx], min_sigma), 2);
    alpha = pow(delta[kdx], 2)/sig_sq;
    beta = delta[kdx]/sig_sq;
    dvec[kdx] ~ gamma(alpha, beta);
  }
}
