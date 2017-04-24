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
  row_vector<lower=0>[Npairs] weight;                 // vector of weights for the distances
                                                      // the distances (scale function)
  
  // Hyperparameters
  real<lower = 0> tau_shape1;                         // shape1 hyperparameter for tau
  real<lower = 0> tau_shape2;                         // shape2 hyperparameter for tau
  real<lower = 0> gamma_bias;                         // Cauchy scale param for shift term
  real<lower = 0> gamma_rho_sq;                       // Cauchy scale param for scale term
  real<lower = 0> gamma_epsilon;                      // Cauchy scale param for error term
  real<lower = 0> min_sigma;                          // Cauchy scale param for error term
}

transformed data {
  row_vector [Npairs] rel_dvec2; 

  for(kdx in 1:Npairs){
    rel_dvec2[kdx] = pow(dvec[kdx], 2);
  }
}

parameters {
  real<lower = 0> bias;                               // global shift for 1D latent distances                                         
  real<lower=0> rho_sq;                               // global scale for 1D latent distance      
  real<lower = 0> meanvar;                            // global shift for 1D latent distances                                         
  real<lower = 0, upper = 1> tau[N];                  // latent 1D coordinate
  //real<lower = 1e-2, upper = gamma_epsilon> meanvar; // global mean variance for distances           
}


model {
  real sig_sq;
  real alpha;
  real beta;
  row_vector[Npairs] delta;                           // latent distances for 1D-coordinates
  
  for(k in 1:Npairs){
    delta[k] = bias + rho_sq * fabs(tau[i_idx[k]] - tau[j_idx[k]]);
  }
  
  bias ~ cauchy(0, gamma_bias);                           
  rho_sq ~ cauchy(0, gamma_rho_sq);                      
  meanvar ~ cauchy(0, gamma_epsilon);  
  
  for(n in 1:N) {
    tau[n] ~ beta(tau_shape1, tau_shape2);
  }
  
  
  for(kdx in 1:Npairs){
    sig_sq = (1e-3 + meanvar) * (1e-2 + rel_dvec2[kdx]) + min_sigma;
    alpha = pow(delta[kdx], 2)/sig_sq;
    beta = delta[kdx]/sig_sq;
    target += weight[kdx] * gamma_lpdf(dvec[kdx] | alpha, beta);
  }
}
