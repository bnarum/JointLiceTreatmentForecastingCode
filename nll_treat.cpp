#include <TMB.hpp>
using namespace density;
using namespace Eigen;

template<class Type>
Type atanh(Type x){
  return (log(Type(1) + x) - log(Type(1) - x)) / Type(2);
}
VECTORIZE1_t(atanh)

template<class Type>
Type logit_d(Type x, Type d){
  return log((x + d) / (Type(1) - x + d));
}

template<class Type>
Type asinh(Type x){
  return log(x + sqrt(pow(x, 2) + 1));
}
VECTORIZE1_t(asinh)

template<class Type>
Type yeo_johnson(Type x, Type pi){
  return (pow(Type(1) + x, pi) - Type(1)) / pi;
}
VECTORIZE1_t(yeo_johnson)

template<class Type>
Type box_cox(Type x, Type pi){
  return (pow(x, pi) - Type(1)) / pi;
}
VECTORIZE1_t(box_cox)
  
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  
  DATA_MATRIX(v); // binary treatment data
  DATA_MATRIX(y_AF);                                     // lice abundance
  DATA_MATRIX(y_MB);                                     // lice abundance
  DATA_MATRIX(y_ST);                                     // lice abundance
  DATA_MATRIX(A);                         // active status
  DATA_MATRIX(low_limit);                         // active status
  DATA_MATRIX(lice_limit);
  DATA_INTEGER(p);
  DATA_INTEGER(q);
  
  // Parameters
  PARAMETER(beta_0); REPORT(beta_0); ADREPORT(beta_0);                                   // intercept
  PARAMETER(beta_low); REPORT(beta_low); ADREPORT(beta_low); // High lice limit effect              
  
  PARAMETER(phi_rate); REPORT(phi_rate); ADREPORT(phi_rate);
  PARAMETER(phi_scale); REPORT(phi_scale); ADREPORT(phi_scale);
  vector<Type> phi(p); for (int l = 1; l <= p; l++) phi(l - 1) = phi_scale * exp(phi_rate * (Type(l) - Type(1))); REPORT(phi); ADREPORT(phi);
  
  PARAMETER(theta_rate); REPORT(theta_rate); ADREPORT(theta_rate);
  PARAMETER(theta_scale); REPORT(theta_scale); ADREPORT(theta_scale);
  vector<Type> theta(q); for (int l = 1; l <= q; l++) theta(l - 1) = theta_scale * exp(theta_rate * (Type(l) - Type(1))); REPORT(theta); ADREPORT(theta);
  
  PARAMETER_VECTOR(gamma_AF); REPORT(gamma_AF); ADREPORT(gamma_AF); // PARAMETER_VECTOR(log_gamma_AF); vector<Type> gamma_AF = exp(log_gamma_AF);                            // autoregressive parameter for lice
  PARAMETER_VECTOR(gamma_MB); REPORT(gamma_MB); ADREPORT(gamma_MB); // PARAMETER_VECTOR(log_gamma_MB); vector<Type> gamma_MB = exp(log_gamma_MB);                            // autoregressive parameter for lice
  PARAMETER_VECTOR(gamma_ST); REPORT(gamma_ST); ADREPORT(gamma_ST); // PARAMETER_VECTOR(log_gamma_ST); vector<Type> gamma_ST = exp(log_gamma_ST);                            // autoregressive parameter for lice
  PARAMETER(gamma_I_AF); REPORT(gamma_I_AF); ADREPORT(gamma_I_AF); // Interaction on lice effect with treatment
  PARAMETER(gamma_I_MB); REPORT(gamma_I_MB); ADREPORT(gamma_I_MB); // Interaction on lice effect with treatment
  PARAMETER(gamma_I_ST); REPORT(gamma_I_ST); ADREPORT(gamma_I_ST); // Interaction on lice effect with treatment
  PARAMETER(pi_AF); REPORT(pi_AF); ADREPORT(pi_AF);
  PARAMETER(pi_MB); REPORT(pi_MB); ADREPORT(pi_MB);
  PARAMETER(pi_ST); REPORT(pi_ST); ADREPORT(pi_ST);
  // PARAMETER(log_d); Type d = exp(log_d); REPORT(d); ADREPORT(d);
  
  // Misc. quantities
  int num_locs = v.rows();                     // Number of locations
  int num_times = v.cols();                    // Timeseries length
  int r_AF = gamma_AF.size();
  int r_MB = gamma_MB.size();
  int r_ST = gamma_ST.size();
  int r_max = std::max({r_AF, r_MB, r_ST});
  
  matrix<Type> prob(num_locs, num_times);                 // Storing treatment probabilities
  matrix<Type> eta_innovation(num_locs, num_times);       // Storing innovations
  matrix<Type> eta_surprise(num_locs, num_times);         // Storing innovations
  matrix<Type> eta(num_locs, num_times);                  // Storing log-odds predictor
  matrix<Type> eta_bars(num_locs, num_times);                  // Storing log-odds predictor
  int max_lag = std::max({p, q, r_AF, r_MB, r_ST});
  
  // Let the initial values of p_0,..p_max(p,q,d,m) equal V_0,... and the others init to zero
  for (int i = 0; i < num_locs; ++i) {
    for (int t = 0; t < max_lag; ++t) {
      prob(i,t) = v(i,t);
      eta_innovation(i,t) = Type(0);
      eta_surprise(i,t) = Type(0);
    }
  }
  
  // Typset negative log-likelihood 
  Type nll = Type(0.0);

  // Contribution of treatments
  for (int i = 0; i < num_locs; ++i) {
    for (int t = r_max; t < num_times; ++t) {
      // Anticipation
      Type eta_bar = beta_0 + beta_low * low_limit(i, t);
      for (int l = 1; l <= r_max; l++){ // Effect of lice
        if (t - l < 0)  break;
        if (isStructuralZero(A(i, t - l))) break;
        // Each of the stages with different orders
        if (l <= r_AF) {
          eta_bar += gamma_AF(l - 1) * pow(y_AF(i, t - l), pi_AF);
        }
        if (l <= r_MB) {
          eta_bar += gamma_MB(l - 1) * pow(y_MB(i, t - l), pi_MB);
        }
        if (l <= r_ST) {
          eta_bar += gamma_ST(l - 1) * pow(y_ST(i, t - l), pi_ST);
        }
        if (l == 1 && !isStructuralZero(v(i, t - l))){ // Interaction with treatment
          eta_bar += gamma_I_AF * pow(y_AF(i, t - l), pi_AF);
          eta_bar += gamma_I_MB * pow(y_MB(i, t - l), pi_MB);
          eta_bar += gamma_I_ST * pow(y_ST(i, t - l), pi_ST);
        }
      }
      
      // Auto-covariance in the surprise. NB: skip those where surprise and innovation are not complete
      eta(i,t) = eta_bar;
      for (int l = 1; l <= p; l++){ // Auto-regressive
        if (t - l < 0)  break;
        if (isStructuralZero(A(i, t - l))) break;
        eta(i, t) += phi(l - 1) * eta_surprise(i, t - l);
      }
      for (int l = 1; l <= q; l++){ // Moving-average
        if (t - l < 0)  break;
        if (isStructuralZero(A(i, t - l))) break;
        eta(i, t) += theta(l - 1) * eta_innovation(i, t - l);
      }
      
      // Store results
      eta_bars(i, t) = eta_bar;
      if (!isStructuralZero(A(i, t))){// Don't use inactive
        prob(i,t) = invlogit(eta(i,t));
        eta_surprise(i, t) = atanh(v(i, t) - invlogit(eta_bar));
        eta_innovation(i, t) = atanh(v(i,t) - prob(i,t));
        // eta_surprise(i, t) = logit_d(v(i, t), d) - eta_bar;
        // eta_innovation(i, t) = logit_d(v(i,t), d) - eta(i,t);
      } else {
        prob(i,t) = Type(0.0);
        eta_innovation(i,t) = Type(0.0);
        eta_surprise(i,t) = Type(0.0);
      }
      
      // Likelihood above cut-off time for comparison
      if (!isStructuralZero(A(i, t)) && (t >= max_lag) && (t < num_times - 1)){ // Ignore warm-up and exclude last data point for out-of-sample evaluation
        nll -= dbinom_robust(
          v(i,t),
          Type(1),
          eta(i,t), // logit(prob) // prob(i,t), // 
          true
        );
      }
    }
  }
  
  REPORT(prob);
  REPORT(eta_innovation);
  REPORT(eta_surprise);
  REPORT(eta);
  REPORT(eta_bars);
  
  return nll;
}
