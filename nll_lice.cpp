#include <TMB.hpp>
#include <numeric>
using namespace tmbutils;

// Lower bound approximation of gamma(1 + x) (https://en.wikipedia.org/wiki/Stirling%27s_approximation). First number is sqrt(2 * pi).
template<class Type>
inline Type gamma1p_lb(Type x)
{
  return Type(2.5066282746310002) * exp(log(x) / Type(2) - x + x * log(x) + Type(1) / (Type(12) * x + Type(1)));
}

// =======================
// Negative Log-Likelihood
// =======================

template<class Type>
Type objective_function<Type>::operator() ()
{
    // ====
    // Data
    // ====
    DATA_MATRIX(y_AF); // Lice adult female
    DATA_MATRIX(y_MB); // Lice mobile
    DATA_MATRIX(y_ST); // Lice stationary
    DATA_MATRIX(y_RN); // Recruits Neighbour
    DATA_MATRIX(y_RS); // Recruits Self
    DATA_MATRIX(A); // Active status
    DATA_MATRIX(temperature); // Temperature
    DATA_MATRIX(dev_time_stage); // 
    DATA_MATRIX(dev_time_ml); // 
    DATA_MATRIX(dev_time_fm); // 
    DATA_MATRIX(dev_time_larval); //
    DATA_MATRIX(new_hatches);
    DATA_MATRIX(v_med); // treatment
    DATA_MATRIX(v_mec); // treatment
    DATA_MATRIX(v_clf); // cleaner fish
    DATA_MATRIX(fish_counted);
    DATA_INTEGER(q);
      
    // ==========
    // Parameters
    // ==========
    
    // Dispersion
    PARAMETER(log_nu_AF); Type nu_AF = exp(log_nu_AF); REPORT(nu_AF); ADREPORT(nu_AF); 
    PARAMETER(log_nu_MB); Type nu_MB = exp(log_nu_MB); REPORT(nu_MB); ADREPORT(nu_MB);
    PARAMETER(log_nu_ST); Type nu_ST = exp(log_nu_ST); REPORT(nu_ST); ADREPORT(nu_ST);
    
    // Source scaling
    PARAMETER(log_iota_I); Type iota_I = exp(log_iota_I); REPORT(iota_I); ADREPORT(iota_I);
    PARAMETER(log_iota_H); Type iota_H = exp(log_iota_H); REPORT(iota_H); ADREPORT(iota_H);
    PARAMETER(log_iota_N); Type iota_N = exp(log_iota_N); REPORT(iota_N); ADREPORT(iota_N);
    
    // Temperature
    PARAMETER(beta_T); REPORT(beta_T); ADREPORT(beta_T);
    
    // Development time variability
    PARAMETER(log_alpha); Type alpha = exp(log_alpha); REPORT(alpha); ADREPORT(alpha); // Shape parameter discrete kernel
    Type gam1p1dalpha = gamma1p_lb(Type(1) / alpha);
    
    // Lice mortality
    PARAMETER(log_neg_rho); Type rho = -exp(log_neg_rho); REPORT(rho); ADREPORT(rho);
    PARAMETER_VECTOR(log_neg_delta_Vmec); vector<Type> delta_Vmec = -exp(log_neg_delta_Vmec); REPORT(delta_Vmec); ADREPORT(delta_Vmec);
    PARAMETER_VECTOR(log_neg_delta_Vmed); vector<Type> delta_Vmed = -exp(log_neg_delta_Vmed); REPORT(delta_Vmed); ADREPORT(delta_Vmed);
    Type delta_Vmed_tot = sum(delta_Vmed); REPORT(delta_Vmed_tot); ADREPORT(delta_Vmed_tot);
    PARAMETER_VECTOR(log_neg_delta_Vclf); vector<Type> delta_Vclf = -exp(log_neg_delta_Vclf); REPORT(delta_Vclf); ADREPORT(delta_Vclf);
    Type delta_Vclf_tot = sum(delta_Vclf); REPORT(delta_Vclf_tot); ADREPORT(delta_Vclf_tot);
    
    // Transformation constant in log(c + x)
    PARAMETER(log_c); Type c = exp(log_c); REPORT(c); ADREPORT(c);
    
    PARAMETER(theta_scale); REPORT(theta_scale); ADREPORT(theta_scale);
    PARAMETER(theta_rate); REPORT(theta_rate); ADREPORT(theta_rate);
    vector<Type> theta(q); for (int k = 0; k < q; k++) theta(k) = theta_scale * exp(theta_rate * k); REPORT(theta);
    
    // PARAMETER(log_zeta_0); Type zeta_0 = exp(log_zeta_0); REPORT(zeta_0); ADREPORT(zeta_0);
    
    int num_locs = y_AF.rows();
    int num_times = y_AF.cols();
    int development_lag = 20;
    int m_mec = delta_Vmec.size();
    int m_med = delta_Vmed.size();
    int m_clf = delta_Vclf.size();
    int max_lag = std::max({q, development_lag});
    
    // ===
    // NLL
    // ===
    
    // Estimates
    matrix<Type> mus_bar_AF(num_locs, num_times); matrix<Type> mus_bar_MB(num_locs, num_times); matrix<Type> mus_bar_ST(num_locs, num_times);
    matrix<Type> mus_AF(num_locs, num_times); matrix<Type> mus_MB(num_locs, num_times); matrix<Type> mus_ST(num_locs, num_times);
    matrix<Type> surprise_AF(num_locs, num_times); matrix<Type> surprise_MB(num_locs, num_times); matrix<Type> surprise_ST(num_locs, num_times);
    matrix<Type> innovation_AF(num_locs, num_times); matrix<Type> innovation_MB(num_locs, num_times); matrix<Type> innovation_ST(num_locs, num_times);
    
    // Stored for in-sample reporting
    matrix<Type> lambdas_I_AF(num_locs, num_times); matrix<Type> lambdas_I_MB(num_locs, num_times); matrix<Type> lambdas_I_ST(num_locs, num_times);
    matrix<Type> lambdas_H_AF(num_locs, num_times); matrix<Type> lambdas_H_MB(num_locs, num_times); matrix<Type> lambdas_H_ST(num_locs, num_times);
    matrix<Type> lambdas_N_AF(num_locs, num_times); matrix<Type> lambdas_N_MB(num_locs, num_times); matrix<Type> lambdas_N_ST(num_locs, num_times);
    
    // Fill in data to init
    for (int i = 0; i < num_locs; ++i) {
      for (int t = 0; t < development_lag; ++t) {
        mus_AF(i,t) = y_AF(i,t); mus_bar_AF(i,t) = y_AF(i,t); surprise_AF(i,t) = Type(0.0); innovation_AF(i,t) = Type(0.0);
        mus_MB(i,t) = y_MB(i,t); mus_bar_MB(i,t) = y_MB(i,t); surprise_MB(i,t) = Type(0.0); innovation_MB(i,t) = Type(0.0);
        mus_ST(i,t) = y_ST(i,t); mus_bar_ST(i,t) = y_ST(i,t); surprise_ST(i,t) = Type(0.0); innovation_ST(i,t) = Type(0.0);
      }
    }
    
    parallel_accumulator<Type> nll(this);
    for (int i = 0; i < num_locs; ++i) {
      for (int t = development_lag; t < num_times; ++t) {
        if (isStructuralZero(A(i, t))){ // Skip if not active, but store values 
          mus_AF(i,t) = Type(0.0); mus_MB(i,t) = Type(0.0); mus_MB(i,t) = Type(0.0); mus_ST(i,t) = Type(0.0);
          mus_bar_AF(i,t) = Type(0.0); mus_bar_MB(i,t) = Type(0.0); mus_bar_ST(i,t) = Type(0.0);
          surprise_AF(i,t) = Type(0.0); surprise_MB(i,t) = Type(0.0); surprise_ST(i,t) = Type(0.0);
          innovation_AF(i,t) = Type(0.0); innovation_MB(i,t) = Type(0.0); innovation_ST(i,t) = Type(0.0);
          
          lambdas_N_AF(i, t) = Type(0.0); lambdas_N_MB(i, t) = Type(0.0); lambdas_N_ST(i, t) = Type(0.0);
          lambdas_H_AF(i, t) = Type(0.0); lambdas_H_MB(i, t) = Type(0.0); lambdas_H_ST(i, t) = Type(0.0);
          lambdas_I_AF(i, t) = Type(0.0); lambdas_I_MB(i, t) = Type(0.0); lambdas_I_ST(i, t) = Type(0.0);
          continue;
        }
        
        // ===============
        // Base prediction 
        // ===============
        
        // Cumulative rates
        Type P_ST_M; Type P_PA_M; Type P_AD_F; Type P_LS_F;
        Type P_ST_F; Type P_PA_F; Type P_AD_M; Type P_LS_M;
        Type r_ST; Type r_MB; Type r_AF;
        Type temp_scale; Type scl_M; Type scl_F;
        
        Type lambda_I_AF = Type(0.0); Type lambda_I_MB = Type(0.0); Type lambda_I_ST = Type(0.0);
        Type lambda_N_AF = Type(0.0); Type lambda_N_MB = Type(0.0); Type lambda_N_ST = Type(0.0);
        Type lambda_H_AF = Type(0.0); Type lambda_H_MB = Type(0.0); Type lambda_H_ST = Type(0.0);
        Type log_kappa = Type(0.0); Type kappa;
        for (int l = 1; l <= development_lag; l++){
          if (isStructuralZero(A(i, t - l))) break; // Fallowing kills all signals
          // =========
          // Mortality
          // =========
          log_kappa += rho;
          if (!isStructuralZero(v_mec(i, t - l))){ // Treatment mechanical
            for (int k = 0; k < std::min(m_mec, l + 1); k++){
              log_kappa += delta_Vmec(k);
            }
          }
          if (!isStructuralZero(v_med(i, t - l))){ // Treatment medical
            for (int k = 0; k < std::min(m_med, l + 1); k++){
              log_kappa += delta_Vmed(k);
            }
          }
          if (!isStructuralZero(v_clf(i, t - l))){ // Treatment cleaner fish
            for (int k = 0; k < std::min(m_clf, l + 1); k++){
              log_kappa += delta_Vclf(k);
            }
          }
          kappa = exp(log_kappa);
          
          // ================
          // Cumulative rates
          // ================
          scl_M = dev_time_ml(i, t - l) / gam1p1dalpha;
          scl_F = dev_time_fm(i, t - l) / gam1p1dalpha;
          // Rates
          P_ST_M = pweibull(Type(l), alpha, Type(1) * scl_M);
          P_ST_F = pweibull(Type(l), alpha, Type(1) * scl_F);
          P_PA_M = pweibull(Type(l), alpha, Type(3) * scl_M);
          P_PA_F = pweibull(Type(l), alpha, Type(3) * scl_F);
          P_AD_M = pweibull(Type(l), alpha, Type(5) * scl_M);
          P_AD_F = pweibull(Type(l), alpha, Type(5) * scl_F);
          P_LS_M = pweibull(Type(l), alpha, Type(10) * scl_M); // Last (dead of old age) stage
          P_LS_F = pweibull(Type(l), alpha, Type(10) * scl_F); // Last (dead of old age) stage
          
          r_AF = kappa * (P_AD_F - P_LS_F);
          r_MB = kappa * ((P_PA_F - P_AD_F) + (P_PA_M - P_AD_M) + (P_AD_M - P_LS_M));
          r_ST = kappa * ((P_ST_F - P_PA_F) + (P_ST_M - P_PA_M));
          
          // =============
          // Contributions
          // =============
          // Intercept (unexplained sources)
          temp_scale = exp(beta_T * (temperature(i, t - l) - Type(10)));
          lambda_I_AF += iota_I * r_AF * temp_scale;
          lambda_I_MB += iota_I * r_MB * temp_scale;
          lambda_I_ST += iota_I * r_ST * temp_scale;
          // Hatches
          lambda_H_AF += iota_H * r_AF * y_RS(i, t - l);
          lambda_H_MB += iota_H * r_MB * y_RS(i, t - l);
          lambda_H_ST += iota_H * r_ST * y_RS(i, t - l);
          // Neighbors (NB: should build hatch rate into y_RN instead)
          lambda_N_AF += iota_N * r_AF * y_RN(i, t - l);
          lambda_N_MB += iota_N * r_MB * y_RN(i, t - l);
          lambda_N_ST += iota_N * r_ST * y_RN(i, t - l);
        }
        
        Type log_kappa_0 = Type(0.0);
        if (m_mec > 0) log_kappa_0 += delta_Vmec(0) * v_mec(i, t);
        if (m_med > 0) log_kappa_0 += delta_Vmed(0) * v_med(i, t);
        if (m_clf > 0) log_kappa_0 += delta_Vclf(0) * v_clf(i, t);
        Type kappa_0 = exp(log_kappa_0);
        
        Type mu_bar_AF = kappa_0 * (lambda_I_AF + lambda_H_AF + lambda_N_AF);
        Type mu_bar_MB = kappa_0 * (lambda_I_MB + lambda_H_MB + lambda_N_MB);
        Type mu_bar_ST = kappa_0 * (lambda_I_ST + lambda_H_ST + lambda_N_ST);
        
        // Auto-regressive correction
        Type log_zeta_AF = Type(0); Type log_zeta_MB = Type(0); Type log_zeta_ST = Type(0);
        for (int l = 1; l <= q; l++){
          if (isStructuralZero(A(i, t - l))) break; // Fallowing kills memory
          // MA
          if (l <= q) {
            log_zeta_AF += theta(l - 1) * innovation_AF(i, t - l);
            log_zeta_MB += theta(l - 1) * innovation_MB(i, t - l);
            log_zeta_ST += theta(l - 1) * innovation_ST(i, t - l);
          }
        }
        Type mu_AF = mu_bar_AF * exp(log_zeta_AF);
        Type mu_MB = mu_bar_MB * exp(log_zeta_MB); 
        Type mu_ST = mu_bar_ST * exp(log_zeta_ST);
        
        // ============
        // Model errors
        // ============
        // Surprise
        surprise_AF(i,t) = log(c + y_AF(i, t)) - log(c + mu_bar_AF);
        surprise_MB(i,t) = log(c + y_MB(i, t)) - log(c + mu_bar_MB);
        surprise_ST(i,t) = log(c + y_ST(i, t)) - log(c + mu_bar_ST);
        // Innovation
        innovation_AF(i,t) = log(c + y_AF(i, t)) - log(c + mu_AF);
        innovation_MB(i,t) = log(c + y_MB(i, t)) - log(c + mu_MB);
        innovation_ST(i,t) = log(c + y_ST(i, t)) - log(c + mu_ST);
        
        // ============
        // Store values
        // ============
        // Estimates
        mus_AF(i,t) = mu_AF; mus_MB(i,t) = mu_MB; mus_ST(i,t) = mu_ST;
        mus_bar_AF(i,t) = mu_bar_AF; mus_bar_MB(i,t) = mu_bar_MB; mus_bar_ST(i,t) = mu_bar_ST;
        // Partial estimates
        lambdas_N_AF(i, t) = lambda_N_AF; lambdas_N_MB(i, t) = lambda_N_MB; lambdas_N_ST(i, t) = lambda_N_ST; 
        lambdas_H_AF(i, t) = lambda_H_AF; lambdas_H_MB(i, t) = lambda_H_MB; lambdas_H_ST(i, t) = lambda_H_ST;
        lambdas_I_AF(i, t) = lambda_I_AF; lambdas_I_MB(i, t) = lambda_I_MB; lambdas_I_ST(i, t) = lambda_I_ST;
        // ==========================
        // Likelihood of observations
        // ==========================
        if ((max_lag <= t) && (t < num_times - 1)){ // We have an initial warm-up. EXCLUDE the last time step, to be used for oos evaluation
          Type log_mu_AF = log(c + mu_AF);
          Type log_mu_MB = log(c + mu_MB);
          Type log_mu_ST = log(c + mu_ST);
          Type log_n = log(fish_counted(i, t));
          
          nll -= dnbinom_robust(
            fish_counted(i, t) * y_AF(i, t),
            log_n + log_mu_AF,
            log_n + Type(2) * log_mu_AF - log_nu_AF,
            true
          );
          nll -= dnbinom_robust(
            fish_counted(i, t) * y_MB(i,t),
            log_n + log_mu_MB,
            log_n + Type(2) * log_mu_MB - log_nu_MB,
            true
          );
          nll -= dnbinom_robust(
            fish_counted(i, t) * y_ST(i,t),
            log_n + log_mu_ST,
            log_n + Type(2) * log_mu_ST - log_nu_ST,
            true
          );
        }
      }
    }
   
    // State reports
    REPORT(mus_AF); REPORT(mus_MB); REPORT(mus_MB); REPORT(mus_ST);
    REPORT(mus_bar_AF); REPORT(mus_bar_MB); REPORT(mus_bar_ST);
    REPORT(surprise_AF); REPORT(surprise_MB); REPORT(surprise_ST);
    REPORT(innovation_AF); REPORT(innovation_MB); REPORT(innovation_ST);
    REPORT(lambdas_I_AF); REPORT(lambdas_I_MB); REPORT(lambdas_I_ST);
    REPORT(lambdas_H_AF); REPORT(lambdas_H_MB); REPORT(lambdas_H_ST);
    REPORT(lambdas_N_AF); REPORT(lambdas_N_MB); REPORT(lambdas_N_ST);
    return nll;
}
