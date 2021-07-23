# include <Rcpp.h>
using namespace Rcpp;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' This is the main function to compute the estimator in C++
//' @keywords internal
//' @export
// [[Rcpp::export]]
SEXP invalidivsacim_est(NumericVector time,
                        NumericVector event,
                        NumericVector stime, 
                        NumericVector Z,
                        NumericVector Z_c, 
                        NumericMatrix D_status,
                        NumericMatrix D_status_c,
                        NumericMatrix eps_2,
                        NumericMatrix Z_c_dot,
                        NumericVector weights){
  
  // variables for estimation
  int n = time.size();
  int k = stime.size();
  NumericVector B_D(k), dB_D(k), b_numer_D(2 * k), b_denom_D(4 * k), B_Z(k), dB_Z(k), risk_cumsum(k), indik_v(k), M_det(k);
  NumericMatrix dN(n, k),  risk_t(n, k), B_intD(n, k + 1), B_intZ(n, k + 1);
  double tmp_denom_D, temp_eps_11, temp_eps_12;
  double beta = 0, tot_risk = 0;
  double del = 0.01;
  
  // variables for variance estimate
  int pdim = Z_c_dot.ncol();
  //int pZdim = Dt_c_dot.ncol();
  NumericVector B_D_se(k), B_Z_se(k);
  //double B_D_se_tmp = 0;
  NumericMatrix b_var_est(k, k);
  NumericMatrix eps_1(n, 2 * k), eps(n, 2 * k);
  NumericMatrix b_D_dot(pdim, 2 * k);
  NumericMatrix H_dot(2 * k, 2 * k);
  NumericMatrix H_numer_dot(2 * k, 2 * k), H_denom_dot(2 * k, k), H_dot_Z(n, 2 * k);
  NumericVector eps_beta(n);
  double tmp, beta_se = 0;
  
  for (int j = 0; j < k; j++) {
    
    for (int i = 0; i < n; i++) {
      // estimation
      dN(i, j) = (time[i] == stime[j]) ? 1:0;
      dN(i, j) *= event[i];
      risk_t(i, j) = (time[i] >= stime[j]) ? 1:0;
      risk_cumsum[j] += risk_t(i, j);
      b_numer_D[2 * j + 0] += weights[i] * Z_c[i] * exp(B_intD(i, j) + B_intZ(i, j)) * dN(i, j);
      b_numer_D[2 * j + 1] += weights[i] * Z_c[i] * D_status_c(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * dN(i, j);
      b_denom_D[4 * j + 0] += weights[i] * Z_c[i] * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * D_status(i, j);
      b_denom_D[4 * j + 1] += weights[i] * Z_c[i] * D_status_c(i, j) * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * D_status(i, j);
      b_denom_D[4 * j + 2] += weights[i] * Z_c[i] * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * Z[i];
      b_denom_D[4 * j + 3] += weights[i] * Z_c[i] * D_status_c(i, j) * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * Z[i];
    }
    
    
    tmp_denom_D = b_denom_D[4 * j + 0] * b_denom_D[4 * j + 3] - b_denom_D[4 * j + 1] * b_denom_D[4 * j + 2];
    M_det[j] = tmp_denom_D;
    tmp_denom_D = (tmp_denom_D < 0) ? -tmp_denom_D:tmp_denom_D;
    indik_v[j] = (tmp_denom_D < del) ? 0:1;
    if (indik_v[j]) {
      dB_D[j] = (b_numer_D[2 * j + 0] * b_denom_D[4 * j + 3] - b_numer_D[2 * j + 1] * b_denom_D[4 * j + 2]) / M_det[j];
      dB_Z[j] = (-b_numer_D[2 * j + 0] * b_denom_D[4 * j + 1] + b_numer_D[2 * j + 1] * b_denom_D[4 * j + 0]) / M_det[j];
    }
    
    
    // estimation
    for (int i = 0; i < n; i++) {
      B_intD(i, j + 1) += D_status(i, j) * dB_D[j];
      B_intD(i, j + 1) += B_intD(i, j);
      eps_1(i, 2 * j + 0) += Z_c[i] * exp(B_intD(i, j) + B_intZ(i, j)) * (dN(i, j) - risk_t(i, j) * (D_status(i, j) * dB_D[j] + Z[i] * dB_Z[j]));
      eps_1(i, 2 * j + 1) += Z_c[i] * D_status_c(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * (dN(i, j) - risk_t(i, j) * (D_status(i, j) * dB_D[j] + Z[i] * dB_Z[j]));
      B_intZ(i, j + 1) += Z[i] * B_Z[j];
    }
    
    
    // variance estimate
    if (indik_v[j]) {
      for (int i = 0; i < n; i++) {
        temp_eps_11 = eps_1(i, 2 * j + 0);
        temp_eps_12 = eps_1(i, 2 * j + 1);
        eps_1(i, 2 * j + 0) = (temp_eps_11 * b_denom_D[4 * j + 3] - temp_eps_12 * b_denom_D[4 * j + 2]) / M_det[j];
        eps_1(i, 2 * j + 1) = (-temp_eps_11 * b_denom_D[4 * j + 1] + temp_eps_12 * b_denom_D[4 * j + 0]) / M_det[j];
      }
    }
    else{
      for (int i = 0; i < n; i++) {
        eps_1(i, 2 * j + 0) = 0;
        eps_1(i, 2 * j + 1) = 0;
      }
    }
    
    
    // estimation part for beta
    beta += risk_cumsum[j] * dB_D[j];
    
    if (j == 0) {
      B_D[j] = 0 + dB_D[j];
      B_Z[j] = 0 + dB_Z[j];
      tot_risk += risk_cumsum[j] * (stime[j] - 0);
    }
    else{
      B_D[j] = B_D[j - 1] + dB_D[j];
      B_Z[j] = B_Z[j - 1] + dB_Z[j];
      tot_risk += risk_cumsum[j] * (stime[j] - stime[j - 1]);
    }
  }
  
  // time invariant intensity
  beta /= tot_risk;
  
  // variance estimate
  for (int j = 0; j < k; j++) {
    for (int l = 0; l < j; l++) {
      for (int i = 0; i < n; i++) {
        H_numer_dot(2 * l + 0, 2 * j + 0) += Z_c[i] * exp(B_intD(i, j) + B_intZ(i, j)) * dN(i, j) * D_status(i, l);
        H_numer_dot(2 * l + 1, 2 * j + 0) += Z_c[i] * exp(B_intD(i, j) + B_intZ(i, j)) * dN(i, j) * Z[i];
        H_numer_dot(2 * l + 0, 2 * j + 1) += Z_c[i] * D_status_c(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * dN(i, j) * D_status(i, l);
        H_numer_dot(2 * l + 1, 2 * j + 1) += Z_c[i] * D_status_c(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * dN(i, j) * Z[i];
        // H_denom_dot(2 * l + 0, j) += (Z_c[i] * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * D_status(i, j) * D_status(i, l)) * b_denom_D[4 * j + 3];
        // H_denom_dot(2 * l + 0, j) += (Z_c[i] * D_status_c(i, j) * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * Z[i] * D_status(i, l)) * b_denom_D[4 * j + 0];
        // H_denom_dot(2 * l + 0, j) -= (Z_c[i] * D_status_c(i, j) * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * D_status(i, j) * D_status(i, l)) * b_denom_D[4 * j + 2];
        // H_denom_dot(2 * l + 0, j) -= (Z_c[i] * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * Z[i] * D_status(i, l)) * b_denom_D[4 * j + 1];
        // H_denom_dot(2 * l + 1, j) += (Z_c[i] * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * D_status(i, j) * Z[i]) * b_denom_D[4 * j + 3];
        // H_denom_dot(2 * l + 1, j) += (Z_c[i] * D_status_c(i, j) * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * Z[i] * Z[i]) * b_denom_D[4 * j + 0];
        // H_denom_dot(2 * l + 1, j) -= (Z_c[i] * D_status_c(i, j) * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * D_status(i, j) * Z[i]) * b_denom_D[4 * j + 2];
        // H_denom_dot(2 * l + 1, j) -= (Z_c[i] * risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * Z[i] * Z[i]) * b_denom_D[4 * j + 1];
      }
      if (indik_v[j]) {
        H_dot(2 * l + 0, 2 * j + 0) = (H_numer_dot(2 * l + 0, 2 * j + 0) * b_denom_D[4 * j + 3] - H_numer_dot(2 * l + 0, 2 * j + 1) * b_denom_D[4 * j + 2]) / M_det[j];
        //H_dot(2 * l + 0, 2 * j + 0) -= dB_D[j] * H_denom_dot(2 * l + 0, j) / M_det[j];
        H_dot(2 * l + 1, 2 * j + 0) = (H_numer_dot(2 * l + 1, 2 * j + 0) * b_denom_D[4 * j + 3] - H_numer_dot(2 * l + 1, 2 * j + 1) * b_denom_D[4 * j + 2]) / M_det[j];
        //H_dot(2 * l + 1, 2 * j + 0) -= dB_D[j] * H_denom_dot(2 * l + 1, j) / M_det[j];
        H_dot(2 * l + 0, 2 * j + 1) = (-H_numer_dot(2 * l + 0, 2 * j + 0) * b_denom_D[4 * j + 1] + H_numer_dot(2 * l + 0, 2 * j + 1) * b_denom_D[4 * j + 0]) / M_det[j];
        //H_dot(2 * l + 0, 2 * j + 1) -= dB_D[j] * H_denom_dot(2 * l + 0, j) / M_det[j];
        H_dot(2 * l + 1, 2 * j + 1) = (-H_numer_dot(2 * l + 1, 2 * j + 0) * b_denom_D[4 * j + 1] + H_numer_dot(2 * l + 1, 2 * j + 1) * b_denom_D[4 * j + 0]) / M_det[j];
        //H_dot(2 * l + 1, 2 * j + 1) -= dB_D[j] * H_denom_dot(2 * l + 1, j) / M_det[j];
      }
    }
  }

  // 
  // 
  // for (int j = 0; j < k; j++) {
  //   for (int j1 = 0; j1 < pdim; j1++){
  //     for (int l = 0; l < j; l++) {
  //       b_D_dot(j1, j) += H_dot(2 * l + 0, 2 * j + 0) * b_D_dot(j1, l);
  //       b_D_dot(j1, j) += H_dot(2 * l + 0, 2 * j + 1) * b_D_dot(j1, l);
  //       b_D_dot(j1, j) += H_dot(2 * l + 1, 2 * j + 0) * b_D_dot(j1, l);
  //       b_D_dot(j1, j) += H_dot(2 * l + 1, 2 * j + 1) * b_D_dot(j1, l);
  //     }
  //     for (int i = 0; i < n; i++) {
  //       if (indik_v[j]) {
  // 
  //         H_dot_Z(i, 2 * j + 0) += (exp(B_intD(i, j) + B_intZ(i, j)) * dN(i, j) * b_denom_D[4 * j + 3] - D_status_c(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * dN(i, j) * b_denom_D[4 * j + 0]) / M_det[j];
  //         H_dot_Z(i, 2 * j + 0) -= risk_t(i, j) * exp(B_intD(i, j) + B_intZ(i, j)) * D_status(i, j) * b_numer_D[j] / M_det[j];
  //       }
  //       b_D_dot(j1, 2 * j + 0) -= H_dot_Z(i, 2 * j + 0) * Z_c_dot(i, j1);
  //       b_D_dot(j1, 2 * j + 1) -= H_dot_Z(i, 2 * j + 1) * Z_c_dot(i, j1);
  //       
  //     }
  //   }
  // }
  // 
  // 
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2 * k; j++) {
      for (int l = 0; l < j; l++) {
        tmp -= H_dot(l, j) * eps(i, l);
      }
      eps(i, j) += eps_1(i, j) - tmp;
      tmp = 0;
    }
  }
  // 
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < 2 * k; j++) {
  //     for (int j1 = 0; j1 < pdim; j1++) {
  //       eps(i, j) += eps_2(i, j1) * b_D_dot(j1, j);
  //     }
  //   }
  // }
  // 
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < k; j++) {
  //     eps_beta[i] += eps(i, 2 * j + 1) * risk_cumsum[j];
  //   }
  //   eps_beta[i] /= tot_risk;
  // }
  // 
  // 
  // for (int j = 0; j < k; j++) {
  //   for (int l = 0; l < k; l++) {
  //     for (int i = 0; i < n; i++) {
  //       b_var_est(j, l) += eps(i, 2 * j + 1) * eps(i, 2 * l + 1);
  //     }
  //   }
  // }
  for (int j = 0; j < k; j++) {
    for (int l = 0; l < k; l++) {
      for (int i = 0; i < n; i++) {
        b_var_est(j, l) += eps(i, 2 * j + 0) * eps(i, 2 * l + 0);
        b_var_est(j, l) += eps(i, 2 * j + 1) * eps(i, 2 * l + 0);
        b_var_est(j, l) += eps(i, 2 * j + 0) * eps(i, 2 * l + 1);
        b_var_est(j, l) += eps(i, 2 * j + 1) * eps(i, 2 * l + 1);
        
      }
    }
  }
  


  for (int j = 0; j < k; j++) {
    for (int l = 0; l < j; l++) {
      B_D_se[j] += b_var_est(l, j);
    }
    for (int l = 0; l < j; l++) {
      B_D_se[j] += b_var_est(j, l);
    }
    B_D_se[j] += b_var_est(j, j);
    if (j > 0) {
      B_D_se[j] += B_D_se[j - 1];
    }
  }

  for (int j = 0; j < k; j++) {
    B_D_se[j] = sqrt(B_D_se[j]);
  }
  
  for (int i = 0; i < n; i++) {
   beta_se += eps_beta[i] * eps_beta[i];
  }
  
  beta_se = sqrt(beta_se);
  
  List by_prod = List::create(
    Named("Z_c") = Z_c,
    Named("dN") = dN,
    Named("risk_t") = risk_t,
    Named("indik_v") = indik_v,
    Named("b_numer_D") = b_numer_D,
    Named("b_denom_D") = b_denom_D,
    Named("B_intD") = B_intD,
    Named("eps_1") = eps_1,
    Named("eps_2") = eps_2,
    Named("eps") = eps,
    Named("H_dot") = H_dot,
    Named("b_D_dot") = b_D_dot,
    Named("b_var_est") = b_var_est,
    Named("eps_beta") = eps_beta);
  
  return List::create(
    Named("stime") = stime,
    Named("dB_D") = dB_D,
    Named("B_D") = B_D,
    Named("dB_Z") = dB_Z,
    Named("B_Z") = B_Z,
    Named("beta") = beta,
    Named("B_D_se") = B_D_se,
    Named("B_Z_se") = B_Z_se,
    Named("beta_se") = beta_se,
    Named("by_prod") = by_prod);
}
