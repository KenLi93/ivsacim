#' This function is to compute the centralized instrument given the covariates L
#' @keywords internal
#' @importFrom stats model.matrix
#' @export
IV_center <- function (Z, L = NULL) {
  
  # Now we obtain the iid-decomp of \hat\theta
  expit <- function (x) {
    exp(x)/(1 + exp(x))
  }
  
  IV_type <- all(Z == 0 | Z == 1)

  # fitting model for expectation of instrument

  if (IV_type) {
    if (is.null(L)) {
      zmod <- glm(Z ~ 1, family = "binomial")
      E_dot <- matrix(rep(1, length(Z)), ncol = 1) * c(fitted(zmod) * (1 - fitted(zmod))) 
    }
    else {
      Z_data <- data.frame(Z, L)
      zmod <- glm(Z ~ ., data = Z_data, family = "binomial")
      E_dot <- model.matrix(zmod) * c(fitted(zmod) * (1 - fitted(zmod))) 
    }
    
  }
  else{
    if (is.null(L)) {
      zmod <- lm(Z ~ 1)
      E_dot <- matrix(rep(1, length(Z)), ncol = 1)
    }
    else {
      Z_data <- data.frame(Z, L)
      zmod <- lm(Z ~ ., data = Z_data)
      E_dot <- model.matrix(zmod)
    }
    
  }
  
  eps.theta <- as.matrix(iid(zmod))
  p.dim <- dim(eps.theta)[2]
  Z.c <- Z - fitted(zmod)
  
  iv_center <- list(Zc = Z.c,
                    epstheta = eps.theta,
                    Edot = E_dot,
                    pdim = p.dim)
  return(iv_center)
}

#' This function is to generate a grid of treatment status
#' @keywords internal
#' @export
treatment_status <- function (n,
                              K,
                              stime, 
                              treatment_init, 
                              treatment_shift_time, 
                              max.time) {

  D_status <- matrix(0, ncol = K, nrow = n)
  for (i in 1:n) {
    if (treatment_shift_time[i] <= 0 || treatment_shift_time[i] >= max.time) {
      D_status[i, ] = treatment_init[i]
      next
    }
    D_status[i, stime < treatment_shift_time[i]] = treatment_init[i]
    D_status[i, stime >= treatment_shift_time[i]] = 1 - treatment_init[i]
  }
  return(D_status)
}



#' This function is to compute the centralized treatment process given the covariates Z, L at each sorted unique event time
#' @keywords internal
#' @importFrom stats model.matrix
#' @export
trt_center <- function (D_status, Z, L = NULL) {
  N = nrow(D_status)
  k = ncol(D_status)
  # Now we obtain the iid-decomp of \hat\theta
  #expit <- function (x) {
  #  exp(x)/(1 + exp(x))
  #}
  
  #IV_type <- all(Z == 0 | Z == 1)
  
  # fitting model for expectation of instrument
  
  #if (IV_type) {
    if (is.null(L)) {
      #E_dot <- matrix(rep(1, 100), ncol = 1) * c(fitted(zmod) * (1 - fitted(zmod))) 
      E_dot = 0
      eps_theta = 0
      Z_status <- matrix(rep(Z, k), ncol = k)
      pred_D_status <- apply(D_status * Z, 2, mean) / mean(Z) * Z_status + apply(D_status * (1 - Z), 2, mean) / mean(1 - Z) * (1 - Z_status)
    }
  #else {
    #  Z_data <- data.frame(Z, L)
    #  zmod <- glm(Z ~ ., data = Z_data, family = "binomial")
    #  E_dot <- model.matrix(zmod) * c(fitted(zmod) * (1 - fitted(zmod))) 
    #}
    
  # }
  # else{
  #   if (is.null(L)) {
  #     zmod <- lm(Z ~ 1)
  #     E_dot <- matrix(rep(1, length(Z)), ncol = 1)
  #   }
  #   else {
  #     Z_data <- data.frame(Z, L)
  #     zmod <- lm(Z ~ ., data = Z_data)
  #     E_dot <- model.matrix(zmod)
  #   }
  #   
  # }
  
  #eps_theta <- as.matrix(iid(zmod))
  #pdim <- dim(eps.theta)[2]
  D_status_c <- D_status - pred_D_status
  
  trt_center <- list(D_status_c = D_status_c,
                    epstheta = eps_theta,
                    Edot = E_dot,
                    pdim = k)
  return(trt_center)
}