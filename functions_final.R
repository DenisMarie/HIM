require(Matrix)
require(huge)
require(igraph)
# require(mSSL)
require(coda)
require(glmnet)
require(bayesreg)
require(mvtnorm)
require(ggplot2)
require(glmnetUtils)
require(BoomSpikeSlab)
require(truncnorm)
require(extraDistr)
require(igraph)
require(FastGGM)


#' lodget
#'
#' @param A - Matrix
#' @param op - Apply cholesky decomposition by specifying "chol" here.
#'
#' @return v
#'
#' @examples
logdet <- function(A, op) {
  require(Matrix)
  if(!missing(op)) {
    if(op != 'chol') stop("The second argument can only be a string 'chol' if it is specified.")
    v <- 2 * sum(log(diag(chol(A))))
  } else {
    # CHANGED (SJ) : Removed 3 calls to Matrix::expand, and created mat.exp dummy
    mat.exp <- Matrix::expand(lu(Y))
    l <- mat.exp$L
    u <- mat.exp$U
    p <- mat.exp$P
    du <- t(t(diag(u)))
    c <- det(p) %*% prod(sign(du))
    v <- log(c) + sum(log(abs(du)))
  }

  return(v)
}


#' multruncn
#'
#' @param mu
#' @param S
#' @param Y
#' @param Z
#'
#' @return Y
#'
#' @examples
#'
multruncn <- function(mu, S, Y, Z) {
  require(stats)
  n <- nrow(Y)
  # TODO: Parallelize for loop
  for (i in 1:n) {
    Y2 <- Y
    Y2 <- Y2[-i, , drop=F]
    mu1 <- mu[i]
    mu2 <- mu
    mu2 <- mu2[-i, , drop=F]
    S11 <- S[i, i]
    S12 <- S[i, , drop=F]
    S12 <- S12[, -i, drop=F]
    S21 <- t(S12)
    S22 <- S
    S22 <- S22[-i, , drop=F]
    S22 <- S22[, -i, drop=F]
    S22inv <- solve(S22)
    S112 <- S11 - S12 %*% S22inv %*% S21
    s <- sqrt(S112)
    m1 <- mu1 + S12 %*% S22inv %*% (Y2 - mu2)
    p <- stats::pnorm(0, m1, s)
    u <- stats::runif(1)
    z <- Z[i]

    if (z == 0) {
      Y1 <- stats::qnorm(u * p, m1, s)
    } else if (z == 1) {
      Y1 <- stats::qnorm(u * (1 - p) + p, m1, s)
    }

    Y[i] <- Y1
  }

  return(Y)
}


#' log_r_gamma_given_G
#'
#' @param gamma
#' @param gamma_prop
#' @param adj
#' @param a
#' @param b
#'
#' @return
#'
#' @examples
log_r_gamma_given_G <- function(gamma, gamma_prop, adj, a, b) {
  # Compute MH ratio for adding or removing one var
  p <- nrow(gamma)

  # +1 if adding, -1 if removing a var
  gamma_diff <- sum(gamma_prop - gamma)

  # Assumption in paper is that adjacency matrix has 0's along the diagonal,
  # while here is has 1's, so need to subtract diag(p)
  adj <- adj - diag(p)

  log_MH <- gamma_diff * a + b * (t(gamma_prop) %*% adj %*% gamma_prop - t(gamma) %*% adj %*% gamma)

  return(log_MH)
}


log_r_y_probit <- function(gamma, gamma_prop, X, Y, Z, h_0, h_alpha, h_beta) {
  # Compute MH ratio p(Y|gamma_prop) / p(Y|gamma) on log scale
  n <- dim(X)[1]
  p <- dim(X)[2]
  X_gamma <- X[, which(gamma != 0)]
  X_gamma_prop <- X[, which(gamma_prop != 0)]

  # Use logdet rather than log(det()) is case det is large/small
  # Similarly, log1p computes log(1 + p) which is accurate for small p
  core_term <- diag(n) +
    h_0 * (matrix(1, n, 1) %*% matrix(1, 1, n)) +
    h_alpha * (Z %*% t(Z)) +
    h_beta * (X_gamma %*% t(X_gamma))

  core_term_prop <- diag(n) +
    h_0 * (matrix(1, n, 1) %*% matrix(1, 1, n)) +
    h_alpha * (Z %*% t(Z)) +
    h_beta * (X_gamma_prop %*% t(X_gamma_prop))

  log_y_mh_ratio <- 0.5 * logdet(core_term, 'chol') +
    (1 / 2 * t(Y) %*% solve(core_term) %*% Y) -
    0.5 * logdet(core_term_prop, 'chol') -
    (1 / 2 * t(Y) %*% solve(core_term_prop) %*% Y)

  return(log_y_mh_ratio)
}


MCMC_non_linear_model_with_graph_selection_Marie <- function(X, O, Z, h_alpha, h_beta, a, b, gamma, burnin, nmc, adj, summary_only, verbose = F) {
  require(pbapply)

  # Standardize data so that S(i,i) = n
  S <- t(X) %*% X
  n <- dim(X)[1]
  p <- dim(X)[2]
  S <- stats::cor(S) * n

  # Initial guess for Sigma, precision matrix, and adjacency matrix

  # Always keep variable selections
  gamma_save <- matrix(0, p, nmc)

  # Record some diagnostic info
  full_gamma_save <- matrix(0, p, burnin + nmc)

  # Keep track of info to compute acceptance rates
  n_gamma_prop <- 0
  n_gamma_accept <- 0
  n_add_prop <- 0
  n_add_accept <- 0
  n_remove_prop <- 0
  n_remove_accept <- 0

  # Allocate storage for MCMC sample, or just for means if only summary is required
  if (summary_only) {
    Omega_save <- matrix(0, p, p)
    adj_save <- Omega_save
    Y_save <- matrix(0, n, 1)
  } else {
    Omega_save <- array(0, dim=c(p, p, nmc))
    adj_save <- Omega_save
    Y_save <- array(0, dim=c(n, 1, nmc))
  }

  # Number of currently included variables
  p_gamma <- apply(gamma, 2, sum)

  # Simulate latent variables
  X_gamma <- X[, which(gamma != 0)]
  V <- diag(n) +
    0 * (matrix(1, n, 1) %*% matrix(1, 1, n)) +
    h_alpha * (Z %*% t(Z)) +
    h_beta * (X_gamma %*% t(X_gamma))
  m <- matrix(0, n, 1)
  Ytmp <- matrix(1, n, 1)
  Y <- multruncn(m, V, Ytmp, O)

  # MCMC sampling
  # CHANGED (SJ): Added progress bar
  if(verbose)
    pb = pbapply::startpb(min = 1, max = (burnin + nmc))
  for (iter in 1:(burnin + nmc)) {
    if(verbose)
      pbapply::setpb(pb, iter)

    # Print out info every 100 iterations
    if(verbose){
      if (iter %% 100 == 0) {
        cat('Y =', Y[1], '\n')
        cat('Y =', Y[2], '\n')
        cat('Y =', Y[3], '\n')
        cat('Iteration =', iter, '\n')
        cat('Number of included variables =', sum(gamma), '\n')
        cat('Number of add variable moves proposed', n_add_prop, 'and accepted', n_add_accept, '\n')
        cat('Number of remove variable moves proposed', n_remove_prop, 'and accepted', n_remove_accept, '\n')
      }
    }

    # Select an entry at random to toggle
    change_index <- sample(1:p, 1)

    # Keep track of number of add vs. remove moves proposed
    if (gamma[change_index] == 0) {
      n_add_prop <- n_add_prop + 1
    } else {
      n_remove_prop <- n_remove_prop + 1
    }
    n_gamma_prop <- n_gamma_prop + 1

    # Toggle value
    gamma_prop <- gamma
    gamma_prop[change_index] <- abs(gamma[change_index] - 1)

    # Compute MH ratio on log scale
    log_r <- log_r_y_probit(gamma, gamma_prop, X, Y, Z, 0, h_alpha, h_beta) +
      log_r_gamma_given_G(gamma, gamma_prop, adj, a, b)

    # Accept proposal with probability r
    if (log(runif(1)) < log_r) {
      if (gamma[change_index] == 0) {
        gamma[change_index] = 1
        p_gamma <- p_gamma + 1
        n_add_accept <- n_add_accept + 1
      } else {
        gamma[change_index] <- 0
        p_gamma <- p_gamma - 1
        n_remove_accept <- n_remove_accept + 1
      }
      n_gamma_accept <- n_gamma_accept + 1
    }

    # resample latent variable
    X_gamma <- X[, which(gamma != 0)]
    V <- diag(n) +
      0 * (matrix(1, n, 1) %*% matrix(1, 1, n)) +
      h_alpha * (Z %*% t(Z)) +
      h_beta * (X_gamma %*% t(X_gamma))
    m <- matrix(0, n, 1)
    Y <- multruncn(m, V, Y, O)

    if (iter > burnin) {
      gamma_save[, iter - burnin] <- gamma

      if (summary_only) {
        adj_save <- adj_save + adj / nmc
        Y_save <- Y_save + Y / nmc
      } else {
        adj_save[, , iter - burnin] <- adj
        Y_save[, , iter - burnin] <- Y
      }
    }
  }

  ar_gamma <- n_gamma_accept / n_gamma_prop

  # Info for diagnostic purposes
  info <- list()
  info['n_add_prop'] <- n_add_prop
  info['n_add_accept'] <- n_add_accept
  info['n_remove_prop'] <- n_remove_prop
  info['n_remove_accept'] <- n_remove_accept
  info['full_gamma'] <- full_gamma_save

  # Add results to return list
  result <- list()
  result[["gamma_save"]] <- gamma_save
  result[["ar_gamma"]] <- ar_gamma
  result[["adj_save"]] <- adj_save
  result[["info"]] <- info
  result[["Y_save"]] <- Y_save

  return(result)
}


Extract_SPVAS <- function(X = nor_miRNA, Y, XY.spvas, YY.spvas ){

  
}
# SS with two moves and one component to change----
SS_prior_bis.new <- function(X, yc, s, eta = eta, niter = niter, c= 10, se.tau = 5,
                             a.tau = 2, b.tau = 1.5, burn = 1000, c.resid = 100){
  require(extraDistr)
  require(stats)

  s <- as.numeric(s)

  p <- ncol(X)
  n <- nrow(X)

  se2_prev <- se2_store <- 1
  beta_prev <- beta_store <- numeric(p)
  g_prev <- g_store <- logical(p)
  omega_prev <- omega_store <- numeric(p)
  tau_prev <- tau_store <- 0
  #cpt <- 0
  #cptM <- acc <- 0
  # Initialization - for i = 1
  g_store[sample(1:p,10, replace =FALSE)] <- TRUE
  g_prev <- g_store
  omega_store[] <- exp(eta)/(1+exp(eta))
  omega_prev <- omega_store
  tau_store <- 0.5
  tau_prev <- tau_store

  for( i in 2:niter){
    r <- 1
    move <- sample(1:2, 1)
    #cptM[i] <- move
    if (sum(g_prev) == 0){
      g_star <- g_prev
      g_star[sample(1:p, 1)] <- TRUE
    }else{
      if (sum(g_prev) == p) {
        g_star <- g_prev
        g_star[sample(1:p, 1)] <- FALSE
      }else{
        # if ( i== 500) browser()
        if (move == 1 & sum(g_prev) > 1){
          g_star <- g_prev
          g_num <- g_star
          g_den <- g_star
          sub <- sample(1:p, r)
          g_num[sub] <- !g_num[sub]

          XtX_num <- crossprod(X[, g_num])
          XtX_den <- crossprod(X[, g_den])

          A_den <- solve(XtX_den + diag(sum(g_den))/c)
          A_num <- solve(XtX_num  + diag(sum(g_num))/c)

          det_A_num <- 0.5*(determinant(A_num, logarithm = T)$modulus-
                              determinant(diag(sum(g_num))*c, logarithm = T)$modulus)
          det_A_den <- 0.5*(determinant(A_den, logarithm = T)$modulus-
                              determinant(diag(sum(g_den))*c, logarithm = T)$modulus)

          a_num <- A_num %*% crossprod(X[,g_num], yc)
          a_den <- A_den %*% crossprod(X[,g_den], yc)

          l_num <- ((n-1)/2)*log(0.5*(crossprod(yc) - crossprod(a_num, (XtX_num  + diag(sum(g_num))/c))%*%a_num))
          l_den <- ((n-1)/2)*log(0.5*(crossprod(yc) - crossprod(a_den,(XtX_den + diag(sum(g_den))/c) )%*%a_den))

          logLikelihood <- det_A_num-l_num-(det_A_den-l_den)

          logPrior <- sum(stats::dbinom(g_num,size = 1, prob = omega_prev, log =TRUE ))-
            sum(stats::dbinom(g_den,size = 1,prob = omega_prev, log =TRUE ))

          R <- exp(logLikelihood+logPrior)
          tt <- stats::rbinom(1, size =1 , prob = min(R,1))
          if (tt==1) {
            g_star <- g_num
            #acc[i] <- 1
          }
        }
        if (move == 2){
          g_star <- g_prev
          g_num <- g_star
          g_den <- g_star
          sub1 <- sample(which(g_star == TRUE), r, replace =F)
          sub2 <- sample(which(g_star == FALSE), r, replace =F)
          g_num[sub1] <- FALSE
          g_num[sub2] <- TRUE
          XtX_num <- crossprod(X[, g_num ])
          XtX_den <- crossprod(X[, g_den])

          A_num <- solve(XtX_num  + diag(sum(g_num))/c)
          A_den <- solve(XtX_den + diag(sum(g_den))/c)

          det_A_num <- 0.5*(determinant(A_num, logarithm = T)$modulus-
                              determinant(diag(sum(g_num))*c, logarithm = T)$modulus)
          det_A_den <- 0.5*(determinant(A_den, logarithm = T)$modulus-
                              determinant(diag(sum(g_den))*c, logarithm = T)$modulus)

          a_num <- A_num %*% crossprod(X[,g_num], yc)
          a_den <- A_den %*% crossprod(X[,g_den], yc)

          l_num <- ((n-1)/2)*log(0.5*(crossprod(yc) - crossprod(a_num, (XtX_num  + diag(sum(g_num))/c)) %*% a_num))
          l_den <- ((n-1)/2)*log(0.5*(crossprod(yc) - crossprod(a_den,(XtX_den + diag(sum(g_den))/c) ) %*% a_den))

          logLikelihood <- det_A_num-l_num-(det_A_den-l_den)

          logPrior <- sum(stats::dbinom(1*g_num,size = 1,prob = omega_prev, log =TRUE ))-
            sum(stats::dbinom(1*g_den,size = 1,prob = omega_prev, log =TRUE ))

          R <- exp(logLikelihood+logPrior)
          tt <- stats::rbinom(1, size =1 , prob = min(R,1))
          if (tt == 1){
            g_star <- g_num
            #acc[i] = 1
          }
        }
      }}

    # Update g
    g_prev <- g_star
    if(i > burn){
      g_store <- as.numeric(g_star) + as.numeric(g_store)

    }else{
      g_store <- as.numeric(g_star)
    }
    #g_store[,i] <- g_star # save new value g_star for sample i

    # update tau_store with a Metropolis step
    tau_star <- tau_prev

    #tau_num <- tau_star
    tau_den <- tau_star
    tau_num <- truncnorm::rtruncnorm(1, a = 0, mean = tau_star, sd = sqrt(se.tau))

    p_num <- exp(eta+tau_num*s) / (1+ exp(eta+tau_num*s))
    p_den <- exp(eta+tau_den*s) / (1+ exp(eta+tau_den*s))

    num <- sum(stats::dbinom(g_star, size = 1, prob = p_num, log = TRUE))
    den <- sum(stats::dbinom(g_star, size = 1, prob = p_den, log = TRUE))

    ratio <- num + stats::dgamma(tau_num,a.tau,b.tau, log = TRUE) +
      extraDistr::dtnorm(tau_den, a = 0, mean = tau_num, sd = sqrt(se.tau), log= TRUE)
    ratio <- ratio - (den + stats::dgamma(tau_den, a.tau,b.tau, log = TRUE) +
                        extraDistr::dtnorm(tau_num, a = 0, mean = tau_num, sd = sqrt(se.tau), log = TRUE))

    u <- runif(1)
    if (u < exp(ratio)){
      #cpt <- cpt+1
      tau_prev <- tau_num
      if(i > burn)
        tau_store <- tau_store + tau_num
      else
        tau_store <- tau_num
      #tau_store[i] <- tau_num # store new tau
    }else{
      tau_prev <- tau_den
      if(i > burn)
        tau_store <- tau_store + tau_den
      else
        tau_store <- tau_den
    }

    # Update omega
    omega_prev[] <- exp(eta+tau_prev*s)/(1+exp(eta+tau_prev*s)) # store new omega
    if(i > burn){
      omega_store[] <- omega_store[] + omega_prev[]
    }else{
      omega_store[] <- omega_prev[]
    }

    #omega_store[,i] <- exp(eta+tau_store[i]*s)/(1+exp(eta+tau_store[i]*s)) # store new omega
    #print(sum(g_store))

    # Calculate residuals
    # if(i > burn){
    #   A <- solve(crossprod(X[,g_prev]) + diag(sum(g_prev))/c)
    #   a <- A %*% crossprod(X[,g_prev], yc)
    #   beta_prev[g_prev] <- mvtnorm::rmvnorm(1,a, se2_prev*A)
    #   beta_store[g_prev] <- beta_store[g_prev] + beta_prev[g_prev]
    #   se2_prev <- extraDistr::rinvgamma(1,(n-1)/2, 0.5*crossprod(yc-as.matrix(X[,g_prev])%*%beta_prev[g_prev]) )
    #   se2_store <- se2_store + se2_prev
    # }
  }

  return(list(g = g_store/(niter-burn+1),
              tau = tau_store/(niter-burn+1),
              omega = omega_store/(niter-burn+1)
              # ,beta_store = beta_store/(niter-burn+1),
              # se2_store = se2_store/(niter-burn+1),
              # epsilon = yc-X%*%(beta_store/(niter-burn+1))
              ))
}

HIM <- function(Y, X1, X2, S, XY.spvas = NULLL,
                eta=-2.5, a.tau=6, b.tau=1, sub.mechanistic =1, 
                niter=5000, burn=1000, burnin=100, nmc=1000, thres=0.2, ggm = F,
                n_threads = parallel::detectCores()-1, quiet = F) {
  require(pbapply)
  require(huge)
  require(Matrix)
  require(igraph)

  if(quiet)
    pbapply::pboptions(type = "none")

  # Then to run the model
  n <- nrow(X2) # the number of observations
  p <- ncol(X2) # the number of miRNAs
  
  if (sub.mechanistic  == 1){
  res <- list()

  # PARALLEL:
  if(!quiet)
    message("Initializing matrix and computing priors.")
  res <- pbapply::pblapply(1:ncol(X1), function(j){
    yc <- c(X1[,j]-mean(X1[,j])) # the outcome (one mRNA)
    n <- length(yc) # the number of observations
    ss <- as.numeric(S[j,]) # the scores based on biological knowledge
    res_j <- SS_prior_bis.new(X = as.matrix(X2), yc = c(yc), s = c(ss), eta = eta , niter =niter, c=100, se.tau = 10,
                              a.tau = a.tau, b.tau = b.tau, burn = burn)
    return(res_j)
  }, cl = 1)

  # to extract the residuals
  # E.ss <- sapply(res, function(c) c$epsilon ) # matrix of residuals
  # colnames(E.ss) <- colnames(X1)

  # to extract the regression coefficients
  # B.ss <-  sapply(res, function(c) c$beta_store) # matrix of regression coefficients
  # colnames(B.ss) <- colnames(X1)
  # rownames(B.ss) <- colnames(X2)

  # to compute the posterior inclusion probabilities
  E.ss <- B.ss <- matrix(0, nrow= nrow(X1), ncol = ncol(X1))
  g.ss <- sapply(res, function(c) (c$g)) #matrix
  colnames(g.ss) <- colnames(X1)
  rownames(g.ss) <- colnames(X2)
  g.ss[g.ss < thres] <- 0 # if the probability is lower than 0.2 the miRNA is not selected
  g.ss[g.ss !=0 ] <- 1
  id.miRNA.ss <- rowSums(g.ss) == 0 # miRNA that are related to none mRNAs
  id.mRNA.ss <- colSums(g.ss) == 0 # mRNA that are related to none miRNAs
  for (i in 1:ncol(X1)){
    if (length(which(g.ss[,i] ==1) ) > 1){
      cvfit <- cv.glmnet(x = X2[,which(g.ss[,i] ==1)], y = X1[,i], alpha = 0 )
      B.ss[,i] <- X2[,which(g.ss[,i] ==1)]%*%coef(cvfit, s = "lambda.min")[-1,1]
      E.ss[, i] <- X1[,i]-B.ss[,i]
    }else{
      if (length(which(g.ss[,i] ==1)) == 1 ) {
        cvfit <- lm( X1[,i]  ~ X2[,which(g.ss[,i] ==1), drop =FALSE])
        B.ss[,i] <- X2[,which(g.ss[,i] ==1), drop =FALSE]%*%cvfit$coefficients[2]
        E.ss[, i] <- X1[,i]-B.ss[,i]
      }else{
        B.ss[,i] <- 0
        E.ss[, i] <- X1[,i]
      }
    }
  }
 
 
    X <- scale(cbind(B.ss, E.ss, X2[,id.miRNA.ss]))
  colnames(X) <- c(paste0(colnames(X1), "_miRNA"), paste0(colnames(X1)),
                   colnames(X2)[id.miRNA.ss])
  X <- X[,!id.mRNA.ss]
}else{
  
  id.miRNA.mRNA <- sort(match(unique(XY.spvas[,1]), colnames(X2)))
  id.miRNA <- (1:ncol(X2))[-id.miRNA.mRNA]
  X <- scale(cbind(X1, X2[,id.miRNA]))
  colnames(X) <- c(colnames(X1),
                   colnames(X2)[id.miRNA])
  B.ss <- E.ss <- NULL
  
}

## CLinical submodel----- 
  if(ggm){
    res.ggm <- FastGGM(E.ss, lambda =10)
    G.mRNA.adj <- res.ggm$partialCor
    G.pvalue <- res.ggm$p_partialCor
    G.mRNA.adj[G.pvalue<0.05] <- 1
    G.mRNA.adj[G.mRNA.adj !=1] <- 0
    # res.huge.ss <- huge::huge(E.ss, method = "glasso")
    # sel.ss <- huge::huge.select(res.huge.ss, criterion = "stars", rep.num = 30)
    # G.huge.mRNA.adj <- sel.ss$refit
    
    colnames(G.mRNA.adj) <- rownames(G.mRNA.adj) <- colnames(X1)
    #myg.mRNA.adj <- igraph::graph_from_adjacency_matrix(as.matrix(G.huge.mRNA.adj), mode="undirected")

  adj <- Matrix::bdiag(diag(ncol(X1)), G.mRNA.adj , diag(length(id.miRNA.ss[id.miRNA.ss])))
  adj <- adj[!id.mRNA.ss,!id.mRNA.ss]
  diag(adj) <- 1

  gamma <- matrix(0, ncol(X), 1)
  gamma_prop <- matrix(0, ncol(X), 1)
  gamma_prop[8, 1] <- 1
  h0 <- 0
  Z <- matrix(0, nrow(X), 5)
  h_0 <- 0
  h_alpha <- 1
  h_beta <- 6
  a <- -3
  b <- 0.5

  result <- MCMC_non_linear_model_with_graph_selection_Marie(X = as.matrix(X),
                                                             O = as.matrix(Y),
                                                             Z = Z,
                                                             h_alpha = h_alpha,
                                                             h_beta = h_beta,
                                                             a = a, b = b,
                                                             gamma = gamma,
                                                             burnin = burnin,
                                                             nmc = nmc, adj = as.matrix(adj), T)
  gg.ising <- t(result[[1]])

  colnames(gg.ising) = colnames(X)
  
  }else{
    G.huge.mRNA.adj <- Matrix::Matrix(diag(ncol(X1))) # check dim here - ensure the same as if ggm = T
    colnames(G.huge.mRNA.adj) <- rownames(G.huge.mRNA.adj) <- colnames(X1)
    cvfit<- cv.glmnet(as.matrix(X) ,as.factor(Y),  family="binomial")
    gg.ising <- colnames(X)[which(coef(cvfit, s="lambda.min")[-1]!=0)]
  }
  return(list(res = res, G = G.huge.mRNA.adj, gg = gg.ising, 
              X=X, E = E.ss, B= B.ss))

}
