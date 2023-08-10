ATE <- function(Y, treat, X, theta = 1, ATT = FALSE,
                verbose = FALSE, max.iter = 100,
                tol = 1e-10, initial.values = NULL,
                backtrack.alpha = 0.3, backtrack.beta = 0.5) {

  # This is the main function from package ATE("asadharis/ATE").
  # This creates an ATE object. The object contains point
  # estimates and variance estimates. It also has
  # generic S3 methods such as print, plot and summary.
  # Args:
  #   Y: Response vector.
  #   X: Design/Covariate matrix.
  #   treat: Vector of treatment indicator. Must be of
  #          the form (0, 1, ...).
  #   ATT: Indicate if treatment effect on the treated is
  #        to be estimated.
  #   verbose: Indicate if extra statments should be printed
  #            to show progress of the function while it runs.
  #   max.iter: Maximum no. of iterations for BFGS algorithms.
  #   tol: Tolerance of the algorithm for stopping conditions.
  #   initial.values: Matrix of initial values for BFGS algorithm.
  #   backtrack.alpha, backtrack.beta: Parameters for backtrack line
  #                                    search algorithm.
  # Returns:
  #   An object of class 'ATE'. This contains point estimates and
  #   and variance covariance matrix of our estimates. Along with other
  #   information about our object such as data used for estimation.

  # J is the number treatment arms
  J <- length(unique(treat))

  # In some cases if we have a data.frame we need
  # to convert it into a numeric matrix
  if (class(X) == "data.frame") {
    X <- as.matrix(X)
  }
  # Take care of the case of a single covariate
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1, nrow = length(X))
    warning("Data matrix 'X' is a vector, will be treated as n x 1 matrix.")
  }
  # The totoal length of the parameter vector Recall:
  #Using the notation of the Package Vignette we have u_K(X) = cbind(1,X).
  K <- ncol(X) + 1

  # Some simple checks before running the main functions
  if (nrow(X) != length(Y)) {
    stop("Dimensions of covariates and response do not match.")
  }
  if (J == 1) {
    stop("There must be atleast two treatment arms")
  }
  if (!all(0:(J - 1) == sort(unique(treat)))) {
    stop("The treatment levels must be labelled 0,1,2,...")
  }
  if (!is.numeric(theta)) {
    stop("'theta' must be a real number.")
  }
  if (J > 2 & ATT) {
    stop("For ATT == TRUE, must have only 2 treatment arms.")
  }

  # If the inital values are not
  # provided we generate simple initial parameters.
  if (is.null(initial.values)) {
    if (ATT) {
      initial.values <- numeric(K)
    } else {
      initial.values <- matrix(0, ncol = K, nrow = J)
    }
  }

  # For ATT there is only one implementation of BFGS.
  if (ATT & !is.vector(initial.values)) {
    stop("For ATT == TRUE, only need one vector of initial values.")
  }
  # In all other cases this has to be a matrix.
  if (!ATT) {
    if (!is.matrix(initial.values)) {
      stop("Initial values must be a matrix")
    }
    if (any(dim(initial.values) != c(J, K))) {
      stop("Matrix of initial values must have dimensions J x K.")
    }
  }

  # Now we specify the category of the problem The simple case with binary treatment.
  gp <- "simple"

  # The case of average treatment effect on the treated.
  if (ATT) {
    gp <- "ATT"
  }
  # The case of Multiple treatment arms.
  if (J > 2) {
    gp <- "MT"
  }

  # Obtain the estimates for whatever the case may be.
  if (gp == "simple") {
    est <- GetPointEstSimple(initial.values[1, ],
                             initial.values[2, ],
                             X = X, Y = Y, treat = treat,
                             theta = theta, max.iter = max.iter,
                             tol = tol,
                             backtrack.alpha = backtrack.alpha,
                             backtrack.beta = backtrack.beta,
                             verbose = verbose)
    if (verbose) {
      cat("\nEstimating Variance")
    }
    cov.mat <- GetCovSimple(est, Y = Y, treat = treat,
                            u.mat = cbind(1, X),
                            theta = theta)

  } else if (gp == "ATT") {
    est <- GetPointEstATT(initial.values, X = X, Y = Y,
                          treat = treat,
                          theta = theta, max.iter = max.iter,
                          tol = tol,
                          backtrack.alpha = backtrack.alpha,
                          backtrack.beta = backtrack.beta,
                          verbose = verbose)
    if (verbose) {
      cat("\nEstimating Variance")
    }
    cov.mat <- GetCovATT(est, Y = Y, treat = treat,
                         u.mat = cbind(1, X),
                         theta = theta)
  } else if (gp == "MT") {
    est <- GetPointEstMultiple(initial.values, X = X, Y = Y,
                               treat = treat,
                               theta = theta, max.iter = max.iter,
                               tol = tol,
                               backtrack.alpha = backtrack.alpha,
                               backtrack.beta = backtrack.beta,
                               verbose = verbose)
    if (verbose) {
      cat("\nEstimating Variance")
    }
    cov.mat <- GetCovMultiple(est, Y = Y, treat = treat,
                              u.mat = cbind(1, X),
                              theta = theta)
  }

  # Begin building the 'ATE' object.
  res <- est
  res$vcov<- cov.mat
  res$X <- X
  res$Y <- Y
  res$treat <- treat
  res$theta <- theta
  res$gp <- gp
  res$J <- J
  res$K <- K
  res$call <- match.call()

  # Rename some of the elements of the list and organize the output object
  if (gp == "simple") {
    estimate <- c(res$tau.one, res$tau.zero, res$tau)
    names(estimate) <- c("E[Y(1)]", "E[Y(0)]", "ATE")
    res$estimate <- estimate
    res$tau.one <- NULL
    res$tau.zero <- NULL
    res$tau <- NULL
  } else if (gp == "ATT") {
    estimate <- c(res$tau.one, res$tau.zero, res$tau)
    names(estimate) <- c("E[Y(1)|T=1]", "E[Y(0)|T=1]", "ATT")
    res$estimate <- estimate
    res$tau.one <- NULL
    res$tau.zero <- NULL
    res$tau <- NULL
  } else {
    estimate <- res$tau.treatment.j
    names(estimate) <- paste("E[Y(", 0:(J - 1), ")]", sep = "")
    res$estimate <- estimate
    res$tau.treatment.j <- NULL
  }
  # Define the class ATE
  class(res) <- "ATE"
  return(res)
}




GetPointEstSimple <- function(initial.one, initial.two, ...) {
  # Function to get point estimates for the simple case of
  # binary treatment and NOT the case of treatement effect on
  # the treated.
  #
  # Args:
  #   initial.one: The initial vector for first BFGS algorithm (treatment arm).
  #   initial.two: The initial vector for second BFGS algorithm (placebo arm).
  #   ...: Other arguments to be passed to the function from parent.
  #
  # Returns:
  #   A list of estimated lambda values where lambda is the parameter
  #   vector which we optiize over. It also contains the weights
  #   obtained from the estimated lambda values used for covariate
  #   balancing estimates. Finally, the list also contains the
  #   point estimates and convergence indicator.

  #Obtain extra arguments
  args<- list(...)
  Y <- args$Y
  X <- args$X
  treat <- args$treat
  theta <- args$theta
  max.iter <- args$max.iter
  tol <- args$tol
  backtrack.alpha <- args$backtrack.alpha
  backtrack.beta <- args$backtrack.beta
  verbose <- args$verbose

  N <- length(Y)  # Sample size.
  u.mat <- cbind(1, X)  # Obtain u matrix of covariates.
  u.bar <- colMeans(u.mat)  # Get column means.

  # Perform some simple checks.
  if (length(initial.one) != length(u.mat[1, ])) {
    stop("Incorrect length of initial vector")
  }
  if (length(initial.two) != length(u.mat[1, ])) {
    stop("Incorrect length of initial vector")
  }
  if (verbose) {
    cat("Running BFGS algorithm for estimating Weights p: \n")
  }

  # Implement the BFGS algorithm to get the optimal lambda
  # and correspoding weights.
  # This corresponds to \lambda_p and \hat{p}_K in vignette.
  treatment.hat <- BFGSAlgorithm(initial.one, u.mat, u.bar, treat,
                                 theta, backtrack.alpha,
                                 backtrack.beta, max.iter,
                                 tol)

  if (verbose) {
    cat("\nRunning BFGS algorithm for estimating Weights q: \n")
  }

  # Implement the BFGS algorithm to get the optimal lambda_q
  # and correspoding weights \hat{q}_K.
  placebo.hat <- BFGSAlgorithm(initial.two, u.mat, u.bar, 1 - treat,
                               theta, backtrack.alpha,
                               backtrack.beta,
                               max.iter, tol)

  # Obtain estimates for tau1 = E[Y(1)], tau2 = E[Y(0)]
  # and tau = E[Y(1)] - E[Y(0)].
  tau.one  <- sum((treatment.hat$weights * Y)[treat == 1])
  tau.zero <- sum((placebo.hat$weights * Y)[treat == 0])
  tau <- tau.one - tau.zero

  # Obtain estimated weights and lambda's for
  # each treatment arm.
  weights.treat <- treatment.hat$weights
  weights.placebo <- placebo.hat$weights

  lambda.treat <- treatment.hat$results
  lambda.placebo <- placebo.hat$results

  # Throw a warning if algorithm did not converge.
  converge = TRUE
  if (!treatment.hat$converge | !placebo.hat$converge) {
    warning("BFGS Algorithm did not converge for atleast one objective function.")
    converge <- FALSE
  }

  # Return list of objects
  return(list(lambda.treat    = lambda.treat,
              lambda.placebo  = lambda.placebo,
              weights.treat   = weights.treat,
              weights.placebo = weights.placebo,
              tau.one = tau.one, tau.zero = tau.zero,
              tau = tau, converge = converge))
}


GetPointEstATT <- function(initial, ...) {
  # Function to get point estimates average treatment effect
  # on the treated. This case also has binary treatment.
  #
  # Args:
  #   initial: The initial vector for the BFGS algorithm.
  #   ...: Other arguments to be passed to the function.
  #
  # Returns:
  #   List of objects similar to the previous function. However, this
  #   function does not return lambda or weights for the treatment arm.
  #   This is because the we do not need a covariate balancing technique
  #   for the treatment arm in this case.

  #Obtain extra arguments
  args<- list(...)
  Y <- args$Y
  X <- args$X
  treat <- args$treat
  theta <- args$theta
  max.iter <- args$max.iter
  tol <- args$tol
  backtrack.alpha <- args$backtrack.alpha
  backtrack.beta <- args$backtrack.beta
  verbose <- args$verbose

  N <- length(Y)  # Sample size.
  u.mat <- cbind(1, X)  # Design matrix u.

  # Main difference here is the definition of u.bar.
  # u.bar is vector of column means ONLY for those
  # who recieved the treatment.
  u.bar <- colMeans(u.mat[treat == 1, ])

  if (length(initial) != length(u.mat[1, ])) {
    stop("Incorrect length of initial vector")
  }
  if (verbose) {
    cat("\nRunning BFGS algorithm Raphson for estimating Weights q: \n")
  }

  # Run the BFGS algorithm for obtaining the parameter and weights
  # used for covariate balancing.
  placebo.hat <- BFGSAlgorithm(initial, u.mat, u.bar,
                                1 - treat, theta,
                                backtrack.alpha, backtrack.beta,
                                max.iter, tol)

  # Note that treatment effect on treated is simple to estimate
  tau.one <- mean(Y[treat == 1])
  # The calibration estimator for E[Y(0)|T = 1]
  tau.zero <- sum((placebo.hat$weights * Y)[treat == 0])
  tau <- tau.one - tau.zero

  # Weights and lambda vectors
  weights.placebo <- placebo.hat$weights
  lambda.placebo <- placebo.hat$results

  # Warning message if the algorithm did not converge
  converge = TRUE
  if (!placebo.hat$converge) {
    warning("\nBFGS algorithm did not converge for the objective function.")
    converge <- FALSE
  }

  # Return list
  return(list(lambda.placebo = lambda.placebo,
              weights.placebo = weights.placebo,
              tau.one = tau.one, tau.zero = tau.zero,
              tau = tau, converge = converge))
}


GetPointEstMultiple <- function(initial.mat, ...) {
  # Function to get point estimates for treatment effect
  # when we have multiple treatment arms.
  #
  # Args:
  #   initial.mat: A matrix of initial values for the different
  #                BFGS algorithm. Each row is an inital vector for
  #                an algorithm. The total number of rows is J: number of
  #                different treatment arms.
  #   ...: Other arguments to be passed to the function.
  #
  # Returns:
  #   A List with the following objects
  #     lam.mat: A matrix of estimated lambda values. This has the same
  #              dimensions as initial.mat.
  #     weights.mat: A matrix estimated weights. The j-th row
  #                  corresponds to the weights for treatment j.
  #     tau.treatment.j: A vector of estimates of [EY(0), EY(1),..., EY(J-1)].
  #     converge: A boolean indicator of convergence status.


  #Obtain extra arguments
  args<- list(...)
  Y <- args$Y
  X <- args$X
  treat <- args$treat
  theta <- args$theta
  max.iter <- args$max.iter
  tol <- args$tol
  backtrack.alpha <- args$backtrack.alpha
  backtrack.beta <- args$backtrack.beta
  verbose <- args$verbose

  N <- length(Y)  # Sample size.
  J <- length(unique(treat))  # Number of treatment arms.
  u.mat <- cbind(1, X)  # Obtain design matrix u.
  u.bar <- colMeans(u.mat)  # Obtain column means of design matrix.

  # A simple verification of dimensions.
  if (ncol(initial.mat) != length(u.mat[1, ])) {
    stop("Incorrect length of initial vector")
  }
  if (verbose) {
    cat("\nRunning BFGS for estimating Weights: \n")
  }

  # Initialize the matrices and vector which will be returned
  # by this function.
  lam.mat <- matrix(0, ncol = ncol(initial.mat), nrow = J)
  weights.mat <- matrix(0, ncol = N, nrow = J)
  tau.treatment.j <- numeric(J)

  # Loop through the different treatment arms.
  for (j in 0:(J - 1)) {
    # Inidicator of treatment arm.
    temp.treat <- 1 * (treat == j)

    # Implement BFGS algorithm
    treatment.j.hat <- BFGSAlgorithm(initial.mat[j + 1, ],
                                     u.mat, u.bar, temp.treat,
                                     theta,
                                     backtrack.alpha,
                                     backtrack.beta, max.iter, tol)

    # Find estimate for E[Y(j)]
    tau.j.hat <- sum((treatment.j.hat$weights * Y)[temp.treat == 1])

    tau.treatment.j[j + 1] <- tau.j.hat
    lam.mat[j + 1, ] <- treatment.j.hat$results
    weights.mat[j + 1, ] <- treatment.j.hat$weights

    # Warning for non-convergence
    converge = TRUE
    if (!treatment.j.hat$converge) {
      warning(paste("BFGS algorithm did not converge for treatment arm",
                    j))
      converge <- FALSE
    }
  }
  # Return result.

  return(list(lam.mat = lam.mat,
              weights.mat = weights.mat,
              tau.treatment.j = tau.treatment.j,
              converge = converge))
}


GetEESimple <- function(u.mat.lambda.treat, u.mat.lambda.placebo,
                        tau.one, tau.zero, Y, treat, u.mat, theta) {
  # This function calculates the estimating equations for
  # the case of simple binary treatments.
  # Args:
  #   u.mat.lambda.treat: The vector u.mat * lambda.treat.
  #   u.mat.lambda.placebo: The vector u.mat * lambda.placebo.
  #   ...: Other parameters to be passed to the function.
  # Returns:
  #   A matrix with N rows, each row is the EE for the
  #   corresponding subject evaluated at the given values of
  #   lambda.treat and lambda.placebo.

  # gk1 is the EE for lambda_treat (lambda_p in vignette).
  temp1 <- treat * CRFamilyDerivative(u.mat.lambda.treat, theta)
  gk1 <- apply(u.mat, 2, "*", temp1) - u.mat

  # The EE for lambda_placebo (lambda_q in vignette).
  temp2 <- (1 - treat) * CRFamilyDerivative(u.mat.lambda.placebo, theta)
  gk2 <- apply(u.mat, 2, "*", temp2) - u.mat

  # The EE for tau.one, the estimate for E[Y(1)]
  gk3 <- (treat * Y) * CRFamilyDerivative(u.mat.lambda.treat, theta) - tau.one
  # The EE for tau.zero, the estimate for E[Y(0)]
  gk4 <- ((1 - treat) * Y) *
         CRFamilyDerivative(u.mat.lambda.placebo, theta) -
         tau.zero

  # Return the full set of estimatin equations
  cbind(gk1, gk2, gk3, gk4)
}

GetCovSimple <- function(fitted.point.est, ...) {
  # This function calculates the estimate of the
  # covariance matrix for the parameters taus.
  # I.e. cov of E[Y(1)] and E[Y(0)]
  # Args:
  #   fitted.point.est: The output list of the function
  #                     GetPointEstSimple.
  #   ...: Other parameters to be passed.
  # Returns:
  #   The large covariance matrix for the taus,
  #   I.e. cov matrix for tau.one and tau.zero.

  #Obtain extra argments
  args <- list(...)
  Y <- args$Y
  treat <- args$treat
  u.mat <- args$u.mat
  theta <- args$theta

  N <- length(Y)  # Sample size.
  n1 <- sum(treat)  # No. of subjects in treatment arm.
  K<- ncol(u.mat)  # No. of covariates.

  # Obtain the lambda and u * lambda objects
  lambda.treat <- (fitted.point.est$lambda.treat)
  lambda.placebo <- (fitted.point.est$lambda.placebo)
  u.mat.lambda.treat <- u.mat %*% lambda.treat
  u.mat.lambda.placebo <- u.mat %*% lambda.placebo

  # Our point estimates tau's
  tau.one <- fitted.point.est$tau.one
  tau.zero <- fitted.point.est$tau.zero

  # Obtain the estimating equations
  gk <- GetEESimple(u.mat.lambda.treat, u.mat.lambda.placebo,
                    tau.one, tau.zero, Y, treat, u.mat, theta)

  # Calculate the MEAT for the Huber-White sandwich estimate
  meat <- crossprod(gk)/N

  # To calculate the bread of the sanwich estimator
  # we follow the notation of the package Vignette.
  # First we create the 2K*2K matrix A
  temp1 <- treat * CRFamilySecondDerivative(u.mat.lambda.treat, theta)
  temp2 <- (1 - treat) * CRFamilySecondDerivative(u.mat.lambda.placebo, theta)

  # We define the matrix A
  tempA1 <- apply(u.mat, 2, "*", temp1)
  A1 <- crossprod(tempA1, u.mat) / N
  tempA2 <- apply(u.mat, 2, "*", temp2)
  A2 <- crossprod(tempA2, u.mat) / N
  A <- bdiag(A1, A2)

  # Now for the 2*2K matrix C
  C1 <- crossprod(u.mat, temp1 * Y)
  C2 <- crossprod(u.mat, temp2 * Y)
  C <- t(bdiag(C1, C2)) / N

  # Calculate the bread of the sandwich estimator
  # using block matrix inverse formula
  Ainv <- bdiag(solve(A1), solve(A2))
  tempMat <- Matrix(0, ncol = 2, nrow = 2 * K)
  bread <- cbind(rbind(Ainv, C %*% Ainv), rbind(tempMat, -diag(rep(1, 2))))

  # Obtain the large covariance matrix
  large.cov.mat <- (bread %*% tcrossprod(meat, bread)) / N

  # Return the submatrix of cov for tau1 and tau0.
  large.cov.mat[-(1:(2 * K)), -(1:(2 * K))]
}

################################################

GetEEATT <- function(u.mat.lambda.placebo, tau.one,
                     tau.zero, Y, treat, u.mat, theta) {
  # This function calculates the estimating equations for
  # the case of estimating the treatment effect on the treated.
  # Args:
  #   u.mat.lambda.placebo: The vector u.mat * lambda.placebo.
  #   ...: Other parameters to be passed to the function.
  # Returns:
  #   A matrix with N rows, each row is the EE for the
  #   corresponding subject evaluated at the given values of
  #   lambda.placebo.

  n1 <- sum(treat)  # No. in treatment arm.
  N <- length(Y)  # Sample size.
  delta <- n1/N  # Parameter delta.

  #The EE for lambda.placebo
  temp1 <- (1 - treat) * CRFamilyDerivative(u.mat.lambda.placebo, theta)

  gk1 <- apply(u.mat, 2, "*", temp1) - apply(u.mat, 2, "*", treat) / delta
  # THe EE for delta = E[T]
  gk2 <- treat - delta
  # The EE for expected value of tau.one = E[Y(1)|T=1]
  gk3 <- treat * (Y - tau.one) / delta
  # THe EE for expected value of tau.zero = E[Y(0)|T=1]
  gk4 <- temp1 * Y - tau.zero

  cbind(gk1, gk2, gk3, gk4)
}


GetCovATT <- function(fitted.point.est, ...) {
  # This function calculates the estimate of the
  # covariance matrix for the parameter vector.
  # Args:
  #   fitted.point.est: The output list of the function
  #                     GetPointEstATT.
  #   ...: Other parameters to be passed.
  # Returns:
  #   The covariance matrix for our parameters
  #   tau.one and tau.zero.

  #Obtain extra argments
  args <- list(...)
  Y <- args$Y
  treat <- args$treat
  u.mat <- args$u.mat
  theta <- args$theta

  # Initialize some variables.
  N <- length(Y)
  n1 <- sum(treat)
  # Recall: delta = E[T], probability of treatment assignment.
  delta <- n1/N
  K<- ncol(u.mat)

  lambda.placebo <- (fitted.point.est$lambda.placebo)
  u.mat.lambda.placebo <- u.mat %*% lambda.placebo

  tau.one <- fitted.point.est$tau.one
  tau.zero <- fitted.point.est$tau.zero

  # Obtain EE and meat for sandwich estimator
  gk <- GetEEATT(u.mat.lambda.placebo = u.mat.lambda.placebo,
                 tau.one = tau.one,
                 tau.zero = tau.zero, Y = Y,
                 treat = treat, u.mat = u.mat,
                 theta = theta)
  meat <- crossprod(gk) / N

  # Obatain the (K+1)*(K+1) matrix A
  temp1 <- (1 - treat) * CRFamilySecondDerivative(u.mat.lambda.placebo, theta)
  tempA <- apply(u.mat, 2, "*", temp1)
  A <- cbind(crossprod(tempA, u.mat), crossprod(u.mat, treat)/(delta ^ 2)) / N
  A <- rbind(A, c(rep(0, K), -1))

  # Evaluate matrix C of dimension 2*(K+1)
  C <- Matrix(0, ncol = K, nrow = 2)
  C[2, ] <- crossprod(u.mat, temp1 * Y)
  C <- cbind(C, c(0, 0)) / N

  # Evaluate the bread matrix and find covariance est
  A.inv <- solve(A)
  temp.mat <- Matrix(0, ncol = 2, nrow = K + 1)
  bread <- cbind(rbind(A.inv, C %*% A.inv), rbind(temp.mat, -diag(rep(1, 2))))

  # Obtain large covariance matrix.
  large.cov.mat <- (bread %*% tcrossprod(meat, bread)) / N
  # Return the submatrix for covariance of tau0 and tau1.
  large.cov.mat[-(1:(K+1)), -(1:(K+1))]
}

################################################

GetEEMultiple <- function(u.mat.lambdas, lambdas, taus,
                          Y, treat, u.mat, theta) {
  # This function calculates the estimating equations for
  # the case of multiple treatment arms.
  # Args:
  #   u.mat.lambdas: The matrix u.mat * lam.mat.
  #   ...: Other parameters to be passed to the function.
  # Returns:
  #   A matrix with N rows, each row is the EE for the
  #   corresponding subject evaluated at the given values of
  #   lambdas.

  # Get the number of treatment arms
  J <- nrow(lambdas)
  # gk1 is a list of EE for lambda_0, lambda_1, ...
  gk1.list <- vector("list", J)
  # gk2 is the list of EE for tau_0, tau_1, ...
  gk2.list <- vector("list", J)

  # For each treatment arm we obtain the EE for
  # lambda_j and tau_j = E[Y(j)]
  for (j in 1:J) {
    u.mat.lambda.j <- u.mat.lambdas[, j]
    temp1 <- 1 * (treat == j - 1) * CRFamilyDerivative(u.mat.lambda.j, theta)
    gk1.list[[j]] <- apply(u.mat, 2, "*", temp1) - u.mat
    gk2.list[[j]] <- 1 * (treat == j - 1) *
                    Y * CRFamilyDerivative(u.mat.lambda.j, theta) -
                    taus[j]
  }

  # Combine the EE for each treatment arm.
  cbind(do.call(cbind, gk1.list), do.call(cbind, gk2.list))
}


GetCovMultiple <- function(fitted.point.est, ...) {
  # This function calculates the estimate of the
  # covariance matrix for the parameter vector with
  # multiple treatments.
  # The function returns the J * J covariance matrix
  # for the parameters tau_0, tau_1, ..., tau_{J-1}
  # Args:
  #   fitted.point.est: The output list of the function
  #                     GetPointEstMultiple.
  #   ...: Other parameters to be passed.
  # Returns:
  #   The covariance matrix for our parameters
  #   tau0, tau1, ..., tau_{J-1}.

  #Obtain extra argments
  args <- list(...)
  Y <- args$Y
  treat <- args$treat
  u.mat <- args$u.mat
  theta <- args$theta

  # Intialize some variables.
  N <- length(Y)
  K <- ncol(fitted.point.est$lam.mat)
  J <- length(unique(treat))

  # Obtain the estimated lambdas and taus.
  lambdas <- fitted.point.est$lam.mat
  taus <- fitted.point.est$tau.treatment.j
  # Obtain the u.mat * lambdas matrix.
  u.mat.lambdas <- tcrossprod(u.mat, lambdas)

  # The meat matrix as usual.
  gk <- GetEEMultiple(u.mat.lambdas = u.mat.lambdas, lambdas = lambdas,
                      taus = taus, Y = Y, treat = treat,
                      u.mat = u.mat,
                      theta = theta)
  meat <- crossprod(gk) / N

  # Alist is the list of length J.
  # Each element is the K * K matrix, a subset of
  # the large A JK * JK matrix.
  A.list <- vector("list", J)
  # This list stores the inverse of each element of A.list.
  A.inv.list <- vector("list", J)

  # And also the C matrix of size J * JK.
  C.list <- vector("list", J)
  for (j in 1:J) {
    u.mat.lambda.j <- u.mat.lambdas[, j]
    temp1 <- 1 * (treat == j - 1) *
             CRFamilySecondDerivative(u.mat.lambda.j, theta)

    #print(dim(u.mat.lambdas))
    # We define the matrix A.
    tempA <- apply(u.mat, 2, "*", temp1)


    A.list[[j]] <- crossprod(tempA, u.mat) / N
    A.inv.list[[j]] <- solve(A.list[[j]])
    C.list[[j]] <- t(crossprod(u.mat, temp1 * Y) / N)
  }

  A.inv <- bdiag(A.inv.list)
  C <- bdiag(C.list)
  temp.mat <- Matrix(0, ncol = J, nrow = J * K)

  # Obtain the bread matrix.
  bread <- cbind(rbind(A.inv, C %*% A.inv),
                 rbind(temp.mat, -diag(rep(1, J))))

  # Obtain full covariance matrix.
  large.cov.mat <- (bread %*% tcrossprod(meat, bread)) / N
  # Return the sub matrix for covariance of Taus.
  large.cov.mat[-(1:(J * K)), -(1:(J * K))]
}
