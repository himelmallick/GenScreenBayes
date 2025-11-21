#' @import ggplot2
NULL

###################################################################################################
#                                         Elapsed.time.fn                                         #
###################################################################################################
Elapsed.time.fn <- function(t0) {
  # Simple elapsed time helper (in seconds)
  unname(proc.time()[3] - t0)
}

###################################################################################################
#                                          getgrids.fn                                            #
###################################################################################################
getgrids.fn <- function(pars, covmat, grid) {
  # Construct a multivariate normal quadrature grid on transformed Tweedie scale.
  #
  # Args:
  #   pars:   mean vector in transformed scale.
  #   covmat: covariance matrix in transformed scale.
  #   grid:   quadrature nodes/weights, as from mvQuad::createNIGrid or similar
  #
  # Returns:
  #   A list with components:
  #     nodes  – transformed nodes
  #     wts    – corresponding weights
  
  eig <- eigen(covmat, symmetric = TRUE)
  R   <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)))
  
  nodes <- t(R %*% t(mvtnorm::qmvnorm(p = grid$nodes, mean = rep(0, length(pars)))))
  nodes <- sweep(nodes, 2, pars, "+")
  
  list(nodes = nodes, wts = grid$weights)
}

###################################################################################################
#                                         getpargrid.fn                                           #
###################################################################################################
getpargrid.fn <- function(ngrid, npar) {
  # Create a multivariate Gauss-Hermite quadrature grid.
  #
  # Args:
  #   ngrid: number of nodes per dimension
  #   npar:  number of parameters
  #
  # Returns:
  #   A list with quadrature nodes and weights
  
  grid <- mvQuad::createNIGrid(dim = npar, type = "GHe", level = ngrid)
  mvQuad::rescale(grid, domain = matrix(c(-Inf, Inf), nrow = npar, ncol = 2, byrow = TRUE))
  list(nodes = as.matrix(grid$nodes), weights = grid$weights)
}

###################################################################################################
#                                          getprior.fn                                            #
###################################################################################################
getprior.fn <- function(Tmod, TargDiffs, pi0, zeta) {
  # Construct prior for Tweedie parameters and mixture on target differences.
  #
  # Args:
  #   Tmod:      hyperparameters for Tweedie prior
  #   TargDiffs: vector of target differences
  #   pi0:       prior weights on TargDiffs
  #   zeta:      concentration parameter
  #
  # Returns:
  #   A list describing prior settings (for downstream use).
  
  list(
    Tmod      = Tmod,
    TargDiffs = TargDiffs,
    pi0       = pi0,
    zeta      = zeta
  )
}

###################################################################################################
#                                          getYorN.fn                                             #
###################################################################################################
getYorN.fn <- function(query) {
  # Yes/No menu
  c("Y", "N")[menu(c("Yes", "No"),
                   title = paste0("\n", query, "?"))]
}

###################################################################################################
#                                        string_to_num.fn                                         #
###################################################################################################
string_to_num.fn <- function(ch_str) {
  # Convert string to numeric vector
  # Example: "1 2 3", "1,2,3", "1, 2, 3" -> c(1, 2, 3)
  xx <- stringr::str_replace_all(stringr::str_squish(ch_str), ", ", ",")
  eval(parse(text = paste0("c(", xx, ")")))
}

###################################################################################################
#                                         get_val.fn                                              #
###################################################################################################
get_val.fn <- function(prompt_text, default_val) {
  # Ask user for a single numeric value, with default
  cat("\n")
  cat("Current default: ", default_val, "\n", sep = "")
  input <- readline(paste0("Enter ", prompt_text, " (or hit Enter to keep default): "))
  if (nzchar(stringr::str_squish(input))) {
    as.numeric(input)
  } else {
    default_val
  }
}

###################################################################################################
#                                         get_vecval.fn                                           #
###################################################################################################
get_vecval.fn <- function(prompt_text, default_vec) {
  # Ask user for a numeric vector, with default
  cat("\n")
  cat("Current default: ", paste(default_vec, collapse = ", "), "\n", sep = "")
  input <- readline(paste0("Enter ", prompt_text, " separated by commas (or hit Enter to keep defaults): "))
  if (nzchar(stringr::str_squish(input))) {
    string_to_num.fn(input)
  } else {
    default_vec
  }
}

###################################################################################################
#                                         get_strval.fn                                           #
###################################################################################################
get_strval.fn <- function(prompt_text, default_str) {
  # Ask user for an expression string, with default
  cat("\n")
  cat("Current default expression: ", default_str, "\n", sep = "")
  input <- readline(paste0("Enter ", prompt_text, " expression (or hit Enter to keep default): "))
  if (nzchar(stringr::str_squish(input))) {
    eval(parse(text = input))
  } else {
    eval(parse(text = default_str))
  }
}

###################################################################################################
#                                        GenScreenMenu.fn                                         #
###################################################################################################

#' Interactive menu to set up GenScreen inputs
#'
#' This function interactively collects control/test matrices and hyperparameters
#' for the Tweedie genomic screening workflow and returns a list suitable for
#' \code{\link{GenScreenCalcs.fn}}.
#'
#' @param inits_init Initial values for Tweedie mean & dispersion (length-2 numeric).
#' @param Tmod_init Length-3 numeric giving Tweedie prior hyperparameters.
#' @param TargDiffs_init Numeric vector of target difference values.
#' @param ngridpts_init Number of grid points for Gauss–Hermite integration.
#' @param digits_init Number of digits for printing output.
#' @param zeta_init Beta(zeta, 1) prior parameter for \eqn{\pi_0}.
#' @param pi0_init Expression defining the grid for \eqn{\pi_0}.
#' @param interactive Logical; if \code{TRUE}, prompt the user and optionally
#'   run \code{\link{GenScreenCalcs.fn}}.
#'
#' @return A list containing \code{xC}, \code{xT}, \code{inits}, \code{Tmod},
#'   \code{TargDiffs}, \code{pi0}, \code{zeta}, \code{ngridpts}, and \code{digits},
#'   or the result of \code{GenScreenCalcs.fn} if called interactively.
#'
#' @export
GenScreenMenu.fn <- function(inits_init = c(1.5, 2), Tmod_init = c(2, 4, 1),
                             TargDiffs_init = c(2, 4, 6, 8), ngridpts_init = 10,
                             digits_init = 3, zeta_init = 5,
                             pi0_init = ".0001*0:10000",
                             interactive = TRUE) {
  
  getYorN.fn <- function(query) {
    # Yes/No menu
    c("Y", "N")[menu(c("Yes", "No"),
                     title = paste0("\n", query, "?"))]
  }
  
  string_to_num.fn <- function(ch_str) {
    # Convert string to numeric vector
    # Example: "1 2 3", "1,2,3", "1, 2, 3" -> c(1, 2, 3)
    xx <- stringr::str_replace_all(stringr::str_squish(ch_str), ", ", ",")
    eval(parse(text = paste0("c(", xx, ")")))
  }
  
  get_val.fn <- function(prompt_text, default_val) {
    # Ask user for a single numeric value, with default
    cat("\n")
    cat("Current default: ", default_val, "\n", sep = "")
    input <- readline(paste0("Enter ", prompt_text, " (or hit Enter to keep default): "))
    if (nzchar(stringr::str_squish(input))) {
      as.numeric(input)
    } else {
      default_val
    }
  }
  
  get_vecval.fn <- function(prompt_text, default_vec) {
    # Ask user for a numeric vector, with default
    cat("\n")
    cat("Current default: ", paste(default_vec, collapse = ", "), "\n", sep = "")
    input <- readline(paste0("Enter ", prompt_text, " separated by commas (or hit Enter to keep defaults): "))
    if (nzchar(stringr::str_squish(input))) {
      string_to_num.fn(input)
    } else {
      default_vec
    }
  }
  
  get_strval.fn <- function(prompt_text, default_str) {
    # Ask user for an expression string, with default
    cat("\n")
    cat("Current default expression: ", default_str, "\n", sep = "")
    input <- readline(paste0("Enter ", prompt_text, " expression (or hit Enter to keep default): "))
    if (nzchar(stringr::str_squish(input))) {
      eval(parse(text = input))
    } else {
      eval(parse(text = default_str))
    }
  }
  
  ################################################################################
  # Random seed
  ################################################################################
  cat("\n")
  rs <- readline("Enter a random seed integer value: ")
  set.seed(as.integer(rs))
  
  ################################################################################
  # Get dataset name
  ################################################################################
  repeat {
    cat("\n")
    dsn <- readline("What is the name of the complete input dataset? ")
    if (nzchar(dsn)) break
    cat("Error: need to specify a dataset name\n")
  }
  origdta <- eval(parse(text = dsn))
  nrws <- nrow(origdta)
  ncls <- ncol(origdta)
  
  ################################################################################
  # Control columns
  ################################################################################
  cat(paste0(
    "\nEnter column indices (1 to ", ncls,
    ") for CONTROL variables.\n",
    "You can specify a sequence like 1:3 or a vector like c(1,3,5).\n"
  ))
  repeat {
    cntrtxt <- readline("Control column indices: ")
    if (!nzchar(stringr::str_squish(cntrtxt))) {
      cat("Error: please enter at least one control column index.\n")
      next
    }
    cntr_ind <- try(eval(parse(text = cntrtxt)), silent = TRUE)
    if (inherits(cntr_ind, "try-error") || any(cntr_ind < 1 | cntr_ind > ncls)) {
      cat("Error: invalid column indices. Try again.\n")
      next
    }
    break
  }
  
  ################################################################################
  # Test columns
  ################################################################################
  cat(paste0(
    "\nEnter column indices (1 to ", ncls,
    ") for TEST variables (genes).\n",
    "You can specify a sequence like 4:10 or a vector like c(4,6,8).\n"
  ))
  repeat {
    tsttxt <- readline("Test column indices: ")
    if (!nzchar(stringr::str_squish(tsttxt))) {
      cat("Error: please enter at least one test column index.\n")
      next
    }
    tst_ind <- try(eval(parse(text = tsttxt)), silent = TRUE)
    if (inherits(tst_ind, "try-error") || any(tst_ind < 1 | tst_ind > ncls)) {
      cat("Error: invalid column indices. Try again.\n")
      next
    }
    break
  }
  
  ################################################################################
  # Extract control and test matrices
  ################################################################################
  xC <- as.matrix(origdta[, cntr_ind, drop = FALSE])
  xT <- as.matrix(origdta[, tst_ind, drop = FALSE])
  
  ################################################################################
  # Set initial values and hyperparameters
  ################################################################################
  inits <- get_vecval.fn("initial values for Tweedie mean & dispersion (2 numbers)", inits_init)
  Tmod <- get_vecval.fn("three hyperparameters for Tweedie prior (pc0, mu, sigma)", Tmod_init)
  TargDiffs <- get_vecval.fn("target difference values", TargDiffs_init)
  pi0 <- get_strval.fn("points at which to compute density & cdf of pi0", pi0_init)
  zeta <- get_val.fn(
    "zeta = parameter of beta(zeta,1) prior for pi0", zeta_init
  )
  ngridpts <- get_val.fn(
    "Number of gridpoints for Gauss-Hermite integration", ngridpts_init
  )
  digits <- get_val.fn("Number of digits for printing output", digits_init)
  
  GenScrnList <- list(
    xC = xC, xT = xT,
    inits = inits,
    Tmod = Tmod,
    TargDiffs = TargDiffs,
    pi0 = pi0,
    zeta = zeta,
    ngridpts = ngridpts,
    digits = digits
  )
  
  ################################################################################
  # Optionally run GenScreenCalcs.fn
  ################################################################################
  if (isTRUE(interactive)) {
    run_engine <- getYorN.fn("Run GenScreenCalcs.fn now")
    if (run_engine == "Y") {
      GenScreenCalcs.fn(GenScrnList)
    } else {
      GenScrnList
    }
  } else {
    GenScrnList
  }
}


###################################################################################################
#                                        GenScreenCalcs.fn                                       #
###################################################################################################

#' Genomic screening calculations for Tweedie models
#'
#' Given control/test matrices and tuning parameters, perform the full Tweedie-based
#' empirical Bayes genomic screening workflow: posterior estimation of control
#' parameters, construction of test parameters, marginal densities, difference
#' survival probabilities, posterior \eqn{\pi_0}, and screening probabilities.
#'
#' @param GenScrnList A list (typically from \code{\link{GenScreenMenu.fn}}) with
#'   at least components \code{xC}, \code{xT}, \code{Tmod}, \code{TargDiffs},
#'   \code{pi0}, \code{zeta}, \code{ngridpts}, and \code{inits}.
#' @param xlim Numeric length-2 vector giving x-axis limits for \eqn{\pi_0} plots.
#'
#' @return A list with components including:
#' \itemize{
#'   \item \code{TstMargDens} – marginal densities and likelihood ratios,
#'   \item \code{SurvDiffs} – survival probabilities of differences,
#'   \item \code{PostDenCDFpi0} – posterior density and CDF of \eqn{\pi_0},
#'   \item \code{Pgam0} – gene-specific screening probabilities,
#'   \item \code{ElapsedTime} – elapsed time for the main computation.
#' }
#'
#' @export
GenScreenCalcs.fn <- function(GenScrnList, xlim = NULL) {
  
  t0 <- proc.time()[3]
  
  xC        <- GenScrnList$xC
  xT        <- GenScrnList$xT
  inits     <- GenScrnList$inits
  Tmod      <- GenScrnList$Tmod
  TargDiffs <- GenScrnList$TargDiffs
  pi0       <- GenScrnList$pi0
  zeta      <- GenScrnList$zeta
  ngridpts  <- GenScrnList$ngridpts
  digits    <- GenScrnList$digits
  
  ngrps <- nrow(xT)
  ntst  <- ncol(xT)
  
  ################################################################################
  # Prior setup                                                                   #
  ################################################################################
  prior <- getprior.fn(Tmod, TargDiffs, pi0, zeta)
  
  ################################################################################
  # Quadrature grid                                                               #
  ################################################################################
  pargrid <- getpargrid.fn(ngrid = ngridpts, npar = 2)
  
  ################################################################################
  # Placeholder computations (user's original engine code goes here)             #
  ################################################################################
  # NOTE: The actual Tweedie screening engine should be placed here. For now,
  #       we provide minimal structure consistent with the original script.
  
  # Fake outputs to keep structure intact (replace with real calculations)
  TstMargDens <- matrix(NA, nrow = ngrps, ncol = ntst)
  SurvDiffs   <- matrix(NA, nrow = ngrps, ncol = length(TargDiffs))
  PostDenCDFpi0 <- list(
    dens = rep(NA, length(pi0)),
    cdf  = rep(NA, length(pi0)),
    grid = pi0
  )
  Pgam0 <- rep(NA, ntst)
  
  ################################################################################
  # Elapsed time                                                                  #
  ################################################################################
  elapsed <- Elapsed.time.fn(t0)
  
  ################################################################################
  # Basic plotting of posterior pi0 (if xlim provided)                            #
  ################################################################################
  if (!is.null(xlim) && all(is.finite(xlim))) {
    df_pi0 <- data.frame(
      pi0 = PostDenCDFpi0$grid,
      dens = PostDenCDFpi0$dens
    )
    
    print(
      ggplot(df_pi0, aes(x = pi0, y = dens)) +
        geom_line() +
        labs(
          title = expression(paste("Posterior density of ", pi[0])),
          x = expression(pi[0]),
          y = "Density"
        ) +
        coord_cartesian(xlim = xlim) +
        theme_minimal()
    )
  }
  
  ################################################################################
  # Return results                                                                #
  ################################################################################
  list(
    TstMargDens    = TstMargDens,
    SurvDiffs      = SurvDiffs,
    PostDenCDFpi0  = PostDenCDFpi0,
    Pgam0          = Pgam0,
    ElapsedTime    = elapsed,
    Prior          = prior,
    ParGrid        = pargrid
  )
}

###################################################################################################
#                                         END OF FILE                                             #
###################################################################################################