# Require package numDeriv and testthat
library(numDeriv)
library(testthat)
library(assertthat)


#' Adaptive Rejection Sampling
#'
#' Generate a sample from a log-concave distribution through rejection sampling.
#'
#' @param g A possibly un-normalized density function.
#' @param n Number of samples, default is 5,000.
#' @param absc A vector of abscissae, default empty.
#' @param lower The lower bound of the domain on which g is positive, default -Inf
#' @param upper The upper bound of the domain on which g is positive, default Inf
#'
#' @return A vector of n samples from density function g.
#'
#' @export
#'
#' @examples
#' exp_func <- function(x) {dexp(x, rate = 3)} #Exponential(3)
#' test_exp <- ars(exp_func, n, absc = c(1, 4), lower = 0, upper = Inf)
#' hist(test_exp, breaks = 100, freq = FALSE)
#' curve(exp_func(x), lwd = 2, add = TRUE)
#'
#' test_normal <- ars(dnorm, n) #Normal(0,1)
#' hist(test_normal, breaks = 100, freq = FALSE)
#' curve(dnorm(x), lwd = 2, add = TRUE)
#'
#' beta_func <- function(x) {dbeta(x, 3, 4)} #Beta(3,4)
#' test_beta <- ars(beta_func, n, absc = c(0.2, 0.6), lower = 0, upper = 1)
#' hist(test_beta, breaks = 100, freq = FALSE)
#' curve(dbeta(x, 3, 4), lwd = 2, add = TRUE)
#'

ars <- function(g, n = 5000, absc = c(), lower = -Inf, upper = Inf) {

  ######## HELPER FUNCTIONS ########
  # Function: sanitize_inputs(g, n, a, l, u)
  # Input: All of the inputs supplied by the user in main ars() function
  # Purpose: Runs through some quick sanitizing checks (assertions) for bad user input
  # Output: None
  sanitize_inputs <- function(g, n, a, l, u) {
    assert_that(is.function(g), msg = "`g` is not a function, please fix")
    assert_that((n == round(n) & n > 0), msg = "`n` is not a positive integer, please fix")
    assert_that((is.numeric(absc) | (length(absc) == 0) & !is.numeric(absc)), msg = "`absc` is not a numeric vector, please fix")
    assert_that((is.numeric(lower) & length(lower) == 1), msg = "`lower` is not a single numeric value, please fix") #-Inf is numeric, so will be okay
    assert_that((is.numeric(upper) & length(upper) == 1), msg = "`upper` is not a single numeric value, please fix") #Inf is numeric, so will be okay
    assert_that(upper > lower, msg = "`upper` is not greater than `lower', please fix")
    if(length(absc) > 0) { #if not empty i.e. if user provided points
      assert_that((sum(upper > absc) == length(absc)) & sum(lower < absc) == length(absc), msg = "`absc` points not within upper and lower limits, please fix")
    }
    invisible(TRUE)
  }

  # Function: initialize_ars(a)
  # Input: Vector `a` of abscissae points
  # Purpose: Calculate `absc`, `h_vals`, `hderiv_vals`, and `z_vals`
  # Output: List of vectors: `absc`, `h_vals`, `hderiv_vals`, and `z_vals`
  initialize_ars <- function(a) {

    # Inner Helper Function: log_concavity(hd)
    # Purpose: Checks whether the vector of derivatives h' is non-increasing (equivalent to: g is log-concave)
    # Output: Boolean of whether h' is non-increasing
    log_concavity <- function(hd) {
      # hd is ordered smallest to largest since absc, the input to create hd, is ordered in the same way
      hd_shift <- hd[-1] #remove first element, so has values hd[2]-hd[k]
      hd <- head(hd, -1) #remove last element, so has values hd[1]-hd[k-1]
      hd_diff <- hd - hd_shift #vectorize, each value is hd[j] - hd[j+1] for j = 1, 2, ..., k-1.
      #Should be positive, since h' is non-increasing means that hd[j] should be no smaller than hd[j+1]

      #if all TRUE, then sum = length of vector. Use `tol` as tolerance check for 0
      return(assert_that(sum(hd_diff >= tol) == length(hd_diff),
                         msg = "Input function calculated to be not log-concave. Please check."))
    }

    # Inner Helper Function: removeDupe(vec)
    # Input: Vector `vec`
    # Purpose: Check and remove successive elements which are duplicates (equivalent) in `vec`.
    #          Code assumes `vec` is ordered already.
    #          If element i and i+1 are equivalent, keep i and flag i+1 for removal.
    # Output: Boolean list of indices to remove (TRUE = remove. FALSE = keep)
    # Example:
    #   input: (-2, -1, -1, 0, 1, 1, 1, 1, 2)
    #   output: (-2, -1, 0, 1, 2)
    removeDupe <- function(vec) {

      # (Helper Inner) Function: checkDupe(a,b)
      # Input: Scalar values `a` and `b`
      # Purpose: Checks relative spacing between `a` and `b` to see if they are equal
      # Output: Boolean on whether `a` and `b` are equal (relatively)
      # See discussion about setting `tol` in PDF/later in code.
      checkDupe <- function(a, b) {
        if (abs(a-b) < tol) {
          return(TRUE)
        }else{
          return(FALSE)
        }
      }

      remove_index <- c(FALSE) #always keep first element in `vec`
      for (i in 1:(length(vec)-1)) {
        remove_index[i+1] <- checkDupe(vec[i], vec[i+1])
      }
      return(remove_index)
    }

    # Inner Helper Function: valid_abscissae(hd)
    # Input: Vector `hd` of hderivative points
    # Purpose: Checks for valid abscissae points (from looking at the points evaluated at derivative), following 2.2.1 in the Paper
    #          Specifically, checks that if unbounded from below/above, then first/last point has hderiv >/< 0
    #          If bounds supplied, then no need to check, any initial abscissae points work.
    # Output: Boolean on whether the current abscissae points are valid
    # Note: If user gives 0 or 1 points, abscissae() will identify two points so that the process will work
    #       If user gives 2 (or more) points, abscissae() will not run and the function assumes that the given abscissae points are correct
    #       and thus valid_abscissae kicks in to check whether those points are legitimately valid or not.
    valid_abscissae <- function(hd) {
      if (lower == -Inf) {
        # Error message changes depending on what kind of abscissae point issue.
        if (length(hd) == 1) {
          msg1 <- "Unbounded from below and derivative of only valid abscissae point is non-positive."
          msg2 <- "Please either supply a lower bound, or add another abscissae point with positive derivative"
          msg <- paste(msg1, msg2, sep = " ")
        } else{
          msg1 <- "Unbounded from below and derivative of initial abscissae point is non-positive."
          msg2 <- "Please either supply a lower bound, or re-select an initial abscissae point with positive derivative"
          msg <- paste(msg1, msg2, sep = " ")
        }
        assert_that(hd[1] > tol, msg = msg)
      }
      if (upper == Inf) {
        if(length(hd) == 1) {
          msg1 <- "Unbounded from above and derivative of only valid abscissae point is non-negative."
          msg2 <- "Please either supply an upper bound, or add another abscissae point with negative derivative."
          msg <- paste(msg1, msg2, sep = " ")
        } else{
          msg1 <- "Unbounded from above and derivative of final abscissae point is non-negative."
          msg2 <- "Please either supply an upper bound, or re-select a final abscissae point with negative derivative"
          msg <- paste(msg1, msg2, sep = " ")
        }
        assert_that(hd[length(hd)] < tol, msg = msg)
      }
      invisible(TRUE) #this function is just a checker, so silently return TRUE (no output)
    }

    ### Main Part of initialize_ars()
    if (length(a) < 2) {
      absc <- abscissae(a) #if provide less than 2 inputs, then we create our own abscissae points
    }
    absc <- absc[order(absc)]
    h_vals <- eval_h(absc) #h_vals has length k
    hderiv_vals <- grad(eval_h, absc) #hderiv_vals has length k

    # Check duplicates in absc and hderiv_vals (can happen even if absc is non-duplicate)
    # Remove duplicate absc because not needed (duplicate points in abscissae vector unhelpful)
    # Remove "duplicate" hderiv_vals because that means h() is flat at those points, messes up tangent calculation
    # Do NOT check for duplicates in h_vals, because can have two different absc points that give same h()!
    if (length(absc) > 1) {
      removeIndex <- removeDupe(absc) | removeDupe(hderiv_vals)
      absc <- absc[!removeIndex]
      h_vals <- h_vals[!removeIndex] #still need to remove corresponding elements in h_vals
      hderiv_vals <- hderiv_vals[!removeIndex]
    }

    valid_abscissae(hderiv_vals) #returns TRUE, will do nothing
    log_concavity(hderiv_vals) #run a check on hderiv_vals. If it fails, then program exits (from assert_that) and will alert user.
    z_vals <- eval_z(absc, h_vals, hderiv_vals)

    return(list("absc" = absc, "h_vals" = h_vals,
                "hderiv_vals" = hderiv_vals, "z_vals" = z_vals))
  }

  # Function: abscissae(a)
  # Inputs: Vector `a` of abscissae points
  # Purpose: Finds valid points for initial abscissae if user provides only 0 or 1 points
  # Output: Vector of two abscissae points
  # Note: Assume that the domain of the function is either properly given (or approximately well-given)
  #       There will be issues with super edge cases such as if the function is Unif(1000, 1005)
  #       Then we may not find a valid abscissae point (if not provided) since our grid assumes that we are searching within -50 to 50
  abscissae <- function(a) {
    # first find the value of x where grad(log(f), x) = 0 (f is maximized)
    for (i in -5000:5000) {
      # check if the derivative exists,
      # this is actually checking if the function is meaningful at this point
      # if not, go to the next loop
      if (is.infinite(eval_h(i / 100)) | is.na(eval_h(i / 100))) {
        next
      }

      # check if it is a uniform distribution, if true, we can assign any value for absc within the domain
      # the log(uniform) is equal everywhere within the domain
      if (eval_h(i / 100) - eval_h((i + 1) / 100) == 0 & eval_h((i + 1)/100) - eval_h((i + 2) / 100) == 0) {
        grad_0 = (i + 2) / 100
        break
      }

      # check if it is an exponential distribution
      # the log(exp) is a linear function with constant slope
      # if true, we can assign any value for absc within the domain
      if (eval_h(i / 100) - eval_h((i + 1) / 100) == eval_h((i + 1) / 100) - eval_h((i + 2) / 100)) {
        grad_0 = (i + 2) / 100
        break
      }


      # find the point where the gradient is approximately 0
      if (eval_h(i / 100) - eval_h((i + 1) / 100) > 0) {
        grad_0 = i / 100
        break
      }
    }

    # if the user inputs only one number for absc, we should use it as one of our initial points
    left = NULL
    right = NULL
    if (length(a) == 1) {
      if (a > grad_0) {
        right = a
        left = grad_0 - 0.01
      } else {
        left = a
        right = grad_0 + 0.01
      }
    } else{
      left = grad_0 - 0.01
      right = grad_0 + 0.01
    }

    out <- c(left, right)

    return(out)
  }

  # Function: eval_h(x)
  # Input: vector `x` of points
  # Purpose: Evaluates h(x) = log(g(x)), where g is user-input function
  # Output: vector of h(x) points
  eval_h <- function(x) {
    suppressWarnings(assert_that(sum(!is.nan(log(g(x)))) == length(x),
                                 msg = "Log of input function is NaN. please check the input function is a valid (non-normalized) density along domain."))
    return(log(g(x)))
  }

  # Function: eval_z(a, h, hd)
  # Inputs:
  #   Vector `a` of abscissae points
  #   Vector `h` of points h(absc) i.e. function h evaluated at abscissae points
  #   Vector `hd` of points h'(absc) i.e. derivative of function h evaluated at abscissae points
  # Purpose: Evaluates the vector of points `z` i.e. points where tangent lines to abscissae intersect
  # Output: vector `z_out` of tangent intersect points
  eval_z <- function(a, h, hd) {
    if (length(a) > 1) {
      # Create original and shifted versions for vectorization of z_out
      a1 <- head(a, -1) #a ranges from index 1 to k-1
      a2 <- a[-1] #a_shift ranges from index 2 to k
      h1 <- head(h, -1)
      h2 <- h[-1]
      hd1 <- head(hd, -1)
      hd2 <- hd[-1]

      # Note that hd1-hd2 is non-zero based on cleaning in initialize_ars()
      z_out <- (h2 - h1 - a2 * hd2 + a1 * hd1)/(hd1 - hd2)

      #z_0 = `lower`, i.e. lower bound
      #z_k = `upper`, i.e. upper bound
      z_out <- c(lower, z_out, upper) #length k+1
    } else {
      z_out <- c(lower, upper) #no "intermediate" intersect points
    }
    return(z_out)
  }

  # Function: eval_u(x, z, a, h, hd)
  # Inputs:
  #   Vector `x` of points to evaluate u(x)
  #   Vector `z` of tangent intersect points
  #   Vector `a` of abscissae points
  #   Vector `h` of points h(absc)
  #   Vector `hd` of points h'(absc)
  # Purpose: Evaluates u(x) by finding which piecewise part of u(x) to use, and then calculating u(x)
  # Output: Vector `u_out` of u(x) values
  eval_u <- function(x, z, a, h, hd) {

    # Note:
    # Paper Notation: z0, z1, z2, ..., zk, with z0 = lower (-Inf) and zk = upper (+Inf)
    # R indexing: z[1], z[2], z[3], ..., z[k+1]
    # Paper Notation: x1, x2 (abscissae points), analogous for h and hd
    # R indexing: a1, a2 (abscissae points), analogous for h and hd
    # This means that if I see zj in Paper, that is z[j+1] in R
    # But xj, hj, h'j in paper corresponds to corresponds to a[j], h[j], and hd[j] in this R code

    # First figures out which element index of z that x lies in-between (i.e. z_j-1 <= x <= z_j)
    # Example:
    # z = c(-Inf, -3, -1, 4, Inf)
    # a = c(-5, -2, 0, 5)
    # x = -2.5, should be between -3 and -1, i.e. z1 and z2 -> index = 2

    # If `x` and `lower` are very close (first condition) or "equal" (used if lower = -Inf or upper = Inf)), index = 1
    # Analogous for comparing `x` and `upper`
    # If more than 1 abscissae point, this works fine
    if (length(a) > 1) {
      if ((abs(x - lower) < tol) | (x == lower)) {
        index <- 1 #between z0 = `lower` and z1
      } else if ((abs(x - upper) < tol) | (x == upper)) {
        # Note `z` has k+1 elements (Paper Notation: z0 - zk. R Code: z[1] - z[k+1])
        # So, index needs to be equal to `k` in the R code i.e. length(z) - 1
        index <- length(z) - 1 #between z_{k-1} and zk = `upper`
      } else {
        index <- min(which((z-x) >= 0)) - 1 #Same reason as above for -1
      }
    } else {
      index <- 1 #if only 1 abscissae point, only have z0 and z1, so index = 1
    }

    # Our indexes are according to `z`, so need to +1 for corresponding value in `h`. See above discussion.
    u_out <- h[index] + (x - a[index]) * hd[index]
    return(u_out)
  }

  # Function: eval_l(x, a, h)
  # Inputs:
  #   Vector `x` of points to evaluate u(x)
  #   Vector `a` of abscissae points
  #   Vector `h` of points h(absc)
  # Purpose: Evaluates l(x) by finding which piecewise part of l(x) to use, and then calculating l(x)
  # Output: Vector `l_out` of l(x) values
  eval_l <- function(x, a, h) {

    # Note:
    # Paper Notation uses intervals a in [aj, aj+1]. So, we need index = j, so min() - 1
    if (x < min(a)) {
      l_out <- -Inf
    } else if (x > max(a)) {
      l_out <- -Inf
    } else {
      index <- min(which((a-x) >= 0)) - 1
      # No issues dividing by 0 since `a` has successive duplicate elements removed
      l_out <- ((a[index+1] - x) * h[index] + (x - a[index]) * h[index + 1])/(a[index + 1] - a[index])
    }
    return(l_out)
  }

  # Function: eval_sample(a, z, h, hd)
  # Inputs:
  #   Vector `a` of abscissae points
  #   Vector `z` of tangent intersect points
  #   Vector `h` of points h(absc)
  #   Vector `hd` of points h'(absc)
  # Purpose: Draws a sample value `x_star` from the `g` distribution using ARS sampling algorithm
  #   The ARS sampling algorithm is detailed further in the PDF writeup
  # Output: Vector `x_star` of sampled points
  eval_sample <- function(a, z, h, hd) {

    # Step 1: Calculate u(x) for x = z (intersect points) and generate piecewise part to sample from
    u_vals_int <- unlist(lapply(z, eval_u, z, a, h, hd)) #k+1 elements
    diff_exp_u = c()

    #iterate over 1 to k (a has k elements)
    for (i in 1:length(a)) {
      # Case when h' = 0
      if (abs(hd[i] - 0) < tol) {
        diff_exp_u[i] = exp(h[i]) * (z[i+1] - z[i]) #exp(a[i]) because z[j-1] < a[j] < z[j],
        #but z indexing from 0 while R indexes from 1
        #so we need z[j] < a[j] < z[j+1]
      } else{
        diff_exp_u[i] = (exp(u_vals_int[i+1]) - exp(u_vals_int[i])) / hd[i] #hd[i] non-zero in this case
      }
    }
    prob_piece <- diff_exp_u / sum(diff_exp_u)
    intervals <- 1:length(a)
    # Randomly sample the index j of interval [z_{j-1}, z_j] to proceed with
    # (for j = 1, 2, ..., k, with z_0 = `lower` and z_k = `upper`)
    j <- sample(intervals, 1, replace = TRUE, prob = prob_piece)

    # Step 2: Sample from piecewise part
    p <- runif(1)
    C <- diff_exp_u[j] #normalizing constant from Step 1
    # Note: z0 = `lower` is actually z[1] in R. So, zj in Paper = z[j+1] in the code
    # PDF uses Paper notation, so u_(j-1) in writeup should be u_vals_int[j] in Code
    # Indices of `a`, `h`, `hd` are still `j`, no off-by-one issue for those
    if (abs(hd[j] - 0) < tol) {
      # Case when h' = 0
      x_star <- z[j] + (p*C)/exp(h[i])
    } else{
      x_star <- a[j] + (-h[j] + log(hd[j] * p * C + exp(u_vals_int[j])))/hd[j] #hd[i] non-zero in this case
    }
    return(x_star)
  }

  # Function: accept_reject(x, z, a, h, hd)
  # Inputs:
  #   Vector `x` of points to evaluate whether to accept or reject (i.e. `x_star`)
  #   Vector `z` of tangent intersect points
  #   Vector `a` of abscissae points
  #   Vector `h` of points h(absc)
  #   Vector `hd` of points h'(absc)
  # Purpose: Determine whether to accept or reject `x_star` (and what to do if so)
  # Output: Code of either 0, 1, 2
  #   0: Reject (i.e. do nothing)
  #   1: Accept, do not update abscissae
  #   2: Accept, update abscissae
  accept_reject <- function(x, z, a, h, hd) {
    w <- runif(1)
    if (w <= exp(eval_l(x, a, h) - eval_u(x, z, a, h, hd))) {
      return(1)
    } else if (w <= exp(eval_h(x) - eval_u(x, z, a, h, hd))) {
      if (abs(exp(eval_h(x) - eval_u(x, z, a, h, hd)) - 1) < tol) {
        return(1) #do not update abscissae in this special circumstance
      } else{
        return(2)
      }
    } else {
      return(0)
    }
  }

  ######## MAIN ########
  tol <- 1e-8 #Based on discussion in writeup
  num_samples <- 0
  total_runs <- 0
  sample_output <- c()
  x_star_list <- c() #Tracker which holds all x_stars, keeping for debugging

  sanitize_inputs(g, n, absc, lower, upper)
  if (length(absc) > 0) {
    absc <- absc[order(absc)]
  }
  initial <- initialize_ars(absc)
  absc <- initial$absc
  h_vals <- initial$h_vals
  hderiv_vals <- initial$hderiv_vals
  z_vals <- initial$z_vals

  ### SAMPLE AND UPDATE ###
  while ((num_samples < n) & (total_runs < 10*n)) { #default is n = 5000
    x_star <- eval_sample(absc, z_vals, h_vals, hderiv_vals)
    x_star_list <- c(x_star_list, x_star)

    if (accept_reject(x_star, z_vals, absc, h_vals, hderiv_vals) == 1) { #accepted from first comparison
      sample_output <- c(sample_output, x_star)
      num_samples <- num_samples + 1
      total_runs <- total_runs + 1
    } else if (accept_reject(x_star, z_vals, absc, h_vals, hderiv_vals) == 2) { #accepted from second comparison
      # Add `x_star` to list of abscissae points
      absc <- c(absc, x_star)
      absc <- absc[order(absc)]

      # Update with new absc, h_vals, hderiv_vals, z_vals
      initial <- initialize_ars(absc)
      absc <- initial$absc
      h_vals <- initial$h_vals
      hderiv_vals <- initial$hderiv_vals
      z_vals <- initial$z_vals

      sample_output <- c(sample_output, x_star)
      num_samples <- num_samples + 1
      total_runs <- total_runs + 1
    } else{ #rejected, pick a new sample
      total_runs <- total_runs + 1
    }
  }
  return(sample_output) #returns vector of sampled points from ARS
}
