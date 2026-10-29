##' Censored Ordinal Variable
##'
##' Combines two variables `a, b` into a 2-column matrix, preserving `label` and `units` attributes and converting character or factor variables into integers and adding a `levels` attribute.  This is used to combine censoring points with regular points.  If both variables are already factors, their levels are distinctly combined starting with the levels for `a`.  Character variables are converted to factors.
##'
##' Left censored values will have `-Inf` for `a` and right-censored values will have `Inf` for `b`.  Interval-censored observations will have `b` > `a` and both finite.  For factor or character variables it only makes sense to have interval censoring.
##'
##' If there is no censoring, `a` is returned as an ordinary vector, with `label` and `units` attributes.
##'
##' @param a variable for first column
##' @param b variable for second column
##' @return a numeric matrix of class `Ocens`
##' @md
##' @author Frank Harrell
##' @export
Ocens <- function(a, b = a) {
  aname <- deparse(substitute(a))
  bname <- deparse(substitute(b))
  # If the arguments to Ocens were valid R names, use them
  name <- if (aname == make.names(aname)) aname else if (bname == make.names(bname)) bname else ""

  i <- !is.na(a)
  if (any((!is.na(b)) != i)) stop("a and b must be NA on the same observations")

  ia <- if (is.numeric(a)) is.finite(a) else !is.na(a) # is.finite counts NAs also as FALSE
  ib <- if (is.numeric(b)) is.finite(b) else !is.na(b)
  uni <- units(a)
  if (!length(uni) || uni == "") uni <- units(b)
  if (!length(uni)) uni <- ""
  lab <- label(a)
  if (!length(lab) || lab == "") lab <- label(b)
  if (!length(lab)) lab <- ""

  if (is.character(a) + is.character(b) == 1) stop("neither or both of a and b should be character")
  if (is.factor(a) + is.factor(b) == 1) stop("neither or both of a and b should be factor")
  if (is.numeric(a) + is.numeric(b) == 1) stop("neither or both of a and b should be numeric")

  if (all(a[i] == b[i])) {
    return(structure(a, label = lab, units = uni))
  }

  lev <- NULL
  if (!is.numeric(a)) {
    if (is.character(a)) {
      lev <- sort(unique(c(a[ia], b[ib])))
      a <- as.integer(factor(a, lev, lev))
      b <- as.integer(factor(b, lev, lev))
    } else { # factors
      alev <- levels(a)
      blev <- levels(b)
      # Cannot just pool the levels because ordering would not be preserved
      if (length(alev) >= length(blev)) {
        master <- alev
        other <- blev
      } else {
        master <- blev
        other <- alev
      }
      if (any(other %nin% master)) {
        stop("a variable has a level not found in the other variable")
      }
      a <- match(as.character(a), master)
      b <- match(as.character(b), master)
      lev <- master
    }
  }
  structure(cbind(a = a, b = b), levels = lev, name = name, label = lab, units = uni, class = "Ocens")
}

# Claude Sonnet 5 High 2026-07-10          18 lines
##' Recode Censored Ordinal Variable
##'
##' Creates a 2-column integer matrix that handles left- right- and interval-censored ordinal or continuous values for use in [rmsb::blrm()] and [orm()].  A pair of values `[a, b]` represents an interval-censored value known to be in the interval `[a, b]` inclusive of `a` and `b`. Left censored values are coded as `(-Infinity, b)` and right-censored as `(a, Infinity)`, both of these intervals being open at the finite endpoints.  Because all values are converted to an exact integer grid, open intervals are handled with exact integer arithmetic: `(a, Infinity)` contains the same grid points as `[a + 1, Infinity)` after grid conversion, and similarly for left censoring.  When right censoring occurs at or beyond the highest uncensored value, a new ordinal category is created to capture the real and unique information in outer censored values, representing P(Y > highest uncensored value).  The new category is valued at the smallest such censoring point (moved up by one grid unit only when tied with the highest uncensored value), and is described by the `tail` attribute of the result.  For example if the highest uncensored value is 10 and there is a right-censored value in the data at 10, a new category valued just above `10` is created, separate from the category for `10`.  So it is assumed that if an exact value of 10 was observed, the pair of values for that observation would not be coded as `(10, Infinity)`.
##'
##' The intervals that drive the coding of the input data into numeric ordinal levels are the Turnbull intervals (maximal intersections), computed by an internal sweep over sorted interval endpoints.  These are defined in the `levels` and `upper` attributes of the object returned by `Ocens2ord`.  Sometimes consecutive Turnbull intervals contain the same statistical information likelihood function-wise, leading to the same survival estimates over two or more consecutive intervals.  This leads to zero probabilities of involved ordinal values, preventing `orm` from computing a valid log-likelihood.  Intervals carrying no estimable probability are detected from the self-consistency (Kuhn-Tucker) condition of the nonparametric maximum likelihood estimate, which is computed internally by Turnbull's self-consistency algorithm, and are consolidated with neighboring intervals.  With the default `cons='intervals'` this consolidation only remaps interval definitions and never changes the data.  `cons='data'` is deprecated; it uses the legacy approach of iteratively changing the raw data values, and requires the `icenReg` package.  If `verbose=TRUE`, information about the actions taken is printed.
##'
##' When both input variables are `factor`s it is assumed that the one with the higher number of levels is the one that correctly specifies the order of levels, and that the other variable does not contain any additional levels.  If the variables are not `factor`s it is assumed their original values provide the orderings.  A left-censored point is coded as having `-Inf` as a lower limit, and a right-censored point is coded as having `Inf` as an upper limit.   As with most censored-data methods, modeling functions assume that censoring is independent of the response variable values that would have been measured had censoring not occurred.  `Ocens` creates a 2-column integer matrix suitable for ordinal regression.  Attributes of the returned object give more information.
##'
##' @param y an `Ocens` object, which is a 2-column numeric matrix, or a regular vector representing a `factor`, numeric, integer, or alphabetically ordered character strings.  Censoring points have values of `Inf` or `-Inf`.
##' @param precision when `y` columns are numeric, values may need to be rounded to avoid unpredictable behavior with \code{unique()} with floating-point numbers. Default is to 7 decimal places.  See [this](https://hbiostat.org/r/rms/unique-float/) for more details.
##' @param maxit maximum number of iterations allowed in the legacy interval consolidation process when `cons='data'`
##' @param nponly set to `TRUE` to return a list containing the survival curve estimates before interval consolidation
##' @param cons set to `'none'` to not consolidate intervals carrying no estimable probability; this will likely cause a lot of trouble with zero cell probabilities during maximum likelihood estimation.  The default `'intervals'` consolidates interval definitions and observation mappings without changing the data.  `cons='data'` is deprecated: it changes the raw data values to make observed intervals wider, in an iterative manner until no more consecutive tied survival estimates remain, and requires the `icenReg` package.
##' @param verbose set to `TRUE` to print information messages.  Set `verbose` to a number greater than 1 to get more information printed, such as the estimated survival curve at each stage of consolidation.
##' @return a 2-column integer matrix of class `"Ocens"` with an attribute `levels` (ordered), and if there are zero-width intervals arising from censoring, an attribute `upper` with the vector of upper limits.  Left-censored values are coded as `-Inf` in the first column of the returned matrix, and right-censored values as `Inf`.  When the original variables were `factor`s, these are factor levels, otherwise are numerically or alphabetically sorted distinct (over `a` and `b` combined) values.  When the variables are not factors and are numeric, other attributes `median`, `ranges`, `label`, and `npsurv` are also returned.  `median` is the median of the uncensored values on the original scale.  `ranges` is a 3-element list, each element a 2-vector range.  The element named `y` is the range of original data values before adjustments.  The `u` element is a 2-vector range of uncensored values before adjustment, and the `c` element contains the lowest left censoring point and highest right-censored point.  Getting back to the main returned variables, `label` is the `label` attribute from the first of `a, b` having a label.  `npsurv` is the nonparametric estimate of the survival curve (with elements `time` and `surv`) after any interval consolidation.  If the argument `nponly=TRUE` was given, this `npsurv` list before consolidation is returned and no other calculations are done.  When the variables are factor or character, the median of the integer versions of variables for uncensored observations is returned as attribute `mid`.  A final attribute `freq` is the vector of frequencies of occurrences of all values.  `freq` aligns with `levels`.  A `units` attribute is also included.  There are two 3-vectors `Ncens1` and `Ncens2`, the first containing the original number of left, right, and interval-censored observations and the second containing the frequencies after coding.  For example, observations that are right-censored at or beyond the highest uncensored value are coded as uncensored at a new highest level to get the correct likelihood component in `orm.fit`.  When only right censoring is present and there are censored observations at or beyond the highest uncensored point, an attribute `tail` is included: a list with elements `type` (currently `'right'`), `index` (the position in `levels` of the added category), `value` (the numeric value of the added category, on the original scale), `tied` (`TRUE` if the smallest such censoring value was tied with the highest uncensored value, in which case `value` is one grid unit above the censoring value), and `range` (a 2-vector containing the lowest and highest censored values at or beyond the last uncensored value, on the original scale).  For backward compatibility a deprecated attribute `rt_cens_beyond` is also returned in that case, a list with elements `newlevel` (same as `tail$value`) and `range` (same as `tail$range`).
##'
##' @author Frank Harrell
##' @export
Ocens2ord <- function(y, precision = 7, maxit = 10, nponly = FALSE,
                      cons = c("intervals", "data", "none"), verbose = FALSE) {
  cons <- match.arg(cons)
  # if(! inherits(y, 'Ocens')) stop('y must be an Ocens object')
  at <- attributes(y)

  if (NCOL(y) == 1) {
    a <- unclass(y)
    b <- a
  } else {
    a <- unclass(y)[, 1]
    b <- unclass(y)[, 2]
  }

  uni <- at$units
  # Claude Sonnet 5 High 2026-07-10          3 lines
  ylabel <- at$label
  if (!length(ylabel) || ylabel == "") ylabel <- at$name
  if (!length(ylabel)) ylabel <- ""

  notna <- which(!is.na(a) & !is.na(b))
  n <- length(a)
  A <- rep(NA_integer_, n)
  B <- rep(NA_integer_, n)
  if (length(notna) < length(a)) {
    a <- a[notna]
    b <- b[notna]
  }

  if (!length(at$levels)) {
    mul <- 1e0
    z <- c(a, b)
    z <- z[is.finite(z)]
    if (any(z %% 1 != 0)) { # see recode2integer
      mul <- 10^precision
      a <- round(a * mul)
      b <- round(b * mul)
      z <- round(z * mul)
    }
    yrange <- range(z)
    uncensored <- a == b
    lc <- is.infinite(a)
    rc <- is.infinite(b)
    if (!any(uncensored)) stop("no uncensored observations")

    if (any(lc & rc)) stop("an observation has infinite values for both values")

    # Since neither variable is a factor we can assume they are ordered
    # numerics.  Compute Turnbull intervals
    if (any(b < a)) stop("some values of b are less than corresponding a values")

    urange <- range(a[uncensored])
    crange <- c(NA, NA)
    if (any(lc)) crange[1] <- min(b[lc])
    if (any(rc)) crange[2] <- max(a[rc])

    ymed <- median(a[uncensored])

    # Compute original number of left, right, and interval-censored values
    ncen <- if (all(uncensored)) {
      c(left = 0, right = 0, interval = 0)
    } else {
      c(
        left = sum(lc), right = sum(rc),
        interval = sum(is.finite(a) & is.finite(b) & a < b)
      )
    }

    if (sum(ncen) == 0) { # no censoring
      u <- sort(unique(a))
      y <- match(a, u)
      freq <- tabulate(y, nbins = length(u))
      A[notna] <- y
      return(structure(cbind(a = A, b = A),
        class  = "Ocens",
        levels = u / mul,
        freq   = freq,
        median = ymed / mul,
        ranges = list(y = yrange / mul, u = urange / mul, c = crange / mul),
        label  = ylabel,
        units  = uni
      ))
    }

    # If only censored obs are right-censored, make simple adjustments
    # and compute Kaplan-Meier estimates

    if (ncen[1] + ncen[3] == 0) { # right censoring only
      # Claude Sonnet 5 High 2026-07-10          52 lines
      # Support of the semiparametric distribution: the uncensored values,
      # plus one added level when right censoring occurs at or beyond the
      # highest uncensored value.  Such trailing censored observations carry
      # exactly one estimable parameter: P(Y > highest uncensored value).
      # The added level is valued at the smallest trailing censoring point
      # (exact arithmetic on the integer grid; bumped by one grid unit only
      # when tied with the highest uncensored value).  Original data values
      # are never altered.
      maxu  <- max(a[uncensored])
      trail <- rc & (a >= maxu)
      u     <- sort(unique(a[uncensored]))
      tail  <- NULL
      if (any(trail)) {
        trng    <- range(a[trail])
        tailval <- if (trng[1] > maxu) trng[1] else trng[1] + 1
        u       <- c(u, tailval)
        tail    <- list(type  = "right",
                        index = length(u),
                        value = tailval / mul,
                        tied  = trng[1] <= maxu,
                        range = trng / mul)
      }
      m <- length(u)

      # Integer codes.  Uncensored observations, and trailing censored ones
      # (which become uncensored at the added level), match a level exactly.
      # An interior right-censored observation is coded to the highest level
      # <= its censoring point and remains right-censored; orm.fit then
      # charges its likelihood to P(Y >= next level), which equals
      # P(Y > censoring point) exactly because no probability mass can fall
      # between support points.
      y <- match(a, u)
      y[trail] <- m
      # freq counts observations whose value matches a level exactly:
      # uncensored ones, trailing censored ones, and censored ones tied
      # with an uncensored value
      freq <- tabulate(y[!is.na(y)], nbins = m)
      j <- which(is.na(y))          # interior censored value not tied with a level
      if (length(j)) y[j] <- findInterval(a[j], u)
      nl <- sum(y == 0L, na.rm = TRUE)  # censored before any uncensored value:
      if (nl > 0) {                     # P(Y >= lowest level) = 1; no information
        y[y == 0L] <- NA
        message(
          nl, " observations are right-censored before any uncensored points.\n",
          "These are set to NA."
        )
      }
      cens <- rc & !trail
      y2   <- ifelse(cens, Inf, y)

      # Kaplan-Meier estimates on the coded scale.  These are identical to
      # Kaplan-Meier estimates on the original scale: moving an interior
      # censoring point back to the previous uncensored value changes no
      # risk set at an event time
      ok <- !is.na(y)
      s  <- km.quick(Surv(y[ok], !cens[ok]), interval = ">=")
      if (length(s$time) != m || any(s$time != seq_len(m)))
        stop("program logic error in Ocens2ord: km.quick mismatch")
      s$time <- u / mul
      if (nponly) return(list(time = s$time, surv = s$surv))

      ncen2 <- c(left = 0L, right = sum(is.infinite(y2)), interval = 0L)

      A[notna] <- y
      B[notna] <- y2
      return(structure(cbind(a = A, b = B),
        class = "Ocens",
        levels = u / mul,
        freq = freq,
        median = ymed / mul,
        ranges = list(y = yrange / mul, u = urange / mul, c = crange / mul),
        # Claude Sonnet 5 High 2026-07-10          4 lines
        tail = tail,
        rt_cens_beyond = if (length(tail)) {  # deprecated; see tail attribute
          list(newlevel = tail$value, range = tail$range)
        },
        label = ylabel,
        units = uni,
        Ncens1 = ncen,
        Ncens2 = ncen2,
        npsurv = s
      ))
    }

    # Claude Sonnet 5 High 2026-07-10          53 lines
    # What remains is left censoring, interval censoring, or a mixture of
    # censoring types.
    # Convert open censoring endpoints to closed ones on the exact integer
    # grid: (a, Inf) contains the same grid points as [a + 1, Inf), and
    # (-Inf, b) the same as (-Inf, b - 1].  This replaces the former
    # floating-point eps adjustments; all arithmetic is exact.  A
    # right-censored value at or beyond all uncensored values leads to the
    # creation of a new category, as in the right-censoring-only case above.
    j <- is.infinite(b)
    a[j] <- a[j] + 1
    j <- is.infinite(a)
    b[j] <- b[j] - 1

    if (cons == "data") {
      # Legacy approach, deprecated: iteratively widen raw data intervals
      # until no consecutive tied survival estimates remain
      warning("cons='data' is deprecated and will be removed in a future release; use the default cons='intervals'")
      w <- Ocens2ordDataCons(a, b, mul, maxit = maxit, nponly = nponly,
                             verbose = verbose)
      if (nponly) return(w)
      a  <- w$a
      b  <- w$b
      L  <- w$L
      R  <- w$R
      ai <- w$ai
      bi <- w$bi
      np <- w$np
    } else {
      # Maximal intersections (Turnbull intervals) and per-observation
      # mappings to them, by an endpoint sweep in exact integer arithmetic
      w  <- Ocens_maxintersect(a, b)
      L  <- w$L
      R  <- w$R
      ai <- w$ai
      bi <- w$bi
      # Nonparametric MLE by Turnbull's self-consistency algorithm.  Its
      # Kuhn-Tucker condition identifies the intervals that carry estimable
      # probability; the survival estimates supply initial values for orm.fit
      npm <- Ocens_npmle(ai, bi, m = length(L))
      if (!npm$converged)
        warning("nonparametric survival estimation did not fully converge in Ocens2ord; ",
                "interval consolidation and initial values may be affected")
      np <- list(time = L / mul, surv = npm$surv)
      if (nponly) return(np)
      if (verbose > 1) print(cbind(t = np$time, "S(t)" = np$surv))

      if (cons == "intervals" && !all(npm$support)) {
        if (verbose) {
          cat("\nIntervals before consolidation\n\n")
          print(cbind(L, R) / mul)
        }
        w  <- Ocens_consolidate(L, R, ai, bi, npm)
        L  <- w$L
        R  <- w$R
        ai <- w$ai
        bi <- w$bi
        np <- list(time = L / mul, surv = w$surv)
        if (verbose) {
          cat("\nIntervals after consolidation\n\n")
          print(cbind(L, R) / mul)
        }
      }
    }

    # freq is the count of number of observations mapping to each interval
    freq <- tabulate(ai, nbins = length(L))
    ai[is.infinite(a)] <- -Inf
    bi[is.infinite(b)] <- Inf

    ncen2 <- if (all(uncensored)) {
      c(left = 0, right = 0, interval = 0)
    } else {
      c(
        left = sum(is.infinite(ai)), right = sum(is.infinite(bi)),
        interval = sum(is.finite(ai) & is.finite(bi) & (ai < bi))
      )
    }

    A[notna] <- ai
    B[notna] <- bi
    y <- cbind(a = A, b = B)
    dimnames(y) <- list(NULL, NULL)

    return(structure(y,
      class   = "Ocens",
      levels  = L / mul,
      upper   = if (any(L != R)) R / mul,
      freq    = freq,
      median  = ymed / mul,
      ranges  = list(y = yrange / mul, u = urange / mul, c = crange / mul),
      label   = ylabel,
      units   = uni,
      Ncens1  = ncen,
      Ncens2  = ncen2,
      npsurv  = np
    ))
  }

  # Categorical variables as integers
  uncensored <- a == b
  if (!any(uncensored)) stop("no uncensored observations")
  if (any(b < a)) stop("some values of b are less than corresponding a values")
  freq <- tabulate(a[uncensored], nbins = length(at$levels))
  mid <- quantile(a[uncensored], probs = .5, type = 1L)
  A[notna] <- a
  B[notna] <- b
  # Categorical variables cannot be infinite, so no left or rt censoring
  ncen <- c(left = 0, right = 0, interval = sum(!uncensored))
  structure(cbind(a = A, b = B),
    class = "Ocens", levels = at$levels, freq = freq, mid = mid,
    label = ylabel,
    units = uni,
    Ncens1 = ncen, Ncens2 = ncen
  )
}

# Claude Sonnet 5 High 2026-07-10          46 lines
##' Maximal Intersections of Censoring Intervals
##'
##' Computes the Turnbull maximal intersections ("innermost intervals") of a
##' set of closed intervals `[a, b]` whose finite endpoints lie on an exact
##' integer grid, along with the mapping of each observation to the
##' intersections it contains.  A maximal intersection is an interval
##' `[L, R]` where `L` is some observation's left endpoint, `R` is some
##' observation's right endpoint, and no other endpoint lies strictly
##' inside.  These are found by a single sweep over the sorted endpoints:
##' with left endpoints ordered before right endpoints at tied values, the
##' maximal intersections are exactly the places where a left endpoint is
##' immediately followed by a right endpoint.  The nonparametric MLE of the
##' distribution can only put probability mass inside maximal intersections
##' (Turnbull 1976), so these determine the estimable intercepts in `orm`.
##' The intersections contained in a given observation's interval always
##' form a consecutive run, computed here with [findInterval()].
##'
##' @param a vector of closed lower interval endpoints, `-Inf` for left-censored values
##' @param b vector of closed upper interval endpoints, `Inf` for right-censored values
##' @return a list with elements `L` and `R` (lower and upper endpoints of
##'   the maximal intersections, in increasing order) and integer vectors
##'   `ai` and `bi` giving for each observation the first and last maximal
##'   intersection contained in `[a, b]`
##' @author Frank Harrell
##' @noRd
Ocens_maxintersect <- function(a, b) {
  n    <- length(a)
  ep   <- c(a, b)
  type <- rep(0:1, each = n)   # 0 = left endpoint, 1 = right endpoint
  o    <- order(ep, type)      # left endpoints sort before right ones at ties
  ep   <- ep[o]
  type <- type[o]
  i    <- which(type[-(2 * n)] == 0 & type[-1] == 1)
  L    <- ep[i]
  R    <- ep[i + 1]
  # First intersection with L >= a, last with R <= b; exact because all
  # finite values are integers.  The a - 1 device implements the strict
  # inequality L < a; it is invalid at a = -Inf (left censoring), where the
  # first intersection always qualifies.
  ai <- findInterval(a - 1, L) + 1L
  ai[is.infinite(a)] <- 1L
  bi <- findInterval(b, R)
  if (any(bi < ai))
    stop("program logic error in Ocens_maxintersect: an observation contains no maximal intersection")
  list(L = L, R = R, ai = ai, bi = bi)
}

# Claude Sonnet 5 High 2026-07-10          82 lines
##' Nonparametric MLE for Censored Data by Self-Consistency
##'
##' Computes the nonparametric maximum likelihood estimate of the
##' distribution of an arbitrarily censored variable using Turnbull's
##' self-consistency (EM) algorithm, operating on the output of
##' `Ocens_maxintersect`.  Because the maximal intersections contained in an
##' observation's interval form a consecutive run `ai : bi`, each iteration
##' is computed in O(n + m) time using cumulative sums, without forming the
##' n x m clique matrix.  At the solution, the self-consistency (Kuhn-Tucker)
##' condition states that the gradient element for intersection `j`,
##' `d_j = sum_i wt_i [ai_i <= j <= bi_i] / P(clique of i)`, equals
##' `sum(wt)` wherever the probability is positive and is smaller wherever
##' the probability must be zero.  Intersections failing to attain the
##' maximum carry no estimable probability and are flagged in `support` for
##' consolidation by `Ocens_consolidate`.
##'
##' @param ai integer vector: first maximal intersection contained in each observation's interval
##' @param bi integer vector: last maximal intersection contained in each observation's interval
##' @param wt vector of case weights
##' @param m number of maximal intersections
##' @param maxiter maximum number of self-consistency iterations
##' @param tol convergence tolerance for the maximum absolute change in probabilities
##' @param ktol relative tolerance for the Kuhn-Tucker support test
##' @param mtol probability threshold below which an intersection is a
##'   candidate for being dropped.  An intersection is flagged unsupported
##'   only when its probability is `< mtol` \emph{and} its KKT gradient
##'   falls short of the maximum; this two-condition rule is deliberately
##'   conservative, retaining borderline intersections that the EM
##'   tolerance cannot cleanly classify from the gradient alone.
##' @return a list with elements `p` (probabilities of the maximal
##'   intersections), `surv` (P(Y >= intersection j), so `surv[1] = 1`),
##'   `support` (logical: intersections carrying estimable probability),
##'   `iter`, and `converged`
##' @author Frank Harrell
##' @noRd
Ocens_npmle <- function(ai, bi, wt = rep(1e0, length(ai)), m,
                        maxiter = 10000L, tol = 1e-8, ktol = 1e-6,
                        mtol = 1e-5) {
  W <- sum(wt)
  # Weighted tabulate: sum of w within values of idx, over 1 : nb
  wtab <- function(idx, w, nb) {
    z <- numeric(nb)
    s <- rowsum(w, idx)
    z[as.integer(rownames(s))] <- s
    z
  }
  # Gradient d_j = sum_i r_i [ai_i <= j <= bi_i] by a difference array
  grad <- function(r)
    cumsum(wtab(ai, r, m + 1L) - wtab(bi + 1L, r, m + 1L))[seq_len(m)]

  p    <- rep(1e0 / m, m)
  iter <- 0L
  conv <- FALSE
  repeat {
    iter <- iter + 1L
    cp   <- c(0e0, cumsum(p))
    den  <- cp[bi + 1L] - cp[ai]   # probability of observation's clique
    if (any(den <= 0e0))
      stop("program logic error in Ocens_npmle: zero clique probability")
    d    <- grad(wt / den)
    pnew <- p * d / W
    conv <- max(abs(pnew - p)) < tol
    p    <- pnew
    if (conv || iter >= maxiter) break
  }
  # Kuhn-Tucker gradient at the solution determines the support.
  # An intersection is unsupported only if BOTH its probability is
  # negligible (< mtol) AND its KKT gradient falls short of W.  The
  # gradient test alone is unreliable near convergence: the EM tolerance
  # only bounds |d_j / W - 1| by tol / p_j, so a genuinely supported
  # intersection with small p_j can spuriously fail d_j >= W(1 - ktol).
  # Keeping such borderline intersections is conservative -- at worst an
  # unneeded near-degenerate intercept is retained, which orm.fit's
  # step-halving safely absorbs; wrongly dropping a supported one would
  # bias the survival estimate.
  cp  <- c(0e0, cumsum(p))
  den <- cp[bi + 1L] - cp[ai]
  d   <- grad(wt / den)
  list(p = p, surv = rev(cumsum(rev(p))),
       support = (p >= mtol) | (d >= W * (1e0 - ktol)),
       iter = iter, converged = conv)
}

# Claude Sonnet 5 High 2026-07-10          40 lines
##' Consolidate Maximal Intersections Carrying No Estimable Probability
##'
##' Removes maximal intersections flagged as unsupported by `Ocens_npmle`,
##' remapping observations and widening the interval labels, without
##' changing any data.  Observation mappings are snapped inward to supported
##' intersections (`ai` up, `bi` down), which is exact because unsupported
##' intersections carry no probability.  For labeling, each unsupported
##' intersection is absorbed into the following supported one (trailing
##' unsupported intersections into the last), widening `[L, R]` to reflect
##' the region over which the estimate is indeterminate; this matches the
##' grouping formerly produced by consolidating consecutive tied survival
##' estimates.
##'
##' @param L,R endpoints of the maximal intersections
##' @param ai,bi observation mappings from `Ocens_maxintersect`
##' @param npm result of `Ocens_npmle`
##' @return a list with consolidated `L`, `R`, remapped `ai`, `bi`, cluster
##'   probabilities `p`, and the consolidated survival curve `surv`
##' @author Frank Harrell
##' @noRd
Ocens_consolidate <- function(L, R, ai, bi, npm) {
  support <- npm$support
  cl <- cumsum(support)                # cluster of last supported intersection <= j
  up <- ifelse(support, cl, cl + 1L)   # cluster of first supported intersection >= j
  m2 <- max(cl)
  if (m2 == 0L)
    stop("program logic error in Ocens_consolidate: no supported intersections")
  # Snap observation mappings inward to supported intersections
  ai2 <- up[ai]
  bi2 <- cl[bi]
  if (any(ai2 > bi2 | ai2 > m2 | bi2 < 1L))
    stop("program logic error in Ocens_consolidate: an observation is incompatible with the estimable intersections")
  # Widen interval labels over absorbed unsupported intersections
  labcl <- factor(pmin(pmax(up, 1L), m2), levels = seq_len(m2))
  Ln <- as.vector(tapply(L, labcl, min))
  Rn <- as.vector(tapply(R, labcl, max))
  pc <- as.vector(tapply(npm$p, labcl, sum))
  list(L = Ln, R = Rn, ai = ai2, bi = bi2, p = pc,
       surv = rev(cumsum(rev(pc))))
}

# Claude Sonnet 5 High 2026-07-10          99 lines
##' Legacy Data-Changing Interval Consolidation
##'
##' Deprecated implementation of interval consolidation by iteratively
##' widening raw data intervals (`cons='data'`), retained for backward
##' compatibility.  Requires the `icenReg` package.  Inputs `a` and `b` are
##' on the exact integer grid with open censoring endpoints already
##' converted to closed ones.
##'
##' @param a,b closed interval endpoints on the integer grid
##' @param mul integer grid multiplier, for verbose output and returned times
##' @param maxit maximum number of consolidation iterations
##' @param nponly if `TRUE` return the survival curve estimate only
##' @param verbose print information messages
##' @return for `nponly=TRUE` a list with `time` and `surv`; otherwise a
##'   list with possibly modified `a`, `b`, the Turnbull intervals `L`, `R`,
##'   observation mappings `ai`, `bi`, and the survival estimates `np`
##' @author Frank Harrell
##' @noRd
Ocens2ordDataCons <- function(a, b, mul, maxit = 10, nponly = FALSE,
                              verbose = FALSE) {
  if (!requireNamespace("icenReg", quietly = TRUE)) {
    stop("The icenReg package must be installed to use Ocens2ord with cons='data'")
  }
  fmi <- utils::getFromNamespace("findMaximalIntersections", "icenReg")

  iter <- 0
  mto <- function(x) diff(range(x)) > 0 # more than one distinct value
  repeat {
    iter <- iter + 1
    if (iter > maxit) stop("exceeded maxit=", maxit, " iterations for pooling intervals")
    it <- fmi(as.double(a), as.double(b))
    L <- it$mi_l
    R <- it$mi_r
    # The integer Y matrix produced by Ocens is the mappings of observations to the (L, R) Turnbull intervals
    # Indexes created by fmi start with 0, we bump them to 1
    ai <- it$l_inds + 1L
    bi <- it$r_inds + 1L
    if (verbose > 1) prn(cbind(a, b, La = L[ai], Lb = L[bi], Ra = R[ai], Rb = R[bi]) / mul)
    dicen <- data.frame(a = a, b = b, grp = rep(1, length(a)))
    g <- icenReg::ic_np(cbind(a, b) ~ grp, data = dicen, B = c(1, 1)) # bug prevents usage of matrix without formula
    # Note: icenReg::getSCurves() will not run
    s <- g$scurves[[1]]$S_curves$baseline
    k <- length(L) - 1L
    np <- list(time = L / mul, surv = s[1:(k + 1)])
    if (nponly) {
      return(np)
    }
    if (length(np$time) != length(np$surv)) {
      warning("vector length mismatch in icenReg::ic_np result from npsurv=TRUE")
    }
    if (verbose > 1) print(cbind(t = np$time, "S(t)" = np$surv))
    s <- 1e-7 * round(np$surv * 1e7)
    su <- unique(s[duplicated(s)])
    if (!length(su)) break

    # Some consecutive intervals had the same information
    # For these code all the raw data as [lower, upper] where lower is the
    # minimum lower limit in the overlapping intervals, upper is the maximum upper limit
    # Compute distinct values of s that have > 1 Turnbull interval with that s value
    # Find original data corresponding to each su
    # Lookup s for each row of data
    S <- s[ai]

    for (ans in su) {
      j <- which(S == ans)
      if (!length(j)) stop("program logic error in Ocens2ordDataCons")
      if (verbose) {
        cat("\nIntervals consolidated to give unique contributions to survival estimates and likelihood\n\nBefore:\n\n")
        print(cbind(a = a[j] / mul, b = b[j] / mul))
      }
      aj <- a[j]
      bj <- b[j]
      l <- is.infinite(aj)
      r <- is.infinite(bj)
      ic <- (!l) & (!r) & (bj > aj)
      # Try only one remedy per group, using else if ...
      if (any(r) && !any(l)) {
        a[j[!l]] <- min(a[j[!l]])
      } else if (any(r)) {
        a[j[r]] <- min(aj)
      } else if (any(l) && all(bj[l] == max(bj[!r]))) {
        b[j[l]] <- min(bj[!r])
      } else if (any(l)) {
        b[r[l]] <- max(bj)
      } else if ((sum(ic) > 1) && (mto(a[j[ic]]) || mto(b[j[ic]]))) {
        a[j[ic]] <- min(aj[!l])
        b[j[ic]] <- max(bj[!r])
      } else if (any(ic)) {
        a[j] <- min(a[j])
        b[j] <- max(b[j])
      }
      if (verbose) {
        cat("\nAfter:\n\n")
        print(cbind(a = a[j] / mul, b = b[j] / mul))
      }
    }
  }
  list(a = a, b = b, L = L, R = R, ai = ai, bi = bi, np = np)
}

##' Convert `Ocens` Object to Data Frame to Facilitate Subset
##'
##' Converts an `Ocens` object to a data frame so that subsetting will preserve all needed attributes
##' @param x an `Ocens` object
##' @param row.names optional vector of row names
##' @param optional set to `TRUE` if needed
##' @param ... ignored
##' @return data frame containing a 2-column integer matrix with attributes
##' @author Frank Harrell
##' @method as.data.frame Ocens
##' @export
as.data.frame.Ocens <- function(x, row.names = NULL, optional = FALSE, ...) {
  deb <- Fdebug("rmsdebug")
  nrows <- NROW(x)
  deb(nrows)
  row.names <- if (optional) character(nrows) else as.character(1:nrows)
  value <- list(x)
  deb(dim(value[[1]]))
  if (!optional) names(value) <- deparse(substitute(x))[[1]]
  deb(dim(value[[1]]))
  structure(value, row.names = row.names, class = "data.frame")
}

##' Subset Method for `Ocens` Objects
##'
##' Subsets an `Ocens` object, preserving its special attributes.  Attributes are not updated.  In the future such updating should be implemented.
##' @title Ocens
##' @param x an `Ocens` object
##' @param ... the usual rows and columns specifiers
##' @param drop set to `FALSE` to not drop unneeded dimensions
##' @return new `Ocens` object or by default an unclassed vector if only one column of `x` is being kept
##' @author Frank Harrell
##' @md
##' @method [ Ocens
##' @export
"[.Ocens" <- function(x, ..., drop) {
  d <- dim(x)
  at <- attributes(x)
  n <- intersect(names(at), c("name", "label", "units", "levels"))
  x <- unclass(x)
  x <- x[..., drop = FALSE]
  if (missing(drop)) drop <- NCOL(x) == 1
  if (drop) x <- drop(x)
  attributes(x) <- c(attributes(x), at[n])
  if (NCOL(x) == 2) class(x) <- "Ocens"
  x
}

##' is.na Method for Ocens Objects
##'
##' @param x an object created by `Ocens`
##'
##' @returns a logical vector whose length is the number of rows in `x`, with `TRUE` designating observations having one or both columns of `x` equal to `NA`
##' @method is.na Ocens
##' @export
##'
##' @md
##'
##' @examples
##' Y <- Ocens(c(1, 2, NA, 4))
##' Y
##' is.na(Y)
is.na.Ocens <- function(x) as.vector(rowSums(is.na(unclass(x))) > 0)

#' Ocens2Surv
#'
#' Converts an `Ocens` object to the simplest `Surv` object that works for the types of censoring that are present in the data.
#'
#' @param Y an `Ocens` object
#'
#' @returns a `Surv` object
#' @export
#' @md
#'
#' @examples
#' Y <- Ocens(1:3, c(1, Inf, 3))
#' Ocens2Surv(Y)
Ocens2Surv <- function(Y) {
  y <- Y[, 1]
  y2 <- Y[, 2]

  su <- survival::Surv
  if (all(y == y2)) {
    return(su(y))
  } # no censoring
  i <- which(is.finite(y) & is.finite(y2))
  w <- 1 * any(is.infinite(y)) + 2 * any(is.infinite(y2)) + 4 * any(y[i] != y2[i])
  if (w == 1) {
    su(y2, event = y == y2, type = "left")
  } else if (w == 2) {
    su(y, event = y == y2, type = "right")
  } else if (w == 4) {
    su(y, event = rep(3, length(y)), time2 = y2, type = "interval")
  } else {
    su(y, time2 = y2, type = "interval2")
  }
}

##' print Method for Ocens Objects
##'
##' @param x an object created by `Ocens`
##' @param ivalues set to `TRUE` to print integer codes instead of character levels when original data were factors or character variables
##' @param digits number of digits to the right of the decimal place used in rounding original levels when `ivalues=FALSE`
##' @param ... ignored
##' @returns nothing
##' @method print Ocens
##' @export
##' @md
##'
##' @examples
##' Y <- Ocens(1:3, c(1, Inf, 3))
##' Y
##' print(Y, ivalues=TRUE)  # doesn't change anything since were numeric
print.Ocens <- function(x, ivalues = FALSE, digits = 5, ...) {
  y <- matrix(NA, nrow(x), ncol(x)) # to drop attributes of x
  y[] <- x
  a <- y[, 1]
  b <- y[, 2]
  nna <- !is.na(a + b)
  ia <- is.finite(a) & nna
  ib <- is.finite(b) & nna
  ifa <- is.infinite(a) & nna
  ifb <- is.infinite(b) & nna
  lev <- attr(x, "levels")
  if (!ivalues && length(lev)) {
    a[ia] <- lev[a[ia]]
    b[ib] <- lev[b[ib]]
  }
  if (!length(lev)) {
    a <- round(a, digits)
    b <- round(b, digits)
  }
  intcens <- ia & ib & (b > a)
  a <- format(a)
  b <- format(b)
  z <- a
  z[ifa] <- paste0(b[ifa], "-")
  z[ifb] <- paste0(a[ifb], "+")
  z[intcens] <- paste0("[", a[intcens], ",", b[intcens], "]")
  print(z, quote = FALSE)
  invisible()
}

extractCodedOcens <- function(x, what = 1, ivalues = FALSE, intcens = c("mid", "low")) {
  intcens <- match.arg(intcens)
  lev <- attr(x, "levels")
  n <- nrow(x)
  a <- b <- integer(n)
  a[] <- x[, 1] # gets rid of attributes
  b[] <- x[, 2]
  ia <- is.infinite(a)
  ib <- is.infinite(b)
  if (ivalues) {
    a <- a - 1
    b <- b - 1
  } else if (length(lev)) {
    a[!ia] <- lev[a[!ia]]
    b[!ib] <- lev[b[!ib]]
  }
  if (what == 2) {
    return(cbind(a = a, b = b))
  }

  ctype <- integer(n)
  ctype[ia] <- 1 # left censoring
  ctype[ib] <- 2 # right
  ctype[ctype == 0 & (a < b)] <- 3 # interval
  l <- ctype == 1
  r <- ctype == 2
  i <- ctype == 3

  y <- numeric(n)
  y[l] <- b[l]
  y[r] <- a[r]
  y[i] <- if (intcens == "mid") 0.5 * (a + b)[i] else a[i]
  if (what == 1) {
    return(y)
  }
  list(a = a, b = b, y = y, ctype = ctype)
}

# Function determining TRUE/FALSE whether Y is known to be >= j
# a and b are results of Ocens2ord
# Returns NA if censoring prevents determining this
# Left censoring
#   Y >= j can be determined if b <= j
#   FALSE in this case
# Right censoring
#   Y >= j can be determined if a >= j
#   TRUE in this case
# Interval censoring
#   Y >= j can be determined if a >= j | b < j
#   TRUE if a >= j, FALSE if b < j
# Assumes that a and b run from 0 to k
geqOcens <- function(a, b, ctype, j) {
  z <- rep(NA, length(a))
  u <- ctype == 0
  l <- ctype == 1
  r <- ctype == 2
  i <- ctype == 3
  z[u] <- a[u] >= j
  z[l & b <= j] <- FALSE
  z[r & a >= j] <- TRUE
  z[i & a >= j] <- TRUE
  z[i & b < j] <- FALSE
  z
}

# g <- function(a, b) {
#   s <- Ocens2Surv(cbind(a, b))
#   print(s)
#   km.quick(s, interval='>=')
# }
# g(1:3, 1:3)
# g(c(-Inf, 2, 3), c(2.5, 2, 3))
# g(1:3, c(1, 2, Inf))
# g(c(1, 4, 7), c(2, 4, 8))
# g(c(-Inf, 2,4, 6), c(3, 3, 4, Inf))
# a <- c(-Inf, 2, 1, 4, 3)
# b <- c(   3, 3, 1, 5, 3)
# g(a, b)
