#' Bayesian Estimation of "Partial Eta-Squared" in Fixed-Effects AN(C)OVA Designs
#'
#' Provides complete Posterior as well as its various summary statistics for
#' "Partial Eta-Squared" in Fixed-Effects AN(C)OVA Designs.
#'
#'
#' @param  f Observed F-value from a main or an interaction effect.
#' @param  N Total design sample size.
#' @param  df1 Degress of freedom for the effect in question aka "Numerator".
#' @param  df2 Degrees of freedom "Error" aka "Denominator", "Within" or "Residual".
#' @param  alpha shape1 parameter of a "Beta Prior".
#' @param  beta shape2 parameter of a "Beta Prior".
#'
#' @return  Provides graphical as well as full textual description of posterior distribution of
#'          "Partial Eta-Squared".
#'
#' @seealso Use function \code{\link{Beta_ID}} provided in this package to find an appropriate "Beta Prior" in your specific domain
#' of resaerch.
#'
#' @seealso Norouzian R., & L. Plonsky (2017) available at: \url{http://journals.sagepub.com/doi/abs/10.1177/0267658316684904}
#'
#' @seealso This supplemetary document also briefly explains when Partial Eta Squared should be used: \url{https://osf.io/7acye/}
#'
#' @references Norouzian, R. & Plonsky, L. (2017). Eta- and partial eta-squared in L2 research:
#'  A cautionary review and guide to more appropriate usage. Second Language Research, pp.1-15.
#'  doi: 10.1177/0267658316684904
#'
#' @author   Reza Norouzian <rnorouzian@gmail.com>
#' @export
#'
#' @examples
#'
#' # P_eta_sq(f = 10, N = 60, df1 = 2, df2 = 57, alpha = 2.55, beta = 3.42)
#'
#' # Note that, the "P_eta_sq(...)" when applied to ONE-WAY designs, is basically
#' # the same as "ETA-SQUARED" but when "P_eta_sq(...)" is applied to factorial
#' # designs it produces "PARTIAL ETA-SQUARED".





P_eta_sq = function(f, N, df1, df2, alpha, beta){


  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))


  posterior <- function(f, N, df1, df2, petasq, alpha, beta) {


    #prior and likelihood
    prior <- function(petasq) dbeta(petasq, alpha, beta) ## try very-vague-priors
    likelihood <- function(petasq) df(f, df1, df2, ncp = (petasq * N) / (1 - petasq) )

    #marginal likelihood
    marginal <- integrate(function(x) prior(x)*likelihood(x), lower = 0, upper = 1)[[1]]

    #posterior
    Density <- function(x) prior(x)*likelihood(x) / marginal
    return(Density(petasq) )

  }


  options(warn = -1)


  ## Decimal display controller:
  decimal <- function(x, k){

    if(is.character(x)){
      return(x)
    }

    format(round(x, k), nsmall = k, scientific =
             ifelse(x >= 1e+05 || x <= -1e+05 || x <= 1e-05 & x >= -1e-05, T, F) )
  }


  ## Possible range that P.eta.Sq. can take:
  petasq  <- seq(0, 1, len = 1e+4)


  ## 3 different definitions of "ncp" for F-dist. in the Math Stats Literature (JOHNSON, KOTZ[1994, VOL.2, P. 495]):
  NCPumvue = ((df1*(df2-2)*f) / df2 ) - df1

  NCPminMSE = ((df1*(df2-4)*f) / df2 ) - ((df1*(df2-4)) / (df2-2))

  NCPcohen = (f*(df1/df2))*N


  ## Observed P.eta.Sq based on different "ncp" definitions:
  obs.petasq.Cohen = NCPcohen / (NCPcohen + N)      # Or directly could be: (f*df1) / ((f*df1) + df2)

  obs.petasq.umvue = NCPumvue / (NCPumvue + N)      # Higher & Less bias accuracy confirmed

  obs.petasq.minMSE = NCPminMSE / (NCPminMSE + N)



  ## Plot of the posterior:
  post <- posterior(f, N, df1, df2, petasq, alpha, beta) # posterior densities

  post.peak = max(post, na.rm = T)

  ylim = c(0, post.peak + .15*post.peak)

  par(family = "serif")
  plot(petasq, post, ty="l", xlab = bquote(bold("Partial Eta-Squared " (eta[p]^2))),
       ylab = "Density", lwd = 3, font.axis = 2, font.lab = 2,
       las = 1, cex.lab = 1.3, cex.sub = 1.3, xaxs = "i", bty = 'n', cex.axis = 1.2
       , xaxt = "n", ylim = ylim )

  axis(1, at = decimal(seq(0, 1, len = 7), 2 ), font.axis = 2, cex.axis = 1.3 )


  legend("topright", c("Beta Prior", "Posterior"), lty = c(3, 1),
         lwd = 3, col = c("green4", 1), text.font = 2, cex = 1.2,
         text.col = c("green4", 1))


  ## Plot of prior:
  Prior = curve( dbeta(x, alpha, beta),
                 add = T, col = "green4", lty = 3, lwd = 3, from = 0, to = 1)

  prior.mode = (alpha - 1) / ( alpha + beta - 2 )

  prior.peak = dbeta(prior.mode, alpha, beta)


  #posterior mode
  post.mode = optimize(function(petasq) posterior(f,N,df1,df2,petasq,alpha,beta),interval=c(0,1),maximum = T, tol=1e-10)[[1]]

  #posterior mean
  post.mean = integrate(function(petasq) petasq* posterior(f,N,df1,df2,petasq,alpha,beta), 0, 1)[[1]]

  #posterior variance
  post.var = integrate(function(petasq) petasq^2*posterior(f,N,df1,df2,petasq,alpha,beta), 0, 1)[[1]] - post.mean^2

  #posterior sd
  post.sd = sqrt(post.var)

  #post median
  #post.median = petasq[which(cumsum(post) > 0.5)[1]] ## gives wrong result

  low.extreme = par('usr')[3]


  ## line segments indicating the mode:

  segments(post.mode, low.extreme, post.mode, post.peak, lty = 3, col = "gray20", xpd = T)


  ## Texts on the line segments:
  text(prior.mode, prior.peak, bquote(bold("Prior")), col = "green4", pos = 3, cex = 1.2)

  text(post.mode, post.peak, bquote(bold("Posterior")), pos = 3, cex = 1.2, xpd = T)

  text(post.mode, post.peak/ 1.5, bquote(bold("Posterior Mode")), srt = 90, col = "gray40", pos = 3, cex = 1.2)

  mtext(side = 3, bquote(bold("Bayesian Estimation of Partial Eta-Squared " (eta[p]^2) ~ "in L2 research")), cex = 1.4, line =  1)

  par(family = "sans")



  # Determine the 95% Posterior High Density Interval:
  HDI <- function(x, y, target=.95) {
    dx <- diff(x)
    areas <- dx * .5 * (head(y,-1) + tail(y, -1))
    peak <- which.max(areas)
    range <- c(peak, peak)
    found <- areas[peak]
    while(found < target) {
      if(areas[range[1]-1] > areas[range[2]+1]) {
        range[1] <- range[1]-1
        found <- found + areas[range[1]-1]
      } else {
        range[2] <- range[2]+1
        found <- found + areas[range[2]+1]
      }
    }
    val<-x[range]
    attr(val, "indexes")<-range
    attr(val, "area")<-found
    return(val)
  }


  CI = HDI(x = petasq, y = post)

  CI.Lower = CI[1]

  CI.Upper = CI[2]

  segments(CI.Lower, 0, CI.Upper, 0, lend = 1, lwd = 30, col = 'red4')

  points(post.mode, 0, pch = 20, col = 0, cex = 3)


  ## List important results:
  vec <- c(prior.mode, post.mode, CI.Lower, CI.Upper,obs.petasq.umvue, obs.petasq.minMSE,
           obs.petasq.Cohen, post.mean, post.sd)

  names(vec) <- c('prior.mode', 'post.mode', 'Lower CI', 'Upper CI' , 'obs.petasq.umvue',
                  'obs.petasq.minMSE', 'obs.petasq.Cohen', 'post.mean', 'post.sd' )

  lapply(vec, function(x) decimal(x, 7) )

}
