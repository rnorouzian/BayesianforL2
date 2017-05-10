#' Bayesian Process for Cohen's d in t-test Designs
#'
#' Examplifies the Bayesian process to estimate Cohen's d in t-test Designs.
#'
#' @param   N1 Group 1 Sample size.
#' @param   N2 Group 2 Sample size (if a two-sample design).
#' @param   d Observed estimate of Cohen's d.
#'
#' @return   Graphically details the process of bayesian estimation of effect size by focusing on
#'          the role of (1) prior, and (2) Likelihood Function to produce the posterior
#'          for Cohen's d.
#'
#'
#' @author   Reza Norouzian <rnorouzian@gmail.com>
#' @export
#'
#' @examples
#'
#'  # Bayesian_Process(N1 = 50, N2 = 50, d = .4)
#'
#'  # Bayesian_Process(N1 = 50, d = .4)



Bayesian_Process = function(d, N1, N2 = NULL) {

  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))

  options(warn = -1)
  m <- matrix( c( 1,0,  1,3,  2,3,  2,0  ), nrow=2, ncol=4 ); layout(m)

  par(mar = c(5.1, 6.1, 4.1, 2.1), family = 'serif', mgp = c(2.5, 1, 0) )

  curve(dnorm(x), -3, 3, n = 1e3, lwd = 3,
        ylab = NA, main = bquote(bolditalic('Prior'~(delta))),
        xlab = bquote(bolditalic(delta)),
        font = 2, cex = 1.2, cex.sub = 1.2, cex.lab = 1.2, axes = F)

  axis(1, cex.axis = 1.2, font = 2)

  low.exterme = par('usr')[3]
  prior.peak = dnorm(0)

  segments(0, low.exterme, 0, prior.peak, lty = 3, lwd = 2, col = 'gray40')

  text(0, prior.peak/2.3, bquote(bold("Neutral Position for"~(delta))), srt = 90,
       pos = 3, cex = 1.3, font = 2 )


  ## Likelihood Function:

  ## Observed data:
  efN = ifelse(is.null(N2), N1, N1*N2/(N1+N2))   # Effective Sample Size
  df  = ifelse(is.null(N2), N1 - 1, N1 + N2 - 2) # Degrees of freedom
  t = d*sqrt(efN)                                # Observed t-value


  Likelihood = curve( dt(t, df, ncp = x*sqrt(efN) ),
                      col = 'red4', lwd = 3, from = -3, to = 3, n = 1e3,
                      ylab = NA, main = bquote(bolditalic('Likelihood'~(delta))),
                      xlab = bquote(bolditalic(delta)),
                      font = 2, cex.lab = 1.2, axes = F)$y


  Like.mode = optimize(function(x) dt(t, df, ncp = x*sqrt(efN) ), interval = c(-3, 3), maximum = T, tol = 1e-10)[[1]]


  axis(1, at = round(seq(-3, 3, len = 7), 2 ), cex.axis = 1.2, font = 2 )


  segments(Like.mode, par('usr')[3], Like.mode, max(Likelihood), lty = 3, lwd = 2, col = 'gray40')

  text(Like.mode, max(Likelihood)/3.4, bquote(bold("Observed"~ (delta))), srt = 90, pos = 3,
       cex = 1.3, font = 2 )


  ## Computing the posterior distribution of effect size

  posterior <- function(t, N1, N2=NULL, delta) {

    efN = ifelse(is.null(N2), N1, N1*N2/(N1+N2))
    df  = ifelse(is.null(N2), N1 - 1, N1 + N2 - 2)

    #prior and likelihood
    prior <- function(delta) dnorm(delta, 0, 1) # dunif(delta, -2, 2) # {(delta^2)*(1/sqrt(2*pi))*exp(-(delta^2)/2)}
    likelihood <- function(delta) dt(t, df, ncp = delta*sqrt(efN))

    #marginal likelihood
    marginal <- integrate(function(x) prior(x)*likelihood(x), -Inf, Inf)[[1]]

    #posterior
    post <- function(x) prior(x)*likelihood(x) / marginal
    return(post(delta))
  }



  post = curve( posterior(t, N1, N2, x), -3, 3, lwd = 3, col = 'blue',
                n = 1e3, ylab = NA, main = bquote(bolditalic('Posterior'~(delta))),
                xlab = bquote(bolditalic(delta)),
                font = 2, cex.lab = 1.2, axes = F)$y

  axis(1, at = seq(-3, 3, by =1 ), cex.axis = 1.2, font = 2, xpd = T)

  low.exterme = par('usr')[3]


  #posterior mode
  post.mode = optimize(function(delta) posterior(t, N1, N2, delta),interval=c(-3,3),maximum = T, tol = 1e-10)[[1]]

  post.peak = posterior(t, N1, N2, post.mode)

  segments(post.mode, low.exterme, post.mode, post.peak, lwd = 2, lty = 3, col = 'gray40')

  text(post.mode, post.peak/4.3, bquote(bold("Posterior Mode"~ (delta))), srt = 90, pos = 3,
       cex = 1.3, font = 2)

  mtext(side = 3, " \u00D7 ", line = 15, col= 'navy', cex = 6, font = 2,xpd = T)


  arrows(0, 1.8*post.peak, 0, 1.25*post.peak, xpd = NA,    ## NA for xpd is a MUST
         code = 2, lwd = 10, length = .2, angle = 20, lend = 1, col = 'red')

}
