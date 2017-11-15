#' Bayesian Equivalence Testing for Cohen's d in t-test Designs
#'
#' Provides complete results of Bayesian Equivalence Testing with
#' an "additional function" showing the result of changing the equivalence
#' bounds to any different sizes.
#'
#' @param  N1 Group 1 sample size.
#' @param  N2 Group 2 sample size (if a two-sample design).
#' @param  t  Observed t-value.
#' @param  wid Width or aka Scale of a Cauchy prior.
#' @param  dL Lower limit of equivalence bound for d.
#' @param  dU Upper limit of equivalence bound for d.
#'
#' @return  Provides complete results of Bayesian Equivalence Testing using High Desity Region (HDI)
#'          with an "additional function" showing the result of changing the equivalence
#'          bounds to any different sizes.
#'
#' @details Equivalnce Testing simply presents another alternative to absolute hypotheis testing (i.e., testing whether an effect is excatly 0).
#'          Specifically, we allow a bound of values which in our view are practically equivalent to zero to be our null hypothesis (i.e., a composite Null).
#'          However, most Bayeian point decisions are better decided upon using a "Decision-Theoritic" approch (see, Berger, 1985).
#'
#' @author  Reza Norouzian <rnorouzian@gmail.com>
#' @export
#'
#' @examples
#' # Suppose a researcher obtains a t-value of 1, from 2 groups each with 20 participants.
#' # The researcher picks a "wi", wide, Cauchy prior. The researcher uses
#' # dL of -0.1 and dU of 0.1 to see the result of equivalene testing:
#'
#'
#'  Eq_Testing(t = 1, N1 = 20, N2 = 20, dL = -.1, dU = .1 )




Eq_Testing = function(t, N1, N2 = NULL, wid = "wi", dL = -.1, dU = .1){


  
  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))


  options(warn = -1)


  m = matrix( c(1, 2), nrow = 2, ncol = 1 ); layout(m)


  if(dL >= dU){stop("Your Upper value must be larger than your Lower value")}


  if(abs(dL) != abs(dU)) { message("\n\tYou have an \"Unequal Equivalence Bound\", thus we can't provide an extra\n\t function showing the effect of choosing various Unequal bounds.") }



  ## decimal display controller

  decimal <- function(x, k){

    if(typeof(x) == "character"){
      return(x)
    }

    as.numeric(format(round(x, k), nsmall = k, scientific = 
           ifelse(x >= 1e+05 || x <= -1e+05 || x <= 1e-05 & x >= -1e-05, TRUE, FALSE) ))
  }



  ## Some taking care of type of t-tests

  wid = c(m = 1/2, wi = sqrt(2)/2, wii = 1)[[wid]]

  rscale = wid



  ## Observed Effect Size:

  effN = ifelse(is.null(N2), N1, N1*N2/(N1+N2))
  observed.d = t / sqrt(effN)



  ## @@@ Start of "HDR" function:

  HDIofICDF = function( ICDFname , credMass=0.95 , tol=1e-8 , ... ) {

    incredMass = 1.0 - credMass
    intervalWidth = function( lowTailPr , ICDFname , credMass , ... ) {
      ICDFname( credMass + lowTailPr , ... ) - ICDFname( lowTailPr , ... )
    }
    optInfo = optimize( intervalWidth , c( 0 , incredMass ) , ICDFname=ICDFname ,
                        credMass=credMass , tol=tol , ... )
    HDIlowTailPr = optInfo$minimum
    return( c( ICDFname( HDIlowTailPr , ... ) ,
               ICDFname( credMass + HDIlowTailPr , ... ) ) )
  }

  ## @@@ End of "HDR" function



  ## ^^^ Start of UTILITY FUNCTIONS: GIVEs pdf, CDF, INVERSE.CDF OF POSTERIOR

  posterior = Vectorize(function(delta, t, N1, N2 = NULL, rscale = wid, log = FALSE){
    effN = ifelse(is.null(N2), N1, N1*N2/(N1+N2))
    df = ifelse(is.null(N2), N1 - 1, N1 + N2 - 2)
    log.post = dt(t, df, delta * sqrt(effN), log = TRUE) + dcauchy(delta, scale = rscale, log = TRUE)
    ifelse(log, log.post, exp(log.post))
  }, "delta")



  pdf.posterior = function(delta, ...){
    norm.const = integrate(posterior, -Inf, Inf, ..., log = FALSE)[[1]]
    posterior(delta, ..., log = FALSE) / norm.const
  }


  cdf.posterior = function(delta, ...){
    norm.const = integrate(posterior, -Inf, Inf, ..., log = FALSE)[[1]]
    x = integrate(posterior, -Inf, delta, ..., log = FALSE)[[1]]
    x / norm.const
  }


  invcdf.posterior = Vectorize(function(p, t, N1, N2 = NULL, rscale = wid ){
    norm.const = integrate(posterior, -Inf, Inf,
                           t = t, N1 = N1, N2 = N2, rscale = wid, log = FALSE)[[1]]
    effN = ifelse(is.null(N2), N1, N1*N2/(N1+N2))
    d.est = t / sqrt(effN)
    d.se = 1/sqrt(effN)
    limits = qnorm(
      pnorm(
        qnorm(p, d.est, d.se
        ) + c(-1,1)*d.se, d.est, d.se),
      d.est, d.se)
    optimize(
      function(delta)(cdf.posterior(delta, t = t, N1 = N1, N2 = N2, rscale = rscale) - p)^2,
      limits)$minimum
  }, "p")

  ## ^^^ END OF UTILITY FUNCTIONS


  # *** Now get HDR of posterior using HDR function:

  CI = HDIofICDF(invcdf.posterior, tol = 1e-8, t = t, N1 = N1, N2 = N2, rscale = rscale )

  # posterior mean
  post.mean = integrate( function(x,...) x * pdf.posterior(x, ...), -Inf, Inf, t=t, N1=N1, N2=N2, rscale = rscale )[[1]]


  # E(delta^2)
  exp.delta2 = integrate( function(x,...) x^2 * pdf.posterior(x, ...), -Inf, Inf, t=t, N1=N1, N2=N2, rscale = rscale  )[[1]]


  # posterior variance
  post.var =  exp.delta2 - post.mean^2


  # posterior sd
  post.sd = sqrt(post.var)



  ## Determine X-axis range:

  x.min.1 = post.mean - 9 * post.sd

  x.max.1 = post.mean + 9 * post.sd


  ## The dL and dU may be different from x.min.1 and x.max.1 respectively, if so, adjust accordingly.

  x.min = if(dL < x.min.1) { dL + (.05*dL) } else { x.min.1 }

  x.max = if(dU > x.max.1) { dU + (.05*dU) } else { x.max.1 }


  ## plot the posterior:

  post.mode = optimize(function(delta) pdf.posterior(delta, t, N1, N2, rscale ),
                       interval = c(x.min, x.max), maximum = TRUE, tol = 1e-12)[[1]]

  post.peak = pdf.posterior(post.mode, t, N1, N2, rscale = rscale )


  par(family = 'serif', mar = c(.1, 4.1, 3.1, 2.1) )


  cc = curve(pdf.posterior(x, t, N1, N2, rscale), from = x.min, to = x.max, type="l",
             col = "navy", xlab = NA, cex.lab=1.5, font.lab=2 ,las=1, ylab ="Density",
             bty = "n", xlim = c(x.min, x.max), ylim = c(0, post.peak+(.1*post.peak)),
             cex = 2, font = 2, xaxt = "n", mgp = c(2.3, .75, 0), n = 1e3 )

  post.x = cc$x
  post.y = cc$y

  XXX <- post.x >= CI[1] &  post.x <= CI[2]

  low.extreme <- par('usr')[3]

  polygon(c(CI[1], post.x[XXX], CI[2]), c(low.extreme, post.y[XXX], low.extreme), col = rgb(1, 1, 0, .4), border = NA )


  segments(post.mode, low.extreme, post.mode, post.peak, lty = 3)

  text(post.mode, post.peak/2, decimal(post.mode, 2), srt = 90, pos = 3, font = 2)


  lines( post.x, post.y, lwd = 3, col = "navy" )


  segments(CI[1], low.extreme, CI[2], low.extreme, col = rgb(1, 0, 0, .8), lend = 1, lwd = 40 )

  segments(dL, low.extreme, dU, low.extreme, col = adjustcolor('green', alpha.f = .5), lend = 1, lwd = 40)


  points(post.mode, low.extreme/5, pch = 21, col = 0, bg = 0, cex = 1.5)


  axis(side = 1, at = decimal(seq(x.min, x.max, len = 7), 2), mgp = c(3, .75, 0), font = 2, cex.axis = 1.3)
  axis(side = 1, at = 0, cex.axis = 2, mgp = c(3, 1.1, 0), font = 2, col = 0, col.axis = "magenta", tcl = F, line = - 1.4 )

  mtext(side = 1, bquote(bold("Population Effect Size"~(delta))), line = 3, cex = 1.5)

  x1 = dL
  y1 = post.peak+.02*post.peak
  x2 = dU
  y2 = y1
  x.text = (dL+dU)/2
  y.text = post.peak+.05*post.peak


  segments( c(dL, dU), rep(low.extreme, 2), c(dL, dU), c(y1, y2), col = 'green2', lend = 1, lty = 2 )

  segments( c(x1, x2), c(y1, y2), rep(x.text, 2), rep(y.text*1.023, 2), lwd = 2, col = 'magenta' )


  text(x.text, y.text, "Practically Equivalent to ZERO", font = 2, pos = 3, col = 'darkgreen', cex = 1.1, xpd = TRUE)


  points( c(dL, dU), c(y1, y2), pch = 21, col = 'green3', bg = 'green3', cex = 1.2 )


  ## How much is it probable that the equivalence be true in population:

  a = cdf.posterior(delta = dL, t, N1, N2, rscale)
  b = cdf.posterior(delta = dU, t, N1, N2, rscale)

  Post.in.ROPE.Y = (b - a)
  Post.in.ROPE.X = (dU - dL) / 2


  BB = decimal( (b - a)*100, 2)


  legend("top", legend = c("There is", paste(BB,"%", sep=""), "probability that TRUE effect size is equivalent to ZERO"),
         bty="n", inset = c(-1, -.2), text.font=2, cex = 1.2, xpd = TRUE, ncol = 3, text.col = c(1, 'red', 1),
         x.intersp = c(-44.5, -1.5, -16)	 )


  legend("topleft", legend = paste("95% HDI: [",decimal(CI[1], 2), ", ", decimal(CI[2], 2),"]"),
         bty="n", inset=c(-.035,.1), text.font=2, text.col = 'red4', cex = 1.2)


  if(CI[1] > dU || CI[2] < dL) {

    legend("topright", "NOT Practically equivalent to \"0\" ", bty = 'n', inset = c(-.01, .1), cex = 1.2, text.font = 4, text.col = 'magenta2' )

  } else

    if(CI[1] > dL & CI[2] < dU) {

      legend("topright", "Practically equivalent to \"0\" ", bty = 'n', inset = c(-.01, .1), cex = 1.2, text.font = 4, text.col = 'magenta2')

    } else  {

      legend("topright", "No decision can be made ", bty = 'n', inset = c(-.01, .1), cex = 1.2, text.font = 4, text.col = 'magenta2')

    }


  ########################################################################################
  ## How choice of ROPE can affect porortion of posterior that ROPE covers (make function):
  ########################################################################################

  par( mar = c(3.1, 4.1, 6.1, 2.1), mgp = c(2.5, .5, 0) )


  eq.low = ifelse(abs(dL) <= .3, 4, 2)*( - ((dU - dL) / 2) )
  eq.upp = ifelse(abs(dL) <= .3, 4, 2)*(   ((dU - dL) / 2) )


  L = seq(eq.low, 0, length.out = 1e2)
  U = seq(eq.upp, 0, length.out = 1e2)


  aa = sapply(L, function(delta) cdf.posterior(delta, t, N1, N2, rscale) )
  bb = sapply(U, function(delta) cdf.posterior(delta, t, N1, N2, rscale) )


  Eq = (bb - aa)  # porortion of posterior that ROPE covers
  half = (U - L)/2


  plot(half, Eq, type = ifelse(abs(dL) == abs(dU), 'l' ,'n'), lwd = 3, col = 'red4', axes = FALSE, font.axis = 2,
       xlab = NA, ylab = paste("%",'Posterior in ROPE', sep = ""), font.lab = 2, cex.lab = 1.15)


  mtext(side = 1, "Half of ROPE", font = 2, line = 1.5)


  axis(1, at = decimal(seq(0, eq.upp[1], length.out = 7), 2), las = 1, font = 2)
  axis(2, at = seq(0, Eq[1], length.out = 5),
       labels = paste(100*round(seq(0, Eq[1], length.out = 5),
                                2), "%", sep = ""), las = 1, font = 2)


  u = par('usr')

  rect(u[1], u[3], u[2], u[4], col = adjustcolor("grey", .1), border = NA )


  rect(u[1], u[3], Post.in.ROPE.X, Post.in.ROPE.Y,
       col = adjustcolor("yellow",
       alpha.f = ifelse(Post.in.ROPE.Y  <= .2, .1,
       ifelse(Post.in.ROPE.Y > .2 & Post.in.ROPE.Y <= .3, .15,
       ifelse(Post.in.ROPE.Y > .3 & Post.in.ROPE.Y <= .4, .2, .3)))),
       lty = 2 )


  points(Post.in.ROPE.X, Post.in.ROPE.Y,
         pch = 21, cex = 2, bg ='green')

  box()

}
