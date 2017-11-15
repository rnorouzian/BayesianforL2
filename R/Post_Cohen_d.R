#' Bayesian Estimation and Model Selection for Cohen's d in t-test Designs
#'
#' Provides complete Summary of Posterior, and Hypothesis Testing Results (i.e., Bayes Factors).
#'
#' @seealso \url{https://rezanorouzian.shinyapps.io/bayes-t/} for an
#'           interactive and more powerful version of this function
#'            which also compares the frequentist and the Bayesian results
#'            with each other.
#'
#'
#' @param  N1 Group 1 sample size.
#' @param  N2 Group 2 sample size (if a two-sample design).
#' @param  ttype Type of t test.
#' @param  tl Type of Hypothesis.
#' @param  t  Observed t value.
#' @param  wid Width or Scale of a Cauchy prior.
#' @param  dexp Expected direction of Cohens d before doing the study.
#'
#' @return  Provides graphical as well as full textual description of posterior distribution of
#'          effect size. Additionally, does Bayesian hypothesis testing producing Bayes Factors.
#'
#'
#' @author  Reza Norouzian <rnorouzian@gmail.com>
#'
#' @import plotrix
#' @import pracma
#'
#' @export
#'
#' @examples
#' # Suppose a researcher obtains a t value of 1, from 2 groups each with 20 participants.
#' # Thus, ttype is 2. The researcher picks a "wi", wide, Cauchy prior. Suppose
#' # researcher has no prefernce for the direction of Cohen d. Thus, tl 2.
#' # when hypothesis is two tailed, dexp is ignored. The researcher can use:
#'
#'
#' # Post_Cohen_d(t = .5, N1 = 20, ttype = 1)




Post_Cohen_d = function(t, N1, N2 = NULL, tl = 2, ttype, dexp = 1, wid = "wi"){


  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))


  options(warn=-1)


  ## decimal display controller:

  decimal <- function(x, k) { as.numeric(format(round(x, k), nsmall = k, scientific = 
           ifelse(x >= 1e+05 || x <= -1e+05 || x <= 1e-05 & x >= -1e-05, TRUE, FALSE) )) }


  par(mgp=c(2.2, .75, 0), mar=c(6.1, 4.1, 5.1, 1.5) )

  ### Data

  wid <- c(m = 1/2, wi = sqrt(2)/2, wii = 1)[[wid]]

  rscale = wid


  ## Integral bounds:

  AAAAA <- if(tl==1 & dexp ==  1) { 0 } else
    if(tl==1 & dexp == -1){   -Inf   }else
    { -Inf }

  BBBBB <- if(tl==1 & dexp ==  1) { Inf } else
    if(tl==1 & dexp == -1){   0   }else
    { Inf }


  ### Start of utility functions

  posterior = Vectorize(function(delta, t, N1, N2 = NULL, rscale = rscale, log = FALSE){
    effN = ifelse(is.null(N2), N1, N1*N2/(N1+N2))
    df = ifelse(is.null(N2), N1 - 1, N1 + N2 - 2)
    log.post = dt(t, df, delta * sqrt(effN), log = TRUE) + dcauchy(delta, scale = rscale, log = TRUE)
    ifelse(log, log.post, exp(log.post))
  }, "delta")


  pdf.posterior = function(delta, ...){
    norm.const = integrate(posterior, AAAAA, BBBBB, ..., log = FALSE)[[1]]
    posterior(delta, ..., log = FALSE) / norm.const
  }


  cdf.posterior = function(delta, ...){
    norm.const = integrate(posterior, AAAAA, BBBBB, ..., log = FALSE)[[1]]
    x = integrate(posterior, AAAAA, delta, ..., log = FALSE)[[1]]
    x / norm.const
  }


  invcdf.posterior = Vectorize(function(p, t, N1, N2 = NULL, rscale = rscale){
    norm.const = integrate(posterior, AAAAA, BBBBB,
                           t = t, N1 = N1, N2 = N2, rscale = rscale, log = FALSE)[[1]]
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


  ### End of utility functions


  # 95% credible interval (two-sided tests)

  cii <- invcdf.posterior(c(.025, .975), t, N1, N2, rscale = rscale)


  ddd = seq(ifelse(AAAAA==-Inf, -6, 0), ifelse(BBBBB==Inf, 6, 0), len = 1e4)

  post <- pdf.posterior(ddd, t, N1, N2, rscale = rscale )


  modd <- ddd[ which.max(post) ]

  Points <- data.frame(ddd, post)

  ## Compute 95% CI for one-sided t-tests using :

  area <- cumtrapz(Points$ddd, Points$post)
  ind1 <- which(area > .025)[1] ## 95% from the left to the right
  ind2 <- which(area > .975)[1] ## "     "        "        "
  ind3 <- which(area > .5)[1]   ## Median

  ci11 <- Points$ddd[ind1]
  ci22 <- Points$ddd[ind2]


  # 95% CI for one- and two-sided t-tests:

  ci <- if (tl==2){ cii } else { c(ci11, ci22) }

  # posterior median
  post.med <- if(ttype==1 & tl==2){invcdf.posterior(.5, t, N1, rscale = rscale)}else
    if(ttype==2 & tl==2){invcdf.posterior(.5, t, N1, N2, rscale = rscale)}else
    { Points$ddd[ind3] }


  # posterior mean
  post.mean <- if(ttype==1) {

    integrate( function(x,...) x * pdf.posterior(x, ...), AAAAA, BBBBB, t=t, N1=N1, rscale = wid )[[1]]

  }else{

    integrate( function(x,...) x * pdf.posterior(x, ...), AAAAA, BBBBB, t=t, N1=N1, N2=N2, rscale = wid )[[1]]

  }


  ZERO.POINT.POST <- if(ttype==1) {pdf.posterior(0, t, N1, rscale = rscale ) }else
  {  pdf.posterior(0, t, N1, N2, rscale = rscale )  }


  ZERO.POINT.PRIOR <- dcauchy(0, 0, rscale)


  # Posterior SD
  exp.delta2 <- if(ttype==1) {

    integrate( function(x,...) x^2 * pdf.posterior(x, ...), AAAAA, BBBBB, t=t, N1=N1, rscale = wid  )[[1]]

  }else{

    integrate( function(x,...) x^2 * pdf.posterior(x, ...), AAAAA, BBBBB, t=t, N1=N1, N2=N2, rscale = wid  )[[1]]

  }

  post.var <-  exp.delta2 - post.mean^2

  post.sd <- sqrt(post.var)


  ######################
  ## Plot
  ######################


  plot(ddd, post, type="l", col="navyblue", xlab = expression(bold("Population Effect Size " (delta))),cex.lab=1.5,
       font.lab=2 ,las=1, ylab ="Density", lwd=3, bty="n", xlim=c(-6, 6), ylim = c(0, max(post)+(.4*max(post))),
       cex=2, font =2)


  XXX <- ddd >= ci[1] &  ddd <= ci[2]

  polygon(c(ci[1], ddd[XXX], ci[2]), c(0, post[XXX], 0), col=rgb(0, 0, 0, .3), density = 20, border = NA)


  rect(ci[1], 0, ci[2], max(post)+.1*max(post), lty=5, col=rgb(1, 1, 1, 0), border=rgb(1, 0, 0, .2) )


  segments(modd, max(post), modd, max(post)+.1*max(post), lty=3, col="gray80")


  arrows(ci[1], max(post)+.1*max(post), ci[2], max(post)+.1*max(post), col="red", length=.1, code = 3,
         lwd=2, angle = 90,xpd=T)


  legend("topright",c("Prior","Posterior"), lty=c(3, 1), lwd=c(3, 3), col=c("cyan2","navyblue"),
         cex=1.3, inset=c(.01,.08), text.font = 2, box.col = "red4" )

  par(family="serif", xpd = T)


  legend("top", legend=paste("95% CI: [",decimal(ci[1], 4), ", ", decimal(ci[2], 4),"]"),
         bty="n", inset=c(0,.047), text.font=2, cex = 1.5)


  legend("topleft", legend=bquote(bold(SD[~(delta)] == .(decimal(post.sd, 4)))),
         bty="n", inset=c(0, .05), text.font=2, cex = 1.5)


  legend("top", legend=bquote(bold(Mode[~(delta)] == .(decimal(modd, 4)))),
         bty="n", inset=c(0, -.02), text.font=2, cex = 1.5)


  legend("topleft", legend=bquote(bold(Mean[~(delta)] == .(decimal(post.mean, 4)))),
         bty="n", inset=c(0, -.02), text.font=2, cex = 1.5)


  legend("topright", legend=bquote(bold(Median[~(delta)] == .(decimal(post.med, 4)))),
         bty="n", inset=c(0, -.02), text.font=2, cex = 1.5)


  par(family="sans", xpd = F)


  if(tl==1){ segments(ifelse(dexp==1, -6, 6), 0, 0, 0, lwd=3, col="navyblue");
    segments(0, 0, 0, Points$post[Points$ddd==0], lwd=3 , col="navyblue") }

  segments(modd, par("usr")[3], modd, max(post), xpd=T, col="navyblue", lty=3 )

  segments(ifelse(0 > modd, 1, -1), ZERO.POINT.PRIOR, ifelse(0 > modd, 1, -1), ZERO.POINT.POST, xpd=T, col="magenta", lty=2 )

  segments(rep(ifelse(0 > modd, 1, -1), 2), c(ZERO.POINT.PRIOR, ZERO.POINT.POST), rep(0, 2), c(ZERO.POINT.PRIOR, ZERO.POINT.POST), xpd=T, col="magenta", lty=2)


  #############
  ##add prior:
  #############


  if(tl== 2) {  curve(dcauchy(x, 0, rscale), -6, 6, add=TRUE, lty=3, lwd=3, col="cyan2")  }

  if(tl== 1 & dexp > 0) { curve(dcauchy(x, 0, rscale),  0, 6, add=TRUE, lty=3, lwd=3, col="cyan2") }

  if(tl== 1 & dexp < 0) { curve(dcauchy(x, 0, rscale), -6, 0 , add=TRUE, lty=3, lwd=3,  col="cyan2") }

  if(tl== 1){ segments(0, 0, 0, ZERO.POINT.PRIOR, lwd=3, lty=3, col="cyan2" ) }

  points(0, ZERO.POINT.PRIOR, pch = 21, cex = 2.2, bg="green", col="red")

  points(0, ZERO.POINT.POST, pch = 21, cex = 2.2, bg="green", col="red")

  BF10 <- ZERO.POINT.PRIOR*ifelse(tl==1, 2, 1) / ZERO.POINT.POST   ## When one-sided multiply prior by 2 as prior is halved


  ## Center the piechart:
  radius <- .6
  AB <- radius^2 * pi
  alpha <- 2/(1/BF10 + 1) * AB/radius^2
  startpos <- pi/2 - alpha/2


  floating.pie(ifelse(t <= -3 & dexp==-1 || t <= -3 & tl==2, 4.5, -4.5), max(post),c(BF10, 1), radius=radius, col=c("red4", 0 ), lwd=2, startpos=startpos)

  if(tl==2) text(ifelse(t <= -3 & dexp==-1 || t <= -3 & tl==2, 4.5, -4.5), max(post)-(.2*max(post)), bquote(bold(BF[10])== .(decimal(BF10, 3) )), cex = 1.2)

  if(tl==1) text(ifelse(t <= -3 & dexp==-1 || t <= -3 & tl==2, 4.5, -4.5), max(post)-(.2*max(post)), bquote(bold(BF['1+'])== .(decimal(BF10, 3) )), cex = 1.2)

  text(ifelse(t <= -3 & dexp==-1 || t <= -3 & tl==2, 4.5, -4.5), max(post)+(.14*max(post)), expression(paste('Support for H'[1])), cex = 1.2, col="red4")

  text(ifelse(t <= -3 & dexp==-1 || t <= -3 & tl==2, 4.5, -4.5), max(post)-(.14*max(post)), expression(paste('Support for H'[0])), cex = 1.2, col="cyan4")


  ###################################
  ## Show how BF calculatons wee done
  ###################################


  txtt = function(n){
    return(as.character( decimal(n, 2) ) )
  }


  G1 = ZERO.POINT.PRIOR
  G2 = ZERO.POINT.POST


  if(tl==2) text(ifelse(0 > modd, 2.5, -2.5), (G1 + G2)/2,
                 bquote(paste(BF[10]," = ", frac(.(txtt(G1)), .(txtt(G2))),
                              .(ifelse(tl == 1, " \u00D7 2", "")),
                              " = ",
                              .(ifelse(tl == 1, txtt(G1/G2 * 2), txtt(G1/G2 ))), sep = "")),
                 col = "red4", cex = 1.4, xpd=T)


  if(tl==1) text(ifelse(0 > modd, 2.5, -2.5), (G1 + G2)/2,
                 bquote(paste(BF['1+']," = ", frac(.(txtt(G1)), .(txtt(G2))),
                              .(ifelse(tl == 1, " \u00D7 2", "")),
                              " = ",
                              .(ifelse(tl == 1, txtt(G1/G2 * 2), txtt(G1/G2 ))), sep = "")),
                 col = "red4", cex = 1.4, xpd=T)


}

