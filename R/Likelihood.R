#' The Structure of a Likelihood Function
#'
#' Details the Likelihood Function as an ingridient of a Bayesian estimation process.
#'
#'
#' @param  n Sample size.
#'
#'
#' @return  Graphically details the process of obtaining Likelihood values
#'          to produce a "Likelihood Function" using Simulation.
#'
#'
#' @author  Reza Norouzian <rnorouzian@gmail.com>
#' @export
#'
#' @examples
#'
#' # Likelihood(n = 3)


Likelihood = function (n){

  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))

  m = 70; SD = 15
  observations =  rnorm(n, m, SD)

  obs.mean = mean(observations)

  decimal <- function(x, k){

    format(round(x, k), nsmall = k, scientific =
             ifelse(x >= 1e+05 || x <= -1e+05 || x <= 1e-05 & x >= -1e-05, T, F) )
  }



  x.min = obs.mean - 5*SD
  x.max = obs.mean + 5*SD

  x.axis.range = seq(obs.mean - 5*SD, obs.mean + 5*SD, len = 8)

  x.left = 1
  x.right = 150


  #Or re-Write the above Likelihood as:

  par(mgp = c(2, 1, 0), mar = c(5.1, 4.1, 4.1, 3.1), family = 'serif' )

  Likelihood <- function(x) sapply(lapply(x, dnorm, x = observations, SD), prod)


  Like.function = curve( Likelihood, ylab =  bquote(bolditalic('Likelihood'~(mu))),
                         from = x.left, to = x.right, lwd = 3
                         ,axes = F, xlab = NA,n = 1e4, bty = 'n', cex.lab = 1.3)$y


  axis(side = 1, at = decimal( seq(x.left, x.right, len = 8), 1), cex.axis = 1, font = 2 )


  like.mode = obs.mean

  like.peak = max(Like.function)

  #axis(side = 2, at = decimal(seq(0, like.peak+.05*like.peak, len = 4), 6), las = 1  )


  mtext(side = 1, bquote((mu)), line = 3, cex = 1.3 )


  low.extreme = par('usr')[3]


  segments(like.mode, low.extreme, like.mode, like.peak, lty = 3 )


  points( observations, rep(low.extreme, n), pch = 20, cex = 2, xpd = T, col = rgb(0, 0, 1, .8))

  index = 1:n

  text(observations, rep(low.extreme, n), paste("Data ", index, sep = ""),
       col = rgb(1, 0, 0, .9), font = 2, pos = 3, xpd = T, cex = .85 )

  text(like.mode, like.peak/2.5, "Observed Estimate", col = 'gray50', font = 2,
       cex = .85, pos = 3, srt = 90)

  points(like.mode, low.extreme, pch = 22, bg = rgb(0, 1, 0, .7), col = NA, cex = 1.8, xpd = T)

  legend('topleft', legend = "Actually Collected Data", pch = 20, col = rgb(0, 0, 1, .8), text.font = 2,
         text.col = 'red', pt.cex = 1.6, bty = 'n', inset = c (-.03, .001), xpd = T, cex = .85,
         y.intersp = 0)

  legend('topright', legend = 'Observed Estimate (Mean)', pch = 15, col = rgb(0, 1, 0, .9), text.font = 2,
         text.col = 'red', pt.cex = 1.3, bty = 'n', inset = c (-.03, .001), xpd = T, cex = .85,
         y.intersp = 0)

  mtext(side = 3, "Method of Maximum Likelihood", cex = 2, font = 2, line = .5)


  par(family = 'sans')

}
