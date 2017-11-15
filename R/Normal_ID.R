#' Normal Prior Distribution Identifier
#'
#' Uses the subject matter researcher's knowledge to generate
#' a corresponding Normal prior distribution.
#'
#' @param Low researcher's LOWEST plausible value for the parameter.
#' @param High researcher's HIGHEST plausible value for the parameter.
#' @param Cover researcher's coverage for the Low and High values provided.
#'
#' @return Provides graphical as well as full textual description of a suitable Normal
#'         distribution for researcers based on their knowledge about how High or Low
#'         the parameter has been found in the literature. Also, helps researcher
#'         to revise their prior by issuing various messages.
#'
#' @details solves to provide graphical and textual information about
#'          an appropriate Normal prior distribution.
#'
#' @author  Reza Norouzian <rnorouzian@gmail.com>
#' @export
#'
#' @examples
#' # Suppose a researcher needs a Normal prior for a Cohen d effect size that in
#' # his/her view can't be less than -3 and more than +3. The researcher believes
#' # these two limit values cover 95% of all possible values that this parameter
#' # can take:
#'
#'
#'  Normal_ID (Low = -3, High = 3, Cover = '95%')
#'
#'
#'
#' # User can also use any value that is between 0 and 1 for the argument
#' # Cover without using percentage sign:
#'
#'
#'
#' Normal_ID (Low = -3, High = 3, Cover = 95)



Normal_ID = function(Low, High, Cover = NULL){


  original_par = par(no.readonly = T)
  on.exit(par(original_par))


  options(warn = -1)


  q <- c(Low, High)


  coverage  <- if (is.character(Cover)) { as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100

  } else if (is.numeric(Cover)) { Cover / 100 } else { .95 }


  Low.percentile = (1 - coverage) / 2


  p1 = Low.percentile
  p2 = Low.percentile + coverage


  alpha <- c(p1, p2)


  equal = function(a, b, sig = 4) { return (round(a, sig) == round(b, sig)) } # Complex, if L and estimated quantiles[1] are Unequal by 4 digits say TRUE


  decimal <- function(x, k){

    if( equal(x, 0) ){ format( round(0, k), nsmall = k ) }else

    { as.numeric(format(round(x, k), nsmall = k, scientific =
           ifelse(x >= 1e+05 || x <= -1e+05 || x <= 1e-05 & x >= -1e-05, TRUE, FALSE) )) }
  }


  if( p1 <= 0 || p2 >= 1 || q[1] >= q[2] || p1 >= p2 ) {

    par(family = 'serif')

    plot(1, axes = FALSE, type = 'n', ann = FALSE)

    text(1, 1, "Unable to find such a prior", cex = 3.5, col = 'red4', font = 2)

  } else {


    beta <- qnorm(alpha)

    mu.sigma = solve(cbind(1, beta), q)

    mean = mu.sigma[1]
    sd = mu.sigma[2]

    x.min = mean - 5*sd
    x.max = mean + 5*sd

    par(mgp = c(3.7, 1, 0), mar = c(5.1, 5.5, 4.1, 1.1) )

    curve ( dnorm(x, mean, sd), lwd = 4, from = x.min,
            to = x.max, xlab = 'Parameter of Interest', ylab = 'Density',
            n = 1e4, xaxt = 'n', las = 1, font.lab = 2, cex.lab = 1.4,
            frame.plot = FALSE, font.axis = 2, cex.axis = 1.1)

    prior.peak = dnorm(mean, mean, sd)
    low.exterme = par('usr')[3]


    segments(mean, low.exterme, mean,
             prior.peak, lty = 3)

    axis(1, at = decimal(seq(x.min, x.max, length.out = 9), 1), font = 2, cex.axis = 1.3 )

    arrows(q[1], 0, q[2], 0, lwd = 2, col = 'red', angle = 90, code = 3, length = .15)

    text(c(q[1], q[2]), rep(0, 2), c(q[1], q[2]), col = 'blue', pos = 3, font = 2, cex = 2)


    if ( equal(mean, 0) ) { text(mean, prior.peak / 2, "Neutral Position", cex = 1.5, pos = 3, srt = 90, font = 2) }


    mtext(side = 3, "This is the \"Normal Prior\" you have in mind", cex = 1.5, bty = 'n', font = 2)
    mtext(side = 3, bquote(bold((mu) == .(decimal (mean, 3)))), line = -4, cex = 1.8, adj = .05, col = 'red4')
    mtext(side = 3, bquote(bold((sigma) == .(decimal (sd, 3)))), line = -6, cex = 1.8, adj = .05, col = 'red4')


  }

  if( p1 <= 0 || p2 >= 1 || q[1] >= q[2] || p1 >= p2 ) {

    message("\n\tUnable to find such a prior, make sure you have selected the correct values.")

  }else

  { message("\nCAUTION: \"ALWAYS\" visually inspect the shape of the prior generated to see \n \t  if it accurately represents your belief and revise if necessary.\n")

    setNames( c(decimal(mean, 7), decimal(sd, 7) ), c("mu","sigma"))

  }

}
