#' Cauchy Prior Distribution Identifier
#'
#' Uses the subject matter researcher's knowledge to generate
#' a corresponding Cauchy prior distribution.
#'
#' @param Low researcher's LOWEST plausible value for the parameter.
#' @param High researcher's HIGHEST plausible value for the parameter.
#' @param Cover researcher's suggested coverage for the Low and High values provided.
#'
#' @return Provides graphical as well as full textual description of a suitable Cauchy
#'         distribution for researcers based on their knowledge about how High or Low
#'         the parameter has been found in the literature. Also, helps researcher
#'         to revise their prior by issuing various messages.
#'
#' @details Uses optimization techniques to provide graphical and textual information about
#'          an appropriate Cauchy prior distribution.
#'
#' @author  Reza Norouzian <rnorouzian@gmail.com>
#' @export
#'
#' @examples
#' # Suppose a researcher needs a Cauchy prior for a Cohen d effect size that in
#' # his/her view can't be less than -6 and more than +6. The researcher believes
#' # these two limit values cover 90% of all possible values that this parameter
#' # can take:
#'
#'
#'  Cauchy_ID (Low = -6, High = 6, Cover = '90%')
#'
#'
#'
#' # User can also use any value that is between 0 and 1 for the argument
#' # Cover without using percentage sign:
#'
#'
#'
#' Cauchy_ID (Low = -6, High = 6, Cover = 90)
#'



Cauchy_ID = function (Low, High, Cover= NULL){


  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))


  options(warn = -1)



  coverage  <- if (is.character(Cover)) { as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100

  } else if (is.numeric(Cover)) { Cover / 100 } else { .90 }


  Low.percentile = (1 - coverage) / 2


  p1 = Low.percentile
  p2 = Low.percentile + coverage



  ## Start Optimization:

  if( p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1 ) {


    par(family = 'serif')

    plot(1, axes = FALSE, type = 'n', ann = FALSE)

    text(1, 1, "Unable to find such a prior", cex = 3.5, col = 'red4', font = 2)

    return( message("\n\tUnable to find such a prior, make sure you have selected the correct values.") )



  } else {



    f <- function(x) {

      y <- c(Low, High) - qcauchy(c(p1, p2), location=x[1],  scale=x[2])

    }



    ## SOLVE:

    AA <- optim(c(1, 1), function(x) sum(f(x)^2), control=list(reltol=(.Machine$double.eps)) )


    parms = unname(AA$par)

  }


  ## CHECK:
  q <- qcauchy( c(p1, p2), parms[1], parms[2] )


  unequal = function(a, b, sig = 4) { return (round(a, sig) != round(b, sig) ) } # Complex, if Low & High and estimated quantiles[1 & 2] are Unequal by 4 digits say TRUE


  if( p1 <= 0 || p2 >= 1 || Low >= High || p1 >= p2 || unequal(Low, q[1]) || unequal(High, q[2]) ) {

    par(family = 'serif')

    plot(1, axes = FALSE, type = 'n', ann = FALSE)

    text(1, 1, "Unable to find such a prior", cex = 3.5, col = 'red4', font = 2)

    message("\n\tUnable to find such a prior, make sure you have selected the correct values")


  } else


  {


    equal = function(a, b, sig = 4) { return (round(a, sig) == round(b, sig)) } # Complex, if L and estimated quantiles[1] are Unequal by 4 digits say TRUE


    decimal <- function(x, k){

      if( equal(x, 0) ){ format( round(0, k), nsmall = k ) } else

      { as.numeric(format(round(x, k), nsmall = k, scientific =
           ifelse(x >= 1e+05 || x <= -1e+05 || x <= 1e-05 & x >= -1e-05, TRUE, FALSE) )) }
    }


    ## call 'location' mean and 'scale' sd fo simplicity:
    mean = parms[1]
    sd = parms[2]

    x.min = mean - 12*sd
    x.max = mean + 12*sd

    par(mgp = c(3.7, 1, 0), mar = c(5.1, 5.5, 4.1, 1.1) )

    curve ( dcauchy(x, mean, sd), lwd = 4, from = x.min,
            to = x.max, xlab = 'Parameter of Interest', ylab = 'Density',
            n = 1e4, xaxt = 'n', las = 1, font.lab = 2, cex.lab = 1.4,
            frame.plot = F, font.axis = 2, cex.axis = 1.1 )


    axis(1, at = decimal(seq(x.min, x.max, length.out = 9), 1), font = 2, cex.axis = 1.3 )

    low.extreme = par('usr')[3]
    prior.peak = dcauchy(mean, mean, sd)


    segments(mean, low.extreme, mean, prior.peak, lty = 3)


    arrows(q[1], 0, q[2], 0, lwd = 2, col = 'red', angle = 90, code = 3, length = .15)

    text(c(q[1],q[2]), rep(0, 2), round(c(q[1], q[2]), 3), col = 'blue', pos = 3, font = 2, cex = 2, xpd = TRUE)

    mtext(side = 3, "This is the \"Cauchy Prior\" you have in mind", cex = 1.5, bty = 'n', font = 2)
    mtext(side = 3, bquote(bold(Mode == .(decimal (mean, 3)))), line = -4, cex = 1.8, adj = .05, col = 'red4')
    mtext(side = 3, bquote(bold(Scale == .(decimal (sd, 3)))), line = -6, cex = 1.8, adj = .05, col = 'red4')

    cat(message("\nCAUTION: \"ALWAYS\" visually inspect the shape of the prior generated to see \n \t  if it accurately represents your belief and revise if necessary.\n"))
    cat(message("\nNOTE: \"Cauchy\" is like a NORMAL distribution but has VERY VERY EXTENDED tails.\n\tThus, using a coverage of \"90%\" for the low and high values is enough .\n"))

    if (all.equal(mean, 0, tol = 1e-4)) { text(mean, prior.peak / 3, "Neutral Position", cex = 1.5, pos = 3, srt = 90, font = 2) }

    setNames(c(decimal(parms[1], 7), decimal(parms[2], 7) ), c("Mode","Scale") )

  }

}
