#' Beta Prior Distribution Identifier
#'
#' Uses the subject matter researcher's knowledge to generate
#' a corresponding Beta prior distribution.
#'
#' @param Low researcher's LOWEST plausible value for the parameter.
#' @param High researcher's HIGHEST plausible value for the parameter.
#' @param Cover researcher's suggested coverage for the Low and High values provided.
#'
#' @return Provides graphical as well as full textual description of a suitable Beta
#'         distribution for researcers based on their knowledge about how High or Low
#'         the parameter has been found in the literature. Also, helps researcher
#'         to revise their prior by issuing various messages.
#'
#' @details Uses non-linear minimization to provide an optimal Beta prior distribution.
#' @author  Reza Norouzian <rnorouzian@gmail.com>
#' @export
#'
#' @examples
#' # Suppose a researcher needs a Beta prior for a parameter that in his/her
#' # view can't be less than 2% and more than 90%. The researcher believes
#' # these two limit values cover 95% of all possible values that the parameter
#' # of interest to him/her can take:
#'
#'
#'  Beta_ID (Low = "2%", High = "90%", Cover = '95%')
#'
#'
#'
#' # User can also use any value that is between 0 and 1 for all arguments
#' # without using percentage sign for some or all arguments:
#'
#'
#' Beta_ID (Low = .02, High = .8, Cover = "90%")
#'
#'
#' Beta_ID (Low = .02, High = .8, Cover = 90)



Beta_ID = function(Low, High, Cover = NULL){

  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))


  options(warn = -1)


  L <- if(is.character(Low)){as.numeric(substr(Low, 1, nchar(Low)-1)) / 100

  }else{ Low }


  U <- if(is.character(High)){as.numeric(substr(High, 1, nchar(High)-1)) / 100

  }else{ High }


  if(L <= 0 || U >= 1){stop("NOTE: The smallest LOWER value that you can choose is \".000001\"AND the largest UPPER value is \".999999\".") }


  if(L >= U){stop("Put the smaller value for Low, and the larger value for High")}


  coverage  <- if (is.character(Cover)) { as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100

  } else if (is.numeric(Cover)) { Cover / 100 } else { .95 }


  Low.percentile = (1 - coverage) / 2


  p1 = Low.percentile
  p2 = Low.percentile + coverage


  if( p1 <= 0 || p2 >= 1 || L > U || p1 > p2 || coverage >= 1 ){

    par(family = 'serif')

    plot(1, axes = FALSE, type = 'n', ann = FALSE)

    text(1, 1, "Unable to find such a prior", cex = 3.5, col = 'red4', font = 2)

    return( message("\n\tUnable to find such a prior, make sure you have selected the correct values.") )

  } else {


    ## Utility functions:

    #Logistic transformation of the Beta CDF.
    #
    f.beta <- function(alpha, beta, x, lower=0, upper=1) {
      p <- pbeta((x-lower)/(upper-lower), alpha, beta)
      log(p/(1-p))
    }
    #
    # Sums of squares.
    #
    delta <- function(fit, actual) sum((fit-actual)^2)
    #
    # The objective function handles the transformed parameters `theta` and
    # uses `f.beta` and `delta` to fit the values and measure their discrepancies.
    #
    objective <- function(theta, x, prob, ...) {
      ab <- exp(theta) # Parameters are the *logs* of alpha and beta
      fit <- f.beta(ab[1], ab[2], x, ...)
      return (delta(fit, prob))
    }
    #
    # Solve the problem.
    #

    decimal <- function(x, k){

      if(is.character(x)) {

        return(x)
      }

as.numeric(format(round(x, k), nsmall = k, scientific =
           ifelse(x >= 1e+05 || x <= -1e+05 || x <= 1e-05 & x >= -1e-05, TRUE, FALSE) ))
    }

    #@@@@@@@ Get the shape parameters of a beta Dist. knowing its quantiles

    # Initialize:
    x <- c(L, U) ## insert the quantiles here

    x.p <- (function(p) log(p/(1-p)))(c(p1, p2))


    # Solve:
    start <- log(c(1e1, 1e1))

    sol <- nlm(objective, start, x=x, prob=x.p, lower=0, upper=1, typsize=c(1,1),
               fscale=1e-12, gradtol=1e-12)

    parms <- exp(sol$estimate)


    # Check:

    quantiles = qbeta(p = c(p1, p2), parms[1], parms[2])


    beta.mode = (parms[1] - 1) / ( parms[1] + parms[2] - 2 )

    beta.peak = dbeta(beta.mode, parms[1], parms[2])

    beta.med = ( parms[1] - (1/3) ) / ( (parms[1] + parms[2]) - (2/3) )

    beta.mean = parms[1] / ( parms[1] + parms[2] )

    beta.var = ( parms[1] * parms[2] ) / ( ( (parms[1] + parms[2] )^2 ) * (parms[1] + parms[2] + 1 ) )

    beta.sd = sqrt( beta.var )


    par(mar = c(6.1, 4.8, 4.1, 2.1) )

    ylim = c(0, beta.peak+.23*beta.peak)

    ylims = if(L < .013 || parms[1] <= 1 || parms[2] <= 1){ NULL }else{ ylim }


    curve( dbeta(x, parms[1], parms[2]), from = 0, to = 1, lwd = 4, n = 1e4,
          ylab = 'Density', font.lab = 2, cex.lab = 1.3, bty = "n", xpd = TRUE,
          xaxt = 'n', xlab = NA, las = 1, font.axis = 2, ylim = ylims, cex.axis = 1.2)


    axis(1, at = round(seq(0, 1, length.out = 9), 2), font = 2, cex.axis = 1.3 )

    axis(1, at = round(seq(0, 1, length.out = 9), 2),
         labels = paste(100* round(seq(0, 1, length.out = 9), 2),"%", sep=""),
         pos = par("usr")[3] - 1 * 0.05 * (par("usr")[4] - par("usr")[3]),
         tick = FALSE, font = 2, cex.axis = 1.3 ) #2nd axis labels


    segments(beta.mode, par('usr')[3], beta.mode, beta.peak ,lty = 3, xpd = TRUE)

    arrows(quantiles[1], 0, quantiles[2], 0, code = 3, lwd = 2, lend = 1, angle = 90, length = .15, col = 'red')

    text(c(quantiles[1], quantiles[2]), rep(0, 2), decimal(c(quantiles[1], quantiles[2]), 3 ), cex = 2, pos = 3, col = 'blue', font = 2, xpd = T )


    mtext(side = 1, 'Parameter of Interest', line = 4.5, font = 2, cex = 1.3)
    mtext(side = 3, bquote(bold((alpha) == .(decimal(parms[1], 3)))), cex = 1.8, line = -2.6, adj = .05, col = 'red4')
    mtext(side = 3, bquote(bold((beta) == .(decimal(parms[2], 3)))), cex = 1.8, line = -4.5, adj = .05, col = 'red4')
    mtext(side = 3, "This is the \"Beta Prior\" you have in mind", cex = 1.5, font = 2)


    unequal = function(a, b, sig = 3) { return (round(a, sig) != round(b, sig)) } # Complex, if L and estimated quantiles[1] are Unequal by 4 digits say TRUE


    if( unequal(L, quantiles[1]) || unequal(U, quantiles[2]) ) {

      par(family = 'serif')

      plot(1, axes = F, ty = 'n', ann = F)

      text(1, 1, "Unable to find such a prior", cex = 3.5, col = 'red4', font = 2)

    }

    if( unequal(L, quantiles[1]) || unequal(U, quantiles[2]) ) {

      message("\n\tUnable to find such a prior, make sure you have selected the correct values.")

    } else

    {

      message("\nNOTE 1: If you don't specify the coverage for your provided values, we pick\n \t a 95% coverage for you.")

      message("\nNOTE 2: \"ALWAYS\" visually inspect the shape of the prior generated to see\n \tif it accurately represents your belief and revise if necessary.\n")

      message("\nCAUTION: This prior generator uses an advanced problem solving technique  \n \t called \"Non-Linear Minimazation (nlm)\", when user defines impossible \n \t values a message will appear signaling the impossiblity to provide \n \t a solution.\n")


      list(alpha = decimal(parms[1], 8), beta = decimal(parms[2], 8), Mean = decimal(beta.mean, 8),
           Median = decimal(beta.med, 8), Mode = decimal(beta.mode, 8), SD = decimal(beta.sd, 8) )

    }

  }

}
