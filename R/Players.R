#' Normality in Nature
#'
#' Uses the example of players movements on a sports field to show the mechanics
#' of how a "Normal Density" comes about.
#'
#' @param Step ith step or ith coin toss.
#' @param n.Players   Number of players on the field.
#' @param n.Steps     Total number of tosses/steps planned to take.
#' @return Provides positions and distribution of positions for any number of players on a sports field and generates
#'         a histogram after ith step.
#'
#' @details Shows the asymptotic behavior of binomial probability mass functions via simulation.
#' @author Reza Norouzian <rnorouzian@gmail.com>
#' @export
#' @examples
#'
#' # see 15 players' final positions and their distributions after taking
#' # 6 steps out of 16 total steps:
#'
#'
#' Players(Step = 16, n.Players = 1000, n.Steps = 16)
#'
#'
#'
#' # see 1000 players' final positions and their distributions after taking
#' # 16 steps out of 16 total steps:
#'
#'
#' Players(Step = 16, n.Players = 15, n.Steps = 16)



Players <- function(Step, n.Players, n.Steps){


  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))


  par(mar = c(2.2, 1.8, 1.5, 1.8) );
  m = matrix( c(1, 1, 1, 1, 1, 1,   2, 2), nrow = 8, ncol = 1 );
  layout(m)



  if(n.Players < 2) { n.Players = 2;
  message("\n\tYou can't have less then \"2 players\" on the field.")}


  if(Step > n.Steps) { Step = n.Steps;
  message("\n\tYou can't take more steps than what you planned to.\n\tWe picked the maximum steps that you planned to.")}


  if(Step < 0 || n.Steps < 1) { Step = 0; n.Steps = 1;
  message("\n\tYou can't have less than \"0\" steps.\n\tAlso, You can't have less than \"1\" as your total steps.")}



  plot(-6:6, -6:6, ty = "n", axes = F, ann = F)


  axis(1, at = c(-6:-1, 1:6), font.axis = 2, cex.axis = 1.5 )
  axis(1, at = 0, font.axis = 2, cex.axis = 1.9, col.axis = 'red' )


  par = par('usr')
  rect(par[1], par[3], par[2], par[4], col = 'darkseagreen1' )


  points( 0, 0, cex = 7, pch = 20, col = 0)
  points( 0, 0, cex = 40, lwd = 5, col = 0)
  abline(v = 0, lwd = 10, col = 0)


  rect(-6, -6, 6, 6, lwd = 5, border = 0)
  rect(-6.5, -2, -5.5, 2, col = 'darkseagreen1', border = 0, lwd = 5)
  rect(rep(-6.5, 2), rep(-2, 2), rep(-5.5, 2), rep(2, 2), border = 0, density = 10, angle = c(180, 90), col = 0)


  rect(6.5, -2, 5.5, 2, col = 'darkseagreen1', border = 0, lwd = 5)
  rect(rep(6.5, 2), rep(-2, 2), rep(5.5, 2), rep(2, 2), border = 0, density = 10, angle = c(180, 90), col = 0)



  box(col = 'darkseagreen1')

  points( c( rep(par[1]+.01, 2), rep(par[2]-.01, 2) ), rep( c(-2, 2), 2 ), cex = 3, pch = 20, col = 0 )



  ## Set up player movement:
  ##########################

  x <- rep(0, n.Players)                       ## Initial position of players
  y <- seq(from = -6, to = 6, len = n.Players) ## y-position for players



  ## Sample movement of players:
  xStepsMx <- matrix(sample(c(-1, 1)*.5, n.Players*n.Steps, replace = TRUE),
                     nrow = n.Players, ncol = n.Steps)


  ## Position of players:
  xPosMx <- t(sapply(1:nrow(xStepsMx), function(ii) cumsum(xStepsMx[ii,]))) + x



  positions = if (Step > 0){ xPosMx[,Step] } else {  x  }


  segments(positions, y, 0, y, lty = 2, col = 'red')
  points(positions, y, cex = 7, lwd = 3, pch = 21, bg = "white")
  text(positions, y, 1:n.Players, font = 2, cex = 1.5)


  if (Step == 0) {

    plot(1, 1, ty = 'n', axes = F, ann = F)

    text(1, 1, "Let players take at least \"1\" step on the field.", cex = 2.5, font = 2, col = 'red4')

  } else

  {
             # c(3.5, 4, 4, 2.1)
    par( mar = c(3.5, 4, 2, 2.1), mgp = c(2, .5, 0) )

    DENS = density(positions, adjust = 3)

    x.DENS = DENS$x
    y.DENS = DENS$y

    plot( DENS, col = 'red', lwd = 3, xaxt = "n", xlab = "Positions",
         ylab = "Probability", font.axis = 2, font.lab = 2, xlim = c(-6, 6), main = NA, bty = 'n',
         zero.line = F)

    low.exterme = par('usr')[3]

    x.DENS.2 = seq(-6, 6)
    y.DENS.2 = approx(x.DENS, y.DENS, xout = x.DENS.2 )

  # points( y.DENS.2, col = 'blue', cex = 2)

    polygon(DENS, col = rgb(0, 1, 1, .1), border = NA )

    axis(1, at = seq(-6, 6, len = 13), xlab = 'positions', font = 2, cex.axis = 1.2 )


    legend("topleft", legend = c(paste("Steps Taken = ", Step), paste("Players =", n.Players) ),
           bty="n", inset = c(-.02, .01), text.font=2, text.col = 'red4', cex = 1.5)


    }


}
