#' Anatomy of Two Prior Probability Density Functions
#'
#' Using advanced graphics, visualizes the behavior of (1) a Normal prior, and
#' (2) a Cauchy priors each presented in two forms. First, as a complete prior over
#' all possible parameter values, and (b) as a halved prior over positive half of parameter
#' values.
#'
#'
#' @param  type "Normal" or "Cauchy".
#' @param  width a numeric value determining the spread of the distribution.
#' @param  half if TRUE provides half of the chosen prior distribution.
#'
#' @return  Graphically details the exact behavior of (1) a Normal prior, and
#'          (2) a Cauchy prior each presented in two form. First, as complete prior over
#'          all possible parameter values, and (b) as a halved prior over positive half of
#'          parameter values.
#'
#'
#' @author  Reza Norouzian <rnorouzian@gmail.com>
#' @export
#'
#' @examples
#'
#' # Prior_Anatomy(type = "Normal", width = 1, half = F)



Prior_Anatomy = function(type, half, width){

  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))

  #INITIAL DATA
  el_x = 0
  el_y = -.5
  n = 11
  multiplier = 1
  el_a = seq(1, 6, length.out = n) #a values of ellipse
  el_b = seq(1/25, 6/25, length.out = n) #b values

  graphics.off()
  require(plotrix)

  width = if(type== "Cauchy" & width=="wide") {
    sqrt(2)/2
  } else if(type== "Cauchy" & width=="medium") {
    1/2
  } else if(type== "Cauchy" & width=="very wide") {
    1
  } else if(type== "Normal" & is.character(width) ){
    stop(message(cat("You must provide a number")))
  } else {
    width
  }

  par(family="serif")

  decimal <- function(x, k){
    if(any(is.character(x))){
      return(x)
    }
    x_decimal = character(0)
    for (x_ind in seq_along(x)){
      if (abs(x[x_ind]) >= 1e+05 || abs(x[x_ind]) <= 1e-05){
        x_decimal[x_ind] = format(x[x_ind], digits = k+1, scientific = TRUE)
      } else {
        x_decimal[x_ind] = format(round(x[x_ind], k), nsmall = k, scientific = FALSE)
      }
    }
    return(x_decimal)
  }

  #Retain only every other 'a' value
  if (n %%2 !=0){ #When 'n' is odd
    el_a_2 = el_a[seq_along(el_a) %% 2 != 0]
  } else { #When 'n' is even
    el_a_2 = el_a[seq_along(el_a) %% 2 == 0]
  }

  #HALF
  if (half == TRUE){
    ff = el_x #The x_value where curve starts
    AA = c(el_x, el_x + el_a_2) #The x-values you need for segments
    AA2 = c(el_x - el_a_2, el_x, el_x + el_a_2) #x_values of points
    pol_ind = c(1,2) #index of x to cover with polygon
    seg_col = c(rep('navy', 3), rep('green3', length(AA) - 3))
    text_cex = c(rep(1.4, 2),
                 rep(1.1, length(AA) - 2) )
    curve_points_col = c('red', rep('blue', 6))
    curve_points_bg = c('red', rep('blue', 6))
  } else{
    ff = min(el_x - el_a) #The x_value where curve starts
    AA = c(el_x - el_a_2, el_x, el_x + el_a_2) #The x-values you need for segments
    AA2 = c(el_x - el_a_2, el_x, el_x + el_a_2)
    pol_ind = c(ceiling(n/2), ceiling(n/2)+2)
    seg_col = c(rep('green3', floor(length(AA)-3)/2),
                rep('navy', 3),
                rep('green3', length(AA) - 3 - floor(length(AA)-3)/2))
    text_cex = c(rep(1.1, floor(length(AA)-3)/2),
                 rep(1.4, 3),
                 rep(1.1, length(AA) - 3 - floor(length(AA)-3)/2) )
    curve_points_col = c(rep('blue', 6), 'red', rep('blue', 6))
    curve_points_bg = c(rep('blue', 6), 'red', rep('blue', 6))
  }

  #Store curve in a variable for now
  cc = if(type == "Cauchy"){
    curve(dcauchy(x, 0, width)*multiplier, ff, max(el_x + el_a), n = 1e4 )
  } else if(type == "Normal") {
    curve(dnorm(x, mean = 0, sd = width)*multiplier, ff, max(el_x + el_a), n = 1e4 )
  } else {
    stop(message(cat("You need to define which prior you want to explore; *Normal* or *Cauchy*")))
  }

  AA = sort(AA) #To plot lines and points as you want

  #Calculate y-values from AA
  BB <- if(type == "Cauchy"){
    dcauchy(AA, 0, width) * multiplier  #The y-values you need for segments
  } else if(type == "Normal") {
    dnorm(AA, mean = 0, sd = width) * multiplier
  } else {
    stop( message( cat("You need to define which prior you want to explore; *Normal* or *Cauchy*")) )
  }

  #Raise the curve to have 0.2 distance (OPTIONAL)
  #Raise_curve = max(el_b + el_y) + 0.8 - (min(BB) - max(el_b + el_y))
  Raise_curve = 0

  BB = BB + Raise_curve #Raise BB

  plot(1, type = 'n', ann = FALSE,
       xlim = c(min(el_x - el_a), max(el_x + el_a)),
       ylim = c(min(el_y - el_b), max(BB)), axes = FALSE)

  axis(side = 2, at = decimal( seq(0, max(BB), length.out = 4), 2), cex.axis = 1, font.axis = 2, las = 1 )

  # mtext("Density", side = 2, font = 2, cex = 1.6, line = 2)
  draw.ellipse(x = rep(el_x, n), y = rep(el_y, n),
               a = el_a, b = el_b,
               lty = 2, border = 'gray20', col = heat.colors(13, alpha = .1))

  xx = seq(from = AA[pol_ind[1]], to = AA[pol_ind[2]], length.out = 1000)

  yy = if(type == "Cauchy") {
    dcauchy(xx, 0, width)*multiplier
  } else {
    dnorm(xx, mean = 0, sd = width)*multiplier
  }

  yy = yy + Raise_curve

  polygon(c(AA[pol_ind[1]], xx, AA[pol_ind[2]]), c(el_y, yy, el_y), border = NA, col= adjustcolor('green', alpha.f = .2) )
  lines(cc$x, cc$y + Raise_curve, lwd = 3, col = "magenta")

  segments(x0 = AA, x1 = AA, y0 = el_y, y1 = BB, lty = 3, lwd = 2, col = seg_col)

  points(AA2, rep(el_y, length(AA2) ), pch = c(rep(22, 6), 21, rep(22, 6)),
         col = c(rep('blue', 6), 'red', rep('blue', 6)), bg = c(rep('blue', 6), 'red', rep('blue', 6)),
         cex = 3)

  text(seq(-6, 6), rep(el_y - .06, length(AA2) ), -6:6, font = 2, cex = 2)

  #print(BB)
  #print(decimal(BB, 2))
  text(AA, BB, labels = decimal(BB, 2), font = 2, cex = text_cex, col = 1, pos = 3, xpd = TRUE )

  points(AA, BB, pch = 21, col = curve_points_col, bg = curve_points_bg, cex = 1)

  legend("topleft", ifelse(type == "Cauchy" & width == sqrt(2)/2,
                           "Recommended Cauchy Prior for L2 Research",
                           ifelse(type == "Cauchy" & width == 1,
                                  "Standard Cauchy Prior",
                                  ifelse(type == "Cauchy" & width != 1 ||
                                           type == "Cauchy" & width != sqrt(2)/2 ||
                                           type == "Cauchy" & width !="wide",
                                         "Cauchy Prior",
                                         ifelse(type == "Normal" & width ==1,
                                                "Standard Normal", "Non-Standard Normal")))),
         lwd = 3, col = "magenta", bty = 'n', text.font = 2, cex = 1.5, inset = -.07,
         xpd = T, x.intersp = .5, y.intersp = 0)

  arrows( .3, (BB[7]+BB[8])/2.3, 1.3, (BB[7]+BB[8])/1.7, code = 2, length = .15, angle = 20, lwd = 2 )

  points( .3, (BB[7]+BB[8])/2.3, pch = 21, bg = 1, cex = 1.1, lwd = 2 )

  text( 3.45, (BB[7]+BB[8])/1.55, "Concentration of Plausible
        Effect Sizes in L2 Research", cex = 1.1, font = 2, xpd = TRUE )

  #text( 3.5, (BB[7]+BB[8])/1.7, "Effect Sizes in L2 Research", cex = 1.1, font = 2, xpd = T)

  text(-7.5, max(BB)/2, bquote(bold(Density)), cex = 1.4, srt = 90, xpd = TRUE )
  mtext(side = 1, "A Range of Reasonable Sizes for a Population Effect in L2 Research", cex = 1.3, font = 2)
  mtext(side = 1, "Informed by (Plonsky & Oswald, 2014)", cex = 1.3, font = 2, line = 1.1)

  par(family="sans")
}

