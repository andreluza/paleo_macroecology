
closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]}



# create the function to get the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# bin age function
# from here
# https://psmits.github.io/paleo_book/managing-and-processing-data-from-the-paleobiology-database.html

#' Break time data up into bins
#' 
#' Have fun with this. Basic rules. Greater than equal to base, less than top.
#' 
#' @param x vector of ages
#' @param by bin width
#' @param age logical bin age returned, not number (default FALSE, return bin number)
#' @return vector of bin memberships
#' @author Peter D Smits <peterdavidsmits@gmail.com>
bin_ages <- function(x, by = NULL, number = NULL, age = FALSE) {
  
  if(is.null(by) & is.null(number)) {
    return('no scheme given. Specify either bin width or number of bins.')
  }
  
  if(!is.null(by) & !is.null(number)) {
    return('too much information. Specify either bin width OR number of bins, not both.')
  }
  
  # range to bin
  top <- ceiling(max(x))
  bot <- floor(min(x))
  
  # create bins
  if(!is.null(by)) {
    unt <- seq(from = bot, to = top, by = by)
  } else if(!is.null(number)) {
    unt <- seq(from = bot, to = top, length.out = number + 1)
  }
  
  # bin top and bottom
  unt1 <- unt[-length(unt)]
  unt2 <- unt[-1]
  
  # assign memberships
  uu <- map2(unt1, unt2, ~ which(between(x, left = .x, right = .y)))
  
  # what if we want the "age" of the bin, not number?
  if(age == TRUE) {
    unt_age <- map2_dbl(unt1, unt2, ~ median(c(.x, .y)))
  }
  
  # create output vector
  y <- x
  for(ii in seq(length(uu))) {
    if(age == FALSE) {
      y[uu[[ii]]] <- ii
    } else if(age == TRUE) {
      y[uu[[ii]]] <- unt_age[ii]
    }
  }
  y
}
