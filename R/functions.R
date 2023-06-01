
closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]}



# create the function to get the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
