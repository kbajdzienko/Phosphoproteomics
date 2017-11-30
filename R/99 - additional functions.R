#Additional functions

# Convert zeros to NAs
zero.to.na <- function(x) {
  x[x==0] <- NA
  return(x)
}
# Convert empty cells to NAs
empty.to.na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}
# Convert NAs to zeros
na.to.zero <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}

# Function returning TRUE for elements repeating more than once
dupl <- function(x) duplicated(x) | duplicated(x, fromLast = T)

# Check if data are log-transformed
# (should work for MS data, for bases >2 for sure)
is.log <- function(intensity) all(intensity < 1000, na.rm = T)

#Function to format numbers on y scale
fmt_scale <- function(x){format(x,nsmall = 1,scientific = T, digits=2)}
#Function to force equal numbers of breaks on the plot scale
equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}
