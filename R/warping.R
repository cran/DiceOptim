warp <- function(y, wparam) {
  if (length(wparam) == 2) {
    m <- wparam[1]
    s <- wparam[2]
    y <- (y - m)/s
  } else {
    la <- (length(wparam)-2)/3
    a <- wparam[1:la]
    b <- wparam[la + 1:la]
    c <- wparam[2*la + 1:la]
    m <- wparam[length(wparam)-1]
    s <- wparam[length(wparam)]
    y <- (y - m)/s
    f <- y
    for (i in 1:length(a)) f <- f + a[i]*tanh(b[i]*(y+c[i]))
    return(f/sum(a))
  }
}

logLw <- function(x, y) {
  la <- length(x)/3
  a <- x[1:la]
  b <- x[la + 1:la]
  c <- x[2*la + 1:la]
  jacob <- rep(1, length(y))
  y <- (y - mean(y))/sd(y)
  
  for (i in 1:length(a)) {
    jacob <- jacob + a[i]*b[i]*(1- tanh(b[i]*(y+c[i]))^2)
  }
  return(-sum(log(jacob)))
}