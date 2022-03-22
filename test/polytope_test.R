n = 1000

muf <- function(d) rep(0, d)
Sigmaf <- function(d) diag(d)
lbf <- function(d) rep(-Inf, 2*d)
ubf <- function(d) c(0, rep(2, 2*d-1))
Af <- function(d) {
  lower_bounds <- -diag(d)
  upper_bounds <- diag(d)
  upper_bounds[1, ] <- c(2, 1, rep(0, d-2))
  A <- rbind(upper_bounds, lower_bounds)
  return(A)
}

d = 2
d = 10

mu = muf(d)
Sigma = Sigmaf(d)
lb = lbf(d)
ub = ubf(d)
A = Af(d)

qr_tA = qr(t(A))

R = qr.R(qr_tA)
U = qr.Q(qr_tA)
L = t(R)
Q = t(U)

qr_tA$qr[!upper.tri(qr_tA$qr)]

upper.tri(qr_tA$qr)
qr_tA$qr[upper.tri(qr_tA$qr)]


idx = which(arr.ind = TRUE, upper.tri(qr_tA$qr))


foo = volesti::gen_cross(d, "H")
qr_tA = qr(t(foo$A))

R = qr.R(qr_tA)
U = qr.Q(qr_tA)
L = t(R)
Q = t(U)
