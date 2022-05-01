# TruncatedNormal example
Sig <- matrix(c(1,0.9,0.9,1), 2, 2)
mu <- c(-3,0) 
l <- c(-Inf,-Inf)

u <- c(-6,Inf)
A <- matrix(c(1,0,-1,1),2,2)
n <- 1e3; # number of sampled vectors
Y <- met::mvrandn(l - A %*% mu, u - A %*% mu, A %*% Sig %*% t(A), n)
X <- rep(mu, n) + solve(A, diag(2)) %*% Y;

plot(t(X))
 
 
n = 1000
Sigma = .5 * diag(2) + .5 * rep(1, 2) %*% t(rep(1, 2))
mu = rep(0, 2)
lb = rep(0, 2)
ub = rep(Inf, 2)
y = met::mvrandn(lb, ub, Sigma, 1000)
plot(t(y))

x = met::rtmvn(n, mu, Sigma, lb, ub)
plot(x)

met::pmvnorm(mu, Sigma, lb, ub)

met::mvNqmc(lb, ub, Sigma)

d = 3
lb = rep(0, d)
ub = rep(Inf, d)
mu = rep(0, d)
Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
met::mvNqmc(lb, ub, Sigma)

met::pmvn(mu, Sigma, lb, ub)
TruncatedNormal::pmvnorm(mu, Sigma, lb, ub)



TruncatedNormal::mvNqmc(l, u, Sig, 1e4)



# identity rectangle ------------------------------------------------------
d = 1000
lb = rep(-3, d)
ub = rep(3, d)
mu = rep(0, d)
Sigma = diag(d)
n = 5000
samples = met::rtmvn(n, mu, Sigma, lb, ub)

summary(samples)

plot(samples[, c(3, 4)])


# high dimensional check --------------------------------------------------
d = 1000
lb = rep(0, d)
ub = rep(Inf, d)
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1,d))
mu = rep(0, d)

tictoc::tic()
met::pmvn(mu, Sigma, lb, ub, n = 1000, n_est = 12)
tictoc::toc()

tictoc::tic()
TruncatedNormal::pmvnorm(mu, Sigma, lb, ub, B = 12000, type = "qmc")
tictoc::toc()

# even higher dimension --------------------------------------------------
d = 2000
lb = rep(0, d)
ub = rep(Inf, d)
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1,d))
mu = rep(0, d)

tictoc::tic()
met::pmvn(mu, Sigma, lb, ub, n = 1000, n_est = 10)
tictoc::toc()



