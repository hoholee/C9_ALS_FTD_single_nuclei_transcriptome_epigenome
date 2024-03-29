# logging which server was used to run this script
system("echo \"Running this on ### $(hostname) ### \\(* w *)/\"")

# load library
library(scDist)
set.seed(1126490984)

simCellType <- function(D, tau, nn = 100, G = 1000, N1 = 5, N2 = 5, J = 30, label = "A") {
  beta_true <- rep(0, G)
  beta_true[sample(1:G, size = nn)] <- uniformly::runif_on_sphere(n = 1, d = nn, r = D)[1, ]
  y <- matrix(0, nrow = (N1 + N2) * J, ncol = G)
  for (i in 1:G) {
    cntr <- 1
    for (j in 1:N1) {
      omega <- rnorm(1, mean = 0, sd = tau)
      y[cntr:(cntr + J - 1), i] <- omega + rnorm(J)
      cntr <- cntr + J
    }
    for (j in 1:N2) {
      omega <- rnorm(1, mean = 0, sd = tau)
      y[cntr:(cntr + J - 1), i] <- beta_true[i] + omega + rnorm(J)
      cntr <- cntr + J
    }
  }

  response <- rep(0, (N1 + N2) * J)
  response[1:(N1 * J)] <- 1

  samples <- c()
  for (i in 1:(N1 + N2)) {
    samples <- c(samples, rep(i, J))
  }

  meta.data <- data.frame(response = response, patient = as.factor(samples), clusters = label)
  out <- list()
  out$Y <- t(y)
  out$meta.data <- meta.data
  return(out)
}

simData <- function(nct = 10, J = 50, N1, N2, G = 1000, nn = 100, tau = 0.5) {
  Y <- matrix(0, nrow = G, ncol = 0)
  meta.data <- data.frame(
    response = NULL,
    patient = NULL,
    clusters = NULL
  )
  D.true <- rep(0, nct)
  for (i in 1:nct) {
    D.true[i] <- rexp(n = 1, rate = 0.05)
    out <- simCellType(D = D.true[i], J = rpois(n = 1, lambda = J), N1 = N1, N2 = N2, label = letters[i], tau = tau)
    Y <- cbind(Y, out$Y)
    meta.data <- rbind(meta.data, out$meta.data)
  }

  out$Y <- Y
  out$meta.data <- meta.data
  out$D.true <- D.true
  return(out)
}

data <- simData(N1 = 5, N2 = 5)

head(data$meta.data)

out <- scDist(data$Y, data$meta.data,
  fixed.effects = "response",
  random.effects = "patient",
  clusters = "clusters"
)

out$results

names(data$D.true) <- letters[1:length(data$D.true)]
data$D.true

DistPlot(out)


## log sessionInfo
sessionInfo()
