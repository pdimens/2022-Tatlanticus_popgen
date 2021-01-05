k <- (1:10)
lhood <- c(
  -0.1479313291,
  -0.1507955207,
  -0.1640825829,
  -0.1666774728,
  -0.1720553811,
  -0.1774462667,
  -0.1767752362,
  -0.1723598764,
  -0.1692276587,
  -0.1661562606
  )
plot(lhood, pch = 19, col = "black", xlab = "K clusters", ylab = "Marginal Likelihood", main = "Marginal Likelihood of K clusters in the sample data")
