context("ecnoise")
library(RCompNoise)

## Model problem from More' and Wild (2012)
set.seed(101)

alpha <- 0.01

n_dim  <- 10
n_fcn  <- 8
h      <- 1e-6
n_repl <- 1e4

eps_rel_true <- 1e-3

xb <- rep(1, n_dim)

v_test_noise <- rep(0, n_repl)
v_test_noise_rel <- rep(0, n_repl)
v_test_result <- rep("", n_repl)

for (ind in 1:n_repl) {
  p <- 2 * runif(n = n_dim) - 1;
  p <- p / norm(p, type = "2")
  test_fval <- rep(0, n_fcn + 1)
  test_R    <- eps_rel_true * rnorm(n = n_fcn + 1)
  mid       <- floor((n_fcn + 2) / 2)

  for (i in 0:n_fcn) {
    s = 2 * (i + 1 - mid) / n_fcn
    x = xb + s * p * h

    test_fval[i + 1] = norm(x, type = "2")^2 * (1 + test_R[i + 1])
  }

  test_res <- ecnoise(test_fval)

  v_test_noise[ind] <- test_res$fnoise
  v_test_noise_rel[ind] <- test_res$fnoise / test_fval[mid]
  v_test_result[ind] <- test_res$inform
}

## More' and Wild state that they expect ECNoise to be accurate only to within a
## factor of 4. They report that 1% of cases fail at 8 function evaluations; we
## will replicate this case, and consider it successful if our 'failure fraction'
## is less than 1%.
ratio_rel <- v_test_noise_rel / eps_rel_true
frac_fail_loose <- 1 - mean( (0.25 <= ratio_rel) & (ratio_rel <= 4) )

test_that(
  "ECNoise failure rate compatible with published values.",
  {
    expect_true(
      frac_fail_loose <= 0.01
    )
  }
)
