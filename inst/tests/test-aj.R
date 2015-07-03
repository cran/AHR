test_that("Aalen-Johansen estimator reduces to Kaplan-Meier estimator for two-state model without recovery", {
    test.data <- function(n) {
        T <- rexp(n, 0.1)
        C <- runif(n, 0, 10)
        X <- pmin(T, C)
        D <- T <= C
        V <- runif(n, 0, X/4) 

        status <- T <= C
        D[D == 0] <- "cens"
        
        data.frame(time=X, from=0, to=D, id=1:n, status=status)
    }
    
    data <- test.data(100)
    tra <- matrix(FALSE, nrow=2, ncol=2)
    tra[1, 2] <- TRUE
              
    times <- seq(0, 5, length.out=10) ##c(1.5, 3)

    data <- data[order(data$Y),]
  
    ##S <- exp(-times/10)

    fs <- survfit(Surv(time, status) ~ 1, data=data)

    ## estimate cumulative incidence function for event type 1
    fit <- aj(sort(data$time), data, list(target="0 1", states=c("0", "1"), transitions=tra,
                                          censoring="cens", s=0, t="last", covariance=TRUE))

    
    f <- approxfun(fs$time, fs$surv, method="constant", yleft=1, rule=2, f=0)
    g <- approxfun(fs$time, 100 * (fs$surv * fs$std.err)^2, method="constant", yleft=0, rule=2, f=0)

    expect_true(all.equal(f(times), fit$S))
    expect_true(all.equal(g(times), fit$V))
})

test_that("variance and covariance calculations match (competing risks)", {
              T <- rexp(100)
              C <- rexp(100)
              r <- rbinom(100, 2, 0.5)
              r[(r == 0) | (T > C)] <- "cens"
              data <- data.frame(id=1:100, time=pmin(T,C), from=rep(0, 100), to=r)
              data <- data[order(data$time),]
              tra <- matrix(FALSE, nrow=3, ncol=3)
              tra[1, 2:3] <- TRUE
              
              ## estimate cumulative incidence function for event type 1
              fit <- aj(sort(data$time), data, list(target="0 1", states=c("0", "1", "2"), transitions=tra,
                                                     censoring="cens", s=0, t="last", covariance=TRUE))

              expect_true(all.equal(fit$S^2 * diag(fit$logCOV), fit$V))
})
