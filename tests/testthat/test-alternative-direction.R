skip_if_not_installed("metadat")
library(metadat)

dat <- dat.hackshaw1998
dat$yi_neg <- -1 * dat$yi
dat$sei <- sqrt(dat$vi)

test_that("step models are consistent when alternative = 'less'.", {

  steps <- c(.025, .50)
  
  pos_gt <- selection_model(
    data = dat, 
    yi = yi, sei = sei,
    steps = steps,
    priors = NULL
  )

  pos_ls <- selection_model(
    data = dat, 
    yi = yi, sei = sei,
    alternative = "less",
    steps = rev(1 - steps),
    priors = NULL
  )

  neg_gt <- selection_model(
    data = dat, 
    yi = yi_neg, sei = sei,
    steps = rev(1 - steps),
    priors = NULL
  )

  neg_ls <- selection_model(
    data = dat, 
    yi = yi_neg, sei = sei,
    alternative = "less",
    steps = steps,
    priors = NULL
  )
  
})

