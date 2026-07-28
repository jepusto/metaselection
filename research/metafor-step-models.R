library(metafor)

dat <- dat.lehmann2018
res <- rma(yi, vi, data=dat, method="ML")
funnel(res, atransf=exp, at=log(c(0.25,0.5,1,2,4,8)), ylim=c(0,0.8))

steps <- c(0.1, 0.5)
control <- list(optimizer = "nlminb", rel.tol = 1e-10)
# step function selection model: one-sided, HA: beta > 0
sel_gr <- selmodel(
  res, 
  type = "stepfun", alternative = "greater", 
  steps = steps,
  control = control
)
plot(sel_gr)
sel_gr

# step function selection model: one-sided, HA: beta < 0
sel_ls <- selmodel(
  res, 
  type = "stepfun", alternative = "less", 
  steps = 1 - steps,
  control = control
)
plot(sel_ls)
sel_ls

c(sel_gr$beta[1], sel_ls$beta[1])
c(sel_gr$tau2, sel_ls$tau2)
data.frame(gr = sel_gr$delta, ls = rev(sel_ls$delta) / sel_ls$delta[sel_ls$deltas])
c(sel_gr$LRT, sel_ls$LRT)
