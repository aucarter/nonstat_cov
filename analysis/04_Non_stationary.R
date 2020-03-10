nu <- 1 
alpha <- nu + 2 / 2
# log(kappa)
logkappa0 <- log(8 * nu) / 2
# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
logtau0 <- logtau0 - logkappa0
# SPDE model
spde <- inla.spde2.matern(mesh, 
                          B.tau = cbind(logtau0, -1, nu, nu * (mesh$loc[,1] - 5) / 10), 
                          B.kappa = cbind(logkappa0, 0, -1, -1 * (mesh$loc[,1] - 5) / 10),
                          theta.prior.mean = rep(0, 3), 
                          theta.prior.prec = rep(1, 3)) 
theta1 <- c(-1, 2, -1)
theta2 <- c(-1, 2, 1)

Q1 <- inla.spde2.precision(spde, theta = theta1)
Q2 <- inla.spde2.precision(spde, theta = theta2)

sample1 <- as.vector(inla.qsample(1, Q1, seed = 1))
sample2 <- as.vector(inla.qsample(1, Q2, seed = 1))

clik <- list(hyper=list(theta=list(initial=20, fixed=TRUE)))
projloc <- inla.mesh.projector(mesh, loc)

x1 <- inla.mesh.project(projloc, sample1)
x2 <- inla.mesh.project(projloc, sample2)

stk1 <- inla.stack(list(y=x1), A=list(projloc$proj$A), tag='d',
                   effects=list(data.frame(i=1:mesh$n)))
stk2 <- inla.stack(list(y=x2), A=list(projloc$proj$A), tag='d',
                   effects=list(data.frame(i=1:mesh$n)))
formula <- y ~ 0 + f(i, model=spde)
fit1 <- inla(formula, control.family=clik,
             data=data.frame(y=sample1, i=1:mesh$n))
fit2 <- inla(formula, control.family=clik,
             data=data.frame(y=sample2, i=1:mesh$n))

x1.mean <- inla.mesh.project(proj, field=res1$summary.ran$i$mean)
x1.var <- inla.mesh.project(proj, field=res1$summary.ran$i$sd^2)
x2.mean <- inla.mesh.project(proj, field=res2$summary.ran$i$mean)
x2.var <- inla.mesh.project(proj, field=res2$summary.ran$i$sd^2)

do.call(function(...) grid.arrange(..., nrow=2),
        lapply(list(inla.mesh.project(proj, sample1), x1.mean, x1.var,
                    inla.mesh.project(proj, sample2), x2.mean, x2.var),
               levelplot, xlab='', ylab='',
               col.regions=topo.colors(100), scale=list(draw=FALSE)))