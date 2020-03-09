rm(list=ls())
indir <- "C:/Users/allorant/OneDrive - UW/Shared with Everyone/UW/3rdYear/Winter/STAT517/Final_project/nonstat_cov/data/"


yr <- 1981
type <- "ppt"
with.infill <- F
fnames<-paste( type,".complete.Y", substring( format( 1001:1103), 2,4), sep="")
fnames<-paste( indir, fnames, sep="")

yr.index<- yr-1894

temp<- scan( fnames[yr.index],what=list( "a", 1,2,3,4,5,6,7,8,9,10,11,12,"a"))

# extract the missing value logical from last 12 ones and zeroes
missing.ind <- unlist(strsplit( temp[[14]],NULL))
missing.ind <- matrix( missing.ind, ncol=12, byrow=T)

# this is complete data with infilled values where real data is missing
temp<- matrix( unlist( temp[2:13]), ncol=12)

# over write with NA's if infills not wanted. 

if( !with.infill) {
  temp[missing.ind==1] <- NA
}


# Set-up problem
require(INLA)

# Read station locations
scan(paste0(indir,"METAinfo"), skip=1, what=list( "a", 1,1,1))-> look
names(look)<-c("station.id", "lon", "lat", "elev")

# Set-up variables needed for INLA
loc <- vector("numeric")
repl <- vector("numeric")
y <- vector("numeric")
elev <- vector("numeric")

# Fill them with data
count = 1
dat <- temp
dat <- rowSums(dat)

# dataframe
totDF <- data.frame(
  y = log(dat),
  lon = look$lon,
  lat = look$lat,
  elev = look$elev/1000,
  idx = ifelse(!is.na(dat), "est", "pred")
)

totDF <- totDF %>%
  mutate(elev = ifelse(elev == -500, 0, elev))



loc <- cbind(totDF$lon[totDF$idx=="est"],totDF$lat[totDF$idx=="est"])

# Transform data

# Make mesh
mesh <- inla.mesh.2d(loc.domain = loc, max.edge = c(0.5, 1))

# Make A matrix
A <- inla.spde.make.A(mesh = mesh, loc = loc)

# Make spde-object
spde <- inla.spde2.matern(mesh)

# Make stack
## Estimation
stk.est <- inla.stack(data = list(y = totDF %>% filter(idx == "est") %>% dplyr::select(y)),
                      A = list(1, A),
                      effects = list(list(intercept=1,
                                          elev = totDF %>% filter(idx == "est") %>% dplyr::select(elev)),
                                     c(inla.spde.make.index("mesh.idx",
                                                            n.spde=mesh$n))), tag="est")
## Prediction for location where log_ann_prec not observed
loc_pred <- cbind(totDF$lon[totDF$idx=="pred"],totDF$lat[totDF$idx=="pred"])

Ap <- inla.spde.make.A(mesh = mesh, loc = loc_pred)
dim(Ap)

stk.p <- inla.stack(data = list(y = NA),
                      A = list(1, Ap),
                      effects = list(list(intercept=1,
                                          elev = totDF %>% filter(idx == "pred") %>% dplyr::select(elev)),
                                     c(inla.spde.make.index("mesh.idx",
                                                            n.spde=mesh$n))), tag="pred")

stk.full <- inla.stack(stk.est, stk.p)
# Formula
formula <- y~intercept+elev+f(mesh.idx, model = spde, replicate=mesh.idx.repl)-1

# Run inla
res <- inla(formula,
            family = "gaussian",
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE,
              A=inla.stack.A(stk.full)),
            verbose = TRUE,
            control.mode=list(theta=c(0.031138, -0.199624, -0.7506),
                              restart=TRUE), num.threads = 2)

save(res, file = paste0(indir, "res.rda"))
load(paste0(indir, "res.rda"))
# 
# library(ggplot2)
# intercept <- res$marginals.fixed[[1]]
# ggplot(data.frame(inla.smarginal(intercept)), aes(x, y)) +
#   geom_line() +
#   theme_bw()
# marg.variance <- inla.tmarginal(function(x) 1/x,
#                                 res$marginals.hyperpar$`Precision for the Gaussian observations`)
# 
# ggplot(data.frame(inla.smarginal(marg.variance)), aes(x, y)) +
#   geom_line() +
#   theme_bw()
# 
# list_marginals <- res$marginals.fitted.values
# 
# marginals <- data.frame(do.call(rbind, list_marginals))
# marginals$h <- rep(names(list_marginals),
#                           times = sapply(list_marginals, nrow))
# res$cpu

summary(res)

index_est <- inla.stack.index(stk.full, tag = "est")$data

est_mean <- res$summary.fitted.values[index_est, "mean"]
est_ll <- res$summary.fitted.values[index_est, "0.025quant"]
est_ul <- res$summary.fitted.values[index_est, "0.975quant"]

index_pred <- inla.stack.index(stk.full, tag = "pred")$data

pred_mean <- res$summary.fitted.values[index_pred, "mean"]
pred_ll <- res$summary.fitted.values[index_pred, "0.025quant"]
pred_ul <- res$summary.fitted.values[index_pred, "0.975quant"]

dpm <- rbind(
  data.frame(
    east = loc[, 1], north = loc[, 2],
    value = est_mean, variable = "est_mean",
    tag = "est"
  ),
  data.frame(
    east = loc[, 1], north = loc[, 2],
    value = est_ll, variable = "est_ll",
    tag = "est"
  ),
  data.frame(
    east = loc[, 1], north = loc[, 2],
    value = est_ul, variable = "est_ul",
    tag = "est"
  ),
  data.frame(
    east = loc_pred[, 1], north = loc_pred[, 2],
    value = pred_mean, variable = "pred_mean",
    tag = "pred"
  ),
  data.frame(
    east = loc_pred[, 1], north = loc_pred[, 2],
    value = pred_ll, variable = "pred_ll",
    tag = "pred"
  ),
  data.frame(
    east = loc_pred[, 1], north = loc_pred[, 2],
    value = pred_ul, variable = "pred_ul",
    tag = "pred"
  )
)

dpm <- dpm %>%
  mutate(variable = as.factor(dpm$variable))


dpm %>% 
  filter( variable == "est_mean" | variable == "pred_mean") %>%
  ggplot(aes(x = east, y = north, color = value)) +
  facet_wrap(~tag, nrow = 1) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_gradient(
    name = "Rainfall",
    low = "blue", high = "orange"
  )

# 
# rang <- apply(mesh$loc[, c(1, 2)], 2, range)
# proj <- inla.mesh.projector(mesh,
#                             xlim = rang[, 1], ylim = rang[, 2],
#                             dims = c(300, 300)
# )
# mean_s <- inla.mesh.project(proj, res$summary.random$s$mean)
# sd_s <- inla.mesh.project(proj, res$summary.random$s$sd)
# 
# df <- expand.grid(x = proj$x, y = proj$y)
# df$mean_s <- as.vector(mean_s)
# df$sd_s <- as.vector(sd_s)
# 
# library(viridis)
# library(cowplot)
# 
# gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
#   geom_raster() +
#   scale_fill_viridis(na.value = "transparent") +
#   coord_fixed(ratio = 1) + theme_bw()
# 
# gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
#   geom_raster() +
#   scale_fill_viridis(na.value = "transparent") +
#   coord_fixed(ratio = 1) + theme_bw()
# 
# plot_grid(gmean, gsd)