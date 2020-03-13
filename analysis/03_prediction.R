rm(list=ls())
indir <- "C:/Users/allorant/OneDrive - UW/Shared with Everyone/UW/3rdYear/Winter/STAT517/Final_project/nonstat_cov/data/"
outdir <- "C:/Users/allorant/OneDrive - UW/Shared with Everyone/UW/3rdYear/Winter/STAT517/Final_project/output/"

libs <- c('data.table','ggplot2','dplyr','INLA')
for(l in libs){
  if(!require(l,character.only = TRUE, quietly = TRUE)){
    message( sprintf('Did not have the required package << %s >> installed. Downloading now ... ',l))
    install.packages(l) 
  }
  library(l, character.only = TRUE, quietly = TRUE)
}

#########################
## Extracting the data ##
#########################

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

# Assembling them together in a dataframe
totDF <- data.frame(
  y = log(dat),
  lon = look$lon,
  lat = look$lat,
  elev = look$elev/1000,
  idx = ifelse(!is.na(dat), "est", "pred") # the location with missing y are given the tag prediction
)

# p1 <- totDF %>% 
#   ggplot(aes(x = lon, y = lat, color = y), xlim = c(20,55), ylim = c(-130, -60)) +
#   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   geom_point(size = .5, alpha = 0.5) +
#   scale_colour_gradient(
#     name = "Rain (mm)",
#     low = "dodgerblue4",
#     # mid = "olivedrab4",
#     high = "gold1"
#   )+
#   xlab("Longitude") +
#   ylab("Latitude") +
#   theme_bw()
# ggsave(p1, filename =  paste0(outdir, "Desc2.png"))

totDF <- totDF %>%
  mutate(elev = ifelse(elev == -500, 0, elev)) # following the authors' approach we force negative elevation to be zero


#############################
## Estimation & Prediction ##
#############################

# We create two stacks of data, one for which both the predictor and the outcome are observed
# another stack for which we only have the elevation and we want to predict precipitation

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

# We combine these two stacks in one object

stk.full <- inla.stack(stk.est, stk.p)

# Formula
## We have a model with a fixed effect and a spatial effect
formula <- y~intercept+elev+f(mesh.idx, model = spde, replicate=mesh.idx.repl)-1

# Run inla
res <- inla(formula,
            family = "gaussian", # we modelled the data as gaussian as 
                                #in the paper even though this theoretically allows
                                #for negative precipitations
                                # other studies have used a gamma distribution
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE,
              A=inla.stack.A(stk.full)),
            verbose = TRUE,
            control.mode=list(theta=c(0.031138, -0.199624, -0.7506),
                              restart=TRUE), num.threads = 2)

# save(res, file = paste0(outdir, "res.rda"))
load(paste0(outdir, "res.rda"))

summary(res)

library(ggplot2)
alpha <- res$marginals.fixed[[1]]
ggplot(data.frame(inla.smarginal(alpha)), aes(x, y)) +
  geom_line() +
  theme_bw()

library(ggplot2)
alpha <- res$marginals.fixed[[2]]
p1 <- ggplot(data.frame(inla.smarginal(alpha)), aes(x, y)) +
  geom_line() +
  xlab("Elevation") +
  theme_bw()
ggsave(p1, filename =  paste0(outdir, "elevation.png"))


marg.variance <- inla.tmarginal(function(x) 1/x,
                                res$marginals.hyperpar$`Precision for the Gaussian observations`)
p1 <- ggplot(data.frame(inla.smarginal(marg.variance)), aes(x, y)) +
  geom_line() +
  xlab("Marginal variance") +
  theme_bw()
ggsave(p1, filename =  paste0(outdir, "Precision.png"))

# We save in separate objects
## The mean, 2.5% and 97.5% estimated values at observed locations
est_mean <- res$summary.fitted.values[index_est, "mean"]
est_ll <- res$summary.fitted.values[index_est, "0.025quant"]
est_ul <- res$summary.fitted.values[index_est, "0.975quant"]

index_pred <- inla.stack.index(stk.full, tag = "pred")$data

## The mean, 2.5% and 97.5% predicted values at unobserved locations
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
    value = pred_mean, variable = "Mean",
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


# p <- dpm %>% 
#   filter( variable == "est_mean" | variable == "pred_mean") %>%
#   ggplot(aes(x = east, y = north, color = value)) +
#   facet_wrap(~tag, nrow = 1) +
#   geom_point(size = 1, alpha = 0.5) +
#   scale_color_gradient(
#     name = "Rain",
#     low = "blue", high = "yellow"
#   ) +
#   theme_bw()
# 
# ggsave(p, filename =  paste0(outdir, "pred.png"))

# Plot of the predictions
levels(dpm$variable)[levels(dpm$variable)=="pred_mean"] <- "Mean"
levels(dpm$variable)[levels(dpm$variable)=="pred_ll"] <- "2.5%"
levels(dpm$variable)[levels(dpm$variable)=="pred_ul"] <- "97.5%"
theme_set(theme_minimal())
p1 <- dpm %>% 
  filter( tag == "pred") %>%
  ggplot(aes(x = east, y = north, color = value), xlim = c(20,55), ylim = c(-130, -60)) +
  facet_wrap(~variable, nrow = 2) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  geom_point(size = .5, alpha = 0.5) +
  scale_colour_gradient(
    name = "Rain (mm)",
    low = "dodgerblue4",
    # mid = "olivedrab4",
    high = "gold1"
  )+
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()
ggsave(p1, filename =  paste0(outdir, "predAndCI.png"))

# p <- dpm %>% 
#   filter( variable == "pred_mean") %>%
#   ggplot(aes(x = east, y = north, color = value)) +
#   facet_wrap(~variable, nrow = 1) +
#   geom_point(size = 1, alpha = 0.5) +
#   scale_color_gradient(
#     name = "Rainfall",
#     low = "blue", high = "orange"
#   )
# ggsave(p, filename =paste0(outdir,"predPlot.png"))

## We may also be interested in predicting the random field alone
rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh,
                            xlim = rang[, 1], ylim = rang[, 2],
                            dims = c(300, 300)
)
mean_s <- inla.mesh.project(proj, res$summary.random$mesh.idx$mean)
sd_s <- inla.mesh.project(proj, res$summary.random$mesh.idx$sd)
# 
df <- expand.grid(x = proj$x, y = proj$y)
df$mean_s <- as.vector(mean_s)
df$sd_s <- as.vector(sd_s)

library(viridis)
library(cowplot)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

spField <- plot_grid(gmean, gsd)
ggsave(gmean, filename =  paste0(outdir, "spatialFieldM.png"))
ggsave(gsd, filename =  paste0(outdir, "spatialFieldSd.png"))
