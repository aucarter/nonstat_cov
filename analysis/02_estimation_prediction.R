rm(list=ls())

if (Sys.info()["sysname"] == "Linux") {
  lib <- "/homes/allorant/rlib"
  j_root <- "/home/j/"
  h_root <- "/homes/allorant/"
} else {
  j_root <- "J:/"
  drive <- "L:/"
  h_root <- "H:/"
}

indir <- "C:/Users/allorant/OneDrive - UW/Shared with Everyone/UW/3rdYear/Winter/STAT517/Final_project/nonstat_cov/data/"
 

libs <- c('data.table','ggplot2','dplyr','tidyr','INLA','elevatr')
for(l in libs){
  if(!require(l,character.only = TRUE, quietly = TRUE)){
    message( sprintf('Did not have the required package << %s >> installed. Downloading now ... ',l))
    install.packages(l) 
  }
  library(l, character.only = TRUE, quietly = TRUE)
}

dt <- read.csv(paste0(indir, "prepped_data.csv"))

dt1981 <- dt %>%
  filter(year == 1981)

d <- dt1981 %>%
  group_by(lat, lon) %>%
  dplyr::select(
    lat, lon, elev, log_ann_prec
  )
head(d)

# simple mapping
library(leaflet)
library(viridis)
library(htmlwidgets)
library(webshot)

pal <- colorNumeric("viridis", domain = d$log_ann_prec)

leaflet(d) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~lon, lat = ~lat, color = ~ pal(log_ann_prec)) %>%
  addLegend("bottomright",
            pal = pal, values = ~log_ann_prec,
            title = "Log precipitation"
  ) %>%
  addScaleBar(position = c("bottomleft")) 


##################
## Estimation ###
#################

# Modeling
## Mesh construction
coo <- cbind(d$lon, d$lat)
mesh <- inla.mesh.2d(loc = coo, offset = c(.5, 1),
                     cutoff = .1, max.edge = c(3, 6)
)

# mesh <- inla.mesh.2d(
#   loc = coo, offset = c(50, 100),
#   cutoff = 1, max.edge = c(30, 60)
# )

plot(mesh)
points(coo, col = "red")
mesh$n

spde1 <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

A <- inla.spde.make.A(mesh = mesh, loc = coo)
dim(A)
nrow(coo)
mesh$n
rowSums(A)

# Make stack
stk.est <- inla.stack(data = list(y = d$log_ann_prec),
                      A = list(1, A),
                      effects = list(list(intercept=1, elev = d$elev),
                                     c(inla.spde.make.index("mesh.idx", n.spde=mesh$n))), tag="est")



##############################
## Estimation & Prediction ###
##############################


# We want the domain to match the boundaries of the country
library(maps)
usa <- map_data("usa")
usa <- map_data("usa") # we already did this, but we can do it again
ggplot() + geom_polygon(data = usa, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3)
# library(raster)
# r <- getData(name = "GADM", country = 'USA', level = 1)
# border <- r %>% 
#   fortify() %>%
#   filter(NAME_1 != c("Alaska","Hawaii")) %>%
#   dplyr::select(long,lat)

border <- usa %>%
  dplyr::select(long,lat)
border <- cbind(border$long,border$lat)
bb<-bbox(border)
x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = 100)
y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = 100)
coop <- as.matrix(expand.grid(x, y))

ind <- point.in.polygon(
  coop[, 1], coop[, 2],
  border[, 1], border[, 2]
)
coop <- coop[which(ind == 1), ]
plot(coop, asp = 1)
plot(border)
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
dim(Ap)



saveRDS(df_elev_epqs, file = paste0(indir,"elevation.RDS"))
if(!file.exists(paste0(indir,"elevation.RDS"))){
  prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  df_elev_epqs <- get_elev_point(data.frame(coop), prj = prj_dd, src = "epqs")
  saveRDS(df_elev_epqs, paste0(indir,"elevation.RDS"))
}
df_elev_epqs <- readRDS(paste0(indir,"elevation.RDS"))
# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(elev = df_elev_epqs@data$elevation),
                 c(inla.spde.make.index("mesh.idx", n.spde=mesh$n)))
)

stk.full <- inla.stack(stk.est, stk.p)
# Formula
formula <- y~intercept+elev+f(mesh.idx, model = spde1)-1

res <- inla(formula,
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE,
              A = inla.stack.A(stk.full)
            )
)

index <- inla.stack.index(stk.full, tag = "pred")$data

pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]

dpm <- rbind(
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_mean, variable = "pred_mean"
  ),
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_ll, variable = "pred_ll"
  ),
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_ul, variable = "pred_ul"
  )
)
dpm$variable <- as.factor(dpm$variable)

ggplot(dpm) + geom_tile(aes(east, north, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Rainfall",
    low = "blue", high = "orange"
  ) +
  theme_bw()

rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh,
                            xlim = rang[, 1], ylim = rang[, 2],
                            dims = c(300, 300)
)

mean_s <- inla.mesh.project(proj, res$summary.random$mesh.idx$mean)
sd_s <- inla.mesh.project(proj, res$summary.random$mesh.idx$sd)

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

plot_grid(gmean, gsd)

# 
# 
# border.ll <- SpatialPolygons(list(Polygons(list(
#   Polygon(border)), '0')), 
#   proj4string = CRS("+proj=longlat +datum=WGS84"))
# border <- spTransform(border.ll, 
#                       CRS("+proj=utm +units=km +zone=10 +south"))
# 
# apply(bbox(border), 1, diff)
# prange <- c(100, 0.5)
# bprior <- list(prior = 'gaussian', param = c(0,1))
# 
# loc.ll <- SpatialPoints(d[,1:2], border.ll@proj4string) 
# loc <- spTransform(loc.ll, border@proj4string)
# mesh <- inla.mesh.2d(loc, max.edge = 200, cutoff = 35,
#                      offset = 150)