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
 

libs <- c('data.table','ggplot2','dplyr','tidyr','INLA')
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

# check to make sure each of the clusters are unique for lat long 
nrow(unique(dplyr::select(totalDF, latnum, longnum)))

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
mesh <- inla.mesh.2d(loc = coo, offset = c(50, 100),
                     cutoff = 1, max.edge = c(30, 60)
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

# Formula
formula <- y~intercept+elev+f(mesh.idx, model = spde1)-1

res <- inla(formula,
            data = inla.stack.data(stk.est),
            control.predictor = list(
              compute = TRUE,
              A = inla.stack.A(stk.est)
            )
)

pred_mean <- res$summary.fitted.values[, "mean"]
pred_ll <- res$summary.fitted.values[, "0.025quant"]
pred_ul <- res$summary.fitted.values[, "0.975quant"]


# Make stack
stk.est <- inla.stack(data = list(y = d$log_ann_prec),
                      A = list(1, A),
                      effects = list(list(intercept=1, elev = d$elev),
                                     c(inla.spde.make.index("mesh.idx", n.spde=mesh$n))), tag="est")

# Run inla
res <- inla(formula, family = "gaussian",
            data = inla.stack.data(stk.est), 
            control.predictor = list(A=inla.stack.A(stk.est)), 
            verbose = TRUE, 
            control.mode=list(theta=c(4.2366, 9.8822, -0.7506),restart=TRUE),
            num.threads = 2)

# We want the domain to match the boundaries of the country

library(raster)
r <- getData(name = "elevation", country = 'USA', level = 0)
border <- r %>% 
  fortify() %>%
  dplyr::select(long,lat)

border <- cbind(border$long,border$lat)
bb<-bbox(border)
x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = 50)
y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = 50)
coop <- as.matrix(expand.grid(x, y))

ind <- point.in.polygon(
  coop[, 1], coop[, 2],
  border[, 1], border[, 2]
)
coop <- coop[which(ind == 1), ]
plot(coop, asp = 1)



border.ll <- SpatialPolygons(list(Polygons(list(
  Polygon(border)), '0')), 
  proj4string = CRS("+proj=longlat +datum=WGS84"))
border <- spTransform(border.ll, 
                      CRS("+proj=utm +units=km +zone=10 +south"))

apply(bbox(border), 1, diff)
prange <- c(100, 0.5)
bprior <- list(prior = 'gaussian', param = c(0,1))

loc.ll <- SpatialPoints(d[,1:2], border.ll@proj4string) 
loc <- spTransform(loc.ll, border@proj4string)
mesh <- inla.mesh.2d(loc, max.edge = 200, cutoff = 35,
                     offset = 150)