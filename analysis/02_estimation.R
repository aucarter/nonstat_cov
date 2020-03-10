indir <- "C:/Users/allorant/OneDrive - UW/Shared with Everyone/UW/3rdYear/Winter/STAT517/Final_project/nonstat_cov/data/"
with.infill <- FALSE
yr <- 1981
type <- "ppt"
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
idx <- !is.na(dat)

	# Remove unobserved
dat <- dat[idx]
lon <- look$lon[idx];
lat <- look$lat[idx];

	# Select data to use for this year
y <- c(y, log(dat))

	# Select locations to use for this ye
loc <- rbind(loc, cbind(lon, lat))

	# Set replication vector
repl <- c(repl, rep(count, length(lon)))

	# Get elevation covariate
tmpEl <- look$elev[idx]/1000
elev <- c(elev, tmpEl)

	
# Transform data

# Make mesh
mesh <- inla.mesh.2d(loc.domain = loc, max.edge = c(0.5, 1))

# Make A matrix
A <- inla.spde.make.A(mesh = mesh, loc = loc)

# Make spde-object
spde <- inla.spde2.matern(mesh)

# Make stack
stk.est <- inla.stack(data = list(y = y),
                      A = list(1, A),
                      effects = list(list(intercept=1,
                                          elev = elev),
                                     c(inla.spde.make.index("mesh.idx",
                                                            n.spde=mesh$n))), tag="est")

# Formula
formula <- y~intercept+elev+f(mesh.idx, model = spde, replicate=mesh.idx.repl)-1

# Run inla
res <- inla(formula,
            family = "gaussian",
            data = inla.stack.data(stk.est),
            control.predictor = list(A=inla.stack.A(stk.est)),
            verbose = TRUE,
            control.mode=list(theta=c(0.031138, -0.199624, -0.7506),
                              restart=TRUE), num.threads = 2)
saveRDS(res, paste0(indir, "resINLA.RDS"))
# Extract mean at desired locations for the random field
predMean <- drop(A%*%res$summary.random$mesh.idx$mean) # +res$summary.fixed[1]
predSd <- drop(A%*%res$summary.random$mesh.idx$sd)

df <- data.frame(x = lon, y = lat)
df$mean_s <- as.vector(predMean)
df$sd_s <- as.vector(predSd)

library(viridis)
library(cowplot)

gmean <- ggplot(df, aes(x = lon, y = lat, fill = mean_s)) +
  geom_tile() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
  geom_tile() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

plot_grid(gmean, gsd)

# Spatial effect
sEffect <- A%*%res$summary.random$mesh.idx$mean

# Fixed effects
fEffect <- res$summary.fixed[1,1]+res$summary.fixed[1,2]*elev

# Predict it all!
Mean <- sEffect@x + fEffect

dpm <- rbind(
  data.frame(
    east = loc[, 1], north = loc[, 2],
    value = Mean, variable = "mean"
  ),
  data.frame(
    east = loc[, 1], north = loc[, 2],
    value = sEffect@x, variable = "spatial effect"
  ),
  data.frame(
    east = loc[, 1], north = loc[, 2],
    value = fEffect, variable = "fixed effect"
  )
)
dpm$variable <- as.factor(dpm$variable)

ggplot(dpm) + geom_tile(aes(east, north, fill = value)) +
  # facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Rainfall",
    low = "blue", high = "orange"
  ) +
  theme_bw()
saveRDS(dpm, paste0(indir, "meanEstEffect.RDS"))

