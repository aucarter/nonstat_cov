### Pull precipitation data from https://www.image.ucar.edu/public/Data/RData.USmonthlyMet.bin
### and clean it for analysis

## Setup
library(data.table); library(ggplot2)

## Read in and inspect
attach("data/RData.USMonthlyMet.bin")

## Prep the precipitation array
USppt
str(USppt)
# Remove observations that are NA for each year
no_na <- function(x) {
    return(!any(is.na(x)))
}
full_idx <- lapply(1:dim(USppt)[1], function(i) {which(apply(USppt[i,,], 2, no_na))})
# Check that 1981 has 7040 full observations
length(full_idx[[1981 - 1895 + 1]])

## Attach annual precip to location and elevation info
dt <- rbindlist(lapply(1:length(full_idx), function(i) {
    ann_prec <- colSums(USppt[i,,full_idx[[i]]])
    cbind(year = 1894 + i, USpinfo[full_idx[[i]], c("elev", "lon", "lat")], ann_prec)
}))
dt[, log_ann_prec := log(ann_prec)]

## Save
write.csv(dt, "data/prepped_data.csv", row.names = F)

## Plot 1981
ggplot(dt[year == 1981], aes(x = lon, y = lat, color = log_ann_prec)) + 
    geom_point(size = 0.1) + coord_equal() + theme_classic() +
    scale_color_distiller(palette = "Spectral", limits = c(min(dt[year == 1981]$log_ann_prec), 
                                                           max(dt[year == 1981]$log_ann_prec))) 
