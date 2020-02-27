### Pull precipitation data from https://www.image.ucar.edu/public/Data/RData.USmonthlyMet.bin
### and clean it for analysis

## Setup
library(data.table)

## Read in and inspect
attach("data/RData.USMonthlyMet.bin")

## Prep the precipitation array
USppt
str(USppt)
# Subset to year 1981
p81 <- USppt[,,]
# Remove observations that are NA for each year
no_na <- function(x) {
    return(!any(is.na(x)))
}
full_idx <- lapply(1:dim(USppt)[1], function(i) {which(apply(USppt[i,,], 2, no_na))})
# Check that 1981 has 7040 full observations
length(full_idx[[1981 - 1895 + 1]])

## Attach log annual precip to location and elevation info
dt <- rbindlist(lapply(1:length(full_idx), function(i) {
    log_ann_prec <-log(colSums(USppt[i,,full_idx[[i]]]))
    cbind(year = 1894 + i, USpinfo[full_idx[[i]], c("elev", "lon", "lat")], log_ann_prec)
}))

write.csv(dt, "data/prepped_data.csv", row.names = F)
