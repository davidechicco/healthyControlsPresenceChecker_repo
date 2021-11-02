setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)


source("healthyControlsPresenceChecker.r")
GEO_code <- "GSE19429"
this_outcome <- healthyControlsCheck(GEO_code, TRUE)

source("healthyControlsPresenceChecker.r")
GEO_code <- "GSE34111"
this_outcome <- healthyControlsCheck(GEO_code, TRUE)

source("healthyControlsPresenceChecker.r")
GEO_code <- "GSE47407"
this_outcome <- healthyControlsCheck(GEO_code, TRUE)

source("healthyControlsPresenceChecker.r")
GEO_code <- "GSE30174"
this_outcome <- healthyControlsCheck(GEO_code, TRUE)
