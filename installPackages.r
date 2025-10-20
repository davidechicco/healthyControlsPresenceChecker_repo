setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

## Script that installs the R packages requested by the getGeneExpressionFromGEO() function

cat(":: Installation of the R packages required ::\n\n")

# Here we install the CRAN missing packages
list.of.packages <- c("xml2", "markdown", "knitr", "rmarkdown", "pacman", "dplyr", "geneExpressionFromGEO", "pacman") # other packages
new_packages_to_install <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new_packages_to_install)) install.packages(new_packages_to_install, repos="https://cran.mirror.garr.it/CRAN/")

# Here we install the Bioconductor missing packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("annotate", "Biobase", "GEOquery")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
if(length(bioCpackagesNotInstalled)) cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("pacman")
if(length(list.of.packages) >= 1) p_load(list.of.packages, character.only = TRUE)
if(length(listOfBiocPackages) >= 1) p_load(listOfBiocPackages, character.only = TRUE)

