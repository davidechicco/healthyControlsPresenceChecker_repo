setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

args = commandArgs(trailingOnly=TRUE)


cat(":: Installing / loading the R packages ::\n:: required by the script ::\n\n")

# Here we install the CRAN missing packages
list.of.packages <- c("easypackages", "xml2", "markdown", "knitr", "rmarkdown", "pacman", "dplyr") # other packages
new_packages_to_install <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new_packages_to_install)) install.packages(new_packages_to_install, repos="https://utstat.toronto.edu/cran/")

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

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)

inputGEOcode <- NULL
inputVerboseFlag <- FALSE

# test if there is at least one argument: if not, return an error
if (length(args) != 2) {
  stop("At least 2 argument must be supplied\n", call.=FALSE)
} else {

  inputGEOcode <- toString(args[1])
  inputVerboseFlag <- toString(args[2])
  
}

#' Function that returns numeric values with 2 decimal numbers.
#' 
#' @param x input numeric value with N decimal numbers.
#' @return a numeric value with 2 decimal numbers.
#' @examples
#' aaa <- dec_two(8.31232)
 dec_two <- function(x) {
  return (format(round(x, 2), nsmall = 2));
}



#' Function that reads in a URL to check and verifies if it exists (function taken from https://stackoverflow.com/a/12195574 )
#' 
#' @param url the URL of a webpage
#' @return the output of a webpage verification check
#' @examples y <- readUrl("http://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html")
readUrl <- function(url) {
    out <- tryCatch(
        {
            # Just to highlight: if you want to use more than one 
            # R expression in the "try" part then you'll have to 
            # use curly brackets.
            # 'tryCatch()' will return the last evaluated expression 
            # in case the "try" part was completed successfully


            readLines(con=url, warn=FALSE) 
            # The return value of `readLines()` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped inside a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("URL does not seem to exist:", url))
            message("Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            return("EMPTY_STRING")
        },
        warning=function(cond) {
            message(paste("URL caused a warning:", url))
            message(cond)
            # Choose a return value in case of warning
            return("EMPTY_STRING")
        },
        finally={
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        # If you want more than one expression to be executed, then you 
        # need to wrap them in curly brackets ({...}); otherwise you could
        # just have written 'finally=<expression>' 
            message(paste("Processed URL:", url))
        }
    )    
    return(out)
}

#' Function that reads in the GEO code of a dataset, and returns true if there's at least a feature
#' containing the healthy controls 
#'
#' @param datasetGeoCode the GEO code of a dataset.
#' @param verbose a boolean flag stating if helping messages should be printed or not
#' @return a boolean value
#' @examples
#' healthyControlsCheckOutcome <- healthyControlsCheck("GSE3268", FALSE)
healthyControlsCheck <- function(datasetGeoCode, verbose = FALSE) 
{

            if(verbose == TRUE) cat("\n\n> > > healthyControlsCheck() started with ", datasetGeoCode, " as input\n\n", sep="")

            GSE_code <- datasetGeoCode
            
            # check   URL
            checked_html_text <- "EMPTY_STRING"
            checked_html_text <- xml2::read_html("https://ftp.ncbi.nlm.nih.gov/geo/series/")
            
            checked_html_text_url <- "EMPTY_STRING"
            url_to_check <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", datasetGeoCode)
            GSE_code_for_url <- GSE_code
            GSE_code_for_url <- substr(GSE_code_for_url,1,nchar(GSE_code_for_url)-3)
            GSE_code_for_url <- paste0(GSE_code_for_url, "nnn")
            complete_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", GSE_code_for_url, "/", GSE_code)
           
           checked_html_text_url <- lapply(complete_url, readUrl)
           
           gset <- NULL
#             
            if(all(checked_html_text == "EMPTY_STRING")) {
         
                    cat("The web url https://ftp.ncbi.nlm.nih.gov/geo/series/ is unavailable right now. Please try again later. The function will stop here\n")
                    return(NULL)
                    
            } else if(all(checked_html_text_url == "EMPTY_STRING" | is.null(checked_html_text_url[[1]]) )) {
         
                    cat("The web url ", complete_url," is unavailable right now (Error 404 webpage not found). The GEO code might be wrong. The function will stop here\n", sep="")
                    return(NULL)        
                    
            } else {

                gset <- GEOquery::getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)
                
                thisGEOplatform <- toString((gset)[[1]]@annotation)

                if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
                gset <- gset[[idx]]
                
                if(verbose == TRUE) cat("=== === === === === ", GSE_code, " === === === === ===  \n", sep="")
                
               healthyControlWordPresent <- grepl("healthy control", (gset@phenoData@data)) %>% any()
	      if(healthyControlWordPresent == TRUE) {
	      
		       if(verbose == TRUE) cat(":: The keyword \"healthy control\" was found in this dataset annotations (", GSE_code, ") ", sep="")
		       healthy_control_indexes <- which(grepl("healthy control", (gset@phenoData@data)))
		       cat("on ", length(healthy_control_indexes), " feature(s)\n", sep="")
		       
		       countFeatures <- 1
		       
		       for(i in healthy_control_indexes){
                    this_feature <- (gset@phenoData@data)[i] %>% colnames()  
                    if(verbose == TRUE)  cat("\n(", countFeatures, ") \"", this_feature, "\" feature\n", sep="")
                    if(verbose == TRUE) (gset@phenoData@data)[i] %>% table() %>% print()
                    
                    thisFeatureGroups <- (gset@phenoData@data)[i] %>% table()
                    thisFeatureGroupsNames <- thisFeatureGroups %>% names()
                    numGroupsInThisFeature <- thisFeatureGroups  %>% nrow()

                    for(k in 1:length(thisFeatureGroups)) {
                            cat(thisFeatureGroupsNames[k], ": ", sep="")
                            thisFeatureGroupPerc <- thisFeatureGroups[[k]] * 100 / (gset@phenoData@data) %>% nrow()
                            cat("\t", dec_two(thisFeatureGroupPerc), "%\n", sep="")
                    }
                    
                    countFeatures <- countFeatures + 1
		       }
	      } else { 
		      if(verbose == TRUE) cat(":: The keyword \"healthy control\" was NOT found among the annotations of this dataset (", GSE_code, ")\n", sep="") 
	      }            
            }
                        
            if(verbose == TRUE)  cat("=== === === === === === === === === === === ===  \n", sep="")
            
            cat("\nhealthyControlsCheck() call output: were healthy controls found in the ", GEO_code, " dataset? ", healthyControlWordPresent, "\n", sep="")
            return(healthyControlWordPresent)
}   

  GEO_code <- inputGEOcode
  
  if(inputVerboseFlag == "TRUE" || inputVerboseFlag == "True" || inputVerboseFlag == "true") { 
      inputVerboseFlag <-  TRUE
    } else if(inputVerboseFlag == "FALSE" || inputVerboseFlag == "False" || inputVerboseFlag == "false") { 
      inputVerboseFlag <- FALSE
    }

  
  this_outcome <- healthyControlsCheck(GEO_code, inputVerboseFlag)

# GEO_code <- "GSE30174"
# this_outcome <- healthyControlsCheck(GEO_code, TRUE)
# cat("Healthy controls found in the ", GEO_code, ": ", this_outcome, "\n", sep="")
# 
# GEO_code <- "GSE30174"
# this_outcome <- healthyControlsCheck(GEO_code, TRUE)
# cat("Healthy controls found in the ", GEO_code, ": ", this_outcome, "\n", sep="")
# 
# GEO_code <- "GSE19429"
# this_outcome <- healthyControlsCheck(GEO_code, TRUE)
# cat("Healthy controls found in the ", GEO_code, ": ", this_outcome, "\n", sep="")
# 
# GEO_code <- "GSE34111"
# this_outcome <- healthyControlsCheck(GEO_code, TRUE)
# cat("Healthy controls found in the ", GEO_code, ": ", this_outcome, "\n", sep="")
# 
# GEO_code <- "GSE47407"
# this_outcome <- healthyControlsCheck(GEO_code, TRUE)
# cat("Healthy controls found in the ", GEO_code, ": ", this_outcome, "\n", sep="")
