setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

list.of.packages <- c("easypackages", "xml2", "markdown", "knitr", "rmarkdown", "pacman", "dplyr", "geneExpressionFromGEO") # other packages
listOfBiocPackages <- c("annotate", "Biobase", "GEOquery")
library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)

#' Function that reads in the GEO code of a dataset, and returns true if there's at least a feature
#' containing the healthy controls 
#'
#' @param datasetGeoCode the GEO code of a dataset.
#' @param verbose a boolean flag stating if helping messages should be printed or not
#' @return a boolean value
#' @examples
#' healthyControlsCheckOutcome1 <- healthyControlsCheck("GSE3268", FALSE)
#' healthyControlsCheckOutcome2 <- healthyControlsCheck("GSE19429", TRUE)
#' healthyControlsCheckOutcome3 <- healthyControlsCheck("GSE34111", FALSE)
#' healthyControlsCheckOutcome4 <- healthyControlsCheck("GSE47407", TRUE)
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
           
           checked_html_text_url <- lapply(complete_url, geneExpressionFromGEO::readUrl)
           
           gset <- NULL
             
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
                
                healthyKeywords <- "healthy|Healthy"
                healthyWordPresent <- grepl(healthyKeywords, (gset@phenoData@data)) %>% any()
                if(healthyWordPresent == TRUE) {
                
                    if(verbose == TRUE) cat(":: The keyword \"healthy\" was found in this dataset annotations (", GSE_code, ")\n", sep="")
                    healthy_indexes <- which(grepl(healthyKeywords, (gset@phenoData@data)))
		            cat("on ", length(healthy_indexes), " feature(s)\n", sep="")
                    
                    countFeatures <- 1
                    for(i in healthy_indexes){
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
                    if(verbose == TRUE) cat(":: The keyword \"healthy\" was NOT found among the annotations of this dataset (", GSE_code, ")\n", sep="") 
                }  
	      
	      controlKeywords <- "control|Control|controls|Controls"
	      healthyControlWordPresent <- grepl(controlKeywords, (gset@phenoData@data)) %>% any()
	      if(healthyControlWordPresent == TRUE) {
	      
		       if(verbose == TRUE) cat(":: The keyword \"control\" was found in this dataset annotations (", GSE_code, ") ", sep="")
		       healthy_control_indexes <- which(grepl(controlKeywords, (gset@phenoData@data)))
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
		      if(verbose == TRUE) cat(":: The keyword \"control\" was NOT found among the annotations of this dataset (", GSE_code, ")\n", sep="") 
	      }            
            }
                        
            if(verbose == TRUE)  cat("=== === === === === === === === === === === ===  \n", sep="")
            
            outcome <- (healthyControlWordPresent | healthyWordPresent)
            if(verbose == TRUE)  cat("\nhealthyControlsCheck() call output: were healthy controls found in the ", GSE_code, " dataset? ", outcome, "\n", sep="")
            return(outcome)
}   

  
