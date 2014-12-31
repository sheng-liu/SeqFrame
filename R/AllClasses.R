
## -----------------------------------------------------------------------------
## Dependencies

##' @import flowCore
##' @import flowViz
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @importMethodsFrom Biobase AnnotatedDataFrame
##' @import flowStats
##' @import rtracklayer

# welcome message
.onLoad <- function(SeqData,Seqata) {
    packageStartupMessage("Loading package: SeqFrame")
}

library(flowCore)
library(flowViz)
library(Biobase)
library(flowStats)
library(rtracklayer)

## -----------------------------------------------------------------------------
## Class SeqFrame


##' @exportClass SeqFrame
setClass(Class="SeqFrame",
         contains="flowFrame",
         representation(annotation="data.frame")) 


# the SeqFrame constructor
# constructs flowFrame to SeqFrame with annotation information
##' @export SeqFrame
SeqFrame=function(exprs=matrix(),
                  parameters=AnnotatedDataFrame(),
                  description=list(),
                  annotation=data.frame()){
    new("SeqFrame",
        exprs=exprs,
        parameters=parameters,
        description=description,
        annotation=annotation)
}

## -----------------------------------------------------------------------------
## Notes:

# Deciding on using GRanges or data.frame as the format of annotation.
# data.frame is better in this case, as it is easier to manipulate and display,
# althought they essentially contain the same informaiton.

# > getClass("SeqFrame")
# Class "SeqFrame" [in ".GlobalEnv"]
# Slots: 
# Name:     annotation      exprs           parameters              description
# Class:    data.frame      NcdfOrMatrix    AnnotatedDataFrame      list                        
# Extends: "flowFrame"


