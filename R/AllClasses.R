# add a new slot called annotation

# for better display on mac for now
# Rstudio graphics not clear
# options(device = "quartz")



##' @import flowCore
##' @import flowViz
##' @importClassesFrom Biobase AnnotatedDataFrame
##' @importMethodsFrom Biobase AnnotatedDataFrame
##' @import flowStats
##' @import rtracklayer


.onLoad <- function(SeqData,Seqata) {
    packageStartupMessage("Loading package: SeqFrame")
}

# .onAttach <- function(SeqData,Seqata) {
#     packageStartupMessage("Done")
# }


library(flowCore)
library(flowViz)
library(Biobase)
library(flowStats)
library(rtracklayer)
# though flowViz load it, somehow it is not loaded when load the package


# for Class AnnotatedDataFrame and its assoicated functions
# its function selectChannels marks one of my function
# it is importaed for flowCore, so presumably it would have imported
# 
# Error in (function (widget, window)  : 
# could not find function "AnnotatedDataFrame"

## -----------------------------------------------------------------------------
## Class SeqFrame

# setClass(Class="SeqFrame",
#          representation(annotation="data.frame"),
#          contains="flowFrame") 

# getClass("SeqFrame")

# guid <- function(len=10){
#     ltrs <- c(LETTERS,letters)
#     paste(c(sample(ltrs,1),sample(c(ltrs,0:9),len-1,replace=TRUE)),collapse="")
# }


# > getClass("flowFrame")
# Class "flowFrame" [package "flowCore"]
# Slots:
# Name:               exprs         parameters        description
# Class:       NcdfOrMatrix AnnotatedDataFrame               list


# setClass(Class="SeqFrame",
#          contains="flowFrame",
#          representation(annotation="data.frame")) 
# change the sequence won't change the sequence slot printed out

##' @exportClass SeqFrame
setClass(Class="SeqFrame",
         contains="flowFrame",
         representation(annotation="data.frame")) 

# data.frame is better in this case, as it is easier to manipulate and display, althought they are essentially the same. 
# when construct, user has to pass in the annotation they want to look at, it may not be the same as the parameters were calculated, such as promoter for histone coverage, while exon counts for the same gene's expression. as long as they are corresponding, they are good. 

# > getClass("SeqFrame")
# Class "SeqFrame" [in ".GlobalEnv"]
# Slots: 
# Name:     annotation      exprs           parameters              description
# Class:    data.frame      NcdfOrMatrix    AnnotatedDataFrame      list                        
# Extends: "flowFrame"

# since I haven't decide whether and what to do with annotation, just use what's working now annotation




# the SeqFrame constructor 
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




# setGeneric(
#     name="SeqFrame",
#     def=function(obj,bamFile=character(0)){
#         standardGeneric("getReadCoverage")
#     })


# 
# # a SeqFrame constructor makes flowFrame to SeqFrame with annotation information

# setMethod(
#     f="SeqFrame",
#     signature="flowFrame",
#     definition=function(obj,annotation){
#        # sf=SeqFrame(exprs=exprs(obj),
#         sf=new("SeqFrame",
#                exprs=exprs(obj),           
#                parameters=parameters(obj),
#                description=description(obj),
#                annotation=annotation
#                     )
#     })





# it turn out the splited SeqFrame remain SeqFrame instead of flowFrame(even though the print method shows it as flowFrame)
## its annotation also gets inherited but not splited


# setMethod(
#     f="SeqFrame",
#     signature="SeqFrame",
#     definition=function(obj,annotation){
#         # sf=SeqFrame(exprs=exprs(obj),
#         sf=new("SeqFrame",
#                exprs=exprs(obj),           
#                parameters=parameters(obj),
#                description=description(obj),
#                annotation=annotation
#         )
#     })


## change annotation to annotation avoid masking by Biobase and BioGenerics


# function(sf,annotation){
#         # sf=SeqFrame(exprs=exprs(obj),
#         sf=new("SeqFrame",
#                exprs=exprs(obj),           
#                parameters=parameters(obj),
#                description=description(obj),
#                annotation=annotation
#         )
#     })









