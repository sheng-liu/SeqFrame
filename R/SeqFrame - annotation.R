# SeqFrame - Accessors
# 
# 
###############################################################################
## so add a id column (the row.numbers), use it subset annotation, transfer to grangs and open genome browser
## todo: remove or hide id columne as row names in flowCore

## generics

##'@exportMethod annotation
setGeneric(
    name="annotation",
    def=function(obj){
        standardGeneric("annotation")
    })

##'@exportMethod annotation<-
setGeneric(
    name="annotation<-",
    def=function(obj,value){
        standardGeneric("annotation<-")
    })


## methods
setReplaceMethod(
    f="annotation",
    signature="SeqFrame",
    definition=function(obj,value){
        initialize(obj,annotation=value)    
    })

setMethod(
    f="annotation",
    signature="SeqFrame",
    definition=function(obj){
        obj@annotation
    })


