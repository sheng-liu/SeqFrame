# SeqFrame - methods
# 
# 
###############################################################################
## so add a id column (the row.numbers), use it subset annotation, transfer to
## grangs and open genome browser todo: remove or hide id columne as row names
## in flowCore

## -----------------------------------------------------------------------------
## SeqFrame annotation method

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
        # methods::initialize(obj,annotation=value) 
        obj@annotation=value
        obj  # this cost me a lot of time to figure out
    })

# > annotation(x)=annotation(x)[i,]
# > x
# SeqFrame object 'anonymous'
# with 0 genomic features and -1 epigenetic observables:
#     [1] name        description range       minRange    maxRange   
# <0 rows> (or 0-length row.names)
# 1 keywords are stored in the 'description' slot

setMethod(
    f="annotation",
    signature="SeqFrame",
    definition=function(obj){
        obj@annotation
    })

## -----------------------------------------------------------------------------
## SeqFrame show method
## this is a version from flowFrame

##' @exportMethod show
setMethod("show",
          signature=signature(object="SeqFrame"),
          definition=function(object)
          {
              dm <- dim(exprs(object))
              cat(paste("SeqFrame object '", identifier(object),
                        "'\nwith ", dm[1], " genomic features and ", 
                        dm[2]-1, " epigenetic observables:\n", sep=""))
              show(pData(parameters(object)))
              cat(paste(length(description(object)), " keywords are stored in the ",
                        "'description' slot\n", sep = ""))
              return(invisible(NULL))
          })

## -----------------------------------------------------------------------------
## SeqFrame subset methods

## "[" function
## by indices 

##' @exportMethod "["
setMethod(f="[",
          signature=c(x="SeqFrame"),
          definition=function(x, i, j, ..., drop=FALSE)
          {
              # based on location of i,j, defines subsetting method
              # sf[1:100,],  missing(j)=T
              # sf[,1:100],  missing(i)=T
              switch(1+missing(i)+2*missing(j),
{
    # subset exprs, and annotation
    exprs(x) = exprs(x)[i, j] 
    annotation(x)=annotation(x)[i , j]
},
{
    x = subsetKeywords(x, j)
    exprs(x) = exprs(x)[ , j]
    annotation(x) = annotation(x)[ , j]
},
{
    exprs(x) = exprs(x)[i,  ]
    #annotation(x) = annotation(x)[i, ]
    annotation(x) = annotation(x)[i,  ]
},
{
    exprs(x) = exprs(x)[ ,  ]
    annotation(x) = annotation(x)[, ]
} )
return(x)

          })


## TODO:
## ADD subset by logical vectors
## ADD subsetting with $
## ADD subset by filter

## -----------------------------------------------------------------------------
## SeqFrame split method

##' @exportMethod split
setMethod(
    f="split",
    signature=c(x="SeqFrame",f="logicalFilterResult"),
    definition=function(x,f,drop=FALSE,...){
        
        annotation=annotation(x)
        class(x)="flowFrame"
        # split as flowFrame
        list=flowCore::split(x,f)
        
        # merge exprs with its own annotation
        sf.list=lapply(list,function(frame){
            df=merge(exprs(frame),annotation,by="id")
            
            cat("Construct SeqFrame\n")
            sf=df2sf(df,keyword(frame))
            
            # re-order exprs(sf) columns to its original sequence
            cln=colnames(exprs(frame))
            exprs(sf)=exprs(sf)[,cln]
            sf
        })
        
        return(sf.list)
    })


