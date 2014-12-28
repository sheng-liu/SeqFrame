## split method on SeqFrame

## dispatch split on SeqFrame
## what it does is split as flowFrame then add annotation and reconstruct

# the splited SeqFrame remain SeqFrame instead of flowFrame (even though the print method shows it as flowFrame)
## its annotation also gets inherited but not splited
## @export split.SeqFrame
# setMethod(
#     f="split",
#     signature=c(x="SeqFrame",f="logicalFilterResult"),
#     definition=function(x,f,drop=FALSE,...){
#         
#         # split
#         list=flowCore::split(x,f)
#         # list=split(x,f)
#         
#         # merge exprs with its own annotation
#         sf.list=lapply(list,function(frame){
#             df=merge(exprs(frame),annotation(frame),by="id")
#             cat("Construct SeqFrame\n")
#             sf=df2sf(df,keyword(frame))
#         })
#         
#         return(sf.list)
#     })

# Error during wrapup: evaluation nested too deeply: infinite recursion / options(expressions=)?
# Error in loadMethod(function (x, f, drop = FALSE, ...)  : 
# error in evaluating the argument 'method' in selecting a method for function 'loadMethod': evaluation nested too deeply: infinite recursion / options(expressions=)?

# S4 method cause infiite loop, as uses split inside of the split 
# use S3 method instead, TODO, define a S4 split disptach on SeqFrame

# Error: segfault from C stack overflow


# split.SeqFrame=function(x,f,drop=FALSE,...){
#     # split
#     list=flowCore::split(x,f)
#     #class(x)="flowFrame"
#     #list=split(x,f)
#     # merge exprs with its own annotation
#     sf.list=lapply(list,function(frame){
#         df=merge(exprs(frame),annotation(frame),by="id")
#         cat("Construct SeqFrame\n")
#         sf=df2sf(df,keyword(frame))
#     })
#     return(sf.list)
# }



# it issue a note, then it is fine. 

# Note: method with signature ‘flowFrame#logicalFilterResult’ chosen for function ‘split’,
# target signature ‘SeqFrame#logicalFilterResult’.
# "SeqFrame#ANY" would also be valid
# 
# Function: split (package base)
# x="ANY", f="ANY"
# x="ANY", f="Vector"
# x="flowFrame", f="ANY"
# x="flowFrame", f="character"
# x="flowFrame", f="factor"
# x="flowFrame", f="filter"
# x="flowFrame", f="filterSet"
# x="flowFrame", f="logicalFilterResult"
# x="flowFrame", f="manyFilterResult"
# x="flowFrame", f="multipleFilterResult"
# x="flowFrame", f="numeric"
# x="flowSet", f="ANY"
# x="flowSet", f="character"
# x="flowSet", f="factor"
# x="flowSet", f="filter"
# x="flowSet", f="filterResult"
# x="flowSet", f="filterResultList"
# x="flowSet", f="list"
# x="flowSet", f="numeric"
# x="list", f="Vector"
# x="ncdfFlowList", f="character"
# x="ncdfFlowList", f="factor"
# x="ncdfFlowSet", f="character"
# x="ncdfFlowSet", f="factor"
# x="ncdfFlowSet", f="filter"
# x="ncdfFlowSet", f="filterResultList"
# x="ncdfFlowSet", f="list"
# x="RangedData", f="ANY"

# x="SeqFrame", f="ANY"
# x="SeqFrame", f="logicalFilterResult"
# (inherited from: x="flowFrame", f="logicalFilterResult")
# x="SeqFrame", f="missing"
# (inherited from: x="SeqFrame", f="ANY")

# x="SummarizedExperiment", f="ANY"
# x="UCSCData", f="ANY"
# x="UCSCData", f="Vector"
# x="Vector", f="ANY"
# x="Vector", f="Vector"

## try purely flowCore then convert to SeqFrame
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

# change flowFrame to SeqFrame
#         sf.list=lapply(list,function(frame){
#             flowFrame2SeqFrame(frame,annotation)
#         })