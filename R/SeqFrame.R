# SeqFrame-methods
# 
# 
###############################################################################
##' @name SeqFrame
##' @aliases SeqFrame
##' @title SeqFrame
##' @rdname SeqFrame-methods
##' @docType methods
##' @description method for coerce spread sheet like data into SeqFrame. data type include data.fame, comma seperated file, tab delimited file and file formats representing annotated genomic intervals, i.e. gff, bed, bed15, bedGraph, wig and bigWig format. Also include a facility for coercing flowFrame into SeqFrame. 


## if data have data have columns, i.e. "chr", "start","end","symbol";"strand" is optional. it will put into annotation slot, one can use annotation() method to access the corresponding genomic annotation, else it is simply converted to seqFrame with empty annotation, for simple visulization. then it needs to be all numbers. 

##' @usage SeqFrame(df=NULL,keyword=NULL,file=character(0),con=character(0),annotation=data.frame(),)
##' 
##' 
##' @param df data.frame to be converted to SeqFrame. if data.frame have data have columns, i.e. "chr", "start","end","symbol";"strand" is optional.it will put into annotation slot, one can use annotation() method to access the corresponding genomic annotation, else it is simply converted to seqFrame with empty annotation, for simple visulization. then it needs to be all numbers. 

##' @param keyword a list including descriptions for the columns. see keyword in flowFrame for detail.
##' 
##' @param file comma seperated file, tab delimited file. Header is required. if data have columns, i.e. "chr", "start","end","symbol";"strand" is optional. it will put into annotation slot, else it is simply converted to seqFrame with empty annotation, for simple visulization. then it needs to be all numbers.
##' 
##' @param con filename to file formats representing annotated genomic intervals, i.e. gff, bed, bed15, bedGraph, wig and bigWig format. File type is automatically detected by file extension. 
##' 
##' @param annotation a data.frame containing annotation information, need to have at least five columns i.e. "chr", "start","end","symbol";"strand" is optional.


##' @details SeqFrame is the anolog of flowFrame in flowCore. It provide an accessor for spread sheet like genomic data to be viewed as events in flow cytometry. 


## SeqSet is optional. as one can simply use flowSet ?
##' @return A SeqFrame containing data and annotation. Users can use all flowCore and flowViz visulization functions to subset, split data based on its geometric property, and use annotation to access the corresponding annotation. 
##' 
##'  @examples 
##' 
##' file="/Users/shengliu/DoScience/DoScience/Projects/Clover/oocytes methylation /DNAme/NGO-GVO-methylation.csv"
##' SeqFrame(file=file)
##' 
##' @exportMethod SeqFrame
setGeneric(
    name="SeqFrame",
    def=function(df=NULL,annotation=data.frame(),keyword=NULL,file=character(0),con=character(0)){
        standardGeneric("SeqFrame")
    })

## -----------------------------------------------------------------------------
## dispatch on data.frame
## making a data.frame into SeqFrame for visualization
## data.frame need to have five columns "chr", "start","end","symbol",stand is optional
## if data.frame have five columns, it will put into annotation slot, else it is
## simply converted to seqFrame with empty annotation, for simple visulization
## then it needs to be all numbers. 

setMethod(
    f="SeqFrame",
    signature=c(df="data.frame"),
    definition=function(df,keyword=NULL){ 
        sf=df2sf(df,keyword=keyword) 
        return(sf)
    })

# .local <- function (df) 
# {
#     nam <- deparse(substitute(df))
#     cat("name:", nam, "\n")
#     df
# }
# .local(df)
# data.name=deparse(substitute(df))
# print("dispatch on data.frame")
# keyword=list(FILENAME=data.name,GUID=data.name)
# print(keyword)


## dispatch on character
## making a csv file into SeqFrame 
setMethod(
    f="SeqFrame",
    signature=c(df=NULL,file="character"),
    definition=function(df,file="character"){ 
    	sf=file2sf(file=file)
        return(sf)
    	})

## dispatch on flowFrame
## making a flowFrame file into SeqFrame 
## sometime when using the flowCore, it is useful to convert a flowFrame to 
## SeqFrame so one can easily access the important annotation information

setMethod(
    f="SeqFrame",
    signature=c(df="flowFrame"),
    definition=function(df,annotation=data.frame()){ 
    	sf=ff2sf(ff=df,annotation=annotation)
        return(sf)
    	})


## dispatch on file formats representing annotated genomic intervals
## GFF, BED, BED15, BEDGRAPH, WIG, BIGWIG
## gff, bed, bed15, bedGraph, wig, bigWig
 ## see help ??UCSCData



setMethod(
    f="SeqFrame",
    signature=c(df="missing",file="missing",con="character"),
    definition=function(con){ 
    

        ## return is a GRanges object
        ## all of the file format have a column "score", 
        ## convert that one into SeqFrame
        ## chr is called seqnames
        gr=import(con=con)
        
        ## 5th column (blocks) in bed format is grlist, need to be removed to output as data.frame. see "test UCSCData.rmd" for detail
        ## TODO: best way would be convert it into character string
        ## for now remove the metacolumn called blocks
        ## which==
        gr=gr[,-5]
        df=as.data.frame(gr)
        ## change seqnames into chr
        colnames(df)[1]="chr"
        
        sf=SeqFrame(df)
        return(sf)
    })


