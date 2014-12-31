# SeqFrame-methods
# 
# 
###############################################################################
##' @name SeqFrame
##' @aliases SeqFrame
##' @title SeqFrame
##' @rdname SeqFrame-methods
##' @docType methods
##' @description method for coerce spread sheet like data into SeqFrame. data
##'   type include data.fame, comma seperated file, tab delimited file and file
##'   formats representing annotated genomic intervals, i.e. gff, bed, bed15,
##'   bedGraph, wig and bigWig format. Also include a facility for coercing
##'   flowFrame into SeqFrame.

##' @usage 
##' SeqFrame(df=NULL,keyword=NULL,file=character(0),con=character(0),annotation=data.frame(),)
##' 
##' @param df data.frame to be converted to SeqFrame. if data.frame have data
##'   have columns, i.e. "chr", "start","end","symbol";"strand" is optional.it
##'   will put into annotation slot, one can use annotation() method to access
##'   the corresponding genomic annotation. One can also construct SeqFrame
##'   without annotation information for simple visulization.

##' @param keyword a list including descriptions for the columns. see keyword in
##'   flowFrame for detail.
##'   
##' @param file comma seperated file, tab delimited file. Header is required. if
##'   data have columns, i.e. "chr", "start","end","symbol";"strand" is
##'   optional. it will put into annotation slot, else it is simply converted to
##'   seqFrame with empty annotation, for simple visulization. then it needs to
##'   be all numbers.
##'   
##' @param con filename to file formats representing annotated genomic
##'   intervals, i.e. gff, bed, bed15, bedGraph, wig and bigWig format. File
##'   type is automatically detected by file extension.
##'   
##' @param annotation a data.frame containing annotation information, need to
##'   have at least five columns i.e. "chr", "start","end","symbol";"strand" is
##'   optional.
##'   
##'   
##' @details SeqFrame is the anolog of flowFrame in flowCore. It provide an 
##'   accessor for spread sheet like genomic data to be viewed as events in flow
##'   cytometry.
##'   
##'   the required field for annotation is "chr", "start","end","symbol". 
##'   Optionally "strand" (if not present, will be repalced with "*").
##'   Information will be feed into annotation slot of SeqFrame. Users can use
##'   annotation method to retrieve the corresponding genomic annotation. One
##'   can create SeqFrame without annotation, for simple visulization of the
##'   data.

##' @return A SeqFrame containing data and annotation. Users can use all
##'   flowCore and flowViz visulization functions to subset, split data based on
##'   its geometric property, and use annotation to access the corresponding
##'   annotation.
##' 
##'  @examples
##'  
##' # csv file 
##' csv.file=system.file("extdata", "NGO-GVO-methylation.csv",package="SeqFrame") 
##' csv=SeqFrame(file=csv.file); csv 
##' test.csv=csv[200:210,];test.csv
##' annotation(test.csv);exprs(test.csv)
##' 
##' # bigWig file
##' bw.file=system.file("tests", "test.bw", package = "rtracklayer")
##' bw=SeqFrame(con=bw.file)
##' exprs(bw);annotation(bw)
##' 
##' # bedGraph
##' bedGraph.file=system.file("tests", "test.bedGraph", package = "rtracklayer")
##' bg=SeqFrame(con=bedGraph.file)
##' exprs(bg); annotation(bg)
##' 
##' # bed.file
##' bed.file=system.file("tests", "test.bed", package = "rtracklayer")
##' bed=SeqFrame(con=bed.file)
##' exprs(bed)

## FIXME
## # gff.file
## gff.file=system.file("tests", "v1.gff", package = "rtracklayer")
## gff=SeqFrame(con=gff.file)
## exprs(gff)
## 

## -----------------------------------------------------------------------------
## SeqFrame method

##' @exportMethod SeqFrame
setGeneric(
    name="SeqFrame",
    def=function(df=NULL,annotation=data.frame(),keyword=NULL,file=character(0),con=character(0)){
        standardGeneric("SeqFrame")
    })

## -----------------------------------------------------------------------------
## dispatches
## dispatch on data.frame
## data.frame need to have five columns "chr", "start","end","symbol",stand is
## optional

setMethod(
    f="SeqFrame",
    signature=c(df="data.frame"),
    definition=function(df,keyword=NULL){ 
        sf=df2sf(df,keyword=keyword) 
        return(sf)
    })

## -----------------------------------------------------------------------------
## dispatch on character
## making a csv file into SeqFrame 
setMethod(
    f="SeqFrame",
    signature=c(df=NULL,file="character"),
    definition=function(df,file="character"){ 
    	sf=file2sf(file=file)
        return(sf)
    	})

## -----------------------------------------------------------------------------
## dispatch on flowFrame
## convert a flowFrame file into SeqFrame so one can easily access the important
## annotation information

setMethod(
    f="SeqFrame",
    signature=c(df="flowFrame"),
    definition=function(df,annotation=data.frame()){ 
    	sf=ff2sf(ff=df,annotation=annotation)
        return(sf)
    	})

## -----------------------------------------------------------------------------
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
        
        ## 5th column (blocks) in bed format is grlist, need to be removed to
        ## output as data.frame. see "test UCSCData.rmd" for detail
        ## TODO: best way would be convert it into character string for now
        ## remove the metacolumn called blocks

        gr=gr[,-5]
        df=as.data.frame(gr)
        ## change seqnames into chr
        colnames(df)[1]="chr"
        
        sf=SeqFrame(df)
        return(sf)
    })


## -----------------------------------------------------------------------------
## Notes
## SeqSet is optional. as one can simply use list and making it a flowSet to do
## multiple data set manipulation.


