## helpers

## df2sf

## one need to select the parameter column to be output
## the annotation column can be done automaticly. 

## "chr","start","end","strand","symbol"
## "readcounts", coverage, percent methylation 


# file="/Users/shengliu/DoScience/DoScience/Projects/Clover/oocytes methylation /DNAme/NGO-GVO-methylation.csv"
# 
# df=read.csv(file=file,as.is=T,header=T)
# keyword=list(FILENAME=data.name,GUID=data.name)
## keyword is a list 
##' @export df2sf
df2sf=function(df,keyword=NULL){
    
    # use df as FILENAME if FILENAME is not supplied
    data.name=deparse(substitute(df))
    # note once df value changed, substitute(df) will produce the new df value instead of its origianl name (promised object)
    
    # add strand if there isn't 
    cln=colnames(df)
    
    # check strand column
    if (length(which(cln=="strand"))==0) {
        strand=rep("*",dim(df)[1])
        df=cbind(df,strand)
    }
    
    # check name colum for UCSC bed file
    if (length(which(cln=="name"))!=0){
        symbol=df$name
        df=cbind(df,symbol)
    }
    
    # check symbol column 
    if (length(which(cln=="symbol"))==0) {
        symbol=rep("UCSCData",dim(df)[1])
        df=cbind(df,symbol)
    }
    
    # add ID column
    id=as.integer(row.names(df))
    df=cbind(df,id)
    cln=colnames(df)
    
    # choose exprs columns to be columns other than annotation
    # this has the risk that some columns may not be numeric
    # a safer solution is only pick one column "score", then one
    # loses calculatable information
    conserved.cln=c("chr","start","end","strand","symbol")
    unconserved.cln=setdiff(cln,conserved.cln)
    #unconserved.cln=unconserved.cln[complete.cases(unconserved.cln)]
    
    
    # another subsetting could be subset based on numeric or character, put all characters into annotation metacolumns, numeric into exprs, then will not lose any information in the conversion
    
    
    
    # subset data.frame directly with a vector 
    # exprs.mx=sapply(unconserved.cln,function(x){df[[x]]})  # matrix
    # df[unconserved.cln] or df[,unconserved.cln] both right 
    # one view df as list of columns, the other as column indexing
    unconserved.df=df[unconserved.cln]
    
    # check column type, remove column if not numeric
    # need to convert to data.frame, otherwise all viewed as character in mx, 
    # as.numeric will remove character including header
    # subset by type, remove columns that are "character"
    type=sapply(unconserved.df,class)    
    exprs.df=unconserved.df[which(type!="character")]
    exprs.mx=as.matrix(exprs.df)
    
    # annotation
    conserved.cln=c(conserved.cln,"id")
    
    
    #     anno.mx=sapply(conserved.cln,function(x){df[[x]]})
    #     anno.df=as.data.frame(anno.mx,stringAsFactors=F) # matrix, and it is still factor
    #     
    # 
    #     # if one of the column is null, delete that column
    #     # currently maintain all column
    #     anno.df=transform(anno.df,
    #                       id=as.integer(as.character(id)),
    #                       start=as.integer(as.character(start)),
    #                       end=as.integer(as.character(end)),
    #                       symbol=as.character(symbol),
    #                       chr=as.character(chr)
    #     )
    
    
    #     # subset data.frame directly with vector
    #     anno.df=df[conserved.cln]
    #     # add all rest character columns to metacolumn
    #     character.df=unconserved.df[which(type=="character")]
    #     anno.df=cbind()
    
    character.cln=unconserved.cln[which(type=="character")]
    anno.cln=c(conserved.cln,character.cln)
    anno.df=df[anno.cln]
    #conserved.cln=conserved.cln[complete.cases(conserved.cln)]
    
    # GRanges form for future use if necessary
    # sf.anno=df2gr(anno.df)  
    
    
    # construction of AnnotatedDataFrame 
    # must contain this five "name, desc, range, minRange and maxRange" 
    # to be able to plot right
    name=colnames(exprs.df)
    desc=colnames(exprs.df)  # add description for each channel
    range=ceiling(apply(exprs.df,2,max))
    minRange=rep(0,dim(exprs.df)[2])
    maxRange=ceiling(apply(exprs.df,2,max))
    parameters.adf=data.frame(name,desc,range,minRange,maxRange,
                              stringsAsFactors = FALSE)
    parameters.adf=AnnotatedDataFrame(data=parameters.adf)
    
    if(is.null(keyword)) {
        # cat("keyword is null, using default data name as keyword.\n")
        keyword=list(FILENAME=data.name,GUID=data.name)
    }
    sf=new("SeqFrame",
           exprs=exprs.mx,
           parameters=parameters.adf,
           description=keyword,
           annotation=anno.df)  #sf.anno
    return(sf)
}


## two important description (keyword)
## access description
# description(sf)$FILENAME=filename
# description(sf)$GUID=filename

# or use keyword method to access description
# keyword(sf)$FILENAME=filename
# keyword(sf)$GUID=filename


# access parameters
# pData(parameters(sf))

##------------------------------------------------------------------------------
##' @export file2sf
file2sf=function(file){
    df=read.csv(file=file,as.is=T,header=T)
    data.name=basename(file)
    keyword=list(FILENAME=data.name,GUID=data.name)
    df2sf(df,keyword)
}


##------------------------------------------------------------------------------
##' @export veggie
veggie=function(){
    veg=c("avocado","melon","olive","pepper","pumpkin","vanilla","tomato","squash","broccoli","chickpea","bean","celery")
    
    
    name=paste(sample(veg,1),sample(0:9,1),sample(0:9,1),sep="")
    return(name)
}




#     num=c(sample(0:9,3,replace=T))
#     num=as.integer(sample(0:9,3,replace=T))
#     sample(0:9,3)


# guid <- function(len=10){
#     ltrs <- letters
#     paste(sample(c(letters,0:9),replace=TRUE),collapse="")
# }

##------------------------------------------------------------------------------
##' @export SeqFrame.table
SeqFrame.table=function(sf){
    merge(annotation(sf),exprs(sf),by="id") 
}




##' @export ff2sf
## flowFrame2SeqFrame=function(ff,annotation){
ff2sf=function(ff,annotation){
    # sf=SeqFrame(exprs=exprs(ff),
    sf=new("SeqFrame",
           exprs=exprs(ff),           
           parameters=parameters(ff),
           description=description(ff),
           annotation=annotation
    )
}


##------------------------------------------------------------------------------


