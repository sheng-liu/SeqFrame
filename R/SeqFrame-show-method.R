
## SeqFrame show method
## this is from flowFrame


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