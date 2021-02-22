#' @title  read the data for expression profile
#'
#' @param data The path and file name for inputing file. File type was .csv.
#' @export
#' @return  a list, includeing gene The gene for the first column, path The input file path, data The expression matrix data and sample The sample name for first row..
#' @example
#'
#' @return
#'
#'



readEXP<-function(data){
  object=list()
  if(!file.exists(data)){
    stop("Please check the filename or file path !")
  }
  object$path=data
  f=read.csv(data,header = F,stringsAsFactors = F)
  object$gene=as.character(f[-1,1])
  object$sample=as.character(f[1,-1])
  objectt$data=as.matrix(f[-1,-1])
  return(object)

}
