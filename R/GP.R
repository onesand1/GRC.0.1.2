#' @title  The corelated analysis of expression profile for EMT score and gene expression
#
#' @param expression  expression profile
#' @param gsva the option for GSVA analysis refer to the package GSVA.
#' @param geneset  a string for gmt file path or a list for one or more gene lists.
#' @param  log   log=T, The expression profile will be normalized by log2(RPKM+1).
#' @param plot A logical parameter. This parameter is defined to show whether the correlation is visualized.default plot=TRUE.

GP<-function(expression,geneset,log=TRUE,gsva="",plot=TRUE){
  library(ggplot2)
  if(plot){
    if(!dir.exists("cor_plot")){
      dir.create("cor_plot")
    }
  }
  result=list()
  samp=as.character(colnames(expression))
  gene=as.character(rownames(expression))
  if(log){
    data=log2(data+1)
    rownames(data)=gene
    colnames(data)=samp

  }


  if(class(geneset)=="character"){
    if(!file.exists(geneset)){
      stop("Please check the filename or file path !")
    }
    gmt<- getGmt(geneset)

  }
  if(class(geneset)=="list"){
    n=length(geneset)
    name=names(geneset)
    for(i in 1:n){
      res=paste(unlist(geneset[i]),collapse = '\t')
      sink('Geneset.gmt')
      cat(name[i])
      cat('\tNA\t')
      cat(res)
      cat('\n')
      sink()


    }
    gmt<- getGmt('Geneset.gmt')
  }
  result$geneSet.name=as.character(names(gmt))
  result$geneSet.gene=as.character(gmt[[1]]@geneIds)

  int_gene=intersect(result$geneSet.gene,gene)
  n=length(int_gene)
  if(n<2){
    stop("The number of intersecting genes between genesets and expression profile no more than 1!")
  }
  result$geneSet.score=as.vector(gsva(data,gmt,min.sz=5,gsva))
  n=length(gene)
  cor_res=c()
  for(i in 1:n){
    cor=cor.test(as.numeric(result$geneSet.score),as.numeric(t(expression[i,])),method ="spearman")
    P=cor$p.value
    rho=cor$estimate
    re=cbind(gene[i],P,rho)
    cor_res=rbind( cor_res,re)
    if(plot){
     # pdf(paste("cor_plot/",gene[i],"_cor.pdf",sep=""))
      dat=cbind(result$geneSet.score,as.numeric(t(expression[i,])))
      colnames(dat)=c("Score","Expression")
      dat=data.frame(dat)
      pp=ggplot(dat, aes(x=Score, y=Expression)) +
        geom_point() +    # Use hollow circles
        geom_smooth(method=lm)+
        labs(title=paste(gene[i],", Cor=",round(rho,3),", P-value=",format(P, scientific = TRUE,digits=4)
,sep=""),x ="Score", y = "Expression")
      ggsave(pp, file=paste("cor_plot/",gene[i],"_cor.pdf",sep=""))

    }

  }
  cor_res=data.frame(cor_res)
  colnames(cor_res)=c('gene','P-value',"Spearman's coefficient")
  result$cor.res=cor_res
  return(result)

}
