#' @title A function for regulation between gane and genes
#'
#' @param target the gene that will be check in expression profile
#' @param data the dataset for gene expression. Consisting of a matrix
#' @param plot a logical value for ploting results.default TRUE, print the figure to window.
#' @param type The ragulation type for target. if type=1,  thie model will caculate the downstream genes for target gene; if type=2, it will caculate the  upstream genes for target gene.
#' @param scale a logical value indicating whether the expression profile needs to normalized by log2(X+1).default scale=TRUE.
#' @param group The group for the target gene expresion profile. default is median and group="med"; it also can set "quar" means  gene expression value less than a quartile and more than three quartile.
#' @param FC  The fold change for difference between target gene and other genes.
#' @param FDR the False discovery rate for results.if type=2,the FDR eual to P-value.
#' @return The result for the regulator analysis.
#'
#' @example
#'
#' @export
#'
#'
regulator<-function(target,data,plot=TRUE,type="1",scale=TRUE,group='med',FC=2,FDR=0.05){
  results=list()

  if(class(data)[1]!="matrix"){
    stop("The data is not a matrix!")

  }
  if(dim(data)[2]<5){
  stop("The number of samples must more than 5!")

  }
  if(target==""){
    stop("The target gene isn't indicated !")
  }
  results$root=target
  gene=as.character(rownames(data))
  twhich=which(gene==as.character(target))
  if(length(twhich)!=1){
    stop("The target gene isn't indicated in your data or more than one target is indicated!")
  }
  sample=colnames(data)
  tdata=as.numeric(t(data[twhich,]))
  data=data[-twhich,]
  gene=gene[-twhich]
  n=dim(data)[1]
  if(type==1){
    results$type="Down-regulation"
    gg=as.numeric(quantile(tdata))
    if(group=="med"){
      high=which(tdata>gg[3])
      low=which(tdata<=gg[3])

    }
    if(group=="quar"){
      high=which(tdata>=gg[4])
      low=which(tdata<=gg[2])

    }
    if(length(high)<3 | length(low)<3){
      warning(paste("The number of samples must more than 3 in every group for gene '",gene[i] ,"' !",sep=""))

    }else{
      res=c()
      for(i in 1:n){
        if(scale){
          p<-t.test(log2(as.numeric(data[i,high])+1),log2(as.numeric(data[i,low])+1),var=TRUE)$p.value
          d<-mean(log2(as.numeric(data[i,high])+1))-mean(log2(as.numeric(data[i,low])+1))
        }else{
          p<-wilcox.test(as.numeric(data[i,high]),as.numeric(data[i,low]))$p.value
          d<-mean(as.numeric(data[i,high]))/mean(as.numeric(data[i,low]))
        }
        re=cbind(gene[i],p,d)
        res=rbind(res,re)
      } }
    if(!length(dim(res)[1])){
      if(as.numeric(res[2])<as.numeric(FDR)){
        results$target=as.character(res[1])
        results$FDR="NA"
        results$P=as.character(res[2])
        results$FC=as.character(res[3])
        if(scale){
          results$reg=as.numeric(as.numeric(results$FC)>0)
        }else{
          results$reg=as.numeric(as.numeric(results$FC)>1)
        }

        results$output=res
        if(plot){
          plot.reg.up(results)

        }
      }else{
        warning("There is not value for the results !")

      }

    }else{
      ord<-order(as.numeric(res[,2]),decreasing=FALSE)
      new_res<-res[ord,]
      fdr<-rep(1,dim(new_res)[1]); for(j in 1:dim(new_res)[1]) fdr[j]<-as.numeric(new_res[j,2])*dim(new_res)[1]/j
      if(scale){

        f_thres=which(abs(as.numeric(new_res[,3]))>log2(as.numeric(FC)))

      }else{
        f_thres=which(abs(as.numeric(new_res[,3]))>as.numeric(FC))

      }
      p_thres=which(fdr<as.numeric(FDR))
      DGE<-cbind(new_res,fdr)
      colnames(DGE)=c("gene","P-value","Fold Change","FDR")
      ss=intersect(p_thres,f_thres)
      if(length(ss)>0){
        res=DGE[ss,]
        results$target=as.character(res[,1])
        results$FDR=as.character(res[,4])
        results$FC=as.character(res[,3])
        if(scale){
          results$reg=as.numeric(as.numeric(results$FC)>0)
        }else{
          results$reg=as.numeric(as.numeric(results$FC)>1)
        }
        results$output=res
        if(plot){
          plot.reg.down(results)

        }
      }else{

        warning("There is not value for the results !")

      }
}
  }


  if(type==2){
    results$type="Up-regulation"
    res=c()
    for(i in 1:n){
      gdata=as.numeric(t(data[i,]))
      gg=as.numeric(quantile(gdata))
      if(group=="med"){
        high=which(gdata>gg[3])
        low=which(gdata<=gg[3])

      }
      if(group=="quar"){
        high=which(gdata>=gg[4])
        low=which(gdata<=gg[2])

      }
      if(length(high)<3 | length(low)<3){
        warning(paste("The number of samples must more than 3 in every group for gene '",gene[i] ,"' !",sep=""))

      }else{
      if(scale){
        p<-t.test(log2(as.numeric(tdata[high])+1),log2(as.numeric(tdata[low])+1),var=TRUE)$p.value
        d<-mean(log2(as.numeric(tdata[high])+1))-mean(log2(as.numeric(tdata[low])+1))
      }else{
        p<-wilcox.test(as.numeric(tdata[high]),as.numeric(tdata[low]))$p.value
        d<-mean(as.numeric(tdata[high]))/mean(as.numeric(tdata[low]))
      }
      re=cbind(gene[i],p,d)
      res=rbind(res,re)
}
    }
    if(!length(dim(res)[1])){
      if(as.numeric(res[2])<as.numeric(FDR)){
      results$target=as.character(res[1])
      results$FDR="NA"
      results$P=as.character(res[2])
      results$FC=as.character(res[3])
      if(scale){
        results$reg=as.numeric(as.numeric(results$FC)>0)
      }else{
        results$reg=as.numeric(as.numeric(results$FC)>1)
      }
      results$output=res
      if(plot){
        plot.reg.up(results)

      }
      }else{
        warning("There is not value for the results !")

      }

    }else{
      ord<-order(as.numeric(res[,2]),decreasing=FALSE)
      new_res<-res[ord,]
      fdr<-rep(1,dim(new_res)[1]); for(j in 1:dim(new_res)[1]) fdr[j]<-as.numeric(new_res[j,2])*dim(new_res)[1]/j
      if(scale){

        f_thres=which(abs(as.numeric(new_res[,3]))>log2(as.numeric(FC)))

      }else{
        f_thres=which(abs(as.numeric(new_res[,3]))>as.numeric(FC))

      }


      p_thres=which(fdr<as.numeric(FDR))
      DGE<-cbind(new_res,fdr)
      colnames(DGE)=c("gene","P-value","Fold Change","FDR")
      ss=intersect(p_thres,f_thres)
      if(length(ss)>0){
        res=DGE[ss,]
        results$target=as.character(res[,1])
        results$FDR=as.character(res[,4])
        results$P=as.character(res[,2])
        results$FC=as.character(res[,3])
        if(scale){
          results$reg=as.numeric(as.numeric(results$FC)>0)
        }else{
          results$reg=as.numeric(as.numeric(results$FC)>1)
        }
        results$output=res
        if(plot){
          plot.reg.up(results)

        }
      }else{

        warning("There is not value for the results !")

      }
}
  }
  return(results)
}


