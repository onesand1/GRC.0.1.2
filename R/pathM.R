#' @title pathway build function
#' @description Based on the expression level, this function constructs a predicted pathway map

pathM<-function(data,FC=2,FDR=0.05,group="med",scale=TRUE){

  results=weight=list()

  if(class(data)[1]!="matrix"){
    stop("The data is not a matrix!")

  }
  if(dim(data)[2]<6){
    stop("The number of samples must more than 6!")

  }
  sample=colnames(data)
  gene=rownames(data)
  gene1=gene2=c()
  m=length(gene)
  results=rep(list(NA),m)
  names(results)=gene
  arrow=color=c()
  for(j in 1:m){
    tdata=as.numeric(t(data[j,]))
    dat=data[-j,]
    gen=gene[-j]
    n=dim(dat)[1]
    gg=as.numeric(quantile(tdata))
    if(length(which(gg>0))>3){
      if(group=="med"){
        high=which(tdata>gg[3])
        low=which(tdata<=gg[3])

      }
      if(group=="quar"){
        high=which(tdata>=gg[4])
        low=which(tdata<=gg[2])

      }
      if(length(high)<3 | length(low)<3){
        results[[j]]=list()
        gene1=c(gene1,gene[j])


      }else{
        res=c()
        for(i in 1:n){
          if(scale){
            p<-t.test(log2(as.numeric(dat[i,high])+1),log2(as.numeric(dat[i,low])+1),var=TRUE)$p.value
            d<-mean(log2(as.numeric(dat[i,high])+1))-mean(log2(as.numeric(dat[i,low])+1))
          }else{
            p<-wilcox.test(as.numeric(dat[i,high]),as.numeric(dat[i,low]))$p.value
            d<-mean(as.numeric(dat[i,high]))/mean(as.numeric(dat[i,low]))
          }
          re=cbind(gen[i],p,d)
          res=rbind(res,re)
        }

    ord<-order(as.numeric(res[,2]),decreasing=FALSE)
    new_res<-res[ord,]
    fdr<-rep(1,dim(new_res)[1]); for(k in 1:dim(new_res)[1]) fdr[k]<-as.numeric(new_res[k,2])*dim(new_res)[1]/k
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
      if(length(ss)==1){
        gene2=c(gene2,as.character(res[1]))
        results[[j]]=list(edges=as.character(res[1]))
        if(scale){
          if(as.numeric(res[3])<0){
            nam=paste0(gene[j],"~",as.character(res[1]))
            aa=c("tee")
            names(aa)=nam
            arrow=c(arrow,aa)
            co=c("red")
            names(co)=nam
            color=c(color,co)

          }else{
            nam=paste0(gene[j],"~",as.character(res[1]))
            aa=c("open")
            names(aa)=nam
            arrow=c(arrow,aa)
            co=c("blue")
            names(co)=nam
            color=c(color,co)
          }
        }else{
          if(as.numeric(res[3])<1){
            nam=paste0(gene[j],"~",as.character(res[1]))
            aa=c("tee")
            names(aa)=nam
            arrow=c(arrow,aa)
            co=c("red")
            names(co)=nam
            color=c(color,co)

          }else{
            nam=paste0(gene[j],"~",as.character(res[1]))
            aa=c("open")
            names(aa)=nam
            arrow=c(arrow,aa)
            co=c("blue")
            names(co)=nam
            color=c(color,co)

          }
        }


      }else{
        gr=as.character(res[,1])
        gene2=c(gene2,gr)
      results[[j]]=list(edges=gr)
      weight=as.numeric(res[,3])

      if(scale){
        ww=which(weight<0)
        wh=which(weight>0)
        if(length(ww)>0)
         { nam=paste0(gene[j],"~",as.character(gr[ww]))
          aa=rep("tee",length(gr[ww]))
          names(aa)=nam
          arrow=c(arrow,aa)
          co=rep("red",length(gr[ww]))
          names(co)=nam
          color=c(color,co)
          }
        if(length(wh)>0)
          {nam=paste0(gene[j],"~",as.character(gr[wh]))
          aa=rep("open",length(gr[wh]))
          names(aa)=nam
          arrow=c(arrow,aa)
          co=rep("blue",length(gr[wh]))
          names(co)=nam
          color=c(color,co)}

        }else{

          ww=which(weight<1)
          wh=which(weight>1)

          if(length(ww)>0)
          { nam=paste0(gene[j],"~",as.character(gr[ww]))
          aa=rep("tee",length(gr[ww]))
          names(aa)=nam
          arrow=c(arrow,aa)
          co=rep("red",length(gr[ww]))
          names(co)=nam
          color=c(color,co)
          }
          if(length(wh)>0)
          {nam=paste0(gene[j],"~",as.character(gr[wh]))
          aa=rep("open",length(gr[wh]))
          names(aa)=nam
          arrow=c(arrow,aa)
          co=rep("blue",length(gr[wh]))
          names(co)=nam
          color=c(color,co)}
        }

      }




    }else{

      results[[j]]=list()
      gene1=c(gene1,gene[j])

    }
      }
    }else{
      results[[j]]=list()
      gene1=c(gene1,gene[j])

    }
}
  nodes=as.character(gene)
  fill <- rep("lightblue", length(nodes))
  names(fill) <- nodes
  edgeList=results
  rem=setdiff(gene1,gene2)
  graph <- new("graphNEL",nodes=nodes, edgeL=edgeList, edgemode="directed")
  if(length(rem)>0)
  {graph<- removeNode(rem, graph)}
  rag=layoutGraph(graph,attrs=list(graph=list(rankdir="LR")))
  nodeRenderInfo(rag) <- list(fontsize=5,cex=1.5,fill=fill,shape="ellipse")
  edgeRenderInfo(rag) <- list(lty=1, lwd=2, col=color,arrowhead=arrow,arrowtail="none")

  #plot
  renderGraph(rag)
  legend("bottomright",c("down-regulation","up-regulation"),y.intersp=1,xpd=T,inset = -0.2, pch=c("|",">"),col=c("red","blue"),bty="n", title= "Edge type", horiz=F,lwd=2)
}
