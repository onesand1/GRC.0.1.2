#' @title Diferential analysis of mRNA expression
#'
#' @description
#' @param group a character of group column name for sample groups.
#' @param data a data frame,including 3 or more columns ('sample','gene name','group') for one gene or many genes.
#' @param exp_type the expression data formula mode. including 'RPKM','FPKM','TPM'.

SL<-function(data,group="group",exp_type="rpkm",plot=TRUE){
  if(plot){
    if(!dir.exists("DE_plot")){
    dir.create("DE_plot")
}
  }
  nn=colnames(data)
  dd=which(nn==group)
  expression=nn[-c(1,dd)]
  result=list(gene=expression,P=c(),FC=c(),q_P=c(),FC_type="",method="")
   if(dim(data)[1]<4){
    stop("The number of samples not more than 3 !")

  }else{
    gg=unique(as.character(data[,group]))
    if(length(gg)!=2){
      stop("The groups need to have 2 levels (2 groups) !")
    }
    mut=as.character(data[,group])
    mut_s=which(mut==gg[1])
    wild_s=which(mut==gg[2])
    se=length(mut_s)+length(wild_s)


    if(se<4 | length(mut)<4 | length(wild_s)<4){
      stop("The number of samples not more tha 3 !")

    }else{
      data=data[c(mut_s,wild_s),]
      m=length(expression)
      if(exp_type=="rpkm" | exp_type=="fpkm"){
        result$FC_type="log2(Fold Change)"
        result$method="t-test"
        m=length(expression)
        if(m==1){
          exp=as.numeric(data[,expression])
          test=t.test(log2(exp[mut_s]+1),log2(exp[wild_s]+1),var=TRUE)
          result$P=test$p.value
          result$FC=mean(log2(exp[mut_s]+1))-mean(log2(exp[wild_s]+1))
          if(plot){
          pdf(paste("plot/",expression,"_box.pdf",sep=""))
          boxplot(exp[mut_s],exp[wild_s],col=c("red","blue"),names=c(gg[1],gg[2]), ylab="Expression level", notch =TRUE,main=expression)
          p=paste("P-value = ",format(test$p.value, scientific = TRUE,digits=4),sep="")
          fold=paste(p,", Log2(Fold change)=",round(result$FC,4),sep="")
          legend("top",fold,bty="n")
          dev.off()
          }

        }else{
          for(i in 1:m){
            exp=as.numeric(data[,expression[i]])
            test=t.test(log2(expression[mut_s]+1),log2(expression[wild_s]+1),var=TRUE)
            result$P=c(result$P,test$p.value)
            result$FC=c(result$FC,mean(log2(expression[mut_s]+1))-mean(log2(expression[wild_s]+1)))
            if(plot){
            pdf(paste("plot/",expression[i],"_box.pdf",sep=""))
            boxplot(exp[mut_s],exp[wild_s],col=c("red","blue"),names=c(gg[1],gg[2]), ylab="Expression level", notch =TRUE,main=expression[i])
            p=paste("P-value = ",format(test$p.value, scientific = TRUE,digits=4),sep="")
            fold=paste(p,", Log2(Fold change)=",round(result$FC,4),sep="")
            legend("top",fold,bty="n")
            dev.off()
            }
          }
		  orderP=order(as.numeric(result$P))
          fdr<-rep(1,length(group));
          for(j in 1:length(group)) {
            fdr[j]<-as.numeric(result$P[j])*length(group)/orderP[j]
}
        result$q_P=fdr
      }}else{
        result$FC_type="Wilcoxon two.sided"
        result$method="Wilcoxon"
        m=length(expression)
        if(m==1){
        test=wilcox.test(expression[mut_s],expression[wild_s],alternative="two.sided")
        result$P=test$p.value
        if(plot){
        pdf(paste("plot/",expression,"_box.pdf",sep=""))
        boxplot(exp[mut_s],exp[wild_s],col=c("red","blue"),names=c(gg[1],gg[2]), ylab="Expression level", notch =TRUE,main=expression)
        p=paste("P-value = ",format(test$p.value, scientific = TRUE,digits=4),sep="")
        fold=paste(p,", Two.side",sep="")
        legend("top",fold,bty="n")
        dev.off()
        }
        }else{

          for(i in 1:m){
            exp=as.numeric(data[,expression[i]])
            test=t.test(log2(expression[mut_s]+1),log2(expression[wild_s]+1),var=TRUE)
            result$P=c(result$P,test$p.value)
            if(plot){
            pdf(paste("plot/",expression[i],"_box.pdf",sep=""))
            boxplot(exp[mut_s],exp[wild_s],col=c("red","blue"),names=c(gg[1],gg[2]), ylab="Expression level", notch =TRUE,main=expression[i])
            p=paste("P-value = ",format(test$p.value, scientific = TRUE,digits=4),sep="")
            fold=paste(p,", Two.side",sep="")
            legend("top",fold,bty="n")
            dev.off()
}
          }
		  orderP=order(as.numeric(result$P))
          fdr<-rep(1,length(group));
          for(j in 1:length(group)) {
            fdr[j]<-as.numeric(result$P[j])*length(group)/orderP[j]
          }
          result$q_P=fdr
}

      }
    result$gene=expression
    return(result)

  }

}
}
