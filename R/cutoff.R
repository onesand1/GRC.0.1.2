#' @title cutoff for survival analysis
#' @description Desgin three method to identify the cutoff for survival analysis, including median,  quartile and Youden index.
#' @param data The data needs three columns: expression value, survival times and status. The status indicator, normally 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death).
#' @param type Type refers to the segmentation threshold of survival analysis, which is divided into three types: median(M), quartile (Q) and Youden index (Y).default type="M".
#' @param time This parameter refers to the units of time, days(d), months(m), and years(y).default time='d'.
#' @param cut This parameter is a threshold used when selecting type= "Y". marker values to use as a cut-off for calculation of sensitivity and specificity.
#' @param plot A logical parameter. This parameter is defined to show whether the survival curve is visualized.default plot=TRUE.
#' @param file The survival curve plot file name. The file type is .pdf.


cutoff<-function(data,type="M",ttype="d",cut=1,plot=TRUE,file="Rplot"){
  data=data.frame(data)
  dd=list()
  expression=as.numeric(data[,1])
group=rep(0,length(expression))
  if(toupper(type)=="M"){
    label="The data is grouped based on the median"
    med=median(expression)
    low=which(expression<=med)
    high=which(expression>med)
    if(ttype=="d"){
      lab="days"

    }
    if(ttype=="m"){

      lab="months"
    }
    if(ttype=="y"){

      lab="years"
    }
  }
  if(toupper(type)=="Q"){
    label="The data is grouped based on the quartile (25% and 75%)"
    med1=quantile(expression)[2]
    med2=quantile(expression)[4]
    high=which(expression>med2)
    low=which(espression<med1)
    if(ttype=="d"){
      lab="days"

    }
    if(ttype=="m"){

      lab="months"
    }
    if(ttype=="y"){

      lab="years"
    }
  }
  if(toupper(type)=="Y"){
    label="The data is grouped based on the Youden index"
    status=as.character(data[,3])
    status[!is.na(status)]=1
    status[is.na(status)]=0
    time=as.numeric(data[,2])
    if(ttype=="d"){
      lab="days"
      time=time
    }
    if(ttype=="m"){
      time=time*30
      lab="months"
    }
    if(ttype=="y"){
      time=time*30*12
      lab="years"
    }
    cut1=cut*365

    ROC=survivalROC(Stime = time,status = status,marker = expression,predict.time = cut1,method= "KM" )
    cut.op= ROC$cut.values[which.max(ROC$TP-ROC$FP)]
    plot(ROC$FP,ROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),
    xlab = paste( "FP","\n", "AUC = ",round(ROC$AUC,3)),
    ylab = "TP",main = paste(cut,"-year survival ROC",sep=""), col="red")

    abline(0,1)

    legend("bottomright",c("Gene expression"),col="red",lty=c(1,1))
    high=which(expression>cut.op)
    low=which(expression<cut.op)

  }
group[high]="High"
group[low]="Low"
data=cbind(data,group)
colnames(data)=c("expression","time","status","group")
na=which(group!=0)
data=data[na,]
fit=survfit(Surv(time,status)~group,data=data)
dd=list(fit=fit,data=data)
if(plot){
  ggsurvplot(fit, data = data, pval = TRUE, pval.method = TRUE,risk.table = F,xlab=paste("Times (",lab,")",sep=""),surv.median.line = "hv")

  ggsave(paste(file,".pdf",sep=""))
  cat("Warning:\n\n")
  cat(paste(label,".\n",sep=""))
  cat(paste("The survival curve saved in '",file,".pdf'.",sep=""))
  }
invisible(dd)
}

