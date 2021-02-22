#'  plot for regulation and correlation function
#' @title  Visual function for analysis
#' @export
#'
#' @param object the regulator function,GL or GP result object.
#' @param type the figure type. line represent the regulated genes in line; circle represent the regulated genes in circle.
#' @param  r radius for every gene backgroud.

plot.reg.down<-function(object,r=2,type="line"){
  gene=object$target
  target=object$root
  reg=object$reg

  n=length(gene)

  x=seq(15,115,length.out = n)
  plot(1,2,xlim=c(0,130),ylim=c(0,110),type="n",axes=F,xlab ="", ylab = "")

  draw.circle(50,95,r,border="green",col="green",lty=1,density=50,angle=30,lwd=10)
  text(50,95,target)
  if(type=="line"){
  for(i in 1:n){
    if(reg[i]>0){
      angle=30
      col="green"
      y1=22+r
    }else{
      angle=90
      col="red"
      y1=23+r
    }

    ##code: 1-end，2-start，3-all；length: end arrow length; angle: arrow angle
    arrows(x0 = 50, y0 = 93, x1 = x[i], y1 = y1, length = 0.1,col = "pink",lwd=2,lty=1,angle=angle)
    draw.circle(x[i],20,r,border=col,col=col,lty=1,density=50,angle=30,lwd=10)
    text(x[i],20,as.character(gene[i]))
  }
  }
  if(type=="circle"){
    k=seq(1.2*pi,1.8*pi,length.out = n)
    y=70*sin(k)+90
    for(i in 1:n){
      if(reg[i]){
        angle=30
        col="green"
        y1=y[i]+r+2
      }else{
        angle=90
        col="red"
        y1=y[i]+3+r
      }

      ##code: 1-end，2-start，3-all；length: end arrow length; angle: arrow angle
      arrows(x0 = 50, y0 = 93, x1 = x[i], y1 =y1, length = 0.1,col = "pink",angle=angle)
      draw.circle(x[i],y[i],r,border=col,col=col,lty=1,density=50,angle=30,lwd=10)
      text(x[i],y[i],as.character(gene[i]))
    }


  }

}


