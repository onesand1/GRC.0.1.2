#'  plot for regulation and correlation function
#' @title  Visual function for analysis
#' @export
#'
#' @param object the regulator function,GL or GP result object.
#' @param type the figure type. line represent the regulated genes in line; circle represent the regulated genes in circle.
#' @param  r radius for every gene backgroud.

plot.reg.up<-function(object,r=2,type="line"){
  gene=object$target
  target=object$root
  reg=object$reg

  n=length(gene)

  x=seq(1,100,length.out = n)
  plot(1,2,xlim=c(0,110),ylim=c(0,110),type="n",axes=F,xlab ="", ylab = "")

  draw.circle(50,20,r,border="green",col="green",lty=1,density=50,angle=30,lwd=10)
  text(50,20,target)
  if(type=="line"){
    k=seq(0.2*pi,0.8*pi,length.out = n)
    y1=r*sin(k)+20
    x1=seq(48,52,length.out = n)
    for(i in 1:n){
      if(reg[i]>0){
        angle=30
        col="green"
        y2=y1[i]+3

      }else{
        angle=90
        col="red"
        y2=y1[i]+5
      }

      ##code: 1-end，2-start，3-all；length: end arrow length; angle: arrow angle
      arrows(x0 = x[i], y0 = 93, x1 =x1[i], y1 =y2, length = 0.1,col = "pink",lwd=2,lty=1,angle=angle)
      draw.circle(x[i],93,r,border=col,col=col,lty=1,density=50,angle=30,lwd=10)
      text(x[i],93,as.character(gene[i]))
    }
  }

  if(type=="circle"){
    k=seq(0.2*pi,0.8*pi,length.out = n)
    y=70*sin(k)+20
    y1=r*sin(k)+20

    x1=seq(48,52,length.out = n)
    for(i in 1:n){
      if(reg[i]){
        angle=30
        col="green"
        y2=y1[i]+3
        y3=y[i]-r-2
      }else{
        angle=90
        col="red"
        y2=y1[i]+5
        y3=y[i]-r-3
      }

      ##code: 1-end，2-start，3-all；length: end arrow length; angle: arrow angle
      arrows(x0 =x[i], y0 = y3, x1 =x1[i], y1 =y2, length = 0.1,col = "pink",angle=angle)
      draw.circle(x[i],y[i],r,border=col,col=col,lty=1,density=50,angle=30,lwd=10)
      text(x[i],y[i],as.character(gene[i]))
    }


  }

}


