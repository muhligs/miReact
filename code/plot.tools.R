# plot tools
#################################
# used to create a heatplot style key legend with colors and values and title
# width, height and relative position is settable. Use values slightly below/above 0/1 for margin position.
# all sizes are relative to plot area.
setKey <- function(cols,keylength=0.25,label=NULL, title="", title.cex=1,position=c(0,0), height=0.05){
  xrange <- par("usr")[2]-par("usr")[1]
  yrange <- par("usr")[4]-par("usr")[3]
  ncols <- length(cols)
  increment <- keylength/ncols
  llcornerx <- par("usr")[1]+xrange*position[1]
  llcornery <- par("usr")[3]+yrange*position[2]
  rect(xleft = llcornerx+xrange*increment*1:ncols,ybottom = llcornery,xright = llcornerx+xrange*increment*2:(ncols+1),ytop = llcornery+yrange*height,col = cols,border = NA,xpd=T)
  text(xpd=T,llcornerx+xrange*(keylength/2),llcornery+yrange*height,title,pos = 3)
  if(!is.null(label))text(xpd=T,llcornerx+xrange*(keylength*(0:(length(label)-1))/(length(label)-1)),llcornery,label,pos = 1)
}
# example usage
#setKey(colorRampPalette(c("blue","gray","red"))(100), title="quantile expression",pos=c(-0.1,-0.1),label=c("0","1/2","1"),height=0.05,keylength = 0.25)
#################################

#################################
colorme <- function(data,ramp=c("darkblue","white","darkred"),numcol=99, quantcol=F){
if(quantcol) data <- ecdf(data)(data)
data <- cut(data,breaks=seq(min(data),max(data),len=numcol+1), include.lowest=T)
colors <- colorRampPalette(ramp)(numcol)[data]
return(colors)
}

#################################

#################################
bardot <- function(data,groups,col=NULL,cex=0.5,pch=20,ylab="expression",main="",ylims=NULL,srt=90,labelpos=0.3,labelcex=1, ylabcex=1, labeladj = 0.25, labfont=1,ltext=NULL){
  set.seed(2)  
  if(is.null(col)) col <- sample(grDevices::rainbow(length(unique(groups))))[match(groups,unique(groups))]
  
  yrange <- range(data)
  yspan <- diff(yrange)
  if(is.null(ylims)) ylims=c(yrange[1]-yspan*0.3,yrange[2])
  plot(match(groups,unique(groups))+runif(length(groups),-0.3,0.3),data,col=col,cex=cex,pch=pch,frame.plot = F, axes=F,ylab=ylab,ylim=ylims,xlab="",main=main,cex.lab=ylabcex)
  axis(2)
  if(is.null(ltext)) ltext <- unique(groups) else ltext <- ltext
  text(ltext,x = 1:length(unique(groups))+labeladj,y = yrange[1]-yspan*labelpos,srt=srt,xpd=NA,col=col[match(unique(groups),groups)], cex=labelcex,pos=2,font=labfont)
}
#################################

#################################
centroid <- function(data,groups){
  res <- t(as.data.frame(lapply(unique(groups),function(x)colMeans(data[groups==x,]))))
  rownames(res) <- unique(groups)
  return(res)
}
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################