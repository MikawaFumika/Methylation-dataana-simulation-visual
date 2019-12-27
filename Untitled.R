
hist(Mcpgcov2i[,4],xlim=c(0,10000),breaks=50,mgp=c(3, 0.5, -0.33),tck=0.03,las=1,xlab="Enhancer length",ylab="Frequency",main="",col="grey")

#2i细胞的
#处理甲基化enhancer位点
setwd("/Users/yeyusong/Desktop/在投/Methylation/data")
methylationsite=read.csv("methypoint.csv",header=T,check.names=FALSE,sep="\t")
#25036
#建立数据矩阵初始向量
methystart=vector()
methyend=vector()
methyscale=vector()
methychr=vector()
for(chr in 1:19)
{
  #制作methylation对应表
  chrmethystart=vector()
  chrmethystart=methylationsite[which(methylationsite[,1]==chr),2]
  chrmethyend=vector()
  chrmethyend=methylationsite[which(methylationsite[,1]==chr),3]
  chrmethyscale=vector()
  chrmethyscale=methylationsite[which(methylationsite[,1]==chr),4]
  chromatin=rep(chr,length.out = length(chrmethyscale))
  methystart=c(methystart,chrmethystart)
  methyend=c(methyend,chrmethyend)
  methyscale=c(methyscale,chrmethyscale)
  methychr=c(methychr,chromatin)
}
chrmethystart=vector()
chrmethystart=methylationsite[which(methylationsite[,1]=="X"),2]
chrmethyend=vector()
chrmethyend=methylationsite[which(methylationsite[,1]=="X"),3]
chrmethyscale=vector()
chrmethyscale=methylationsite[which(methylationsite[,1]=="X"),4]
chromatin=rep("X",length.out = length(chrmethyscale))
methystart=c(methystart,chrmethystart)
methyend=c(methyend,chrmethyend)
methyscale=c(methyscale,chrmethyscale)
methychr=c(methychr,chromatin)

#建立矩阵
Mcpgcov2i <- matrix(methychr)
Mcpgcov2i=data.frame(Mcpgcov2i,methystart)
Mcpgcov2i=data.frame(Mcpgcov2i,methyend)
Mcpgcov2i=data.frame(Mcpgcov2i,methyscale)


#处理2i 数据
setwd("/Users/yeyusong/Desktop/在投/Methylation/data/2i")
filename <- list.files()


for (txt in 1:1)
{
  cell2i=read.table(filename[txt],header=F,check.names=FALSE,sep="\t")
  
  enhancerfraction=vector()
  enhancercpgnum=vector()
  #染色体循环
  for(chr in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X"))
  {
    #制作细胞染色质甲基化对应表 注意45 56
    cell2ichrmethysite=vector()
    cell2ichrmethysite=as.numeric(as.character(cell2i[which(cell2i[,1]==chr),2]))
    cell2ichrmethyfraction=vector()
    cell2ichrmethyfraction=as.numeric(as.character(cell2i[which(cell2i[,1]==chr),4]))
    
    #排序矩阵
    SortMcpgcov2i<- matrix(cell2ichrmethysite)
    SortMcpgcov2i=data.frame(SortMcpgcov2i,cell2ichrmethyfraction)
    SortMcpgcov2i=SortMcpgcov2i[order(SortMcpgcov2i[,1]),]    
    cell2ichrmethysite=SortMcpgcov2i[,1]
    cell2ichrmethyfraction=SortMcpgcov2i[,2]
    
    lengthofchr=length(which(Mcpgcov2i[,1]==chr))
    #匹配计算
    for(i in 1:lengthofchr)
    {
      enhancerpointstart=min(which(cell2ichrmethysite>=Mcpgcov2i[i,2]))
      enhancerpointend=max(which(cell2ichrmethysite<=Mcpgcov2i[i,3]))
      if(enhancerpointstart<=enhancerpointend & enhancerpointstart<100000000 & enhancerpointend<100000000)
      {cpgfraction=cell2ichrmethyfraction[enhancerpointstart:enhancerpointend]
      cpgnumber=length(cpgfraction)
      enhancercpgnum=c(enhancercpgnum,cpgnumber)
      enhancerfraction=c(enhancerfraction,sum(cpgfraction)/cpgnumber)
      }
      else
      {
      enhancercpgnum=c(enhancercpgnum,0)
      enhancerfraction=c(enhancerfraction,"NA")
      }
    }
  }
  
  Mcpgcov2i=data.frame(Mcpgcov2i,enhancercpgnum)
  Mcpgcov2i=data.frame(Mcpgcov2i,enhancerfraction)

}

#存储数据
setwd("/Users/yeyusong/Desktop/在投/Methylation/data/Processeddata")
write.table(Mcpgcov2i,"Mcpgcov2i1.txt",quote = FALSE,row.names = FALSE,col.names = T,sep="\t")

Mcpgcov2i=read.csv("Mcpgcov2i.txt",header=T,check.names=FALSE,sep="\t")


#ser细胞的

#处理甲基化enhancer位点
setwd("/Users/yeyusong/Desktop/在投/Methylation/data")
methylationsite=read.csv("methypoint.csv",header=T,check.names=FALSE,sep="\t")
#25036
#建立数据矩阵初始向量
methystart=vector()
methyend=vector()
methyscale=vector()
methychr=vector()
for(chr in 1:19)
{
  #制作methylation对应表
  chrmethystart=vector()
  chrmethystart=methylationsite[which(methylationsite[,1]==chr),2]
  chrmethyend=vector()
  chrmethyend=methylationsite[which(methylationsite[,1]==chr),3]
  chrmethyscale=vector()
  chrmethyscale=methylationsite[which(methylationsite[,1]==chr),4]
  chromatin=rep(chr,length.out = length(chrmethyscale))
  methystart=c(methystart,chrmethystart)
  methyend=c(methyend,chrmethyend)
  methyscale=c(methyscale,chrmethyscale)
  methychr=c(methychr,chromatin)
}
chrmethystart=vector()
chrmethystart=methylationsite[which(methylationsite[,1]=="X"),2]
chrmethyend=vector()
chrmethyend=methylationsite[which(methylationsite[,1]=="X"),3]
chrmethyscale=vector()
chrmethyscale=methylationsite[which(methylationsite[,1]=="X"),4]
chromatin=rep("X",length.out = length(chrmethyscale))
methystart=c(methystart,chrmethystart)
methyend=c(methyend,chrmethyend)
methyscale=c(methyscale,chrmethyscale)
methychr=c(methychr,chromatin)

#建立矩阵
Mcpgcovser<- matrix(methychr)
Mcpgcovser=data.frame(Mcpgcovser,methystart)
Mcpgcovser=data.frame(Mcpgcovser,methyend)
Mcpgcovser=data.frame(Mcpgcovser,methyscale)



#处理ser 数据
setwd("/Users/yeyusong/Desktop/在投/Methylation/data/ser")
filename <- list.files()


for (txt in 1:20)
{
  cellser=read.table(filename[txt],header=F,check.names=FALSE,sep="\t")
  
  enhancerfraction=vector()
  enhancercpgnum=vector()
  #染色体循环
  for(chr in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X"))
  {
    #制作细胞染色质甲基化对应表 注意45 56
    cellserchrmethysite=vector()
    cellserchrmethysite=as.numeric(as.character(cellser[which(cellser[,1]==chr),2]))
    cellserchrmethyfraction=vector()
    cellserchrmethyfraction=as.numeric(as.character(cellser[which(cellser[,1]==chr),4]))
    
    #排序矩阵
    SortMcpgcovser<- matrix(cellserchrmethysite)
    SortMcpgcovser=data.frame(SortMcpgcovser,cellserchrmethyfraction)
    SortMcpgcovser=SortMcpgcovser[order(SortMcpgcovser[,1]),]    
    cellserchrmethysite=SortMcpgcovser[,1]
    cellserchrmethyfraction=SortMcpgcovser[,2]
    
    lengthofchr=length(which(Mcpgcovser[,1]==chr))
    #匹配计算
    for(i in 1:lengthofchr)
    {
      enhancerpointstart=min(which(cellserchrmethysite>=Mcpgcovser[i,2]))
      enhancerpointend=max(which(cellserchrmethysite<=Mcpgcovser[i,3]))
      if(enhancerpointstart<=enhancerpointend & enhancerpointstart<1000000000 & enhancerpointend<1000000000)
      {cpgfraction=cellserchrmethyfraction[enhancerpointstart:enhancerpointend]
      cpgnumber=length(cpgfraction)
      enhancercpgnum=c(enhancercpgnum,cpgnumber)
      enhancerfraction=c(enhancerfraction,sum(cpgfraction)/cpgnumber)
      }
      else
      {
        enhancercpgnum=c(enhancercpgnum,0)
        enhancerfraction=c(enhancerfraction,"NA")
      }
    }
  }
  
  Mcpgcovser=data.frame(Mcpgcovser,enhancercpgnum)
  Mcpgcovser=data.frame(Mcpgcovser,enhancerfraction)
  
}
cellcpgnumberser=vector()
for(cells in 1:20)
{
  cellser=read.table(filename[cells],header=F,check.names=FALSE,sep="\t")
  cellcpgnumberser[cells]=length(cellser[,1])
}


#存储数据
setwd("/Users/yeyusong/Desktop/在投/Methylation/data/Processeddata")
write.table(Mcpgcovser,"Mcpgcovser.txt",quote = FALSE,row.names = FALSE,col.names = T,sep="\t")

Mcpgcovser=read.csv("Mcpgcovser.txt",header=T,check.names=FALSE,sep="\t")




#2i release

setwd("/Users/yeyusong/Desktop/在投/Methylation/data/")
methylationsite=read.csv("methypoint.csv",header=T,check.names=FALSE,sep="\t")

Mcpgcovrelease=vector()

for(i in 1:1667)
{
  Mcpgcovrelease[i]=1
}
Mcpgcovrelease=data.frame(Mcpgcovrelease)


setwd("/Users/yeyusong/Desktop/在投/Methylation/data/release/")
filename <- list.files()


for (txt in 1:14)
{
  release=read.csv(filename[txt],header=F,check.names=FALSE,sep="\t")
  

enhancerfraction=vector()
enhancercpgnum=vector()

  #染色体循环
    #制作细胞染色质甲基化对应表 注意45 56
    releasechrmethysite=vector()
    releasechrmethysite=as.numeric(as.character(release[which(release[,1]==1),2]))
    releasechrmethyfraction=vector()
    releasechrmethyfraction1=vector()
    releasechrmethyfraction2=vector()
    releasechrmethyfraction1=as.numeric(as.character(release[which(release[,1]==1),4]))
    releasechrmethyfraction2=as.numeric(as.character(release[which(release[,1]==1),5]))
    releasechrmethyfraction=releasechrmethyfraction1/(releasechrmethyfraction2+releasechrmethyfraction1)
    
    #排序矩阵
    SortMcpgcovrelease<- matrix(releasechrmethysite)
    SortMcpgcovrelease=data.frame(SortMcpgcovrelease,releasechrmethyfraction)
    SortMcpgcovrelease=SortMcpgcovrelease[order(SortMcpgcovrelease[,1]),]    
    releasechrmethysite=SortMcpgcovrelease[,1]
    releasechrmethyfraction=SortMcpgcovrelease[,2]
    
    lengthofchr=which(as.character(methylationsite[,1])==1)
    #匹配计算
    
    for(i in lengthofchr)
    {
      enhancerpointstart=min(which(releasechrmethysite>=methylationsite[i,2]))
      enhancerpointend=max(which(releasechrmethysite<=methylationsite[i,3]))
      if(enhancerpointstart<=enhancerpointend & enhancerpointstart<100000000 & enhancerpointend<100000000)
      {cpgfraction=releasechrmethyfraction[enhancerpointstart:enhancerpointend]
      cpgnumber=length(which(cpgfraction!="NaN"))
      {
      if(cpgnumber==0)
        {enhancercpgnum=c(enhancercpgnum,0)
        enhancerfraction=c(enhancerfraction,"NA")}
      else
      {
      enhancercpgnum=c(enhancercpgnum,cpgnumber)
      enhancerfraction=c(enhancerfraction,sum(cpgfraction[which(cpgfraction!="NaN")])/cpgnumber)
      }
      }
      }
      else
      {
        enhancercpgnum=c(enhancercpgnum,0)
        enhancerfraction=c(enhancerfraction,"NA")
      }
     print(i)
    }
  
  Mcpgcovrelease=data.frame(Mcpgcovrelease,enhancercpgnum)
  Mcpgcovrelease=data.frame(Mcpgcovrelease,enhancerfraction)
}
write.table(Mcpgcovrelease,"Mcpgcovrelease.txt",quote = FALSE,row.names = FALSE,col.names = T,sep="\t")

releasemethymean=vector()
releasemethyvar=vector()
for(i in 1:14)
{releasemethymean[i]=mean(as.numeric(as.character(Mcpgcovrelease[which(Mcpgcovrelease[,2*i+3]!='NA'),(2*i+3)])))
releasemethyvar[i]=var(as.numeric(as.character(Mcpgcovrelease[which(Mcpgcovrelease[,2*i+3]!='NA'),(2*i+3)])))}
xreleasemethymean=c(0,0.75,1,1.5,2,2.5,3,3.5,4,4.75,5.5,6.5,7.5,8.5)


## 做图


## fig1a

hold=3
allenhancer2i=vector()
for(i in 1:12)
{
  eachcell=as.numeric(as.character(Mcpgcov2i[,4+2*i])[which(as.character(Mcpgcov2i[,4+2*i])!="NA" & as.numeric(as.character(Mcpgcov2i[,3+2*i]))>hold)])/100
  allenhancer2i=c(allenhancer2i,eachcell)
}

meancell2i=vector()
for(i in 1:12)
{
  eachcell=as.numeric(as.character(Mcpgcov2i[,4+2*i])[which(as.character(Mcpgcov2i[,4+2*i])!="NA" & as.numeric(as.character(Mcpgcov2i[,3+2*i]))>hold)])/100
  meancell2i=c(meancell2i,mean(eachcell))
}

allenhancerser=vector()
for(i in 1:20)
{
  eachcell=as.numeric(as.character(Mcpgcovser[,4+2*i])[which(as.character(Mcpgcovser[,4+2*i])!="NA" & as.numeric(as.character(Mcpgcovser[,3+2*i]))>hold)])/100
  allenhancerser=c(allenhancerser,eachcell)
}

meancellser=vector()
for(i in 1:20)
{
  eachcell=as.numeric(as.character(Mcpgcovser[,4+2*i])[which(as.character(Mcpgcovser[,4+2*i])!="NA" & as.numeric(as.character(Mcpgcovser[,3+2*i]))>hold)])/100
  meancellser=c(meancellser,mean(eachcell))
}

sd(allenhancer2i[which(allenhancer2i!=0 & allenhancer2i!=1)])
sd(allenhancerser[which(allenhancerser!=0 & allenhancerser!=1)])
mean(allenhancer2i)
mean(allenhancerser)
cell2idis=hist(allenhancer2i,breaks=10)
cellserdis=hist(allenhancerser,breaks=10)

cell2ix=vector()
for(i in 1:12)
{cell2ix[i]=rnorm(1,3,0.11)}

cellserx=vector()
for(i in 1:20)
{cellserx[i]=rnorm(1,4,0.11)}

#col=rgb(30/255,144/255,255/255) col=rgb(238/255,44/255,44/255)
par(mar=c(1.5,3.5,1,1))
plot(cell2ix,meancell2i,type="p",tck=0.03,las=1,xlab="",col=rgb(30/255,144/255,255/255),pch=19, ylab="", main="",xlim=c(0.5,4.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)
par(new=TRUE)
plot(cellserx,meancellser,type="p",tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),pch=19, ylab="", main="",xlim=c(0.5,4.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)

#小提琴小提琴
#初始化 8,25 rgb(238/255,44/255,44/255)  he  6,9  rgb(30/255,144/255,255/255)   & 205 85 85
yinzi=0.06
par(new=T)
polygon(c(1,(1-yinzi*(cell2idis$density)),1,(1+yinzi*rev(cell2idis$density))),c(0,seq(0.05,0.95,0.1),1,seq(0.95,0.05,-0.1)),density = NULL, border = F, col = rgb(30/255,144/255,255/255),bty='L')
par(new=T)
polygon(c(0.7,0.7,1.3,1.3),c(0.33,0.34,0.34,0.33),density = NULL, border = F, col ='black', bty='L')
par(new=T)
polygon(c(2,(2-yinzi*(cellserdis$density)),2,(2+yinzi*rev(cellserdis$density))),c(0,seq(0.05,0.95,0.1),1,seq(0.95,0.05,-0.1)),density = NULL, border = F, col = rgb(238/255,44/255,44/255),bty='L')
par(new=T)
polygon(c(1.7,1.7,2.3,2.3),c(0.70,0.71,0.71,0.70),density = NULL, border = F, col ='black', bty='L')
mtext("DNA Methylation(%)",side=2,line=2.0,cex=1.0)
axis(side=1,las=1,at=c(1.5,3.5),mgp=c(0,0.5,0),tck=0.0,las=1,labels=c("Enhancers","Cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(0.05,0.5,0.95),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c('0','50','100'),cex.axis=1.0)


## fig1b


#E455565
setwd("/Users/yeyusong/Desktop/在投/Methylation/data")
E45methylationsite=read.csv("E45_methy_fraction.csv",header=T,check.names=FALSE,sep=",")
E55methylationsite=read.csv("E55_methy_fraction.csv",header=T,check.names=FALSE,sep=",")
E65methylationsite=read.csv("E65_methy_fraction.csv",header=T,check.names=FALSE,sep=",")

hold=3
#E45
E45methyleveldistri=vector()
meanE45methylevel=vector()

for(cells in 1:91)
{
  x=(as.numeric(as.character(E45methylationsite[,4+2*cells]))[which(as.numeric(as.character(E45methylationsite[,(4+2*cells-1)]))>hold)])
  E45methyleveldistri=c(E45methyleveldistri,x)
  meanE45methylevel=c(meanE45methylevel,mean(x))
}

#E55
E55methyleveldistri=vector()
meanE55methylevel=vector()

for(cells in 1:80)
{
  x=(as.numeric(as.character(E55methylationsite[,4+2*cells]))[which(as.numeric(as.character(E55methylationsite[,(4+2*cells-1)]))>hold)])
  E55methyleveldistri=c(E55methyleveldistri,x)
  meanE55methylevel=c(meanE55methylevel,mean(x))
}

#E65
E65methyleveldistri=vector()
meanE65methylevel=vector()

for(cells in 1:96)
{
  x=(as.numeric(as.character(E65methylationsite[,4+2*cells]))[which(as.numeric(as.character(E65methylationsite[,(4+2*cells-1)]))>hold)])
  E65methyleveldistri=c(E65methyleveldistri,x)
  meanE65methylevel=c(meanE65methylevel,mean(x))
}

cellE45dis=hist(E45methyleveldistri,breaks=10)
cellE55dis=hist(E55methyleveldistri,breaks=10)
cellE65dis=hist(E65methyleveldistri,breaks=10)
mean(E45methyleveldistri)
mean(E55methyleveldistri)
mean(E65methyleveldistri)

cellE45x=vector()
for(i in 1:91)
{cellE45x[i]=rnorm(1,4,0.11)}

cellE55x=vector()
for(i in 1:80)
{cellE55x[i]=rnorm(1,5,0.11)}

cellE65x=vector()
for(i in 1:96)
{cellE65x[i]=rnorm(1,6,0.11)}

#rgb(30/255,144/255,255/255) rgb(238/255,44/255,44/255) rgb(0/255,139/255,69/255)
par(mar=c(1.5,3.5,1,1))
plot(cellE45x,meanE45methylevel/100,type="p",tck=0.03,las=1,xlab="",col=rgb(30/255,144/255,255/255),pch=19, ylab="", main="",xlim=c(0.5,6.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)
par(new=TRUE)
plot(cellE55x,meanE55methylevel/100,type="p",tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),pch=19, ylab="", main="",xlim=c(0.5,6.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)
par(new=TRUE)
plot(cellE65x,meanE65methylevel/100,type="p",tck=0.03,las=1,xlab="",col=rgb(0/255,139/255,69/255),pch=19, ylab="", main="",xlim=c(0.5,6.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)

#小提琴小提琴

yinzi=6
par(new=T)
polygon(c(1,(1-yinzi*(cellE45dis$density)),1,(1+yinzi*rev(cellE45dis$density))),c(0,seq(0.05,0.95,0.1),1,seq(0.95,0.05,-0.1)),density = NULL, border = F, col = rgb(30/255,144/255,255/255),bty='L')
par(new=T)
polygon(c(0.7,0.7,1.3,1.3),c(0.17,0.18,0.18,0.17),density = NULL, border = F, col ='black', bty='L')
par(new=T)
polygon(c(2,(2-yinzi*(cellE55dis$density)),2,(2+yinzi*rev(cellE55dis$density))),c(0,seq(0.05,0.95,0.1),1,seq(0.95,0.05,-0.1)),density = NULL, border = F, col = rgb(238/255,44/255,44/255),bty='L')
par(new=T)
polygon(c(1.7,1.7,2.3,2.3),c(0.66,0.67,0.67,0.66),density = NULL, border = F, col ='black', bty='L')
par(new=T)
polygon(c(3,(3-yinzi*(cellE65dis$density)),3,(3+yinzi*rev(cellE65dis$density))),c(0,seq(0.05,0.95,0.1),1,seq(0.95,0.05,-0.1)),density = NULL, border = F, col = rgb(0/255,139/255,69/255),bty='L')
par(new=T)
polygon(c(2.7,2.7,3.3,3.3),c(0.67,0.68,0.68,0.67),density = NULL, border = F, col ='black', bty='L')
mtext("DNA Methylation(%)",side=2,line=2.0,cex=1.0)
axis(side=1,las=1,at=c(2,5),mgp=c(0,0.5,0),tck=0.0,las=1,labels=c("Enhancers","Cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(0.05,0.5,0.95),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c('0','50','100'),cex.axis=1.0)

## fig1c
hold=3
enhancerdisnaive1=vector()
i=1
enhancerdisnaive1=as.numeric(as.character(Mcpgcov2i[,4+2*i])[which(as.character(Mcpgcov2i[,4+2*i])!="NA" & as.numeric(as.character(Mcpgcov2i[,3+2*i]))>hold)])/100

enhancerdisnaive2=vector()
i=6
enhancerdisnaive2=as.numeric(as.character(Mcpgcov2i[,4+2*i])[which(as.character(Mcpgcov2i[,4+2*i])!="NA" & as.numeric(as.character(Mcpgcov2i[,3+2*i]))>hold)])/100
#col=rgb(30/255,144/255,255/255) col=rgb(238/255,44/255,44/255)

x=hist(enhancerdisnaive2,breaks=10,tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255), ylab="", main="",xlim=c(0,1),xaxt="n",bty='L',border= T)
enhancerdisnaive2=x$counts/sum(x$counts)
par(mar=c(2.5,3.2,1,1))
plot(1,1,type='l',xlim=c(0.05,1.05),ylim=c(0,0.6),col='white',tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",bty='L',xaxt='n')
for(i in 1:10)
  {
    par(new=T)
    polygon(c(i/10-0.04,i/10-0.04,i/10+0.04,i/10+0.04),c(0,enhancerdisnaive2[i],enhancerdisnaive2[i],0), density = NULL, border = F, col = rgb(30/255,144/255,255/255),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="")
  }
axis(side=1,las=1,at=c(0.05,0.55,1.05),mgp=c(0,0.5,0),tck=0.0,las=1,labels=c(0,50,100),cex.axis=1.0)
mtext("Methylation(%)",side=1,line=1.2,cex=1.0,at=0.55)
mtext("Probability",side=2,line=1.7,cex=1.0)


x=hist(enhancerdisnaive1,breaks=10,tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255), ylab="", main="",xlim=c(0,1),xaxt="n",bty='L',border= T,xaxt='n')
enhancerdisnaive1=x$counts/sum(x$counts)
par(mar=c(2.5,3.2,1,1))
plot(1,1,type='l',xlim=c(0.05,1.05),ylim=c(0,0.6),col='white',tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",bty='L',xaxt="n")
for(i in 1:10)
{
  par(new=T)
  polygon(c(i/10-0.04,i/10-0.04,i/10+0.04,i/10+0.04),c(0,enhancerdisnaive1[i],enhancerdisnaive1[i],0), density = NULL, border = F, col = rgb(30/255,144/255,255/255),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="")
}
axis(side=1,las=1,at=c(0.05,0.55,1.05),mgp=c(0,0.5,0),tck=0.0,las=1,labels=c(0,50,100),cex.axis=1.0)
mtext("Methylation(%)",side=1,line=1.2,cex=1.0,at=0.55)
mtext("Probability",side=2,line=1.7,cex=1.0)



## fig1d
## 1 3 6
hold=3
enhancerdisnaive1=vector()
i=6
enhancerdisnaive1=as.numeric(as.character(Mcpgcovser[,4+2*i])[which(as.character(Mcpgcovser[,4+2*i])!="NA" & as.numeric(as.character(Mcpgcovser[,3+2*i]))>hold)])/100
x=hist(enhancerdisnaive1,breaks=10,tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255), ylab="", main="",xlim=c(0,1),xaxt="n",bty='L',border= T)
enhancerdisnaive1=x$counts/sum(x$counts)
par(mar=c(2.5,3.2,1,1))
plot(1,1,type='l',xlim=c(0.05,1.05),ylim=c(0,0.7),col='white',tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",bty='L',xaxt='n')
for(i in 1:10)
{
  par(new=T)
  polygon(c(i/10-0.04,i/10-0.04,i/10+0.04,i/10+0.04),c(0,enhancerdisnaive1[i],enhancerdisnaive1[i],0), density = NULL, border = F, col =rgb(238/255,44/255,44/255),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="")
}
axis(side=1,las=1,at=c(0.05,0.55,1.05),mgp=c(0,0.5,0),tck=0.0,las=1,labels=c(0,50,100),cex.axis=1.0)
mtext("Methylation(%)",side=1,line=1.2,cex=1.0,at=0.55)
mtext("Probability",side=2,line=1.7,cex=1.0)


##fig1e

interval=5
hold=3
#E45

cellsdis=matrix()
for(cells in 1:91)
{
  x=(as.numeric(as.character(E45methylationsite[,4+2*cells]))[which(as.numeric(as.character(E45methylationsite[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=c(0,seq(5,95,interval),(100)))
    x=x$counts/(sum(x$counts))
    cellsdis=data.frame(cellsdis,x)
  }
  else
    return
}
indexE45=vector()
for(cells in 2:87)
{
  x=cellsdis[,cells]
  for(cellsvs in 2:87)
  {
    y=cellsdis[,cellsvs]
    indexE45=c(indexE45,sum((y)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexE45)


cellsdisE55=matrix()
for(cells in 1:80)
{
  x=(as.numeric(as.character(E55methylationsite[,4+2*cells]))[which(as.numeric(as.character(E55methylationsite[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=c(0,seq(5,95,interval),(100)))
    x=x$counts/(sum(x$counts))
    cellsdisE55=data.frame(cellsdisE55,x)
  }
  else
    return
}
indexE55=vector()
for(cells in 2:78)
{
  x=cellsdisE55[,cells]
  for(cellsvs in 2:78)
  {
    y=cellsdisE55[,cellsvs]
    indexE55=c(indexE55,sum((y)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexE55)


cellsdisE65=matrix()
for(cells in 1:96)
{
  x=(as.numeric(as.character(E65methylationsite[,4+2*cells]))[which(as.numeric(as.character(E65methylationsite[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=c(0,seq(5,95,interval),(100)))
    x=x$counts/(sum(x$counts))
    cellsdisE65=data.frame(cellsdisE65,x)
  }
  else
    return
}
indexE65=vector()
for(cells in 2:96)
{
  x=cellsdisE65[,cells]
  for(cellsvs in 2:96)
  {
    y=cellsdisE65[,cellsvs]
    indexE65=c(indexE65,sum((y)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexE65)

# naive 和 prime 的 index
hold=3
#naive
cellsdisnaive=matrix()
for(cells in 1:12)
{
  x=(as.numeric(as.character(Mcpgcov2i[,4+2*cells]))[which(as.numeric(as.character(Mcpgcov2i[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=seq(0,100,5))
    x=x$counts/(sum(x$counts))
    cellsdisnaive=data.frame(cellsdisnaive,x)
  }
  else
    return
}
indexnaive=vector()
for(cells in 2:13)
{
  x=cellsdisnaive[,cells]
  for(cellsvs in 2:13)
  {
    y=cellsdisnaive[,cellsvs]
    indexnaive=c(indexnaive,sum((y+0.001)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexnaive)
hist(indexnaive)

#prime
cellsdisprime=matrix()
for(cells in 1:20)
{
  x=(as.numeric(as.character(Mcpgcovser[,4+2*cells]))[which(as.numeric(as.character(Mcpgcovser[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=seq(0,100,5))
    x=x$counts/(sum(x$counts))
    cellsdisprime=data.frame(cellsdisprime,x)
  }
  else
    return
}
indexprime=vector()
for(cells in 2:21)
{
  x=cellsdisprime[,cells]
  for(cellsvs in 2:21)
  {
    y=cellsdisprime[,cellsvs]
    indexprime=c(indexprime,sum((y+0.001)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexprime)
hist(indexprime)

#drawing
par(mar=c(1.7,3,0.5,0))
indexx=c(3,4,5)
indexy=c(mean(indexE45),mean(indexE55),mean(indexE65))
plot(indexx,indexy,type='l',col='grey',tck=0.03,las=1,cex=0.3,xlim=c(2.5,5.5),ylim=c(0,0.4),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexx,indexy,type='p',pch=19,col='black',bg=rgb(102/255,205/255,0/255),tck=0.03,las=1,cex=1,xlim=c(2.5,5.5),ylim=c(0,0.4),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
axis(side=1,las=1,at=c(3,4,5),mgp=c(0,0.5,0),tck=0.0,las=1,labels=c("E4.5","E5.5","E6.5"),cex.axis=1.0)
axis(side=2,las=1,at=c(0,0.2,0.4),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.2","0.4"),cex.axis=1.0)
mtext("HI",side=2,line=1.5,cex=1.2)


par(mar=c(1.7,3,0.5,0))
indexx=c(1,2)
par(new=T)
indexy=c(mean(indexnaive),mean(indexprime)+0.02)
plot(indexx,indexy,type='l',col='grey',tck=0.03,las=1,cex=0.3,xlim=c(0.5,2.5),ylim=c(0,0.4),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexx,indexy,type='p',pch=19,col='black',bg=rgb(102/255,205/255,0/255),tck=0.03,las=1,cex=1,xlim=c(0.5,2.5),ylim=c(0,0.4),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
axis(side=1,las=1,at=c(1,2),mgp=c(0,0.5,0),tck=0.0,las=1,labels=c('Naive','Prime'),cex.axis=1.0)
axis(side=2,las=1,at=c(0,0.2,0.4),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.2","0.4"),cex.axis=1.0)
mtext("HI",side=2,line=1.5,cex=1.2)

#fig1f
par(mar=c(3,3,0.5,0.5))
plot(xreleasemethymean,releasemethymean,type='l',xlim=c(0,8.5),ylim=c(0.08,0.18),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",xaxt='n',yaxt='n',bty="l")
par(new=T)
polygon(c(xreleasemethymean,rev(xreleasemethymean)),c(releasemethymean+0.5*releasemethyvar,rev(releasemethymean-0.5*releasemethyvar)), density = NULL, border = F, col = rgb(240/255,240/255,240/255),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="Time(h)", ylab="Methylation(%)", main="",xaxt='n',yaxt='n')
par(new=T)
plot(xreleasemethymean,releasemethymean,xlim=c(0,8.5),ylim=c(0.08,0.18),cex=1,pch=16,tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",xaxt='n',yaxt='n',bty="l")
par(new=T)
plot(xreleasemethymean,releasemethymean,type='l',xlim=c(0,8.5),ylim=c(0.08,0.18),col='grey',tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",xaxt='n',yaxt='n',bty="l")
axis(side=1,las=1,at=c(0,1,2,3,4,5,6,7,8),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(0,8,1),cex.axis=1.0)
axis(side=2,las=1,at=c(0.08,0.13,0.18),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(8,13,18),cex.axis=1.0)
mtext("Methylation(%)",side=2,line=1.5,cex=1)
mtext("Time(h)",side=1,line=1.5,cex=1)


#fig2a
#2i
enhancercontaincpgnumberx=seq(1,49,2)
enhancercontaincpgnumber=vector()
cells=1
for(cells in 1:12)
{
  enhancercontaincpgnumbery=hist(as.numeric(as.character(Mcpgcov2i[,3+2*cells]))[as.numeric(as.character(Mcpgcov2i[,3+2*cells]))>0 & as.numeric(as.character(Mcpgcov2i[,3+2*cells]))<50],breaks=20,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="enhancer cpg sites",ylab="Frequency",main="",col="grey")
  if(length(enhancercontaincpgnumbery$density)<25)
  {
    bluk=25-length(enhancercontaincpgnumbery$density)
    enhancercontaincpgnumber=c(enhancercontaincpgnumber,enhancercontaincpgnumbery$counts/sum(enhancercontaincpgnumbery$counts),seq(length=bluk,from=0,to=0))
  }
  else
    enhancercontaincpgnumber=c(enhancercontaincpgnumber,enhancercontaincpgnumbery$density/sum(enhancercontaincpgnumbery$density))
}

#ser
enhancercontaincpgnumberser=vector()
for(cells in 1:20)
{
  enhancercontaincpgnumbery=hist(as.numeric(as.character(Mcpgcovser[,3+2*cells]))[as.numeric(as.character(Mcpgcovser[,3+2*cells]))>0 & as.numeric(as.character(Mcpgcovser[,3+2*cells]))<50],breaks=20,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="enhancer cpg sites",ylab="Frequency",main="",col="grey")
  if(length(enhancercontaincpgnumbery$density)<25)
  {
    bluk=25-length(enhancercontaincpgnumbery$density)
    enhancercontaincpgnumberser=c(enhancercontaincpgnumberser,enhancercontaincpgnumbery$counts/(sum(enhancercontaincpgnumbery$counts)),seq(length=bluk,from=0,to=0))
  }
  else
    enhancercontaincpgnumberser=c(enhancercontaincpgnumberser,enhancercontaincpgnumbery$counts/(sum(enhancercontaincpgnumbery$counts)))
}

##yinyinghanshu
par(mar=c(3,4,1,1))
x2ienhancercontain=vector()
y2ienhancercontain=vector()
for(i in 1:25)
  x2ienhancercontain[i]=1
for(i in 1:25)
  y2ienhancercontain[i]=0
for(j in 1:25)
{
  for(i in 1:12)
  {
    if(enhancercontaincpgnumber[j+25*i-25]>y2ienhancercontain[j])
    {y2ienhancercontain[j]=enhancercontaincpgnumber[j+25*i-25]}
  }
} 
for(j in 1:25)
{
  for(i in 1:12)
  {
    if(enhancercontaincpgnumber[j+25*i-25]<x2ienhancercontain[j])
    {x2ienhancercontain[j]=enhancercontaincpgnumber[j+25*i-25]}
  }
} 
x2ienhancercontain=x2ienhancercontain+0.0001

x2ienhancercontain=vector()
y2ienhancercontain=vector()
for(i in 1:25)
  x2ienhancercontain[i]=1
for(i in 1:25)
  y2ienhancercontain[i]=0
for(j in 1:25)
{
  for(i in 1:20)
  {
    if(enhancercontaincpgnumberser[j+25*i-25]>y2ienhancercontain[j])
    {y2ienhancercontain[j]=enhancercontaincpgnumberser[j+25*i-25]}
  }
} 
for(j in 1:25)
{
  for(i in 1:20)
  {
    if(enhancercontaincpgnumberser[j+25*i-25]<x2ienhancercontain[j])
    {x2ienhancercontain[j]=enhancercontaincpgnumberser[j+25*i-25]}
  }
} 
x2ienhancercontain=x2ienhancercontain+0.0001

par(mar=c(2.5,3,0.5,0.5))
fitcpgsitesinenhancerx=seq(1,50,1)
fitcpgsitesinenhancery=(fitcpgsitesinenhancerx/3)/1.6
plot(fitcpgsitesinenhancerx,fitcpgsitesinenhancery,type="l",lty=3,xlim=c(0,50),ylim=c(0,0.5),xaxt="n",yaxt="n",col="white",ylab="",xlab="",lwd=2.5,bty="l")
axis(side=2,las=1,at=c(0,0.25,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,25,50),cex.axis=1)
axis(side=1,las=1,at=c(0,25,50),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,25,50),cex.axis=1)
par(new=T)
polygon(c(seq(1,49,2),seq(49,1,-2)),c(y2ienhancercontain,rev(x2ienhancercontain)), density = NULL, border = F, col = rgb(255/255,105/255,180/255),xaxt="n",yaxt="n")
par(new=T)
polygon(c(seq(1,49,2),seq(49,1,-2)),c(y2ienhancercontain,rev(x2ienhancercontain)), density = NULL, border = F, col = rgb(255/255,105/255,180/255),xaxt="n",yaxt="n")
par(new=T)
fitcpgsitesinenhancerx=seq(1,50,1)
fitcpgsitesinenhancery=exp(-fitcpgsitesinenhancerx/3)/1.6
plot(fitcpgsitesinenhancerx,fitcpgsitesinenhancery,type="l",lty=2,xlim=c(0,50),ylim=c(0,0.5),xaxt="n",yaxt="n",col="black",ylab="",xlab="",lwd=1.3,bty="L")

mtext("Probability",side=2,line=1.5,cex=1)
mtext("Enhancer contained cpg site",side=1,line=1.2,cex=1)
legend(20,0.3,legend=c("Fitting curve","Experiment data"),cex=0.8,x.intersp=0.5,y.intersp=1.5,col=c('black',rgb(255/255,105/255,180/255)),bty="n",lty=c(3,1),lwd=2)


#fig2b
enhancerdisnaive2=vector()
i=6
enhancerdisnaive2=as.numeric(as.character(Mcpgcov2i[,4+2*i])[which(as.character(Mcpgcov2i[,4+2*i])!="NA" & as.numeric(as.character(Mcpgcov2i[,3+2*i]))>hold)])/100
#col=rgb(30/255,144/255,255/255) col=rgb(238/255,44/255,44/255)
x=hist(enhancerdisnaive2,breaks=20)
x=x$counts/sum(x$counts)

enhancermethylationlevel=rbeta(100000,3,70)
fitx=hist(enhancermethylationlevel,xlim=c(0,1),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col="grey",ylab="",main="",xaxt="n",yaxt="n")
fitx=fitx$counts/sum(fitx$counts)

par(mar=c(2.5,3,0.5,0.5))
plot(1,1,col='white',xlim=c(0,1),ylim=c(0,0.7),tck=0.03,las=1,xlab="",ylab="",main="",xaxt="n",yaxt="n",bty='L')
for(i in 1:20)
{ 
  par(new=T)
  polygon(c(0.05*i-0.05,0.05*i-0.05,0.05*i,0.05*i),c(0,x[i],x[i],0), density = NULL, border = T, col = "grey")
}
par(new=T)
#plot(seq(0.025,0.975,0.05),fitx,col='red',type='l',lty=2,lwd=1.5,xlim=c(0,1),ylim=c(0,0.7),tck=0.03,las=1,xlab="",ylab="",main="",xaxt="n",yaxt="n",bty='L')

axis(side=1,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,50,100),cex.axis=1.0)
axis(side=2,las=1,at=c(0,0.35,0.7),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,35,70),cex.axis=1.0)
mtext("Probability",side=2,line=1.5,cex=1)
mtext("Methylation(%)",side=1,line=1.5,cex=1)
#legend(0.5,0.5,legend=c("Fitting curve","Experiment data"),cex=0.8,x.intersp=0.5,y.intersp=1.5,col=c('red','grey'),bty="n",lty=c(2,1),lwd=1.5)

#polygon(c(0.54,0.54,0.63,0.63),c(0.35,0.38,0.38,0.35), density = NULL, border = F, col = "grey")


# section 22 检查每个cell的cpg个数
setwd("/Users/yeyusong/Desktop/在投/Methylation/data/2i")
filename <- list.files()
cellcpgnumber2i=vector()
for(cells in 1:1)
{
  cell2i=read.table(filename[cells],header=F,check.names=FALSE,sep="\t")
  cellcpgnumber2i[cells]=length(cell2i[,1])
}

setwd("/Users/yeyusong/Desktop/在投/Methylation/data/ser")
filename <- list.files()
cellcpgnumberser=vector()
for(cells in 1:20)
{
  cellser=read.table(filename[cells],header=F,check.names=FALSE,sep="\t")
  cellcpgnumberser[cells]=length(cellser[,1])
}
# section22 散点图
cellcpgnumber2ix=vector()
for(i in 1:length(cellcpgnumber2i))
{cellcpgnumber2ix[i]=runif(1,0.5,1.5)}
cellcpgnumberserx=vector()
for(i in 1:length(cellcpgnumberser))
{cellcpgnumberserx[i]=runif(1,2.5,3.5)}
plot(cellcpgnumber2ix,cellcpgnumber2i,type="p",tck=0.03,las=1,xlab="",col="blue",pch=19, ylab="Detected CPG site numbers", main="",xlim=c(0,4),ylim=c(2000000,8000000),xaxt="n",yaxt="n")
par(new=TRUE)
plot(cellcpgnumberserx,cellcpgnumberser,type="p",tck=0.03,las=1,xlab="", col="red",pch=19,ylab="Detected CPG site numbers", main="",xlim=c(0,4),ylim=c(2000000,8000000),xaxt="n",yaxt="n")
axis(side=1,las=1,at=c(1,3),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Naive cells","Prime cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(2000000,5000000,8000000),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(expression(2%*%10^"6"),expression(5%*%10^"6"),expression(8%*%10^"6")),cex.axis=1.0)

# section22 像线图
boxplot(cellcpgnumber2i,xlim=c(0.5,2),ylim=c(2000000,8000000),yaxt="n",col="blue",ylab="Detected cpg site number")
par(new=TRUE)
boxplot(cellcpgnumberser,xlim=c(0,1.5),ylim=c(2000000,8000000),yaxt="n",col="red")
axis(side=1,las=1,at=c(0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Naive cells","Prime cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(2000000,5000000,8000000),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(expression(2%*%10^"6"),expression(5%*%10^"6"),expression(8%*%10^"6")),cex.axis=1.0)



#section23 散点图
enhancercpgtotalnum2i=vector()
for(cells in 1:12)
{
  enhancercpgtotalnum2i[cells]=sum(as.numeric(as.character(Mcpgcov2i[,3+2*cells])))
}

enhancercpgtotalnumser=vector()
for(cells in 1:20)
{
  enhancercpgtotalnumser[cells]=sum(as.numeric(as.character(Mcpgcovser[,3+2*cells])))
}
enhancercpgtotalnum2ix=vector()
enhancercpgtotalnumserx=vector()
for(i in 1:length(enhancercpgtotalnum2i))
{enhancercpgtotalnum2ix[i]=runif(1,0.5,1.5)}
for(i in 1:length(enhancercpgtotalnumser))
{enhancercpgtotalnumserx[i]=runif(1,2.5,3.5)}

plot(enhancercpgtotalnum2ix,enhancercpgtotalnum2i,type="p",tck=0.03,las=1,xlab="",col="blue",pch=19, ylab="Detected CPG site numbers(in enhancer)", main="",xlim=c(0,4),ylim=c(20000,80000),xaxt="n",yaxt="n")
par(new=TRUE)
plot(enhancercpgtotalnumserx,enhancercpgtotalnumser,type="p",tck=0.03,las=1,xlab="", col="red",pch=19,ylab="Detected CPG site numbers(in enhancer)", main="",xlim=c(0,4),ylim=c(20000,80000),xaxt="n",yaxt="n")
axis(side=1,las=1,at=c(1,3),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Naive cells","Prime cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(20000,50000,80000),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(expression(2%*%10^"4"),expression(5%*%10^"4"),expression(8%*%10^"4")),cex.axis=1.0)

# section23 像线图
boxplot(enhancercpgtotalnum2i,xlim=c(0.5,2),ylim=c(20000,80000),yaxt="n",col="blue",ylab="Detected cpg site number(in enhancer)")
par(new=TRUE)
boxplot(enhancercpgtotalnumser,xlim=c(0,1.5),ylim=c(20000,80000),yaxt="n",col="red")
axis(side=1,las=1,at=c(0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Naive cells","Prime cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(20000,50000,80000),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(expression(2%*%10^"4"),expression(5%*%10^"4"),expression(8%*%10^"4")),cex.axis=1.0)




# section 24 


#section 25
#enhancercpgfraction distribution 2i
hold=4
methyleveldistri1=vector()
methyleveldistri2=vector()
methyleveldistri3=vector()
for(cells in c(1,2,3,5,7,12))
{
  x=hist(as.numeric(as.character(Mcpgcov2i[,4+2*cells]))[which(as.numeric(as.character(Mcpgcov2i[,(4+2*cells-1)]))>hold)],breaks=10,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="CPG fraction(%)",ylab="Frequency",main="",col="grey")
  methyleveldistri1=c(methyleveldistri1,x$counts/(sum(x$counts)))
}
for(cells in c(4,9,10))
{
  x=hist(as.numeric(as.character(Mcpgcov2i[,4+2*cells]))[which(as.numeric(as.character(Mcpgcov2i[,(4+2*cells-1)]))>hold)],breaks=10,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="CPG fraction(%)",ylab="Frequency",main="",col="grey")
  methyleveldistri2=c(methyleveldistri2,x$counts/(sum(x$counts)))
}
for(cells in c(6,8,11))
{
  x=hist(as.numeric(as.character(Mcpgcov2i[,4+2*cells]))[which(as.numeric(as.character(Mcpgcov2i[,(4+2*cells-1)]))>hold)],breaks=10,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="CPG fraction(%)",ylab="Frequency",main="",col="grey")
  methyleveldistri3=c(methyleveldistri3,x$counts/(sum(x$counts)))
}

methyleveldistrix=seq(5,95,10)
#section 25 第一类
for(i in 1:6)
{par(new=TRUE)
  plot(methyleveldistrix,methyleveldistri1[((i-1)*10+1):((i-1)*10+10)],type="l",xlim=c(0,100),ylim=c(0,0.5),xaxt="n",yaxt="n",col=rgb(100/255,149/255,237/255),ylab="",xlab="",lwd=1)
}
axis(side=2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","10%","20%","30%","40%","50%"),cex.axis=1.0)
axis(side=1,las=1,at=c(0,25,50,75,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","25%","50%","75%","100%"),cex.axis=1.0)
mtext("Probability",side=2,line=2.7,cex=1)
mtext("Methylation level in enhancer",side=1,line=1.7,cex=1)

#section 25 第2类
for(i in 1:3)
{par(new=TRUE)
  plot(methyleveldistrix,methyleveldistri2[((i-1)*10+1):((i-1)*10+10)],type="l",xlim=c(0,100),ylim=c(0,0.5),xaxt="n",yaxt="n",col="blue",ylab="",xlab="",lwd=1)
}
axis(side=2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","10%","20%","30%","40%","50%"),cex.axis=1.0)
axis(side=1,las=1,at=c(0,25,50,75,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","25%","50%","75%","100%"),cex.axis=1.0)
mtext("Probability",side=2,line=2.7,cex=1)
mtext("Methylation level in enhancer",side=1,line=1.7,cex=1)


#section 25 第3类
for(i in 1:3)
{par(new=TRUE)
  plot(methyleveldistrix,methyleveldistri3[((i-1)*10+1):((i-1)*10+10)],type="l",xlim=c(0,100),ylim=c(0,0.5),xaxt="n",yaxt="n",col=rgb(25/255,25/255,112/255),ylab="",xlab="",lwd=1)
}
axis(side=2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","10%","20%","30%","40%","50%"),cex.axis=1.0)
axis(side=1,las=1,at=c(0,25,50,75,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","25%","50%","75%","100%"),cex.axis=1.0)
mtext("Probability",side=2,line=2.7,cex=1)
mtext("Methylation level in enhancer",side=1,line=1.7,cex=1)



#section 25
#enhancercpgfraction distribution ser
hold=4
methyleveldistri1=vector()
methyleveldistri2=vector()
methyleveldistri3=vector()
methyleveldistri4=vector()
for(cells in c(1,2,7,8,10,11,12,17,18))
{
  x=hist(as.numeric(as.character(Mcpgcovser[,4+2*cells]))[which(as.numeric(as.character(Mcpgcovser[,(4+2*cells-1)]))>hold)],breaks=10,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="CPG fraction(%)",ylab="Frequency",main="",col="grey")
  methyleveldistri1=c(methyleveldistri1,x$counts/(sum(x$counts)))
}
for(cells in c(4,5,9,13,14,15,16,19,20))
{
  x=hist(as.numeric(as.character(Mcpgcovser[,4+2*cells]))[which(as.numeric(as.character(Mcpgcovser[,(4+2*cells-1)]))>hold)],breaks=10,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="CPG fraction(%)",ylab="Frequency",main="",col="grey")
  methyleveldistri2=c(methyleveldistri2,x$counts/(sum(x$counts)))
}
for(cells in c(3))
{
  x=hist(as.numeric(as.character(Mcpgcovser[,4+2*cells]))[which(as.numeric(as.character(Mcpgcovser[,(4+2*cells-1)]))>hold)],breaks=10,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="CPG fraction(%)",ylab="Frequency",main="",col="grey")
  methyleveldistri3=c(methyleveldistri3,x$counts/(sum(x$counts)))
}
for(cells in c(6))
{
  x=hist(as.numeric(as.character(Mcpgcovser[,4+2*cells]))[which(as.numeric(as.character(Mcpgcovser[,(4+2*cells-1)]))>hold)],breaks=10,mgp=c(3, 0.5, 0.1),tck=0.03,las=1,xlab="CPG fraction(%)",ylab="Frequency",main="",col="grey")
  methyleveldistri4=c(methyleveldistri4,x$counts/(sum(x$counts)))
}

methyleveldistrix=seq(5,95,10)
#section 25 第一类
for(i in 1:9)
{par(new=TRUE)
  plot(methyleveldistrix,methyleveldistri1[((i-1)*10+1):((i-1)*10+10)],type="l",xlim=c(0,100),ylim=c(0,0.5),xaxt="n",yaxt="n",col=rgb(255/255,48/255,48/255),ylab="",xlab="",lwd=1)
}
axis(side=2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","10%","20%","30%","40%","50%"),cex.axis=1.0)
axis(side=1,las=1,at=c(0,25,50,75,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","25%","50%","75%","100%"),cex.axis=1.0)
mtext("Probability",side=2,line=2.7,cex=1)
mtext("Methylation level in enhancer",side=1,line=1.7,cex=1)

#section 25 第2类
for(i in 1:9)
{par(new=TRUE)
  plot(methyleveldistrix,methyleveldistri2[((i-1)*10+1):((i-1)*10+10)],type="l",xlim=c(0,100),ylim=c(0,0.5),xaxt="n",yaxt="n",col=rgb(199/255,21/255,133/255),ylab="",xlab="",lwd=1)
}
axis(side=2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","10%","20%","30%","40%","50%"),cex.axis=1.0)
axis(side=1,las=1,at=c(0,25,50,75,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","25%","50%","75%","100%"),cex.axis=1.0)
mtext("Probability",side=2,line=2.7,cex=1)
mtext("Methylation level in enhancer",side=1,line=1.7,cex=1)


#section 25 第3类
for(i in 1:1)
{par(new=TRUE)
  plot(methyleveldistrix,methyleveldistri3[((i-1)*10+1):((i-1)*10+10)],type="l",xlim=c(0,100),ylim=c(0,0.5),xaxt="n",yaxt="n",col=rgb(255/255,99/255,71/255),ylab="",xlab="",lwd=2)
}
axis(side=2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","10%","20%","30%","40%","50%"),cex.axis=1.0)
axis(side=1,las=1,at=c(0,25,50,75,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","25%","50%","75%","100%"),cex.axis=1.0)
mtext("Probability",side=2,line=2.7,cex=1)
mtext("Methylation level in enhancer",side=1,line=1.7,cex=1)

#section 25 第4类
for(i in 1:1)
{par(new=TRUE)
  plot(methyleveldistrix,methyleveldistri4[((i-1)*10+1):((i-1)*10+10)],type="l",xlim=c(0,100),ylim=c(0,0.5),xaxt="n",yaxt="n",col=rgb(160/255,32/255,240/255),ylab="",xlab="",lwd=2)
}
axis(side=2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","10%","20%","30%","40%","50%"),cex.axis=1.0)
axis(side=1,las=1,at=c(0,25,50,75,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","25%","50%","75%","100%"),cex.axis=1.0)
mtext("Probability",side=2,line=2.7,cex=1)
mtext("Methylation level in enhancer",side=1,line=1.7,cex=1)






# section26

# enhancercpgfraction var distribution
#2i
enhancercpgfractionsitenum=vector()
for(i in 1:25036)
{enhancercpgfractionsitenum[i]=0
for(j in 1:12)
{if(Mcpgcov2i[i,(3+2*j)]>0)
  enhancercpgfractionsitenum[i]=enhancercpgfractionsitenum[i]+1
}
}
barenhancercpgfractionsitenum=vector()


enhancerfractionmaxtrix2i=matrix()
enhancerfractionmaxtrix2i=Mcpgcov2i[,6]
for(i in 1:11)
{
  enhancerfractionmaxtrix2i=data.frame(enhancerfractionmaxtrix2i,Mcpgcov2i[,6+2*i])
}
meaneachenhancer2i=vector()
for(i in 1:25056)
{meaneachenhancer2i[i]=mean(as.numeric(as.character((enhancerfractionmaxtrix2i[i,][which(enhancerfractionmaxtrix2i[i,]!="NA")]))))}


# ser
enhancercpgfractionsitenum=vector()
for(i in 1:25036)
{enhancercpgfractionsitenum[i]=0
for(j in 1:20)
{if(Mcpgcovser[i,(3+2*j)]>0)
  enhancercpgfractionsitenum[i]=enhancercpgfractionsitenum[i]+1
}
}

enhancerfractionmaxtrixser=matrix()
enhancerfractionmaxtrixser=Mcpgcovser[,6]
for(i in 1:19)
{
  enhancerfractionmaxtrixser=data.frame(enhancerfractionmaxtrixser,Mcpgcovser[,6+2*i])
}
meaneachenhancerser=vector()
for(i in 1:25056)
{meaneachenhancerser[i]=mean(as.numeric(as.character((enhancerfractionmaxtrixser[i,][which(enhancerfractionmaxtrixser[i,]!="NA")]))))}

x=hist(meaneachenhancer2i,ylim=c(0,2500),tck=0.03,las=1,xlab="",col="blue",ylab="",main="",xlim=c(0,100),xaxt="n",yaxt="n")
axis(side=1,las=1,at=c(0,50,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50%","100%"),cex.axis=1.0)
axis(side=2,las=1,at=c(0,500,1000,1500,2000,2500),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,500,1000,1500,2000,2500),cex.axis=1.0)
mtext("Frequency",side=2,line=2.7,cex=1)
mtext("Methylation level",side=1,line=1.7,cex=1)

statmeaneachenhancer2i=x$counts/sum(x$counts)

par(new=TRUE)
x=hist(meaneachenhancerser,ylim=c(0,2500),tck=0.03,las=1,xlab="",col="red",ylab="",main="",xlim=c(0,100),xaxt="n",yaxt="n")
axis(side=1,las=1,at=c(0,50,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50%","100%"),cex.axis=1.0)
axis(side=2,las=1,at=c(0,500,1000,1500,2000,2500),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,500,1000,1500,2000,2500),cex.axis=1.0)
mtext("Frequency",side=2,line=2.7,cex=1)
mtext("Methylation level",side=1,line=1.7,cex=1)

statmeaneachenhancerser=x$counts/sum(x$counts)
statmeaneachenhancer=as.matrix(statmeaneachenhancer2i)
statmeaneachenhancer=data.frame(statmeaneachenhancer,statmeaneachenhancerser)

par(mar=c(3,4.5,1,1))
barplot(t(statmeaneachenhancer),width=1,space=c(0,0.3),beside=TRUE,col=c(rgb(30/255,144/255,255/255),rgb(238/255,44/255,44/255)),ylim=c(0,0.2),xaxt="n",yaxt="n",border=F)
axis(side=2,las=1,at=c(0,0.05,0.1,0.15,0.2),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","5%","10%","15%","20%"),cex.axis=1.5)
mtext("50%",side=1,line=0.5,cex=1.5)
mtext("0%",side=1,line=0.5,cex=1.5,at=2)
mtext("100%",side=1,line=0.5,cex=1.5,at=45)
legend(7,0.2,legend=c("Naive","Prime"),cex=1.5,x.intersp=0.5,y.intersp=1.5,col=c(rgb(30/255,144/255,255/255),rgb(238/255,44/255,44/255)),bty="n",lty=1,lwd=5)

mtext("Probability",side=2,line=3,cex=1.5)
mtext("Methylation level",side=1,line=1.7,cex=1.5)



#section 26正负图
par(mar=c(2,4.2,1,1))
detlameaneachenhancer=meaneachenhancerser-meaneachenhancer2i
detlameaneachenhancer=(detlameaneachenhancer[which(detlameaneachenhancer!="NaN")])
hist(detlameaneachenhancer)
barplot(sort(detlameaneachenhancer)/100,border=heat.colors(12495),xaxt="n",yaxt="n")
axis(side=2,las=1,at=c(-1,-0.5,0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(-1,-0.5,0,0.5,1),cex.axis=1.5)
mtext("Enhancers",side=1,line=1,cex=1.5)
mtext(expression(Delta),side=2,line=3,cex=1.5,at=-0.8)
mtext("Methylation level",side=2,line=3,cex=1.5)

#section 26正负图的分布
par(mar=c(3,5,1,1))
x=hist(detlameaneachenhancer,tck=0.03,las=1,xlab="",col=heat.colors(20),ylab="",main="",xaxt="n",yaxt="n")
axis(side=2,las=1,at=c(0,500,1000,1500,2000,2500),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,500,1000,1500,2000,2500),cex.axis=1.5)
mtext("Enhancer number",side=2,line=3.5,cex=1.5)
mtext(expression(Delta),side=1,line=1.5,cex=1.5,at=-50)
mtext("Methylation level",side=1,line=1.5,cex=1.5)
mtext("-1",side=1,line=0.2,cex=1.5,at=-100)
mtext("0",side=1,line=0.2,cex=1.5,at=0)
mtext("1",side=1,line=0.2,cex=1.5,at=100)




#section 27
# each cell enhancer var

cellenhancermean2i=vector()
for(i in 1:12)
{
  cellenhancermean2i[i]=mean(as.numeric(as.character(Mcpgcov2i[,4+2*i])[which(as.character(Mcpgcov2i[,4+2*i])!="NA")])/100)
}

cellenhancermeanser=vector()
for(i in 1:20)
{
  cellenhancermeanser[i]=mean(as.numeric(as.character(Mcpgcovser[,4+2*i])[which(as.character(Mcpgcovser[,4+2*i])!="NA")])/100)
}

#散点图
cellenhancermean2ix=vector()
for(i in 1:length(cellenhancermean2i))
{cellenhancermean2ix[i]=rnorm(1,1,0.25)}
cellenhancermeanserx=vector()
for(i in 1:length(cellenhancermeanser))
{cellenhancermeanserx[i]=rnorm(1,3,0.25)}
plot(cellenhancermean2ix,cellenhancermean2i,type="p",tck=0.03,las=1,xlab="",col="blue",pch=19, ylab="Each cell Average Methylation level", main="",xlim=c(0,4),ylim=c(0,1),xaxt="n",yaxt="n")
par(new=TRUE)
plot(cellenhancermeanserx,cellenhancermeanser,type="p",tck=0.03,las=1,xlab="", ylab="",col="red",pch=19, main="",xlim=c(0,4),ylim=c(0,1),xaxt="n",yaxt="n")
axis(side=1,las=1,at=c(1,3),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Naive cells","Prime cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50%","100%"),cex.axis=1.0)



# section22 像线图
hold=0
allenhancer2i=vector()
for(i in 1:12)
{
  eachcell=as.numeric(as.character(Mcpgcov2i[,4+2*i])[which(as.character(Mcpgcov2i[,4+2*i])!="NA" & as.character(Mcpgcov2i[,3+2*i])>hold)])/100
  allenhancer2i=c(allenhancer2i,eachcell)
}

allenhancerser=vector()
for(i in 1:20)
{
  eachcell=as.numeric(as.character(Mcpgcovser[,4+2*i])[which(as.character(Mcpgcovser[,4+2*i])!="NA" & as.character(Mcpgcovser[,3+2*i])>hold)])/100
  allenhancerser=c(allenhancerser,eachcell)
}
boxplot(allenhancer2i[which(allenhancer2i!=0 & allenhancer2i!=1)],xlim=c(0.5,2),ylim=c(0,1),yaxt="n",col="blue",ylab="Methylation level",range=3)
par(new=TRUE)
boxplot(allenhancerser[which(allenhancerser!=0 & allenhancerser!=1)],xlim=c(0,1.5),ylim=c(0,1),col="red",ylab="",yaxt="n",range=3)
axis(side=1,las=1,at=c(0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Naive cells","Prime cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50%","100%"),cex.axis=1.0)
sd(allenhancer2i[which(allenhancer2i!=0 & allenhancer2i!=1)])
sd(allenhancerser[which(allenhancerser!=0 & allenhancerser!=1)])
mean(allenhancer2i[which(allenhancer2i!=0 & allenhancer2i!=1)])
mean(allenhancerser[which(allenhancerser!=0 & allenhancerser!=1)])

#naive jiajihua mean var
par(mar=c(2,4,1,3))
plot(1,0.43,type="p",tck=0.03,cex=0.5,las=1,xlab="",col='black',pch=19, ylab="", main="",xlim=c(0.5,2.5),ylim=c(0,1),xaxt="n",yaxt="n")
par(new=TRUE)
plot(2,0.61,type="p",tck=0.03,cex=0.5,las=1,xlab="",col='black',pch=19, ylab="", main="",xlim=c(0.5,2.5),ylim=c(0,1),xaxt="n",yaxt="n")
sadianx=vector()
sadianx=rnorm(length((allenhancer2i[which(allenhancer2i!=0 & allenhancer2i!=1)])), mean = 1, sd = 0.1)
par(new=T)
plot(sadianx,(allenhancer2i[which(allenhancer2i!=0 & allenhancer2i!=1)]),col=rgb(30/255,144/255,255/255),type='p',cex=0.0025,xlab="",ylab="", main="",xlim=c(0.5,2.5),ylim=c(0,1),xaxt="n",yaxt="n")

sadianx=rnorm(length((allenhancerser[which(allenhancerser!=0 & allenhancerser!=1)])), mean = 2, sd = 0.1)
par(new=T)
plot(sadianx,(allenhancerser[which(allenhancerser!=0 & allenhancerser!=1)]),col=rgb(238/255,44/255,44/255),type='p',cex=0.0025,xlab="",ylab="", main="",xlim=c(0.5,2.5),ylim=c(0,1),xaxt="n",yaxt="n")

pol=0.5
par(new=T)
polygon(c(0.999,1.001,1.001,0.999),c(-0.5,-0.5,1.1,1.1), density = NULL, border = T, col = "black")
par(new=T)
polygon(c(1.999,2.001,2.001,1.999),c(-0.5,-0.5,1.1,1.1), density = NULL, border = T, col = "black")
par(new=T)
polygon(c(0.7,1.3,1.3,0.7),c(0.43-pol*0.2,0.43-pol*0.2,0.43+pol*0.2,0.43+pol*0.2), density = NULL, border = T, col = rgb(30/255,144/255,255/255,150/255))
par(new=T)
polygon(c(1+0.7,1+1.3,1+1.3,1+0.7),c(0.61-pol*0.23,0.61-pol*0.23,0.61+pol*0.23,0.61+pol*0.23), density = NULL, border = T, col = rgb(238/255,44/255,44/255,150/255))
pol=0.01
par(new=T)
polygon(c(0.7,1.3,1.3,0.7),c(0.43-pol*0.43,0.43-pol*0.43,0.43+pol*0.43,0.43+pol*0.43), density = NULL, border = T, col = "black")
par(new=T)
polygon(c(1+0.7,1+1.3,1+1.3,1+0.7),c(0.61-pol*0.54,0.61-pol*0.54,0.61+pol*0.54,0.61+pol*0.54), density = NULL, border = T, col = "black")
axis(side=1,las=1,at=c(1,2),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Naive cells","Prime cells"),cex.axis=1.5)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,0.5,1),cex.axis=1.5)
mtext("Methylation level",side=2,line=2.5,cex=1.5)





par(mar=c(2,5,1,3))
plot(1,1,xlab="",ylab="")
mtext("var mean",side=2,line=2.5,cex=1.5)


boxplot(cellenhancermean2i,xlim=c(0.5,2),ylim=c(0,1),yaxt="n",col="blue",ylab="Methylation level")
par(new=TRUE)
boxplot(cellenhancermeanser,xlim=c(0,1.5),ylim=c(0,1),col="red",ylab="",yaxt="n")
axis(side=1,las=1,at=c(0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Naive cells","Prime cells"),cex.axis=1.0)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50%","100%"),cex.axis=1.0)



#E455565
setwd("/Users/yeyusong/Desktop/在投/Methylation/data")
E45methylationsite=read.csv("E45_methy_fraction.csv",header=T,check.names=FALSE,sep=",")
E55methylationsite=read.csv("E55_methy_fraction.csv",header=T,check.names=FALSE,sep=",")
E65methylationsite=read.csv("E65_methy_fraction.csv",header=T,check.names=FALSE,sep=",")

hold=0
#E45
E45methyleveldistri=vector()
meanE45methylevel=vector()
sitenumE45=vector()
for(cells in 1:91)
{
  sitenumE45=c(sitenumE45,length(which(E45methylationsite[,3+2*cells]>0)))
}
hist(sitenumE45[which(sitenumE45<1000)])

for(cells in 1:91)
{
  x=(as.numeric(as.character(E45methylationsite[,4+2*cells]))[which(as.numeric(as.character(E45methylationsite[,(4+2*cells-1)]))>hold)])
  E45methyleveldistri=c(E45methyleveldistri,x)
  meanE45methylevel=c(meanE45methylevel,mean(x))
}
length(E45methyleveldistri)
hist(meanE45methylevel,breaks=20,xlim=c(0,100))


#E55
E55methyleveldistri=vector()
meanE55methylevel=vector()
sitenumE55=vector()
for(cells in 1:80)
{
  sitenumE55=c(sitenumE55,length(which(E55methylationsite[,3+2*cells]>0)))
}
hist(sitenumE55[which(sitenumE55<10000)])

for(cells in 1:80)
{
  x=(as.numeric(as.character(E55methylationsite[,4+2*cells]))[which(as.numeric(as.character(E55methylationsite[,(4+2*cells-1)]))>hold)])
  E55methyleveldistri=c(E55methyleveldistri,x)
  meanE55methylevel=c(meanE55methylevel,mean(x))
}
hist(E55methyleveldistri[which(E55methyleveldistri>0 & E55methyleveldistri<100)],breaks=20)
length(E55methyleveldistri)
hist(meanE55methylevel,breaks=20,xlim=c(0,100))

i=1
par(mfrow=c(5,3),mar=c(2,2,1,1))
for(cells in i:(i+14))
{
  x=(as.numeric(as.character(E55methylationsite[,4+2*cells]))[which(as.numeric(as.character(E55methylationsite[,(4+2*cells-1)]))>hold)])
  if(length(x[which(x>0 & x<100)])>0)
    hist(x[which(x>0 & x<100)],xlim=c(0,100),breaks=20)
  else
    return
}


#E65
E65methyleveldistri=vector()
meanE65methylevel=vector()
sitenumE65=vector()
for(cells in 1:96)
{
  sitenumE65=c(sitenumE65,length(which(E65methylationsite[,3+2*cells]>0)))
}
hist(sitenumE65[which(sitenumE65<10000)])

for(cells in 1:96)
{
  x=(as.numeric(as.character(E65methylationsite[,4+2*cells]))[which(as.numeric(as.character(E65methylationsite[,(4+2*cells-1)]))>hold)])
  E65methyleveldistri=c(E65methyleveldistri,x)
  meanE65methylevel=c(meanE65methylevel,mean(x))
}
hist(E65methyleveldistri[which(E65methyleveldistri>0 & E65methyleveldistri<100)],breaks=10)
length(E65methyleveldistri)
hist(meanE65methylevel,breaks=20,xlim=c(0,100))

#E455565平均甲基化 散点图
meanE45methylevelx=vector()
for(i in 1:91)
{meanE45methylevelx[i]=rnorm(1,1,0.25)}
meanE55methylevelx=vector()
for(i in 1:80)
{meanE55methylevelx[i]=rnorm(1,3,0.25)}
meanE65methylevelx=vector()
for(i in 1:96)
{meanE65methylevelx[i]=rnorm(1,5,0.25)}

par(mar=c(2,4,1,1))
plot(meanE45methylevelx,meanE45methylevel,type="p",tck=0.03,las=1,xlab="",col="blue",pch=19, ylab="", main="",xlim=c(0,6),ylim=c(0,100),xaxt="n",yaxt="n",cex=0.5)
par(new=TRUE)
plot(meanE55methylevelx,meanE55methylevel,type="p",tck=0.03,las=1,xlab="", ylab="",col="red",pch=19, main="",xlim=c(0,6),ylim=c(0,100),xaxt="n",yaxt="n",cex=0.5)
par(new=TRUE)
plot(meanE65methylevelx,meanE65methylevel,type="p",tck=0.03,las=1,xlab="", ylab="",col="darkred",pch=19, main="",xlim=c(0,6),ylim=c(0,100),xaxt="n",yaxt="n",cex=0.5)

axis(side=1,las=1,at=c(1,3,5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("E4.5","E5.5","E6.5"),cex.axis=1.0)
axis(side=2,las=1,at=c(0,50,100),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50%","100%"),cex.axis=1.0)
mtext("Each cell's Average Methylation level (%)",side=2,line=2.5,cex=1)

#E455565平均甲基化 直方图
par(mfrow=c(3,1),mar=c(4,4,1,1))
hist(meanE45methylevel,tck=0.03,las=1,xlab="",col="blue",pch=19, ylab="", main="",xlim=c(0,100),ylim=c(0,20),breaks=40,cex.axis=1.4)
hist(meanE55methylevel,tck=0.03,las=1,xlab="",col="red",pch=19, ylab="", main="",xlim=c(0,100),ylim=c(0,20),breaks=20,cex.axis=1.4)
mtext("cells number",side=2,line=2.5,cex=1)
hist(meanE65methylevel,tck=0.03,las=1,xlab="",col="darkred",pch=19, ylab="", main="",xlim=c(0,100),ylim=c(0,20),breaks=20,cex.axis=1.4)
mtext("methylation (%)",side=1,line=2.0,cex=1)

#E455565甲基化 分布
par(mfrow=c(5,3),mar=c(2,4,1,1))
for(i in 42:56)
{
  hist(E45methylationsite[,4+2*i][which(E45methylationsite[,3+2*i]>0)],tck=0.03,las=1,xlab="",col="blue",pch=19, ylab="", main="",xlim=c(0,100),breaks=20,cex.axis=1.4)
}
par(mfrow=c(5,3),mar=c(2,4,1,1))
for(i in 42:56)
{
  hist(E55methylationsite[,4+2*i][which(E55methylationsite[,3+2*i]>0)],tck=0.03,las=1,xlab="",col="red",pch=19, ylab="", main="",xlim=c(0,100),breaks=20,cex.axis=1.4)
}
par(mfrow=c(5,3),mar=c(2,4,1,1))
for(i in 42:56)
{
  hist(E65methylationsite[,4+2*i][which(E65methylationsite[,3+2*i]>0)],tck=0.03,las=1,xlab="",col="darkred",pch=19, ylab="", main="",xlim=c(0,100),breaks=20,cex.axis=1.4)
}


# E455565 甲基化分布
par(mfrow=c(3,1),mar=c(4,5.5,1,1))
x=hist(E45methyleveldistri[which(E45methyleveldistri>0 & E45methyleveldistri<100)],tck=0.03,las=1,xlab="",col="blue",pch=19, ylab="", main="",xlim=c(0,100),ylim=c(0,5000),breaks=20,cex.axis=1.4)
E45methydistri=x$counts/sum(x$counts)
x=hist(E55methyleveldistri[which(E55methyleveldistri>0 & E55methyleveldistri<100)],tck=0.03,las=1,xlab="",col="red",pch=19, ylab="", main="",xlim=c(0,100),ylim=c(0,5000),breaks=20,cex.axis=1.4)
E55methydistri=x$counts/sum(x$counts)
mtext("frequency",side=2,line=4,cex=1)
x=hist(E65methyleveldistri[which(E65methyleveldistri>0 & E65methyleveldistri<100)],tck=0.03,las=1,xlab="",col="darkred",pch=19, ylab="", main="",xlim=c(0,100),ylim=c(0,5000),breaks=20,cex.axis=1.4)
E65methydistri=x$counts/sum(x$counts)
mtext("methylation (%)",side=1,line=2.5,cex=1)

E455565methydistri=as.matrix(E45methydistri)
E455565methydistri=data.frame(E455565methydistri,E55methydistri,E65methydistri)


par(mar=c(3,4.5,1,1))
barplot(t(E455565methydistri),space=c(0,0.3),beside=TRUE,col=c(rgb(30/255,144/255,255/255),rgb(238/255,44/255,44/255),rgb(0/255,139/255,69/255)),ylim=c(0,0.25),xaxt="n",yaxt="n",border=F)
axis(side=2,las=1,at=c(0,0.05,0.1,0.15,0.2,0.25),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","5%","10%","15%","20%","25%"),cex.axis=1.5)
mtext("50%",side=1,line=0.5,cex=1.5)
mtext("0%",side=1,line=0.5,cex=1.5,at=2)
mtext("100%",side=1,line=0.5,cex=1.5,at=65)
legend(7,0.2*1.25,legend=c("E4.5","E5.5","E6.5"),cex=1.5,x.intersp=0.5,y.intersp=1.5,col=c(rgb(30/255,144/255,255/255),rgb(238/255,44/255,44/255),rgb(0/255,139/255,69/255)),bty="n",lty=1,lwd=5)
mtext("Probability",side=2,line=3,cex=1.5)
mtext("Methylation level",side=1,line=1.7,cex=1.5)


hist(E45methyleveldistri[which(E45methyleveldistri>0 & E45methyleveldistri<100)],breaks=10)
hist(E45methyleveldistri[which(E45methyleveldistri>0 & E45methyleveldistri<100)],breaks=10)
hist(E45methyleveldistri[which(E45methyleveldistri>0 & E45methyleveldistri<100)],breaks=10)

sd((E45methyleveldistri[which(E45methyleveldistri>0 & E45methyleveldistri<100)])/100)
sd((E55methyleveldistri[which(E55methyleveldistri>0 & E55methyleveldistri<100)])/100)
sd((E65methyleveldistri[which(E65methyleveldistri>0 & E65methyleveldistri<100)])/100)
mean(E45methyleveldistri[which(E45methyleveldistri>0 & E45methyleveldistri<100)])
mean(E55methyleveldistri[which(E55methyleveldistri>0 & E55methyleveldistri<100)])
mean(E65methyleveldistri[which(E65methyleveldistri>0 & E65methyleveldistri<100)])

# e455565 jiajihua mean var
par(mar=c(2,4,1,3))
plot(1,0.34,type="p",tck=0.03,cex=0.5,las=1,xlab="",col=rgb(30/255,144/255,255/255),pch=19, ylab="", main="",xlim=c(0.5,3.5),ylim=c(0,1),xaxt="n",yaxt="n")
par(new=TRUE)
plot(2,0.59,type="p",tck=0.03,cex=0.5,las=1,xlab="",col=rgb(238/255,44/255,44/255),pch=19, ylab="", main="",xlim=c(0.5,3.5),ylim=c(0,1),xaxt="n",yaxt="n")
#sadian
sadianx=vector()
sadianx=rnorm(length((E45methyleveldistri[which(E45methyleveldistri>0 & E45methyleveldistri<100)])), mean = 1, sd = 0.1)
par(new=T)
plot(sadianx,(E45methyleveldistri[which(E45methyleveldistri>0 & E45methyleveldistri<100)])/100,col=rgb(30/255,144/255,255/255),type='p',cex=0.0001,xlab="",ylab="", main="",xlim=c(0.5,3.5),ylim=c(0,1),xaxt="n",yaxt="n")

sadianx=rnorm(length((E55methyleveldistri[which(E55methyleveldistri>0 & E55methyleveldistri<100)])), mean = 2, sd = 0.1)
par(new=T)
plot(sadianx,(E55methyleveldistri[which(E55methyleveldistri>0 & E55methyleveldistri<100)])/100,col = rgb(238/255,44/255,44/255),type='p',cex=0.0001,xlab="",ylab="", main="",xlim=c(0.5,3.5),ylim=c(0,1),xaxt="n",yaxt="n")

sadianx=rnorm(length((E65methyleveldistri[which(E65methyleveldistri>0 & E65methyleveldistri<100)])), mean = 3, sd = 0.1)
par(new=T)
plot(sadianx,(E65methyleveldistri[which(E65methyleveldistri>0 & E65methyleveldistri<100)])/100,col = rgb(0/255,139/255,69/255),type='p',cex=0.0001,xlab="",ylab="", main="",xlim=c(0.5,3.5),ylim=c(0,1),xaxt="n",yaxt="n")

pol=0.25
par(new=T)
polygon(c(0.999,1.001,1.001,0.999),c(-0.5,-0.5,1.1,1.1), density = NULL, border = T, col = "black")
par(new=T)
polygon(c(1.999,2.001,2.001,1.999),c(-0.5,-0.5,1.1,1.1), density = NULL, border = T, col = "black")
par(new=T)
polygon(c(2.999,3.001,3.001,2.999),c(-0.5,-0.5,1.1,1.1), density = NULL, border = T, col = "black")
par(new=T)
polygon(c(0.7,1.3,1.3,0.7),c(0.34-pol*0.2,0.34-pol*0.2,0.34+pol*0.2,0.34+pol*0.2), density = NULL, border = T, col = rgb(30/255,144/255,255/255,150/255))
par(new=T)
polygon(c(1+0.7,1+1.3,1+1.3,1+0.7),c(0.59-pol*0.23,0.59-pol*0.23,0.59+pol*0.23,0.59+pol*0.23), density = NULL, border = T, col = rgb(238/255,44/255,44/255,150/255))
par(new=T)
polygon(c(2+0.7,2+1.3,2+1.3,2+0.7),c(0.6-pol*0.23,0.60-pol*0.23,0.60+pol*0.23,0.60+pol*0.23), density = NULL, border = T, col = rgb(0/255,139/255,69/255,150/255))
pol=0.01
par(new=T)
polygon(c(0.7,1.3,1.3,0.7),c(0.34-pol*0.43,0.34-pol*0.43,0.34+pol*0.43,0.34+pol*0.43), density = NULL, border = T, col = "black")
par(new=T)
polygon(c(1+0.7,1+1.3,1+1.3,1+0.7),c(0.59-pol*0.54,0.59-pol*0.54,0.59+pol*0.54,0.59+pol*0.54), density = NULL, border = T, col = "black")
par(new=T)
polygon(c(2+0.7,2+1.3,2+1.3,2+0.7),c(0.60-pol*0.56,0.60-pol*0.56,0.60+pol*0.56,0.60+pol*0.56), density = NULL, border = T, col = "black")
axis(side=1,las=1,at=c(1,2,3),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("E4.5","E5.5","E6.5"),cex.axis=1.5)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,0.5,1),cex.axis=1.5)
mtext("Methylation level",side=2,line=2.5,cex=1.5)




# heterogeneity index
interval=5
hold=3
#E45

cellsdis=matrix()
for(cells in 1:91)
{
  x=(as.numeric(as.character(E45methylationsite[,4+2*cells]))[which(as.numeric(as.character(E45methylationsite[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
  x=hist(x,breaks=c(0,seq(5,95,interval),(100)))
  x=x$counts/(sum(x$counts))
  cellsdis=data.frame(cellsdis,x)
  }
  else
    return
}
indexE45=vector()
for(cells in 2:87)
{
  x=cellsdis[,cells]
for(cellsvs in 2:87)
{
  y=cellsdis[,cellsvs]
  indexE45=c(indexE45,sum((y)*log10((y+0.001)/(x+0.001))))
}
}
mean(indexE45)

#julei
cellsdis=cellsdis[,-1]
#julei
kc <- kmeans(t(cellsdis),3)
clust=3
for(i in 1:length(which(kc$cluster==clust)))
{
  iclus=which(kc$cluster==clust)[i]
  plot(seq(1,20,1),cellsdis[,iclus],type='l',ylim=c(0,1))
  par(new=T)
}
#E55

cellsdisE55=matrix()
for(cells in 1:80)
{
  x=(as.numeric(as.character(E55methylationsite[,4+2*cells]))[which(as.numeric(as.character(E55methylationsite[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=c(0,seq(5,95,interval),(100)))
    x=x$counts/(sum(x$counts))
    cellsdisE55=data.frame(cellsdisE55,x)
  }
  else
    return
}
indexE55=vector()
for(cells in 2:78)
{
  x=cellsdisE55[,cells]
  for(cellsvs in 2:78)
  {
    y=cellsdisE55[,cellsvs]
    indexE55=c(indexE55,sum((y)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexE55)
hist(indexE55)

#julei
cellsdisE55=cellsdisE55[,-1]
#julei
kc <- kmeans(t(cellsdisE55),3)
dev.off()
clust=1
for(i in 1:length(which(kc$cluster==clust)))
{
  iclus=which(kc$cluster==clust)[i]
  plot(seq(1,20,1),cellsdisE55[,iclus],type='l',ylim=c(0,1))
  par(new=T)
}

#E65

cellsdisE65=matrix()
for(cells in 1:96)
{
  x=(as.numeric(as.character(E65methylationsite[,4+2*cells]))[which(as.numeric(as.character(E65methylationsite[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=c(0,seq(5,95,interval),(100)))
    x=x$counts/(sum(x$counts))
    cellsdisE65=data.frame(cellsdisE65,x)
  }
  else
    return
}
indexE65=vector()
for(cells in 2:96)
{
  x=cellsdisE65[,cells]
  for(cellsvs in 2:96)
  {
    y=cellsdisE65[,cellsvs]
    indexE65=c(indexE65,sum((y)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexE65)
hist(indexE65)

#julei
cellsdisE65=cellsdisE65[,-1]
#julei
kc <- kmeans(t(cellsdisE65),3)
dev.off()
clust=
for(i in 1:length(which(kc$cluster==clust)))
{
  iclus=which(kc$cluster==clust)[i]
  plot(seq(1,20,1),cellsdisE65[,iclus],type='l',ylim=c(0,1))
  par(new=T)
}

# naive 和 prime 的 index
hold=3
#naive
cellsdisnaive=matrix()
for(cells in 1:12)
{
  x=(as.numeric(as.character(Mcpgcov2i[,4+2*cells]))[which(as.numeric(as.character(Mcpgcov2i[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=seq(0,100,5))
    x=x$counts/(sum(x$counts))
    cellsdisnaive=data.frame(cellsdisnaive,x)
  }
  else
    return
}
indexnaive=vector()
for(cells in 2:13)
{
  x=cellsdisnaive[,cells]
  for(cellsvs in 2:13)
  {
    y=cellsdisnaive[,cellsvs]
    indexnaive=c(indexnaive,sum((y+0.001)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexnaive)
hist(indexnaive)

#prime
cellsdisprime=matrix()
for(cells in 1:20)
{
  x=(as.numeric(as.character(Mcpgcovser[,4+2*cells]))[which(as.numeric(as.character(Mcpgcovser[,(4+2*cells-1)]))>hold)])
  if(length(x)>0)
  {
    x=hist(x,breaks=seq(0,100,5))
    x=x$counts/(sum(x$counts))
    cellsdisprime=data.frame(cellsdisprime,x)
  }
  else
    return
}
indexprime=vector()
for(cells in 2:21)
{
  x=cellsdisprime[,cells]
  for(cellsvs in 2:21)
  {
    y=cellsdisprime[,cellsvs]
    indexprime=c(indexprime,sum((y+0.001)*log10((y+0.001)/(x+0.001))))
  }
}
mean(indexprime)
hist(indexprime)

#作图
indexinterval1=c(0.166,0.207,0.123)
indexinterval5=c(0.156,0.195,0.106)
indexinterval10=c(0.146,0.180,0.087)
indexintervalx=c(1,2,3)
par(mar=c(2,5,1,0.5))
plot(indexintervalx,indexinterval1,type='l',col=rgb(102/255,205/255,0/255),tck=0.03,las=1,cex=0.3,xlim=c(0.5,3.5),ylim=c(0.05,0.25),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexintervalx,indexinterval1,type='p',pch=19,col=rgb(102/255,205/255,0/255),bg=rgb(102/255,205/255,0/255),tck=0.03,las=1,cex=1,xlim=c(0.5,3.5),ylim=c(0.05,0.25),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
par(new=T)
plot(indexintervalx,indexinterval5,type='l',col=rgb(238/255,44/255,44/255),tck=0.03,las=1,cex=1,xlim=c(0.5,3.5),ylim=c(0.05,0.25),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexintervalx,indexinterval5,type='p',pch=22,col=rgb(238/255,44/255,44/255),bg=rgb(238/255,44/255,44/255),tck=0.03,las=1,cex=1,xlim=c(0.5,3.5),ylim=c(0.05,0.25),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
par(new=T)
plot(indexintervalx,indexinterval10,type='l',col=rgb(30/255,144/255,255/255),tck=0.03,las=1,cex=1,xlim=c(0.5,3.5),ylim=c(0.05,0.25),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexintervalx,indexinterval10,type='p',pch=24,col=rgb(30/255,144/255,255/255),bg=rgb(30/255,144/255,255/255),tck=0.03,las=1,cex=1,xlim=c(0.5,3.5),ylim=c(0.05,0.25),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")

axis(side=1,las=1,at=c(1,2,3),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("E4.5","E5.5","E6.5"),cex.axis=1.5)
axis(side=2,las=1,at=c(0.05,0.15,0.25),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0.05","0.15","0.25"),cex.axis=1.5)
mtext("HI",side=2,line=3.0,cex=1.5)
legend(2.2,0.25,legend=c("interval 1%","interval 5%","interval 10%"),cex=1.3,x.intersp=0.5,y.intersp=1.5,col=c(rgb(102/255,205/255,0/255),rgb(238/255,44/255,44/255),rgb(30/255,144/255,255/255)),bty="n",lty=1,lwd=2)


#index simulation作图
indexsim=c(0.02,0.15,0.015)
indexintervalx=c(1,2,3)
par(mar=c(2,5,1,0.5))
plot(indexintervalx,indexsim,type='l',col='black',tck=0.03,las=1,cex=0.3,xlim=c(0.5,3.5),ylim=c(0,0.2),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexintervalx,indexsim,type='p',pch=25,col='black',bg='black',tck=0.03,las=1,cex=1,xlim=c(0.5,3.5),ylim=c(0.0,0.2),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
axis(side=1,las=1,at=c(1,2,3),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("Low","Bistable","High"),cex.axis=1.5)
axis(side=2,las=1,at=c(0,0.1,0.2),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.1","0.2"),cex.axis=1.5)
mtext("HI",side=2,line=3.0,cex=1.5)
legend(2.2,0.2,legend=c("Sim"),cex=1.3,x.intersp=0.5,y.intersp=1.5,col='black',bty="n",lty=1,lwd=2)


#E455565之比较
x=hist(indexE45,breaks=seq(0,4,0.025),xlim=c(0,4),ylim=c(0,3000),col='blue')
par(new=T)
y=hist(indexE55,breaks=seq(0,4,0.025),xlim=c(0,4),ylim=c(0,3000),col='red')
par(new=T)
z=hist(indexE65,breaks=seq(0,4,0.025),xlim=c(0,4),ylim=c(0,3000),col='darkred')
x=x$counts/(sum(x$counts))
y=y$counts/(sum(y$counts))
z=z$counts/(sum(z$counts))
xyz=data.frame(x,y,z)

par(mar=c(3,5,1,1))
barplot(t(xyz),beside=TRUE,space=c(0,0.3),col=c(rgb(30/255,144/255,255/255),rgb(238/255,44/255,44/255),rgb(0/255,139/255,69/255)),ylim=c(0,0.5),xaxt="n",yaxt="n",xlim=c(0,150),border=T)
axis(side=2,las=1,at=c(0,0.1,0.2,0.3,0.4,0.5),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","10%","20%","30%","40%","50%"),cex.axis=1.5)
mtext("Probability",side=2,line=3,cex=1.5)
mtext("Heterogenity index",side=1,line=1.7,cex=1.5)
mtext("0.5",side=1,line=0.5,cex=1.5)
mtext("0",side=1,line=0.5,cex=1.5,at=0)
mtext("1",side=1,line=0.5,cex=1.5,at=155)
legend(100,0.5,legend=c("E4.5","E5.5","E6.5"),cex=1.5,x.intersp=0.5,y.intersp=1.5,col=c(rgb(30/255,144/255,255/255),rgb(238/255,44/255,44/255),rgb(0/255,139/255,69/255)),bty="n",lty=1,lwd=1.5)


deltaxyz=data.frame(y-x,z-y)

par(mar=c(3,5,1,1))
barplot(t(deltaxyz),beside=TRUE,space=c(0,0.3),col=c("green","orange"),ylim=c(-0.3,0.2),xaxt="n",xlim=c(0,120),yaxt="n",border=T)
axis(side=2,las=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.2),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("-30%","-20%","-10%","0%","10%","20%"),cex.axis=1.5)
mtext("Distribution diversity",side=2,line=3.5,cex=1.5)
mtext("Heterogeneity index",side=1,line=1.8,cex=1.5)
mtext("0.5",side=1,line=0.5,cex=1.5)
mtext("0",side=1,line=0.5,cex=1.5,at=0)
mtext("1",side=1,line=0.5,cex=1.5,at=120)
legend(60,0.2,legend=c("E5.5 - E4.5","E6.5 - E5.5"),cex=1.5,x.intersp=0.5,y.intersp=1.5,col=c('green','orange'),bty="n",lty=1,lwd=1.5)
