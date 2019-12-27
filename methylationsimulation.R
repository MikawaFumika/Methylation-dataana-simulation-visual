#初始化 8,25 rgb(238/255,44/255,44/255)  he  6,9  rgb(30/255,144/255,255/255)

4,12
9,3

enhancermethylationlevel=vector()
for(i in 1:1000)
{enhancermethylationlevel[i]=runif(1,0,0.9)
}
hist(enhancermethylationlevel)

enhancermethylationlevel=rbeta(1000,8,25)
x=hist(enhancermethylationlevel,xlim=c(0,1),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col="grey",ylab="",main="",xaxt="n",yaxt="n",probability = T)
x=hist(enhancermethylationlevel,xlim=c(0,1),ylim=c(0,0.2/max(x$counts/(sum(x$counts)))*(max(x$density))),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col="grey",ylab="",main="",xaxt="n",yaxt="n",probability = T)
initialfitnaive=x$counts/(sum(x$counts))
x$counts/(sum(x$counts))
axis(side=1,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,0.5,1),cex.axis=1.0)
axis(side=2,las=1,at=c(0,max(x$density)),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,1),cex.axis=1.0)
mtext("Density",side=2,line=1.5,cex=1)
mtext("Methylation level",side=1,line=1.7,cex=1)

#prime fit
par(mar=c(3,4,1,1))
plot(10,0.1,type="p",tck=0.03,cex=0.5,las=1,xlab="",col='white',pch=19, ylab="", main="",xlim=c(0,20),ylim=c(0,0.2),xaxt="n",yaxt="n",bty="l")
par(new=T)
polygon(c(0,as.vector(rbind(seq(1,19,1),seq(1,19,1))),20,20,rev(as.vector(rbind(seq(1,19,1),seq(1,19,1)))),0),c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,rev(as.vector(rbind(statmeaneachenhancer[,2],statmeaneachenhancer[,2])))), density = NULL, border = T, col = rgb(238/255,44/255,44/255))
par(new=T)
plot(seq(1,20,1),statmeaneachenhancer[,2],type="h",tck=0.03,cex=0.5,las=1,xlab="",col='black', ylab="", main="",xlim=c(0,20),ylim=c(0,0.2),xaxt="n",yaxt="n",bty="l")

par(new=T)
plot(seq(0.5,19.5,1),initialfitnaive,type="l",lty=2,lwd=4,tck=0.03,cex=0.5,las=1,xlab="",col='black', ylab="", main="",xlim=c(0,20),ylim=c(0,0.2),xaxt="n",yaxt="n",bty="l")
axis(side=1,las=1,at=c(0,10,20),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","50%","100%"),cex.axis=1.5)
axis(side=2,las=1,at=c(0,0.1,0.2),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.1","0.2"),cex.axis=1.5)
mtext("Methylation level",side=1,line=1.8,cex=1.5)
mtext("Probability",side=2,line=2.4,cex=1.5)
legend(1,0.2,,legend=c("Prime fitting"),cex=1.5,x.intersp=0.5,y.intersp=1.5,col='black',bty="n",lty=3,lwd=4)


#naive fit
par(mar=c(3,4,1,1))
plot(10,0.1,type="p",tck=0.03,cex=0.5,las=1,xlab="",col='white',pch=19, ylab="", main="",xlim=c(0,20),ylim=c(0,0.2),xaxt="n",yaxt="n",bty="l")
par(new=T)
polygon(c(0,as.vector(rbind(seq(1,19,1),seq(1,19,1))),20,20,rev(as.vector(rbind(seq(1,19,1),seq(1,19,1)))),0),c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,rev(as.vector(rbind(statmeaneachenhancer[,1],statmeaneachenhancer[,1])))), density = NULL, border = T, col = rgb(30/255,144/255,255/255))
par(new=T)
plot(seq(1,20,1),statmeaneachenhancer[,1],type="h",tck=0.03,cex=0.5,las=1,xlab="",col='black', ylab="", main="",xlim=c(0,20),ylim=c(0,0.2),xaxt="n",yaxt="n",bty="l")

par(new=T)
plot(seq(0.5,19.5,1),initialfitnaive,type="l",lty=2,lwd=4,tck=0.03,cex=0.5,las=1,xlab="",col='black', ylab="", main="",xlim=c(0,20),ylim=c(0,0.2),xaxt="n",yaxt="n",bty="l")
axis(side=1,las=1,at=c(0,10,20),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0%","50%","100%"),cex.axis=1.5)
axis(side=2,las=1,at=c(0,0.1,0.2),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.1","0.2"),cex.axis=1.5)
mtext("Methylation level",side=1,line=1.8,cex=1.5)
mtext("Probability",side=2,line=2.4,cex=1.5)
legend(10,0.2,,legend=c("Naive fitting"),cex=1.5,x.intersp=0.5,y.intersp=1.5,col='black',bty="n",lty=3,lwd=4)


#canshu
u=0.009
k=3
v=0.2

u=0.49
k=3
v=1

u=0.005
k=3
v=0.06

u=0.004
k=5
v=0.005

u=0.0001
k=3
v=0.001

u=0.01
k=3
v=0.06

#canshu redefine
k=3

zeta=0.01
zeta=0.12
zeta=0.42
u=zeta
v=0.2+(2.1*zeta)^4

u=0.05
u1=0.9
v=0.1
v=0.01
v=0.04

m=vector()
for(i in 1:1000)
{
  mi=sample(1:50,1)
  while(runif(1,0,1)>(exp(-mi/3)/1.6))
  {mi=sample(1:50,1)}
  m[i]=mi
}
hist(m)
m=sample(c(4,5,6),1000,replace=T)

hold=5

a=0.5

#阈值
correlativenumber=5
correlativesites=vector()
Mmethylationlevel=vector()

#迭代
enhancermethylationlevel=rbeta(1000,3,70)
hist(enhancermethylationlevel)

cellsdis=matrix()
par(mfrow=c(5,3),mar=c(1,1,1,1))
for(step in 1:30)
{
  for(i in 1:1000)
  {
    correlativesites=seq(i-correlativenumber,i+correlativenumber,1)
    correlativesites=correlativesites[correlativesites>0 & correlativesites<1001]
    Mmethylationlevel[i]=mean(enhancermethylationlevel[correlativesites])
  }
# Mmethylationlevel=mean(enhancermethylationlevel)


for(i in 1:1000)
{
Ebetaplus=u+u1*(enhancermethylationlevel[i])^k/((enhancermethylationlevel[i])^k+v)+(Mmethylationlevel[i]-enhancermethylationlevel[i])*a
if(Ebetaplus>1)
  Ebetaplus=1
if(Ebetaplus<0)
  Ebetaplus=0
c=m[i]*(Ebetaplus)
d=m[i]*(1-Ebetaplus)
enhancermethylationlevel[i]=rbeta(1,c,d)
}
  x=hist(enhancermethylationlevel[which(m>hold)],xlim=c(0,1),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),ylab="",main="",xaxt="n",yaxt="n",probability = T,bty="l")
  mtext(c("0","0.5","1"),at=c(0,0.5,1),side=1,line=-0.3,cex=0.5)
  
  x=x$counts/(sum(x$counts))
  cellsdis=data.frame(cellsdis,x)
}

index=vector()
for(cells in 2:16)
{
  x=cellsdis[,cells]
  for(cellsvs in 2:16)
  {
    y=cellsdis[,cellsvs]
    index=c(index,sum((y)*log10((y+0.001)/(x+0.001))))
  }
}
mean(index)
#0.01
#0.15
#0.016
hist(index)

y=cellsdis[,2]
for(cells in 3:16)
{y=y+cellsdis[,cells]}
y=y/15
for(cells in 2:16)
{
  x=cellsdis[,cells]
  index=c(index,sum((y)*log10((y+0.001)/(x+0.001))))
}
mean(index)

#fig3a
#小提琴小提琴
#初始化  80,144,200  144,220,200    220,85,85
yinzi=0.8
z=2
par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5))
plot(10,0.1,type="p",tck=0.03,cex=0.5,las=1,xlab="",col='white',pch=19, ylab="", main="",xlim=c(0,16),ylim=c(0,1),xaxt="n",yaxt="n",bty="l")
for(i in 2:16)
{
par(new=T)
polygon(c(i-1,(i-1-yinzi*rev(cellsdis[,i+z])),i-1,(i-1+yinzi*cellsdis[,i+z])),c(1,seq(0.975,0.025,-0.05),0,seq(0.025,0.975,0.05)),density = NULL, border = T, col = rgb((220+10-i*5)/255,(160)/255,(85-10+2*i)/255))
}
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,50,100),cex.axis=1)
axis(side=1,las=1,at=c(seq(1,15,1)),mgp=c(0,0.5,0),tck=0,las=1,labels=seq(0,14,1),cex.axis=1)

mtext("Methylation(%)",side=2,line=1.5,cex=1)
mtext("Cycles",side=1,line=1.5,cex=1)

cellsdis=cellsdis[,-6]

#fig3d
#1 7 15
j=2

par(mar=c(2.5,3.2,1,1))
plot(1,1,type='l',xlim=c(0.05,1.05),ylim=c(0,0.3),col='white',tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",bty='L',xaxt='n',yaxt='n')
for(i in 1:20)
{
  par(new=T)
  polygon(c(i/20-0.02,i/20-0.02,i/20+0.02,i/20+0.02),c(0,cellsdis[i,j],cellsdis[i,j],0), density = NULL, border = F, col=rgb(238/255,44/255,44/255),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="")
}
axis(side=1,las=1,at=c(0.05,0.55,1.05),mgp=c(0,0.5,0),tck=0.0,las=1,labels=c(0,50,100),cex.axis=1.0)
axis(side=2,las=1,at=c(0,0.15,0.3),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,15,30),cex.axis=1.0)

mtext("Probability(%)",side=2,line=1.7,cex=1.0)

mtext("Methylation(%)",side=1,line=1.2,cex=1.0,at=0.55)


#fig3e 
meanmethylationsimulation=vector()
for(i in 2:16)
meanmethylationsimulation[i-1]=sum(cellsdis[,i]*seq(0.025,0.975,0.05))
varmethylationsimulation=vector()
for(i in 2:16)
varmethylationsimulation[i-1]=sum(cellsdis[,i]*(seq(0.025,0.975,0.05)-meanmethylationsimulation[i-1])^2)

methylationsimulationx=seq(1,15,1)
par(mar=c(2.5,2.5,0.5,0.5))
plot(methylationsimulationx,meanmethylationsimulation,type='l',xlim=c(0.5,15.5),ylim=c(0.2,0.8),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",xaxt='n',yaxt='n',bty="L")
par(new=T)
polygon(c(methylationsimulationx,rev(methylationsimulationx)),c(meanmethylationsimulation+0.5*varmethylationsimulation,rev(meanmethylationsimulation-0.5*varmethylationsimulation)), density = NULL, border = F, col = rgb(240/255,240/255,240/255),tck=0.03,mgp=c(0,0.5,0),las=1,xlab="Time(h)", ylab="Methylation(%)", main="",xaxt='n',yaxt='n')
par(new=T)
plot(methylationsimulationx,meanmethylationsimulation,xlim=c(0.5,15.5),ylim=c(0.2,0.8),cex=0.8,pch=16,tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",xaxt='n',yaxt='n',bty="l")
par(new=T)
plot(methylationsimulationx,meanmethylationsimulation,type='l',xlim=c(0.5,15.5),ylim=c(0.2,0.8),col='grey',tck=0.03,mgp=c(0,0.5,0),las=1,xlab="", ylab="", main="",xaxt='n',yaxt='n',bty="l")
axis(side=1,las=1,at=seq(1,15,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(0,14,1),cex.axis=0.9)
axis(side=2,las=1,at=c(0.2,0.5,0.8),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(20,50,80),cex.axis=1.0)
mtext("Methylation(%)",side=2,line=1.5,cex=1)
mtext("Cycles",side=1,line=1.5,cex=1,at=8)


#fig3e final

meanmethylationsimulation=vector()
for(i in 1:15)
  meanmethylationsimulation[i]=0
meanmethylationsimulation=data.frame(meanmethylationsimulation)


for(i in 1:10)
{
#yi xia以下can shu以下参以下参数调整
v=0.01*i
a=0.5
u=0.05
u1=0.9
k=3
hold=5
#阈值
correlativenumber=10
correlativesites=vector()
Mmethylationlevel=vector()

#迭代
enhancermethylationlevel=rbeta(1000,3,70)
hist(enhancermethylationlevel)

cellsdis=matrix()
par(mfrow=c(5,3),mar=c(1,1,1,1))
for(step in 1:30)
{
  for(i in 1:1000)
  {
    correlativesites=seq(i-correlativenumber,i+correlativenumber,1)
    correlativesites=correlativesites[correlativesites>0 & correlativesites<1001]
    Mmethylationlevel[i]=mean(enhancermethylationlevel[correlativesites])
  }
  # Mmethylationlevel=mean(enhancermethylationlevel)
  
  
  for(i in 1:1000)
  {
    Ebetaplus=u+u1*(enhancermethylationlevel[i])^k/((enhancermethylationlevel[i])^k+v)+(Mmethylationlevel[i]-enhancermethylationlevel[i])*a
    if(Ebetaplus>1)
      Ebetaplus=1
    if(Ebetaplus<0)
      Ebetaplus=0
    c=m[i]*(Ebetaplus)
    d=m[i]*(1-Ebetaplus)
    enhancermethylationlevel[i]=rbeta(1,c,d)
  }
  x=hist(enhancermethylationlevel[which(m>hold)],xlim=c(0,1),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),ylab="",main="",xaxt="n",yaxt="n",probability = T,bty="l")
  mtext(c("0","0.5","1"),at=c(0,0.5,1),side=1,line=-0.3,cex=0.5)
  
  x=x$counts/(sum(x$counts))
  cellsdis=data.frame(cellsdis,x)
}

meanvector=vector()
for(i in 3:17)
  meanvector[i-2]=sum(cellsdis[,i]*seq(0.025,0.975,0.05))

meanmethylationsimulation=data.frame(meanmethylationsimulation,meanvector)
}



#col=rgb(30/255,144/255,255/255) col=rgb(238/255,44/255,44/255)
par(mar=c(1.5,3,1,1))
for(j in 1:5)
{
vmeanx=vector()
for(i in 1:15)
{vmeanx[i]=rnorm(1,j,0.1)}
plot(vmeanx,meanmethylationsimulation[,j+1],type="p",tck=0.03,las=1,xlab="",col=rgb((144+j*7)/255,(200-j*11)/255,(20+j*6)/255),pch=19, ylab="", main="",xlim=c(0.5,10.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)
par(new=T)
}
for(j in 6:10)
{
  vmeanx=vector()
  for(i in 1:15)
  {vmeanx[i]=rnorm(1,j,0.1)}
  plot(vmeanx,meanmethylationsimulation[,j+1],type="p",tck=0.03,las=1,xlab="",col=rgb((220-j*14)/255,(85+j*6)/255,(85+j*12)/255),pch=19, ylab="", main="",xlim=c(0.5,10.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)
  par(new=T)
}
axis(side=1,las=1,at=seq(1,10,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(0.01,0.1,0.01),cex.axis=0.9)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50","100"),cex.axis=1.0)
mtext("Methylation(%)",side=2,line=2,cex=1.0)
mtext("V",side=1,line=1.5,cex=1.0)


#col=rgb(30/255,144/255,255/255) col=rgb(238/255,44/255,44/255)
par(mar=c(2.5,3,1,1))
for(j in 1:10)
{
  vmeanx=vector()
  for(i in 1:15)
  {vmeanx[i]=rnorm(1,j,0.1)}
  plot(vmeanx,meanmethylationsimulation[,j+1],type="p",tck=0.03,las=1,xlab="",col=rgb((0)/255,(0)/255,(0)/255),pch=19, ylab="", main="",xlim=c(0.5,10.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)
  par(new=T)
}
axis(side=1,las=1,at=seq(1,10,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(0.01,0.1,0.01),cex.axis=0.8)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50","100"),cex.axis=1.0)
mtext("Methylation(%)",side=2,line=2,cex=1.0)
mtext(expression(V),side=1,line=1.5,cex=1.0)

#fig3f


v=0.1
v=0.01
v=0.04

a=0.5

correlativenumber=10

indexcycle=vector()
for(i in 1:15)
  indexcycle[i]=1


for(j in seq(1,15,1))
{

  cellsdis=matrix()
  
  for(cycleloop in 1:10)
  {
  enhancermethylationlevel=rbeta(1000,3,70)
  par(mfrow=c(5,3),mar=c(1,1,1,1))
  for(step in 1:(j+2))
  {
    for(i in 1:1000)
    {
      correlativesites=seq(i-correlativenumber,i+correlativenumber,1)
      correlativesites=correlativesites[correlativesites>0 & correlativesites<1001]
      Mmethylationlevel[i]=mean(enhancermethylationlevel[correlativesites])
    }
    # Mmethylationlevel=mean(enhancermethylationlevel)
    
    
    for(i in 1:1000)
    {
      Ebetaplus=u+u1*(enhancermethylationlevel[i])^k/((enhancermethylationlevel[i])^k+v)+(Mmethylationlevel[i]-enhancermethylationlevel[i])*a
      if(Ebetaplus>1)
        Ebetaplus=1
      if(Ebetaplus<0)
        Ebetaplus=0
      c=m[i]*(Ebetaplus)
      d=m[i]*(1-Ebetaplus)
      enhancermethylationlevel[i]=rbeta(1,c,d)
    }
    x=hist(enhancermethylationlevel[which(m>hold)],xlim=c(0,1),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),ylab="",main="",xaxt="n",yaxt="n",probability = T,bty="l")
    mtext(c("0","0.5","1"),at=c(0,0.5,1),side=1,line=-0.3,cex=0.5)
    
    x=x$counts/(sum(x$counts))
    cellsdis=data.frame(cellsdis,x)
  }
  }
  
  index=vector()
  for(cells in seq((j+2),(j+2)*10,j+2))
  {
    x=cellsdis[,cells]
    for(cellsvs in seq((j+2),(j+2)*10,j+2))
    {
      y=cellsdis[,cellsvs]
      index=c(index,sum((y)*log10((y+0.001)/(x+0.001))))
    }
  }
  indexcycle[j]=mean(index)
  print(j)
}

par(mar=c(2.5,3.2,0.5,0.5))
indexx=seq(1,15,1)
indexy=indexcycle
plot(indexx,indexy,type='l',col='grey',tck=0.03,las=1,cex=0.3,xlim=c(0.5,15.5),ylim=c(0,0.08),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexx,indexy,type='p',pch=19,col='black',bg=rgb(102/255,205/255,0/255),tck=0.03,las=1,cex=0.7,xlim=c(0.5,15.5),ylim=c(0,0.08),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")

axis(side=1,las=1,at=seq(1,15,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(0,14,1),cex.axis=0.9)
axis(side=2,las=1,at=c(0,0.02,0.04,0.06,0.08),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.02","0.04","0.06","0.08"),cex.axis=1.0)
mtext("HI",side=2,line=2.2,cex=1.0)
mtext("Cycles",side=1,line=1.5,cex=1.0)





##fig3f gai


v=0.2
v=0.02
v=0.06
v=0.04

a=0.5

correlativenumber=10

indexv=vector()
for(i in 1:10)
  indexv[i]=1


for(v in seq(1,10,1))
{
    v=v*0.01
    enhancermethylationlevel=rbeta(1000,3,70)
    
    cellsdis=matrix()
    par(mfrow=c(5,3),mar=c(1,1,1,1))
    for(step in 1:30)
    {
      for(i in 1:1000)
      {
        correlativesites=seq(i-correlativenumber,i+correlativenumber,1)
        correlativesites=correlativesites[correlativesites>0 & correlativesites<1001]
        Mmethylationlevel[i]=mean(enhancermethylationlevel[correlativesites])
      }
      # Mmethylationlevel=mean(enhancermethylationlevel)
      
      
      for(i in 1:1000)
      {
        Ebetaplus=u+u1*(enhancermethylationlevel[i])^k/((enhancermethylationlevel[i])^k+v)+(Mmethylationlevel[i]-enhancermethylationlevel[i])*a
        if(Ebetaplus>1)
          Ebetaplus=1
        if(Ebetaplus<0)
          Ebetaplus=0
        c=m[i]*(Ebetaplus)
        d=m[i]*(1-Ebetaplus)
        enhancermethylationlevel[i]=rbeta(1,c,d)
      }
      x=hist(enhancermethylationlevel[which(m>hold)],xlim=c(0,1),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),ylab="",main="",xaxt="n",yaxt="n",probability = T,bty="l")
      mtext(c("0","0.5","1"),at=c(0,0.5,1),side=1,line=-0.3,cex=0.5)
      
      x=x$counts/(sum(x$counts))
      cellsdis=data.frame(cellsdis,x)
    }
    index=vector()
    for(cells in 6:20)
    {
      x=cellsdis[,cells]
      for(cellsvs in 6:20)
      {
        y=cellsdis[,cellsvs]
        index=c(index,sum((y)*log10((y+0.001)/(x+0.001))))
      }
    }
    indexv[v/0.01]=mean(index)
    print(v)
}



par(mar=c(2.5,3,0.5,0.5))
indexx=seq(1,10,1)
indexy=indexv
plot(indexx,indexy,type='l',col='grey',tck=0.03,las=1,cex=0.3,xlim=c(0.5,10.5),ylim=c(0,0.6),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexx,indexy,type='p',pch=19,col='black',bg=rgb(102/255,205/255,0/255),tck=0.03,las=1,cex=1.2,xlim=c(0.5,10.5),ylim=c(0,0.6),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")

axis(side=1,las=1,at=seq(1,10,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(0.01,0.1,0.01),cex.axis=0.9)
axis(side=2,las=1,at=c(0,0.2,0.4,0.6),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.2","0.4","0.6"),cex.axis=1.0)
mtext("HI",side=2,line=1.5,cex=1.0)
mtext("V",side=1,line=1.5,cex=1.0)

#函数图像
u=0.05
u1=0.9
k=3
v=0.2
x=seq(0,1,0.01)
y=u+u1*x^k/(x^k+v)
par(mar=c(3,4,1,1))
plot(x,y,xlim=c(0,1),ylim=c(0,1),type="l",tck=0.03,las=1,xlab="",col=rgb(30/255,144/255,255/255),ylab="",main="",xaxt="n",yaxt="n",bty="l",lwd=2)
par(new=TRUE)
y=seq(0,1,0.01)
plot(x,y,xlim=c(0,1),ylim=c(0,1),type="l",tck=0.03,las=1,xlab="",col="black",ylab="",main="",xaxt="n",yaxt="n",bty="l")
par(new=TRUE)
k=3
v=0.02
y=u+u1*x^k/(x^k+v)
plot(x,y,xlim=c(0,1),ylim=c(0,1),type="l",tck=0.03,las=1,xlab="",col=rgb(205/255,85/255,85/255),ylab="",main="",xaxt="n",yaxt="n",bty="l",lwd=2)
par(new=T)
k=3
v=0.06
x=seq(0,1,0.01)
y=u+u1*x^k/(x^k+v)
plot(x,y,xlim=c(0,1),ylim=c(0,1),type="l",tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),ylab="",main="",xaxt="n",yaxt="n",bty="l",lwd=2)
axis(side=1,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,0.5,1),cex.axis=1.5)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c(0,0.5,1),cex.axis=1.5)
mtext(expression(phi(x)),side=2,line=2.2,cex=1.5)
mtext(expression(x),side=1,line=1.7,cex=1.5)
legend(0.5,0.4,legend=c("high methy","bistability","low methy"),cex=1.3,x.intersp=0.5,y.intersp=1.5,col=c(rgb(205/255,85/255,85/255),rgb(238/255,44/255,44/255),rgb(30/255,144/255,255/255)),bty="n",lty=1,lwd=1.5)


#fig4a

#index matrix
v=0.2
v=0.02
v=0.06
v=0.04
a=0.5

correlativenumber=10

indexmatrix=vector()
for(i in 1:10)
  indexmatrix[i]=1
indexmatrix=data.frame(indexmatrix)

for(a in seq(0,0.5,0.05))
{
  indexvector=vector()
  for(correlativenumber in seq(10,1,-1))
  {
    enhancermethylationlevel=rbeta(1000,3,70)
    
    cellsdis=matrix()
    par(mfrow=c(5,3),mar=c(1,1,1,1))
    for(step in 1:30)
    {
      for(i in 1:1000)
      {
        correlativesites=seq(i-correlativenumber,i+correlativenumber,1)
        correlativesites=correlativesites[correlativesites>0 & correlativesites<1001]
        Mmethylationlevel[i]=mean(enhancermethylationlevel[correlativesites])
      }
      # Mmethylationlevel=mean(enhancermethylationlevel)
      
      
      for(i in 1:1000)
      {
        Ebetaplus=u+u1*(enhancermethylationlevel[i])^k/((enhancermethylationlevel[i])^k+v)+(Mmethylationlevel[i]-enhancermethylationlevel[i])*a
        if(Ebetaplus>1)
          Ebetaplus=1
        if(Ebetaplus<0)
          Ebetaplus=0
        c=m[i]*(Ebetaplus)
        d=m[i]*(1-Ebetaplus)
        enhancermethylationlevel[i]=rbeta(1,c,d)
      }
      x=hist(enhancermethylationlevel[which(m>hold)],xlim=c(0,1),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),ylab="",main="",xaxt="n",yaxt="n",probability = T,bty="l")
      mtext(c("0","0.5","1"),at=c(0,0.5,1),side=1,line=-0.3,cex=0.5)
      
      x=x$counts/(sum(x$counts))
      cellsdis=data.frame(cellsdis,x)
    }
    index=vector()
    for(cells in 3:17)
    {
      x=cellsdis[,cells]
      for(cellsvs in 3:17)
      {
        y=cellsdis[,cellsvs]
        index=c(index,sum((y)*log10((y+0.001)/(x+0.001))))
      }
    }
    indexvector[11-correlativenumber]=mean(index)
    print(a/0.05)
  }
  indexmatrix=data.frame(indexmatrix,indexvector)
}

par(mar=c(2.5,3,0.5,0.5))
plot(1,1,col='white',xlim=c(0,0.55),ylim=c(0,0.4),tck=0.03,las=1,xlab="",ylab="",main="",xaxt="n",yaxt="n",bty='L')
indexmatrixx=seq(0,0.5,0.05)
par(new=T)
for(i in 1:5)
{
plot(indexmatrixx,as.numeric(as.character(indexmatrix[(i*2-1),]))[-1],xlim=c(0,0.5),ylim=c(0.05,0.75),type='p',col=rgb((295-i*40)/255,0,0),mgp=c(0,0.2,0),tck=0.03,las=1,cex=1,xlab="",ylab="",yaxt="n",xaxt='n',bty='L',pch=16)
par(new=TRUE)
}
for(i in 1:5)
{
  plot(indexmatrixx,as.numeric(as.character(indexmatrix[(i*2-1),]))[-1],xlim=c(0,0.5),ylim=c(0.05,0.75),type='l',col=rgb((295-i*40)/255,0,0),mgp=c(0,0.2,0),tck=0.03,las=1,lwd=1,xlab="",ylab="",yaxt="n",xaxt='n',bty='L')
  par(new=TRUE)
}
mtext("Heterogeneity",side=2,line=2,cex=1)
mtext(expression(alpha),side=1,line=1.3,cex=1)
axis(side=1,las=1,at=seq(0,0.5,0.1),labels=seq(0,0.5,0.1),mgp=c(0,0.2,0),tck=0.01,las=1,cex.axis=1.0)
axis(side=2,las=1,at=c(0,0.25,0.5,0.75),labels=c(0,0.25,0.5,0.75),mgp=c(0,0.2,0),tck=0.01,las=1,cex.axis=1.0)
legend("topleft",legend=c("L=10","L=8","L=6","L=4","L=2"),cex=1.0,x.intersp=0.2,y.intersp=1.5,bty="n",col=c(rgb((6*50-50)/255,0,0),rgb((5*50-50)/255,0,0),rgb((4*50-50)/255,0,0),rgb((3*50-50)/255,0,0),rgb((2*50-50)/255,0,0),rgb((1*50-50)/255,0,0)),lty=1,lwd=1.5)


#subgroup simulation

#canshu
v=0.04
u1=0.9
k=3
u=0.05
m=vector()
for(i in 1:1000)
{
  mi=sample(1:50,1)
  while(runif(1,0,1)>(exp(-mi/3)/1.6))
  {mi=sample(1:50,1)}
  m[i]=mi
}
hist(m)

hold=3

a=0.5

#阈值
correlativenumber=3
correlativesites=vector()
Mmethylationlevel=vector()

#divideparameter
divide=0.3



#ifdivide
ifdivide=vector()

#divideparameter
divide=0.3

#loop time
timelast=3

#chuzhi

averagecellmethylation1=vector()
for(i in 1:15)
  averagecellmethylation1[i]=0
averagecellmethylation1=data.frame(averagecellmethylation1)

index1=vector()

enhancermethylationlevel=rbeta(1000,3,70)
matrix_enhancermethylationlevel=data.frame(enhancermethylationlevel)
for(i in 1:14)
{
  enhancermethylationlevel=rbeta(1000,3,70)
  matrix_enhancermethylationlevel=data.frame(matrix_enhancermethylationlevel,enhancermethylationlevel)
}
#matrix die daidiedai
enhancermethylationlevel=rbeta(1000,3,70)
matrix_enhancermethylationlevel_compare=data.frame(enhancermethylationlevel)
for(i in 1:14)
{
  enhancermethylationlevel=rbeta(1000,3,70)
  matrix_enhancermethylationlevel_compare=data.frame(matrix_enhancermethylationlevel_compare,enhancermethylationlevel)
}

for(simtimes in 1:16)
{
for(timesteps in 1:timelast)
{

#iteration one step

for(cells in 1:15)
{ 
  if(runif(1,0,1)<divide)
  {
    ifdivide[cells]=1
  enhancermethylationlevel=matrix_enhancermethylationlevel[,cells]
  for(i in 1:1000)
  {
    correlativesites=seq(i-correlativenumber,i+correlativenumber,1)
    correlativesites=correlativesites[correlativesites>0 & correlativesites<1001]
    Mmethylationlevel[i]=mean(enhancermethylationlevel[correlativesites])
  }

  for(i in 1:1000)
  {
    Ebetaplus=u+u1*(enhancermethylationlevel[i])^k/((enhancermethylationlevel[i])^k+v)+(Mmethylationlevel[i]-enhancermethylationlevel[i])*a
    if(Ebetaplus>1)
      Ebetaplus=1
    if(Ebetaplus<0)
      Ebetaplus=0
    c=m[i]*(Ebetaplus)
    d=m[i]*(1-Ebetaplus)
    enhancermethylationlevel[i]=rbeta(1,c,d)
  }
  matrix_enhancermethylationlevel[,cells]=enhancermethylationlevel
  }
  else
  {
    ifdivide[cells]=0
  }
}

#subgroup

samplecells=seq(1,(15+sum(ifdivide)),1)
samplecells=sample(samplecells,15)
for(samplecellloop in 1:15)
{
  if(samplecells[samplecellloop]<16)
    matrix_enhancermethylationlevel_compare[,samplecellloop]=matrix_enhancermethylationlevel[,(samplecells[samplecellloop])]
  else
    matrix_enhancermethylationlevel_compare[,samplecellloop]=matrix_enhancermethylationlevel[,which(ifdivide==1)[((samplecells[samplecellloop])-15)]]
}
matrix_enhancermethylationlevel=matrix_enhancermethylationlevel_compare


#fig
cellsdis=matrix()
par(mfrow=c(5,3),mar=c(1,1,1,1))
for(cells in 1:15)
{
enhancermethylationlevel=matrix_enhancermethylationlevel[,cells]
x=hist(enhancermethylationlevel[which(m>hold)],xlim=c(0,1),breaks=seq(0,1,0.05),tck=0.03,las=1,xlab="",col=rgb(238/255,44/255,44/255),ylab="",main="",xaxt="n",yaxt="n",probability = T,bty="l")
mtext(c("0","0.5","1"),at=c(0,0.5,1),side=1,line=-0.3,cex=0.5)
x=x$counts/(sum(x$counts))
cellsdis=data.frame(cellsdis,x)
}
print(timesteps)
}
#jilushuju bao han zai xun h包含在循环nei包含在循环内()
averagecellmethylation=vector()
methylationvec=seq(0.025,0.975,0.05)
for(i in 1:15)
  averagecellmethylation[i]=sum(cellsdis[,(i+1)]*methylationvec)
averagecellmethylation1=data.frame(averagecellmethylation1,averagecellmethylation)

index=vector()
for(cells in seq(2,16,1))
{
  x=cellsdis[,cells]
  for(cellsvs in seq(2,16,1))
  {
    y=cellsdis[,cellsvs]
    index=c(index,sum((y)*log10((y+0.001)/(x+0.001))))
  }
}
index1=c(index1,mean(index))
}



#fig time var
#col=rgb(30/255,144/255,255/255) col=rgb(238/255,44/255,44/255)
par(mfrow=c(1,1),mar=c(3,3,1,1))
for(simtime in 0:5)
{
for(j in 1:16)
{
  vmeanx=vector()
  vmeany=vector()
  for(i in 1:15)
  {vmeanx[i]=rnorm(1,j,0.1)}                    
  plot(vmeanx,averagecellmethylation1[,(simtime*16+j+1)],type="p",tck=0.03,las=1,xlab="",col='black',pch=19, ylab="", main="",xlim=c(0.5,16.5),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',cex=0.3)
  par(new=T)
}
}
axis(side=1,las=1,at=seq(1,16,3),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(3,48,9),cex.axis=0.9)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50","100"),cex.axis=1.0)
mtext("Methylation(%)",side=2,line=1.7,cex=1.0)
mtext("Time(h)",side=1,line=1.5,cex=1.0)

#fig time mean var
par(mar=c(2.5,3.2,0.5,0.5))

indexx=seq(1,16,1)
indexy=vector()
for(i in 1:16)
  indexy[i]=0
for(j in 1:6)
{
for(i in 1:16)
  indexy[i]=indexy[i]+mean(as.numeric(as.character(averagecellmethylation1[,j*16-16+i+1])))
}
plot(indexx,indexy/6,type='l',col='grey',tck=0.03,las=1,cex=0.3,xlim=c(0.5,16.5),ylim=c(0,1),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexx,indexy/6,type='p',pch=19,col='black',bg=rgb(102/255,205/255,0/255),tck=0.03,las=1,cex=0.7,xlim=c(0.5,16.5),ylim=c(0,1),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
axis(side=1,las=1,at=seq(1,16,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(3,48,3),cex.axis=0.9)
axis(side=2,las=1,at=c(0,0.5,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","50","100"),cex.axis=1.0)
mtext("Methylation(%)",side=2,line=1.7,cex=1.0)
mtext("Time(h)",side=1,line=1.5,cex=1.0)

#fig time  var hi
par(mar=c(2.5,3,0.5,0.5))
indexx=seq(1,16,1)
indexy=index1[indexx]+index1[16+indexx]+index1[16*2+indexx]+index1[16*3+indexx]+index1[16*4+indexx]+index1[16*5+indexx]
indexy=indexy/6
plot(indexx,indexy,type='l',col='grey',tck=0.03,las=1,cex=0.3,xlim=c(0.5,16.5),ylim=c(0,0.15),xlab="",ylab="",xaxt='n',yaxt='n',bty="l",lwd=2)
par(new=T)
plot(indexx,indexy,type='p',pch=19,col='black',bg=rgb(102/255,205/255,0/255),tck=0.03,las=1,cex=0.7,xlim=c(0.5,16.5),ylim=c(0,0.15),xlab="",ylab="",xaxt='n',yaxt='n',bty="l")
axis(side=1,las=1,at=seq(1,16,1),mgp=c(0,0.5,0),tck=0.03,las=1,labels=seq(3,48,3),cex.axis=0.9)
axis(side=2,las=1,at=c(0,0.05,0.1,0.15),mgp=c(0,0.5,0),tck=0.03,las=1,labels=c("0","0.05","0.1","0.15"),cex.axis=1.0)
mtext("Heterogeneity",side=2,line=2.1,cex=1.0)
mtext("Time(h)",side=1,line=1.5,cex=1.0)

