##Figures and statistics for antibiotic persistance resistance paper

#FIGURE 1, mutation rate and fitness costs of for rif and kan

#package is needed for the std.error function

require(plotrix)

#setting the directory for the files for this figure and statistics
setwd("C:/Users/michael/Dropbox/Idaho Lab Files/AB_res_pers/figures/fig1_mut_rate_fit/")

#first testing whether the replicative mutation rate differs

############## figure new, graphical display of experimental design

dev.new(width=5, height=4)

plot(1, type="n", axes=F, xlab="Time (days)", ylab="",xlim=c(-10,85),ylim=c(0,10))
axis(side=1,at=c(0,15,30,45,60,75),labels=c("0","15","30","45","60","75"),col="black")
text(x=7.5,y=.25,labels="No Antibiotic")
rect(0,-.3,15.1,.1,xpd=NA,col="#999999",border=FALSE)
text(x=22.5,y=.25,labels="+ Antibiotic")
rect(14.8,-.3,30.2,.1,xpd=NA,col="#191919",border=FALSE)
text(x=52.5,y=.25,labels="No Antibiotic")
rect(30,-.3,75,.1,xpd=NA,col="#999999",border=FALSE)

rect(-5,.75,5,1.25,col="#999999",border=NA)

###########figure 1 part a

#reading in the file
repl_rate <- read.table("rep_mut.csv",header=TRUE,sep=",")

#STATISTICS,simple 1 factor aov with 3 replicates on each side, F=23.28, p=.00849
repl_rate.aov <- aov(data=repl_rate,formula=repl_mut_rate~Antibiotic)
summary(repl_rate.aov)

#STATISTICS STILL, calculating the 95% confidence interval for each value
#calculating for rif
rif_repl <- subset(repl_rate, repl_rate$Antibiotic=="Rif")
rif_repl.mean <- mean(rif_repl$repl_mut_rate)
rif_repl.stderr <-sd(rif_repl$repl_mut_rate)/sqrt(3)
rif_repl.u95 <- rif_repl.mean+rif_repl.stderr*1.96
rif_repl.l95 <- rif_repl.mean-rif_repl.stderr*1.96


kan_repl <- subset(repl_rate, repl_rate$Antibiotic=="Kan")
kan_repl.mean <- mean(kan_repl$repl_mut_rate)
kan_repl.stderr <-sd(kan_repl$repl_mut_rate)/sqrt(3)
kan_repl.u95 <- kan_repl.mean+kan_repl.stderr*1.96
kan_repl.l95 <- kan_repl.mean-kan_repl.stderr*1.96


#creating dataframe with plotting values
names <- c("Kanamycin", "Rifampicin")
means <- c(kan_repl.mean,rif_repl.mean)
stderr <- c(kan_repl.stderr,rif_repl.stderr)
u95 <- c(kan_repl.u95,rif_repl.u95)
l95 <- c(kan_repl.l95,rif_repl.l95)

#creating the dataset for plotting
repl_plot <- as.data.frame(cbind(names,means,stderr,u95,l95))

#adding in location information, this determines how the points are correlated

#deleting objects used in the creation of the dataframe so that they don't create mess later
remove(names,means,stderr,u95,l95,kan_repl.l95,kan_repl.mean,kan_repl.u95,kan_repl.stderr,rif_repl.l95,rif_repl.mean,rif_repl.u95,rif_repl.stderr, kan_repl,rif_repl)

mp <- barplot(as.numeric(as.vector(repl_plot$means)),names.arg=as.vector(repl_plot$names),xlab="Resistance Phenotype",ylab="Mutation rate (per cell per generation)",col="darkgrey",border=NA)
par(mar=c(3,5,3,5))
barplot(as.numeric(as.vector(repl_plot$means)),names.arg=as.vector(stat_plot$names),xlab="Resistance Phenotype",ylab="Mutations per CFU per gen",col="darkgrey",border=NA,ylim=c(0,1.6*10^-7),cex.lab=1.5,cex.names=1.5,axes=FALSE)
fig1_axis <- c(0,5.0*10^-8,1.0*10^-7,1.5*10^-7)
axis(side=2,at=fig1_axis,labels=fig1_axis,cex.axis=1.5)
#preliminary plot, will be combined in actual plot with other graphs
# Plot the vertical lines of the error bars
# The vertical bars are plotted at the midpoints
segments(mp,as.numeric(as.vector(repl_plot$u95)),mp,as.numeric(as.vector(repl_plot$l95)),lwd=4,col="black")
# Now plot the horizontal bounds for the error bars
# 1. The lower bar

segments(mp[1],1.45*10^-7,mp[2],1.45*10^-7,lwd=2,col="black")
segments(mp[1],1.45*10^-7,mp[1],1.3*10^-7,lwd=2,col="black")
segments(mp[2],1.45*10^-7,mp[2],1.3*10^-7,lwd=2,col="black")
text(x=((mp[1]+mp[2]))/2,y=1.5*10^-7,labels="**",cex=3)
mtext(side=3,text="A",at=.3,cex=4,line=-2.5)


###########figure 1 part b
############stationary phase mutation rate KAN
stat_mut <- read.table("stat_mut.csv",header=TRUE,sep=",")
stat_mut$mut_rate <- stat_mut$mut_rate *24
rif_stat <- subset(stat_mut,Antibiotic=="Rif")
kan_stat <- subset(stat_mut,Antibiotic=="Kan")
#paired wilcox rank test used because of difference in sample variance and other complications
#paired is acceptable because values really did come from the same populations
#p=.003906

stat_rate.aov <- aov(data=stat_mut,formula=mut_rate~Antibiotic)
summary(stat_rate.aov)

#prepping for the plot
rif_stat.mean <- mean(rif_stat$mut_rate)
rif_stat.stderr <-std.error(rif_stat$mut_rate)
rif_stat.u95 <- rif_stat.mean+rif_stat.stderr*1.96
rif_stat.l95 <- rif_stat.mean-rif_stat.stderr*1.96

kan_stat.mean <- mean(kan_stat$mut_rate)
kan_stat.stderr <-std.error(kan_stat$mut_rate)
kan_stat.u95 <- kan_stat.mean+kan_stat.stderr*1.96
kan_stat.l95 <- kan_stat.mean-kan_stat.stderr*1.96

names <- c("Kanamycin", "Rifampicin")
means <- c(kan_stat.mean,rif_stat.mean)
stderr <- c(kan_stat.stderr,rif_stat.stderr)
u95 <- c(kan_stat.u95,rif_stat.u95)
l95 <- c(kan_stat.l95,rif_stat.l95)

#creating the dataset for plotting
stat_plot <- as.data.frame(cbind(names,means,stderr,u95,l95))


#deleting objects used in the creation of the dataframe so that they don't create mess later
remove(names,means,stderr,u95,l95,kan_stat.l95,kan_stat.mean,kan_stat.u95,kan_stat.stderr,rif_stat.l95,rif_stat.mean,rif_stat.u95,rif_stat.stderr, kan.slopes,rif.slopes)



#preliminary plot, will be combined in actual plot with other graphs for this 3 part figure 
layout(1)
mp <- barplot(as.numeric(as.vector(stat_plot$means)),names.arg=as.vector(stat_plot$names),xlab="Resistance Phenotype",ylab="Mutation rate (per cell per day)",col="darkgrey",border=NA)
par(mar=c(3,5,3,5))
barplot(as.numeric(as.vector(stat_plot$means)),names.arg=as.vector(stat_plot$names),ylab="Mutations per CFU per day",col="darkgrey",border=NA,ylim=c(1*10^-8,1*10^-6),cex.lab=1.5,cex.names=1.5,axes=FALSE)
fig1b_axis <- c(1*10^-8,5*10^-7,1.0*10^-6)
axis(side=2,at=fig1b_axis,labels=fig1b_axis,cex.axis=1.5)



segments(mp,as.numeric(as.vector(stat_plot$u95)),mp,as.numeric(as.vector(stat_plot$l95)),lwd=4,col="black")

#signifigance display
segments(mp[1],9*10^-7,mp[2],9*10^-7,lwd=2,col="black")
segments(mp[1],9*10^-7,mp[1],8.6*10^-7,lwd=2,col="black")
segments(mp[2],9*10^-7,mp[2],1.6*10^-7,lwd=2,col="black")
text(x=((mp[1]+mp[2]))/2,y=9.3*10^-7,labels="**",cex=3)


#############Figure 1 part c, fitness assay

fitness_0 <- read.table("rif_kan_mutant_fitness.csv",sep=",",header=T)

#calculated means and std deviation

fitness_means <- aggregate(as.data.frame(fitness_0$Fitness), list(Anti=fitness_0$Antibiotic,Number=fitness_0$Number), FUN="mean")
std.error <- aggregate(as.data.frame(fitness_0$Fitness), list(Anti=fitness_0$Antibiotic,Number=fitness_0$Number), FUN="std.error")

std.error$Anti <- NULL
std.error$Number <- NULL

#piecing the means back together
fit_ave_0 <- cbind(fitness_means,std.error)
#renaming columns
colnames(fit_ave_0) <- c("Antibiotic","Number", "fit_mean","stderror")
#calculating the 95% confidence intervals
fit_ave_0$u95 <- fit_ave_0$fit_mean+(1.96*fit_ave_0$stderror)
fit_ave_0$l95 <- fit_ave_0$fit_mean-(1.96*fit_ave_0$stderror)

plotCI(x=fit_ave_0$Number, y=fit_ave_0$fit_mean,ui=fit_ave_0$u95,li=fit_ave_0$l95, pch=16,xlab="Clone Number",ylab="Relative Fitness")


fit_0_stats <- aov(fit_ave_0$fit_mean ~fit_ave_0$Antibiotic)
summary(fit_0_stats)
wilcox.test(formula=fit_ave_0$fit_mean~fit_ave_0$Antibiotic,alternative="two.sided",exact=FALSE)

par(mar=c(3,5,3,5))
boxplot(fit_ave_0$fit_mean~fit_ave_0$Antibiotic,ylab="Relative Fitness",cex.lab=1.5,cex.axis=1.5,pch=16,frame=F,xaxt='n',ylim=range(0.0:2.0))

segments(1,1.8,2,1.8,lwd=2,col="black")
segments(1,1.8,1,1.7,lwd=2,col="black")
segments(2,1.8,2,1.7,lwd=2,col="black")
text(x=1.5,y=1.9,labels="***",cex=3)
mtext(side=3,text="C",at=.8,cex=4,line=-2.5)
#mtext(side=1,text="Resistance Phenotype",at=1.5,cex=1.5,padj=2)
mtext(side=1,text="Kanamycin",at=1,cex=1.5,padj=.6)
mtext(side=1,text="Rifampicin",at=2,cex=1.5,padj=.6)


#setting up the screen for figure 1
split.screen(c(1,2))
screen(1)
split.screen(c(2,1))
screen(3)
screen(4)
screen(2)

###########FIGURE 2
########### biofilm concentrations over time

#setting the directory
setwd("C:/Users/michael/Dropbox/Idaho Lab Files/AB_res_pers/figures/fig2_biofilm_tracking/")


###########kan
kan_biofilm <- read.table("kan_biofilm.csv",sep=",",header=TRUE)
kan_biofilm$log_res <- log(kan_biofilm$Fraction_Res,base=10)

kan_biofilms_means <- aggregate(as.data.frame(kan_biofilm$log_res), list(Time = kan_biofilm$Day), FUN="mean")
kan_biofilm_stder <- aggregate(as.data.frame(kan_biofilm$log_res), list(Time = kan_biofilm$Day), FUN="std.error")

kan_biofilm_stder$Time <- NULL
kan_biofilm_plot <- cbind(kan_biofilms_means,kan_biofilm_stder)

colnames(kan_biofilm_plot) <- c("Time","mean","stderr")

kan_biofilm_plot$u95 <- kan_biofilm_plot$mean+kan_biofilm_plot$stderr*1.96
kan_biofilm_plot$l95 <- kan_biofilm_plot$mean-kan_biofilm_plot$stderr*1.96


#reading in rif data
layout(1)

rif_biofilm <- read.table("rif_biofilm.csv",sep=",",header=TRUE)
rif_biofilm$log_res <- log(rif_biofilm$Res_Frac,base=10)

rif_biofilms_means <- aggregate(as.data.frame(rif_biofilm$log_res), list(Time = rif_biofilm$Time), FUN="mean")
rif_biofilm_stder <- aggregate(as.data.frame(rif_biofilm$log_res), list(Time = rif_biofilm$Time), FUN="std.error")

rif_biofilm_stder$Time <- NULL
rif_biofilm_plot <- cbind(rif_biofilms_means,rif_biofilm_stder)

colnames(rif_biofilm_plot) <- c("Time","mean","stderr")

rif_biofilm_plot$u95 <- rif_biofilm_plot$mean+rif_biofilm_plot$stderr*1.96
rif_biofilm_plot$l95 <- rif_biofilm_plot$mean-rif_biofilm_plot$stderr*1.96

#plotting

layout(matrix(1:3,3,1),widths=c(.5,.5,.5),heights=c(1,5,5))
par(mar=c(0,6,0,2))
plot(1, type="n", axes=F, xlab="", ylab="",xlim=c(0,75),ylim=c(0,1))
text(x=7,y=.2,labels="No Antibiotic",cex=1.5)
rect(-1,-.5,15.1,0,xpd=NA,col="gray68",border=FALSE)
text(x=22.5,y=.2,labels="+ Antibiotic",cex=1.5)
rect(14.8,-.5,30.2,0,xpd=NA,col="green4",border=FALSE)
text(x=52.5,y=.2,labels="No Antibiotic",cex=1.5)
rect(30,-.5,76,0,xpd=NA,col="gray68",border=FALSE)


par(mar=c(5,6,3,2))
plotCI(x=kan_biofilm_plot$Time,y=kan_biofilm_plot$mean,ui=kan_biofilm_plot$u95,li=kan_biofilm_plot$l95,pch=16,axes=FALSE,col=c("gray68","gray68","green4","green4","green4","green4","gray68","gray68","gray68"),ylim=c(-8,0),xlim=c(0,75),cex=1.5,cex.lab=1.2,xlab="",ylab="")
points(x=kan_biofilm_plot$Time,y=kan_biofilm_plot$mean,type="b",col=c("gray68","gray68","green4","green4","green4","green4","gray68","gray68","gray68"),lwd=2,cex=2,pch=16)
segments(x0=kan_biofilm_plot$Time[2],x1=kan_biofilm_plot$Time[3],y0=kan_biofilm_plot$mean[2],y1=kan_biofilm_plot$mean[3],col="green4",lwd=4)
segments(x0=kan_biofilm_plot$Time[3],x1=kan_biofilm_plot$Time[4],y0=kan_biofilm_plot$mean[3],y1=kan_biofilm_plot$mean[4],col="green4",lwd=4)
segments(x0=kan_biofilm_plot$Time[4],x1=kan_biofilm_plot$Time[5],y0=kan_biofilm_plot$mean[4],y1=kan_biofilm_plot$mean[5],col="green4",lwd=4)
segments(x0=kan_biofilm_plot$Time[5],x1=kan_biofilm_plot$Time[6],y0=kan_biofilm_plot$mean[5],y1=kan_biofilm_plot$mean[6],col="green4",lwd=4)



box(bty="L",lwd=2)
#axis(side=1,at=c(0,15,30,45,60,75),)
labelsY=parse(text=c(parse(text='10^-8'),parse(text='10^-6'),parse(text='10^-4'),parse(text='10^-2'),parse(text='1')))
axis(side=2,at=c(-8,-6,-4,-2,0),labels=labelsY,las=1,cex.axis=1.5)
text(x=5,y=-1,xpd=NA,labels="A",cex=4)

par(mar=c(5,6,3,2))
plotCI(x=rif_biofilm_plot$Time,y=rif_biofilm_plot$mean,ui=rif_biofilm_plot$u95,li=rif_biofilm_plot$l95,pch=16,axes=FALSE,xlab="Time (days)",ylab="",col=c("grey48","grey48","green4","green4","green4","grey48","grey48","grey48"),ylim=c(-8,0),xlim=c(0,75),cex=1.5,cex.lab=1.5,cex.axis=1.5)
points(x=rif_biofilm_plot$Time,y=rif_biofilm_plot$mean,type="b",col=c("grey48","grey48","green4","green4","green4","grey48","grey48","grey48"),lwd=2,cex=2,pch=16)
segments(x0=rif_biofilm_plot$Time[2],x1=rif_biofilm_plot$Time[3],y0=rif_biofilm_plot$mean[2],y1=rif_biofilm_plot$mean[3],col="green4",lwd=4)
segments(x0=rif_biofilm_plot$Time[3],x1=rif_biofilm_plot$Time[4],y0=rif_biofilm_plot$mean[3],y1=rif_biofilm_plot$mean[4],col="green4",lwd=4)
segments(x0=rif_biofilm_plot$Time[4],x1=rif_biofilm_plot$Time[5],y0=rif_biofilm_plot$mean[4],y1=rif_biofilm_plot$mean[5],col="green4",lwd=4)
box(bty="L",lwd=2)
axis(side=1,at=c(0,15,30,45,60,75),cex.axis=1.5)
text(x=5,y=-1,xpd=NA,labels="B",cex=4)

labelsY=parse(text=c(parse(text='10^-8'),parse(text='10^-6'),parse(text='10^-4'),parse(text='10^-2'),parse(text='1')))
axis(side=2,at=c(-8,-6,-4,-2,0),labels=labelsY,las=1,cex.axis=1.5)
text(x=-10.5,y=1.5,srt=90,labels="Fraction of the biofilm resistant",xpd=NA,cex=1.5)

#######Figure 3, comparison of rate of accumulation to mutation rates

## converting estimates of mutation rates in generations to mutation rates per unit time



########################FIGURE 4 PLANKTONIC
setwd("~/Dropbox/Idaho Lab Files/AB_res_pers/figures/fig3_plank_tracking/")


###########reading data sets
kan_plank <- read.table("kan_plank.csv",sep=",",header=TRUE)
rif_plank <- read.table("rif_plank.csv",sep=",",header=TRUE)

#######making log resistance
kan_plank$log_ratio <- log(kan_plank$Fraction,base=10)
rif_plank$log_ratio <- log(rif_plank$Fraction,base=10)


#######subsetting data

kan_plank.1 <- subset(kan_plank, kan_plank$Biofilm==1)
kan_plank.1 <- kan_plank.1[with(kan_plank.1, order(Generations)), ]

kan_plank.2 <- subset(kan_plank, kan_plank$Biofilm==2)
kan_plank.2 <- kan_plank.2[with(kan_plank.2, order(Generations)), ]

kan_plank.3 <- subset(kan_plank, kan_plank$Biofilm==4)
kan_plank.3 <- kan_plank.3[with(kan_plank.3, order(Generations)), ]



rif_plank.1 <- subset(rif_plank, rif_plank$Biofilm==1)
rif_plank.1 <- rif_plank.1[with(rif_plank.1, order(Generations)), ]

rif_plank.2 <- subset(rif_plank, rif_plank$Biofilm==2)
rif_plank.2 <- rif_plank.2[with(rif_plank.2, order(Generations)), ]

rif_plank.3 <- subset(rif_plank, rif_plank$Biofilm==3)
rif_plank.3 <- rif_plank.3[with(rif_plank.3, order(Generations)), ]











#1 point from b2 and 1 point from b4 were trimmed out, mistakes in plating resulted in no data
#b2 day 14, b4 day 24
layout(matrix(1:2,2,1),widths=c(1,1),heights=c(1,1))
par(mar=c(1,6,4,2))
plot(x=kan_plank.1$Generations,y=kan_plank.1$log_ratio,type="l",lwd=3,ylim=range(0:-8),ylab="",xlab="",xlim=c(0,250),bty="n",axes=FALSE)

labelsY=parse(text=c(parse(text='10^-8'),parse(text='10^-6'),parse(text='10^-4'),parse(text='10^-2'),parse(text='1')))
axis(side=2,at=c(-8,-6,-4,-2,0),labels=labelsY,las=1)

#axis(side=1,at=c(0,50,100,150,200,250))
box(bty="L",lwd=2)
text(x=245,y=.8,labels="A",xpd=NA,cex=2)
#rectangle for overnight fraction of mutants and 95% confidence interval
rect(0,-8.0639,260,-7.555,col="lightgray",border=NA)

points(x=kan_plank.2$Generations,y=kan_plank.2$log_ratio,type="l",lwd=3)
points(x=kan_plank.3$Generations,y=kan_plank.3$log_ratio,type="l",lwd=3)

### plotting rif
par(mar=c(4,6,1,2))
plot(x=rif_plank.1$Generations,y=rif_plank.1$log_ratio,type="l",lwd=3,ylim=range(0:-8),ylab="",xlab="Generations",xlim=c(0,250),bty="n",axes=FALSE)
labelsY=parse(text=c(parse(text='10^-8'),parse(text='10^-6'),parse(text='10^-4'),parse(text='10^-2'),parse(text='1')))
axis(side=2,at=c(-8,-6,-4,-2,0),labels=labelsY,las=1)

axis(side=1,at=c(0,50,100,150,200,250))
box(bty="L",lwd=2)
text(x=245,y=.8,labels="B",xpd=NA,cex=2)
#rectangle for overnight fraction of mutants and 95% confidence interval
rect(0,-6.888175,260,-6.677918,col="lightgray",border=NA)

points(x=rif_plank.2$Generations,y=rif_plank.2$log_ratio,type="l",lwd=3)
points(x=rif_plank.3$Generations,y=rif_plank.3$log_ratio,type="l",lwd=3)
text(x=-50,y=2,srt=90,labels="Fraction Resistant",xpd=NA,cex=1.2)




#mutation rate comparison

setwd("C:/Users/Michael/Dropbox/Idaho Lab Files/AB_res_pers/figures/mutation_rate_accumulation/")

mut_rate_comp <- read.table("mutation_rate_counts.csv",header=TRUE,sep=",")
mut_rate_comp_kan <- subset(mut_rate_comp,antibiotic=='kan')
mut_rate_comp_rif <- subset(mut_rate_comp,antibiotic=='rif')

par(mar=c(5,5,3,3))
plot(x=mut_rate_comp_kan$time[2:11],y=mut_rate_comp_kan$mutation_rate[2:11],type='l',lwd=3,xlim=range(-3:12),xaxt="n",yaxt="n",ylim=c(0,1.05*10^-6),xlab="",ylab="",cex.lab=1.5)
points(x=mut_rate_comp_kan$time[1],y=mut_rate_comp_kan$mutation_rate[1],pch=16,cex=2)
axis(side = 1,at = c(-2,1,4,7,10),labels = c("","1","4","7","10"),cex.axis=1.5)
mtext(side = 1,line = 1,at = -2,text="Stationary",cex=1.5)
mtext(side = 1,line = 2.5,at = -2,text="Phase",cex=1.5)
mtext(side=1,line=2.5,at=5.7,text = "Number of generations per day",cex=1.5)
segments(x0 = -2,x1=11,y0=6.79*10^-8,y1=6.79*10^-8,lwd=3,col='ivory4')
text(x=6,y=1*10^-7,labels = "Observed rate in biofilms",cex=1.5,col="ivory4")

axis(side = 2,at=c(1*10^-8,2.5*10^-7,5*10^-7,7.5*10^-7,1*10^-6),cex.axis=1.5)
mtext(side=2,line=2.75,at=5*10^-7,text="Per day mutation rate",cex=1.5)
text(x=-2,y=9.8*10^-7,labels="A",cex=4)


par(mar=c(5,5,3,3))
plot(x=mut_rate_comp_rif$time[2:11],y=mut_rate_comp_rif$mutation_rate[2:11],type='l',lwd=3,xlim=range(-3:12),xaxt="n",yaxt="n",ylim=c(0,4.5*10^-7),xlab="",ylab="",cex.lab=1.5)
points(x=mut_rate_comp_rif$time[1],y=mut_rate_comp_rif$mutation_rate[1],pch=16,cex=2)
axis(side = 1,at = c(-2,1,4,7,10),labels = c("","1","4","7","10"),cex.axis=1.5)
mtext(side = 1,line = 1,at = -2,text="Stationary",cex=1.5)
mtext(side = 1,line = 2.5,at = -2,text="Phase",cex=1.5)
mtext(side=1,line=2.5,at=5.7,text = "Number of generations per day",cex=1.5)

segments(x0 = -2,x1=11,y0=1.60*10^-7,y1=1.60*10^-7,lwd=3,col='ivory4')
text(x=8,y=1.35*10^-7,labels = "Observed rate in biofilms",cex=1.5,col="ivory4")

axis(side = 2,at=c(1*10^-8,1*10^-7,2*10^-7,3*10^-7,4*10^-7),cex.axis=1.5)
mtext(side=2,line=2.75,at=2*10^-7,text="Per day mutation rate",cex=1.5)
text(x=-2,y=4.25*10^-7,labels="B",cex=4)
