
library(gdata)
library(plyr)
library(ggplot2)
library(reshape2)
library(AICcmodavg)
library(dismo)
library(gbm)
library(pROC)
library(grid)

## set working directory (must change to run code from GitHub)
setwd(dir="~/Dropbox/Rasmussen and Harper brach database/Biology.Letters")

## read in processed dataset, as output by "Process.Data.R".

data <- read.csv("OS.brachs.processed.csv",header=TRUE)

data$stage_top <- data$Interval_Top

#make a column indicating whether each genus is a Lazarus (range-through) genus in the suceeding timeslice
data$stage <- factor(data$stage_top, levels=c(1:23))
data$stage2 <- factor(data$stage_top, levels=c(1:23))
data <- data[order(data$Genus, -data$stage_top),]
nexts <- data[c("Interval","Genus","stage_top")]
rows <- nrow(nexts)
nexts <- nexts[(2:nrow(nexts)),]
nexts[rows,] <- c(NA,NA,NA)
rownames(nexts) <- seq(nrow(nexts))
data$Gap.Next <- ifelse(data$Genus == nexts$Genus & data$Interval != (nexts$Interval -1),1,0)

## calculate max occupancy for each interval

max.occ <- function(df) max(df$occupancy)
maxoccs <- ddply(data,.(Interval),each(max.occ))
data <- merge(data,maxoccs)

## round predictors to reduce false precision and so that number of categories is similar in all cases

round.fac <- 10

data$log.richness <- round(log(data$richness),0)
data$occupancy <- data$occupancy/data$max.occ
data$occupancy <- round((data$occupancy*round.fac),0)*round.fac
data$Prop.Cratonic <- round((data$Prop.Cratonic*round.fac),0)*round.fac
data$PropTrunc <- round((data$PropTrunc*round.fac),0)*round.fac
data$MaxLat <- ceiling(data$MaxLat/round.fac)*round.fac
data$MinLat <- floor(data$MinLat/round.fac)*round.fac
data$AbsLatRange <- round((data$AbsLatRange/round.fac),0)*round.fac
data$MeanAbsLat <- round(data$MeanAbsLat /round.fac,0) *round.fac
data$gcd <- round(data$gcd/1000,0)*1000
data$log.richness <- round(data$log.richness,1)
data$MeanDepth <- data$MeanDepth
data$DepthRange <- round(data$DepthRange,1)
data$MinDepth <- round(data$MinDepth,0)
data$MaxDepth <- round(data$MaxDepth,0)
data$Prop.North <- round((data$Prop.North*round.fac),0)*round.fac
data$cooling.latrange <- round(data$cooling.latrange/round.fac,0)*round.fac
data$cooling.latrange <- ifelse(data$cooling.latrange==91,90,data$cooling.latrange)
data$cooling.start.latrange <- round(data$cooling.start.latrange/round.fac,0)*round.fac
data$warming.latrange <- round(data$warming.latrange/round.fac,0)*round.fac
data$warming.latrange <- ifelse(data$warming.latrange==91,90,data$warming.latrange)
data$warming.start.latrange <- round(data$warming.start.latrange/round.fac,0)*round.fac
data$GH.IH.extinct <- ifelse(data$cooling.latrange==0,1,0)
data$IH.GH.extinct <- ifelse(data$warming.latrange==0,1,0)
data$Age <- round(data$Age,0)

## tabulate number of extinctions per timeslice
numex <- function(df) sum(df$Ex)
means <- ddply(data,.(stage_top),each(numex))
data <- merge(data,means)

#reduce to Katian and latest Katian
# Katian  cutoff: 456
data <- data[(data$stage_top > 443 & data$stage_top < 456 ),]
data$interval <- data$Interval_Top

## examine extinctions per interval
num.ex <- function(df) sum(df$Ex)
exes <- ddply(data,.(interval),each(num.ex ))

##exclude intervals with few extinctions
drops <- c(448.5,446.8,446.3,446.0)

data <- data[ ! data$interval %in% drops, ]

## make data frame with only variables to include in models
data2 <- data[c("Order","Genus","interval","Ex","log.richness","gcd","Prop.Cratonic","PropTrunc","MaxDepth","MinDepth", "DepthRange","MaxLat","MinLat","AbsLatRange","occupancy","continents","Age","Pedicle.state")]

LateK   <- data2[(data2$interval == 445.6),]
EarlyK <- data2[(data2$interval > 446),]
Hirn <- data2[(data2$interval == 443.7),]

### remove weird duplicates
Hirn$cull.factor <- paste(Hirn$Genus,Hirn$gcd,Hirn$PropTrunc)
hirn.cull <- c("Clorilamnulella 0 100","Plectothyrella 0 100")
Hirn <- Hirn[ ! Hirn$cull.factor %in% hirn.cull, ]

LateK$cull.factor <- paste(LateK$Genus,LateK$gcd,LateK$PropTrunc)
LateK.cull <- c("Christiania 0 100","Diambonia 0 100","Kiaeromena 0 0")
LateK <- LateK[ ! LateK$cull.factor %in% LateK.cull, ]

## Now make model for each interval, contrast predictions
df <- EarlyK

##make weights proportional to class frequency
obs <- length(df$Ex)
exes <- sum(df$Ex)
ExFreq <- (exes/obs)
SurFreq <- (1-(exes/obs))
MaxFreq <- max(ExFreq,SurFreq)
ExWeight <- 1/(ExFreq/MaxFreq)
SurWeight <- 1/(SurFreq/MaxFreq)
train.weights <- ifelse(df$Ex==1,ExWeight,SurWeight)
null.weights <- rep(1,length(train.weights))

## model for all analyses
EarlyK.mod <- gbm.step(data=df , gbm.x = 5:18, gbm.y = 4,family = "bernoulli", tree.complexity = 5,learning.rate = 0.001, bag.fraction = .5,site.weights=train.weights,n.folds=20,n.trees=100,step.size=100,max.trees=10000,silent=FALSE)


rel.infs <- data.frame(summary(EarlyK.mod))
colnames(rel.infs) <- c("predictor","E.Katian.rel.inf.")

rel.infs$predictor <- revalue(rel.infs$predictor,c("AbsLatRange"="Paleolatitude Range","continents"="Num. Paleocontinents","DepthRange"="Depth Range","gcd" = "Great Circle Dist.","log.richness"="Log Richness","MaxDepth"="Maximum Depth","MaxLat"="Max. Paleolatitude","MinDepth"="Minimum Depth","MinLat"="Min. Paleolatitude","occupancy"="Grid Cell Occupancy","Prop.Cratonic"="% Cratonic","PropTrunc"="% Discontinuous","Age" =  "Genus Age","Pedicle.state" = "Pedunculate"))

rel.infs$predictor = factor(rel.infs$predictor, levels=rel.infs$predictor[order(rel.infs$E.Katian.rel.inf.)], ordered=TRUE)

pdf("EarlyK.Relative.Influence.plot.pdf", width = 4.5, height =4.4) 
p <- ggplot(rel.infs,aes(E.Katian.rel.inf.,predictor))
p + geom_point(size=4) + theme_bw() + xlab("Relative Influence (%)") + ylab("Predictor")
dev.off()

## make predictions

prediction <- predict(EarlyK.mod,newdata=EarlyK,n.trees=EarlyK.mod$n.trees,type="response")
LateK.pred <- predict(EarlyK.mod,newdata=LateK,n.trees=EarlyK.mod$n.trees,type="response")
Hirn.pred <- predict(EarlyK.mod,newdata=Hirn,n.trees=EarlyK.mod$n.trees,type="response")

EarlyK$dataset <- rep("EarlyK",nrow(EarlyK))
EarlyK.pred <- predict(EarlyK.mod,newdata=EarlyK,n.trees=EarlyK.mod$n.trees,type="response")
EarlyK <- cbind(EarlyK,EarlyK.pred)
EarlyK <- EarlyK[order(EarlyK.pred),] 
EarlyK$risk.seq <- seq(1,nrow(EarlyK))

# make ROC plots
makeROC <- function(df){
  df <- unique(df)
  ROC <- roc(df$Ex ~ df$EarlyK.pred,ci=TRUE, of="auc")
  sensitivity <- ROC$sensitivities
  specificity <- ROC$specificities
  taxa <- rep(length(df$Ex),length(sensitivity))
  ex.prop <- rep(mean(df$Ex),length(sensitivity))
  ex.num <- rep(sum(df$Ex),length(sensitivity))
  AUC <- rep(ROC$auc,length(sensitivity))
  AUC.975 <- rep(ROC$ci[3],length(sensitivity))
  AUC.025 <- rep(ROC$ci[1],length(sensitivity))
  sequence <- seq(1,length(sensitivity))
  output <- data.frame(sensitivity,specificity,taxa,ex.prop,ex.num,AUC,AUC.975,AUC.025,sequence)
  return(output)}

eK.ROC <- ddply(EarlyK,.(dataset),each(makeROC))

eK.ROC.single <- eK.ROC[(eK.ROC$sequence==1),]

eK.ROC.single$AUC <- paste(rep("AUC :",nrow(eK.ROC.single)),round(eK.ROC.single$AUC,2))
eK.ROC.single$taxa <- paste(rep("Gen. :",nrow(eK.ROC.single)),eK.ROC.single$taxa)
eK.ROC.single$ex.num <- paste(rep("N ext:",nrow(eK.ROC.single)),eK.ROC.single$ex.num)
eK.ROC.single$x <- rep(0,nrow(eK.ROC.single))
eK.ROC.single$y <- rep(.99,nrow(eK.ROC.single))

pdf("EarlyK.ROC.plot.pdf", width = 4, height =4) 
p <- ggplot(eK.ROC,aes(1-specificity,sensitivity))
p + geom_abline(linetype=3) + geom_path(size=1) + coord_equal() + theme_bw() + xlab("False extinction prediction rate")+ ylab("True extinction prediction rate")  + geom_text(data=eK.ROC.single,aes(x=.55,y=.35,label=AUC),size=5,colour="black",hjust=0) + geom_text(data=eK.ROC.single,aes(x=.55,y=.22,label=taxa),size=5,colour="black",hjust=0) + geom_text(data=eK.ROC.single,aes(x=.55,y=.09,label=ex.num),size=5,colour="darkgray",,hjust=0) 
dev.off()

LateK$dataset <- rep("LateK",nrow(LateK))
EarlyK.pred <- predict(EarlyK.mod,newdata=LateK,n.trees=EarlyK.mod$n.trees,type="response")
LateK <- cbind(LateK,EarlyK.pred)
LateK <- LateK[order(EarlyK.pred),] 
LateK$risk.seq <- seq(1,nrow(LateK))

Hirn$dataset <- rep("Hirn",nrow(Hirn))
EarlyK.pred <- predict(EarlyK.mod,newdata=Hirn,n.trees=EarlyK.mod$n.trees,type="response")
Hirn <- cbind(Hirn,EarlyK.pred)
Hirn <- Hirn[order(EarlyK.pred),] 
Hirn$risk.seq <- seq(1,nrow(Hirn))

AllData <- rbind(LateK,Hirn)

### make bar charts
AllData$Outcome <- ifelse(AllData$Ex == 1,"Extinction","Survivor")
AllData$Interval <- as.factor(ifelse(AllData$dataset == "LateK","A. Late Katian","B. Late Hirnantian"))
AllData$Interval <- relevel(AllData$Interval, "A. Late Katian")

AllData.single <- AllData[(AllData$risk.seq == 1),]
AllData.single$x <- rep(.45,nrow(AllData.single))
AllData.single$y <- rep(5,nrow(AllData.single))

Foliomena <- c("Foliomena","Christiania","Kassinella","Dedzetina","Kozlowskites","Cyclospira","Leptestiina")

Foliomena_fauna <- AllData[(AllData$Genus %in% Foliomena & AllData$Interval == "A. Late Katian"),]
Foliomena_fauna <- Foliomena_fauna[1:4,]

#pdf("Genus.Risk.Barplot.pdf", width = 10, height =6) 
#p <- ggplot(AllData,aes(risk.seq,EarlyK.pred,colour=Outcome,label=Genus))
#p + geom_segment(aes(yend = 0,xend=risk.seq),size=1.2) + facet_wrap(~Interval,ncol=1) + geom_text(size=1.3,colour="black",angle = 90,hjust = -.05) + coord_cartesian(ylim=c(0,1.2),xlim=c(0,203)) + theme_bw() + theme(legend.justification=c(1.8,0), legend.position=c(1,.1)) + theme(legend.title = element_text(colour="black", size=16,face="plain"),legend.text = element_text(colour="black", size = 16)) + ylab("Relative risk prediction from early Katian model") +  theme(strip.text.x = element_text(size = 16)) + theme(axis.title=element_text(size=16))+ theme(strip.background = element_rect(fill = NA,colour=NA)) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + scale_colour_manual(values=c("black","gray")) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + theme(axis.title.x = element_blank()) + geom_segment(data=Foliomena_fauna,aes(y=EarlyK.pred + .4,yend =EarlyK.pred + .2,x=risk.seq,xend=risk.seq),arrow = arrow(length = unit(0.1,"cm"))) 
#dev.off()

pdf("Genus.Risk.Barplot.Vertical.pdf", width = 6.32, height =6.87) 
p <- ggplot(AllData,aes(EarlyK.pred,risk.seq,colour=Outcome,label=Genus))
p + geom_segment(aes(xend = 0,yend=risk.seq),size=.9) + facet_wrap(~Interval,ncol=2) + geom_text(size=1,colour="black",hjust = -.09,face="bold") + coord_cartesian(xlim=c(0,1.1),ylim=c(0,203)) + theme_bw() + theme(legend.justification=c(1.8,0), legend.position=c(1.05,.15)) + theme(legend.title = element_text(colour="black", size=14,face="plain"),legend.text = element_text(colour="black", size = 14)) + xlab("Relative risk prediction from early Katian model") + theme(strip.background = element_blank(),strip.text.x = element_blank()) + theme(axis.title=element_text(size=12)) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + scale_colour_manual(values=c("black","gray")) + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank()) + theme(axis.title.y = element_blank()) + geom_segment(data=Foliomena_fauna,aes(x=EarlyK.pred + .3,xend =EarlyK.pred + .1,y=risk.seq,yend=risk.seq),arrow = arrow(length = unit(0.1,"cm"))) + scale_y_reverse() + geom_text(data=AllData.single,aes(x,y,label=Interval),colour="black",hjust = 0)
dev.off()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

ROCs <- ddply(AllData,.(dataset),each(makeROC))

ROCs$Interval <- as.factor(ifelse(ROCs$dataset == "LateK","C. Late Katian","D. Late Hirnantian"))
ROCs$Interval <- relevel(ROCs$Interval, "C. Late Katian")

ROCs.single <- ROCs[(ROCs$sequence==1),]
ROCs.single$AUC <- paste(rep("AUC :",nrow(ROCs.single)),round(ROCs.single$AUC,2))
ROCs.single$taxa <- paste(rep("Gen. :",nrow(ROCs.single)),ROCs.single$taxa)
ROCs.single$ex.num <- paste(rep("N ext:",nrow(ROCs.single)),ROCs.single$ex.num)
ROCs.single$x <- rep(0,nrow(ROCs.single))
ROCs.single$y <- rep(.99,nrow(ROCs.single))

pdf("Interval.ROC.plots.pdf", width = 6.25, height =3.5) 
p <- ggplot(ROCs,aes(1-specificity,sensitivity))
p + geom_abline(linetype=3) + geom_path(size=1) + facet_wrap(~Interval,ncol=2) + coord_equal() + theme_bw() + xlab("False extinction prediction rate")+ ylab("True extinction prediction rate") + theme(strip.background = element_blank(),strip.text.x = element_blank()) + geom_text(data=ROCs.single,aes(x=.55,y=.35,label=AUC),size=5,colour="black",hjust=0) + geom_text(data=ROCs.single,aes(x=.55,y=.22,label=taxa),size=5,colour="black",hjust=0) + geom_text(data=ROCs.single,aes(x=.55,y=.09,label=ex.num),size=5,colour="darkgray",,hjust=0) + theme(panel.margin = unit(.4, "lines")) + geom_text(data=ROCs.single,aes(x,y,label=Interval),colour="black",hjust = 0)
dev.off()

