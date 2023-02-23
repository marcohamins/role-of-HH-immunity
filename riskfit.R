########  Details ########
### Author: Marco Hamins-Puertolas 
### Code to output figures and risk estimates from manuscript

###### Read in necessary libraries #####
library(patchwork)
library(readxl)
library(ggplot2)
library(scales)
library(plyr)
library(lubridate)
library(magrittr)
library(tidymv)
library(glmmTMB)
library(ggpubr)
library(stringr)

##### Import data and set up analyses ######
curDir <- getwd() # set working directory to /role-of-HH-immunity on local comp
figFolder <- paste0(curDir,"/Figures/")

# Read in results from manuscript
load(file = paste0(curDir,"/Data/importanceMat.rData"))
load(file = paste0(curDir,"/Data/housedfCohortData.rData"))

# Read in results from github file modelfit.R
# load(file = paste0(curDir,"/Data/importanceMat_update.rData"))
# load(file = paste0(curDir,"/Data/housedfCohortData_update.rData"))

df <- householdCohortData

impMat <- importanceMat

# rename features
impMat$Feature <- rep(c("Gender",
                        "JE Pre HAI", "DENV-1 Pre HAI","DENV-2 Pre HAI","DENV-3 Pre HAI","DENV-4 Pre HAI",
                        "Max. DENV Pre HAI","Min. DENV Pre HAI","Avg. DENV Pre HAI","Var. DENV Pre HAI",
                        "JE Post HAI", "DENV-1 Post HAI","DENV-2 Post HAI","DENV-3 Post HAI","DENV-4 Post HAI",
                        "Max. DENV Post HAI","Min. DENV Post HAI","Var. DENV Post HAI","Avg. DENV Post HAI",
                        "JE HAI Diff.", "DENV-1 HAI Diff.","DENV-2 HAI Diff.","DENV-3 HAI Diff.","DENV-4 HAI Diff.",
                        "Max. DENV HAI Diff.","Min. DENV HAI Diff.","Avg. DENV HAI Diff.","Var. DENV HAI Diff.",
                        "JE HAI Ratio", "DENV-1 HAI Ratio","DENV-2 HAI Ratio","DENV-3 HAI Ratio","DENV-4 HAI Ratio",
                        "Max. DENV HAI Ratio","Min. DENV HAI Ratio","Avg. DENV HAI Ratio","Var. DENV HAI Ratio",
                        "Time b/w samples", "Pre Sample Date","Post Sample Date", "Age at post sample", "Age at enrollment"),100)


###### define functions for future use ######

#### function that creates a dataframe on which to perform FoI analysis ####
df2prop <- function(df,yearVec,ageVec,serobound){
  df <- df[df$sampDate_post_year %in% yearVec,]
  
  p = rep(0, (length(ageVec)-1))
  n = rep(0, (length(ageVec)-1))
  inf = rep(0, (length(ageVec)-1))
  mid = rep(0, (length(ageVec)-1))
  for(i in 1:(length(ageVec)-1)){
    dfhold <- df[((df$ageAtFollow) >= ageVec[i])&(df$ageAtFollow) < ageVec[i+1],]
    dfhold$seropos <- ifelse((dfhold$d1_post>=serobound)|(dfhold$d2_post>=serobound)|(dfhold$d3_post>=serobound)|(dfhold$d4_post>=serobound),1,0)
    
    p[i] = mean(round(dfhold$seropos)) 
    n[i] = length(dfhold$seropos)
    inf[i] = sum(round(dfhold$seropos))
    mid[i] = mean(c(ageVec[i],ageVec[i+1]))
  }
  newdf <- data.frame (a=ageVec[2:length(ageVec)],
                       p=p,
                       n=n,
                       inf=inf,
                       mid=mid
  )
  return(newdf)
}

####### figure 1 -- descriptive #####
# run constant FOI model
# 2017, run fit on individuals from 1-30
ageVec <- c(1:30)
pvec17 <- df2prop(df,2017,ageVec,20)

fit=glm(cbind(pvec17$inf,pvec17$n-pvec17$inf)~offset(log(pvec17$a)),family=binomial(link="cloglog"), data=pvec17)
phi=exp(coef(summary(fit))[1])
pvec17$fit = 1-exp(-phi*pvec17$a)
pvec17$foi = exp(fit$coef)

ageVec <- c(1:93) # add so we can plot entire age range
pvec17 <- df2prop(df,2017,ageVec,20)
phi=exp(coef(summary(fit))[1])
pvec17$fit = 1-exp(-phi*pvec17$a)
pvec17$foi = exp(fit$coef)

print("FOI fit and CI")
print(exp(coef(summary(fit))[1] + c(-1,0,1)*coef(summary(fit))[2]*1.96))
print(1-exp(-exp(coef(summary(fit))[1] + c(-1,0,1)*coef(summary(fit))[2]*1.96)))

# make actual plot
quantLocs <- as.vector(quantile(df[(df$sampDate_pre_year < 2018)&(df$ageAtEnr>0),]$ageAtEnr,na.rm=T,probs = seq(0, 1, by = 0.1)))
p1 <- ggplot()
p1 <- p1 + geom_line(pvec17,mapping=aes(x=a,y=fit),size=1.2)
p1 <- p1 + xlab("Age (years)") + ylab("Seroprevalence")
p1 <- p1 + coord_cartesian(xlim = c(0,93),ylim=c(0,1))
p1 <- p1 + theme_pubr(base_size = 18,legend="right")
p1 <- p1 + stat_summary_bin(df[(df$sampDate_pre_year< 2018)&(df$ageAtEnr>0),],mapping=aes(x=ageAtEnr,y=as.numeric(!((d1_pre<20)&(d2_pre<20)&(d3_pre<20)&(d4_pre<20)))),geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot")

quantLocs <- as.vector(quantile((df$ageAtEnr),na.rm=T,probs = seq(0, 1, by = 0.1)))
p2 <- ggplot(df[df$ageAtEnr>0,],aes(y=HAI_DEN_ENR_AVG,x=ageAtEnr)) + theme_pubr(base_size = 18,legend = "right")
p2 <- p2 + geom_smooth(data=df[df$ageAtEnr>0,],mapping=aes(y=HAI_DEN_ENR_AVG,x=ageAtEnr)
                       ,method = "glm", formula = y ~ splines::ns(x, 4),col="black")
p2 <- p2 + xlab("Age (years)") + ylab("Avg. DENV (HAI)")
p2 <- p2 + scale_y_continuous(trans=pseudo_log_trans(2),breaks = c(10,20,40,80,160),labels = c("<10","20","40","80","160"))
p2 <- p2 + stat_summary_bin(geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot")
p2 <- p2 + coord_cartesian(ylim=c(10,160))

print((p1 + p2)  + plot_annotation(tag_levels = 'A'))

ggsave(filename="figure1.png",path = figFolder,
       width = 16, height = 10, device='png', dpi=800)

####### figure 2 model predictions #######
## hai pre v post 
preDenv <- df %>% dplyr::select(c("d1_pre","d2_pre","d3_pre","d4_pre")) #df[,c(322:325)]
postDenv <-df %>% dplyr::select(c("d1_post","d2_post","d3_post","d4_post")) 
divDenv <- df %>% dplyr::select(c("div_d1","div_d2","div_d3","div_d4")) 

# find the serotype that has highest pre or post HAIs
maxPre <- unlist(lapply(X = 1:length(df$max_div),FUN = function(x){
  hold <- which(divDenv[x,] == df$max_div[x])
  if(length(hold)>1){
    hold <- hold[which(postDenv[x,hold] == max(postDenv[x,hold]))]
    if(length(hold)>1){
      hold <- hold[1]
    }
  }
  return(preDenv[x,hold])
}))

maxPost <- unlist(lapply(X = 1:length(df$max_div),FUN = function(x){
  hold <- which(divDenv[x,] == df$max_div[x])
  if(length(hold)>1){
    hold <- hold[which(postDenv[x,hold] == max(postDenv[x,hold]))]
    if(length(hold)>1){
      hold <- hold[1]
    }
  }
  return(postDenv[x,hold])
}))

df$max_preSero <- maxPre
df$max_postSero <- maxPost

x <- c(0:5120)
y <- 4*x
dfline<- data.frame(x = x, y=y)

df$trainInf_forplotting <- ifelse(df$trainInf ==0,ifelse(df$prediction_round,1,0),2)
p1 <- ggplot(df) 
p1 <- p1 + geom_smooth(dfline,mapping=aes(x=x,y=y),col="black",lwd=.65)
p1 <- p1 + geom_jitter(data = df[df$trainInf_forplotting==0,],aes(x=max_preSero,y=max_postSero,fill=as.factor(trainInf_forplotting),col=as.factor(trainInf_forplotting)),pch=21,alpha=.5,height=.1,width=.25)
p1 <- p1 + geom_jitter(data = df[df$trainInf_forplotting==1,],aes(x=max_preSero,y=max_postSero,fill=as.factor(trainInf_forplotting),col=as.factor(trainInf_forplotting)),pch=21,alpha=.5,height=.1,width=.25)
p1 <- p1 + geom_jitter(data=df[df$trainInf_forplotting==2,],aes(x=max_preSero,y=max_postSero,fill=as.factor(trainInf_forplotting)),pch=21,height=.1,width=.25)
p1 <- p1 + theme_pubr(base_size = 18,legend = "top") + facet_wrap(~ageCatPreFollow,nrow=1)
p1 <- p1 + scale_x_continuous(name = "Pre interval HAI",trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p1 <- p1 + scale_y_continuous(name = "Post interval HAI", trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p1 <- p1 + scale_fill_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p1 <- p1 + scale_color_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p1 <- p1 + guides(fill=guide_legend("Predicted Infection"),color=FALSE)

# feature importance 
impMatSum <- impMat %>% group_by(Feature) %>% dplyr::summarise(mean = mean(Gain),
                                                               qslow = quantile(Gain,c(.05)),
                                                               qsup = quantile(Gain,c(.95)),
                                                               sd = sd(Gain))
impMatSum <- impMatSum %>% dplyr::arrange(desc(mean))

p3 <- ggplot(impMatSum[1:10,],aes(y=reorder(Feature, mean),x=mean)) + geom_bar(stat="identity", position=position_dodge(),orientation = "y",fill="#009E73") 
p3 <- p3 + geom_errorbar(aes(xmin=qslow,xmax=qsup),col="black",orientation = "y")
p3 <- p3 + xlab("Gain") + ylab("Feature")
p3 <- p3 + theme_pubr(base_size = 18)

# combine into one plot
p <- (p3/(p1) + plot_annotation(tag_levels = list(c('B','C'))))
print(p)
ggsave(filename="figure2.png",path = figFolder,
       width = 20, height = 10, device='png', dpi=700)

####### figure 3 age and # infections by multiple factors #######
bootstrapnum <- 100
houses <- unique(df$houseAnon)

quantLocs <- round(as.vector(quantile((df$ageAtFollow),na.rm=T,probs = seq(0, 1, by = 0.05))))
quantLocs <- c(0,1,5,6,11,16,21,26,31,36,41,51,61,71,105)
quantLocs2 <- c(0,1,6,11,16,21,26,31,36,41,51,105)
df$personTime <- df$timeDiff/365.25
pd = position_dodge(1)

# print total ratio
dfInfRatio= df %>%
  dplyr::group_by(infection,.drop = FALSE) %>%
  dplyr::summarise(time = sum(personTime),
                   n = dplyr::n()) %>%
  mutate(totalTime = sum(time)) %>%
  mutate(inc = n/totalTime)

# create a single df to get sizing information
dfInfTypeYear = df %>%
  dplyr::group_by(sampDate_post_year, infection,.drop = FALSE) %>%
  dplyr::summarise(time = sum(personTime),
                   n = dplyr::n()) %>%
  mutate(totalTime = sum(time)) %>%
  mutate(inc = n/totalTime)

dfInfTypeAge = df %>%
  dplyr::group_by(cut(df$ageAtFollow, breaks = quantLocs2), infection,.drop = FALSE) %>%
  dplyr::summarise(time = sum(personTime),
                   n = dplyr::n()) %>%
  mutate(totalTime = sum(time)) %>%
  mutate(inc = n/totalTime)

dfInfTypeAgePrimary = df[df$prediction_round,] %>%
  dplyr::group_by(cut(df[df$prediction_round,]$ageAtFollow, breaks = quantLocs2),as.factor((d1_pre < 20)&(d2_pre < 20)&(d3_pre < 20)&(d4_pre < 20)),.drop = FALSE) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  mutate(freq = n / sum(n),
         totaln = sum(n))

names(dfInfTypeAge) <- c("ageAtFollow","infection","time","n","totalTime","inc")

names(dfInfTypeAgePrimary) <- c("ageAtFollow","sus","n","freq","totaln")

estMatYear <- matrix(0,bootstrapnum,dim(dfInfTypeYear)[1])
estMatAge <- matrix(0,bootstrapnum,dim(dfInfTypeAge)[1])
estMatAgePrimary <- matrix(0,bootstrapnum,dim(dfInfTypeAgePrimary)[1])


estMatYearRatio <- matrix(0,bootstrapnum,length(unique(dfInfTypeYear$sampDate_post_year)))
estMatAgeRatio <- matrix(0,bootstrapnum,length(unique(dfInfTypeAge$ageAtFollow)))
estMatAllRatio <- matrix(0,bootstrapnum)

dfInfTypeYearRatio <- data.frame(year = unique(dfInfTypeYear$sampDate_post_year),
                                 ratio = rep(0,length(unique(dfInfTypeYear$sampDate_post_year))))
dfInfTypeAgeRatio <- data.frame(ageAtFollow = unique(dfInfTypeAge$ageAtFollow),
                                ratio = rep(0,length(unique(dfInfTypeAge$ageAtFollow))))

dfInfType <- df %>%
  dplyr::group_by(infection,.drop = FALSE) %>%
  dplyr::summarise(time = sum(personTime),
                   n = dplyr::n()) %>%
  mutate(totalTime = sum(time)) %>%
  mutate(inc = n/totalTime)

popRatio <- dfInfType$inc[2]/dfInfType$inc[1]

for(j in 1:bootstrapnum){
  houseSamp <- sample(houses,size = round(length(houses)),replace = T)
  tempdf <- df[df$houseAnon %in% houseSamp,]
  
  tempdf2 <- tempdf %>% dplyr::slice(rep(1:dplyr::n(), n=as.numeric(table(houseSamp)[tempdf$houseAnon])))
  
  dfInfTypeYear = tempdf2 %>%
    dplyr::group_by(sampDate_post_year, infection,.drop = FALSE) %>%
    dplyr::summarise(time = sum(personTime),
                     n = dplyr::n()) %>%
    mutate(totalTime = sum(time)) %>%
    mutate(inc = n/totalTime)
  
  dfInfTypeAll <- tempdf2 %>%
    dplyr::group_by(infection,.drop = FALSE) %>%
    dplyr::summarise(time = sum(personTime),
                     n = dplyr::n()) %>%
    mutate(totalTime = sum(time)) %>%
    mutate(inc = n/totalTime)
  
  dfInfTypeAge = tempdf2 %>%
    dplyr::group_by(cut(tempdf2$ageAtFollow, breaks = quantLocs2), infection,.drop = FALSE) %>%
    dplyr::summarise(time = sum(personTime),
                     n = dplyr::n()) %>%
    mutate(totalTime = sum(time)) %>%
    mutate(inc = n/totalTime)
  
  dfInfTypeAgePrimary = tempdf2[tempdf2$prediction_round,] %>%
    dplyr::group_by(cut(tempdf2[tempdf2$prediction_round,]$ageAtFollow, breaks = quantLocs2),as.factor((d1_pre < 20)&(d2_pre < 20)&(d3_pre < 20)&(d4_pre < 20)),.drop = FALSE) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    mutate(freq = n / sum(n),
           totaln = sum(n))
  
  estMatYear[j,] <- dfInfTypeYear$inc
  estMatAge[j,] <- dfInfTypeAge$inc
  estMatAgePrimary[j,] <- dfInfTypeAgePrimary$freq
  
  
  estMatYearRatio[j,] <- dfInfTypeYear[(dfInfTypeYear$infection == "Asymp. infection"),]$inc/dfInfTypeYear[(dfInfTypeYear$infection == "Symp. infection"),]$inc
  estMatAgeRatio[j,] <- dfInfTypeAge[(dfInfTypeAge$infection == "Asymp. infection"),]$inc/dfInfTypeAge[(dfInfTypeAge$infection == "Symp. infection"),]$inc
  estMatAllRatio[j] <- dfInfTypeAll[(dfInfTypeAll$infection == "Asymp. infection"),]$inc/dfInfTypeAll[(dfInfTypeAll$infection == "Symp. infection"),]$inc
}


estYear <- rep(.1,dim(dfInfTypeYear)[1])
estupYear <- rep(.1,dim(dfInfTypeYear)[1])
estlowYear <- rep(.1,dim(dfInfTypeYear)[1])

estAge <- rep(.1,dim(dfInfTypeAge)[1])
estupAge <- rep(.1,dim(dfInfTypeAge)[1])
estlowAge <- rep(.1,dim(dfInfTypeAge)[1])

estYearRatio <- rep(.1,dim(estMatYearRatio)[2])
estupYearRatio <- rep(.1,dim(estMatYearRatio)[2])
estlowYearRatio <- rep(.1,dim(estMatYearRatio)[2])

estAgeRatio <- rep(.1,dim(estMatAgeRatio)[2])
estupAgeRatio <- rep(.1,dim(estMatAgeRatio)[2])
estlowAgeRatio <- rep(.1,dim(estMatAgeRatio)[2])

estAgePrimary <- rep(.1,dim(dfInfTypeAgePrimary)[1])
estupAgePrimary <- rep(.1,dim(dfInfTypeAgePrimary)[1])
estlowAgePrimary <- rep(.1,dim(dfInfTypeAgePrimary)[1])


for(i in 1:(dim(dfInfTypeYear)[1])){
  estupYear[i] <- quantile(estMatYear[,i],.975,na.rm=T)
  estlowYear[i] <- quantile(estMatYear[,i],.025)
  estYear[i] <- median(estMatYear[,i])
}
for(i in 1:(dim(dfInfTypeAge)[1])){
  estupAge[i] <- quantile(estMatAge[,i],.975,na.rm=T)
  estlowAge[i] <- quantile(estMatAge[,i],.025,na.rm=T)
  estAge[i] <- median(estMatAge[,i],na.rm=T)
}

for(i in 1:(dim(estMatYearRatio)[2])){
  estupYearRatio[i] <- quantile(estMatYearRatio[,i],.975,na.rm=T)
  estlowYearRatio[i] <- quantile(estMatYearRatio[,i],.025)
  estYearRatio[i] <- median(estMatYearRatio[,i])
}

estupYearRatio[is.infinite(estupYearRatio)] <- 200
estlowYearRatio[is.infinite(estlowYearRatio)] <- 200
estYearRatio[is.infinite(estYearRatio)] <- 200

for(i in 1:(dim(estMatAgeRatio)[2])){
  estupAgeRatio[i] <- quantile(estMatAgeRatio[,i],.975,na.rm=T)
  estlowAgeRatio[i] <- quantile(estMatAgeRatio[,i],.025,na.rm=T)
  estAgeRatio[i] <- median(estMatAgeRatio[,i],na.rm=T)
}
estupAgeRatio[is.infinite(estupAgeRatio)] <- 200
estlowAgeRatio[is.infinite(estlowAgeRatio)] <- 200
estAgeRatio[is.infinite(estAgeRatio)] <- 200

for(i in 1:(dim(dfInfTypeAgePrimary)[1])){
  estupAgePrimary[i] <- quantile(estMatAgePrimary[,i],.975,na.rm=T)
  estlowAgePrimary[i] <- quantile(estMatAgePrimary[,i],.025,na.rm=T)
  estAgePrimary[i] <- median(estMatAgePrimary[,i],na.rm=T)
}

estupAllRatio<- quantile(estMatAllRatio,.975,na.rm=T)
estlowAllRatio<- quantile(estMatAllRatio,.025,na.rm=T)
estAllRatio <- median(estMatAllRatio,na.rm=T)


names(dfInfTypeAge) <- c("ageAtFollow","infection","time","n","totalTime","inc")
names(dfInfTypeAgePrimary) <- c("ageAtFollow","sus","n","freq","totaln")

dfInfTypeYear$inc <- estYear
dfInfTypeYear$incU <- estupYear
dfInfTypeYear$incL <- estlowYear

dfInfTypeAge$inc <- estAge
dfInfTypeAge$incU <- estupAge
dfInfTypeAge$incL <- estlowAge

dfInfTypeYearRatio$ratio <- estYearRatio
dfInfTypeYearRatio$ratioU <- estupYearRatio
dfInfTypeYearRatio$ratioL <- estlowYearRatio

dfInfTypeAgeRatio$ratio <- estAgeRatio
dfInfTypeAgeRatio$ratioU <- estupAgeRatio
dfInfTypeAgeRatio$ratioL <- estlowAgeRatio

dfInfTypeAgePrimary$freq <- estAgePrimary
dfInfTypeAgePrimary$freqU <- estupAgePrimary
dfInfTypeAgePrimary$freqL <- estlowAgePrimary

upperBvals <- function(df){
  upperB <- sapply(1:length(df$ageAtFollow),
                   FUN = function(x){
                     hold = str_split(str_split(df$ageAtFollow[x],",")[[1]][2],"]")[[1]][1]
                     return(as.numeric(hold))
                   })
  return(upperB)
}

lowerBvals <- function(df){
  lowerB <- sapply(1:length(df$ageAtFollow),
                   FUN = function(x){
                     hold = str_split(str_split(df$ageAtFollow[x],",")[[1]][1],"\\(")[[1]][2]
                     return(as.numeric(hold))
                   })
  return(lowerB)
}

centralBvals <- function(df){
  centralB <- sapply(1:length(df$ageAtFollow),
                     FUN = function(x){
                       hold1 = str_split(str_split(df$ageAtFollow[x],",")[[1]][2],"]")[[1]][1]
                       hold2 = str_split(str_split(df$ageAtFollow[x],",")[[1]][1],"\\(")[[1]][2]
                       return((as.numeric(hold1) + as.numeric(hold2))/2)
                     })
  return(centralB)
}

dfInfTypeAge$upperB <- upperBvals(dfInfTypeAge)
dfInfTypeAge$lowerB <- lowerBvals(dfInfTypeAge)
dfInfTypeAge$centralB <- centralBvals(dfInfTypeAge)

dfInfTypeAgeRatio$upperB <- upperBvals(dfInfTypeAgeRatio)
dfInfTypeAgeRatio$lowerB <- lowerBvals(dfInfTypeAgeRatio)
dfInfTypeAgeRatio$centralB <- centralBvals(dfInfTypeAgeRatio)

dfInfTypeAgePrimary$upperB <- upperBvals(dfInfTypeAgePrimary)
dfInfTypeAgePrimary$lowerB <- lowerBvals(dfInfTypeAgePrimary)
dfInfTypeAgePrimary$centralB <- centralBvals(dfInfTypeAgePrimary)

p1 <- ggplot(dfInfTypeYear[(dfInfTypeYear$infection != "No infection"),],aes(x=sampDate_post_year-1,y=inc,group=as.factor(infection),col=as.factor(infection)))
p1 <- p1 + geom_point(stat="identity",size = 2,position = pd)
p1 <- p1 + geom_errorbar(aes(ymin=incL,ymax=incU),width=.25,position = pd)
p1 <- p1 + xlab("Year") 
p1 <- p1 + ylab("Incidence")
p1 <- p1 + labs(col="")
p1 <- p1 + theme_pubr(base_size = 18,legend="top",x.text.angle = 0)
p1 <- p1 + scale_color_manual(labels = c("Symp. infection" = "Clinical infection", "Asymp. infection" = "Subclinical infection"),values =rev(c("#00AFBB", "#E7B800", "#FC4E07")))
p1 <- p1 + coord_cartesian(ylim=c(0,.2))

p2 <- ggplot(dfInfTypeAge[(dfInfTypeAge$infection != "No infection"),],aes(x=centralB,y=inc,group=as.factor(infection),col=as.factor(infection)))
p2 <- p2 + geom_point(stat="identity",size = 2,position = pd)
p2 <- p2 + geom_errorbar(aes(ymin=incL,ymax=incU),position = pd)
p2 <- p2 + xlab("Age") 
p2 <- p2 + ylab("Incidence")
p2 <- p2 + labs(col="")
p2 <- p2 + theme_pubr(base_size = 18,legend="top",x.text.angle = 0)
p2 <- p2 + scale_color_manual(labels = c("Symp. infection" = "Clinical infection", "Asymp. infection" = "Subclinical infection"),values =rev(c("#00AFBB", "#E7B800", "#FC4E07")))
p2 <- p2 + coord_cartesian(ylim=c(0,.2))

p3 <- ggplot(dfInfTypeAgePrimary[dfInfTypeAgePrimary$sus==T,],aes(x=centralB,y=freq,group=sus))
p3 <- p3 + geom_point(stat="identity",size = 2)
p3 <- p3 + geom_errorbar(aes(ymin=freqL,ymax=freqU))
p3 <- p3 + xlab("Age") 
p3 <- p3 + ylab("Proportion primary infections")
p3 <- p3 + labs(col="")
p3 <- p3 + theme_pubr(base_size = 18,legend="right",x.text.angle = 0)
p3 <- p3 + guides(fill=guide_legend(reverse=TRUE))

p4 <- ggplot(dfInfTypeYearRatio,
             aes(x=year-1,
                 y=ratio))
p4 <- p4 + geom_rect(aes(xmin=min(year)-2,xmax=max(year)+2,ymin=estlowAllRatio,ymax=estupAllRatio),fill="light grey")
p4 <- p4 + geom_hline(yintercept=estAllRatio,linetype="dotted")
p4 <- p4 + geom_errorbar(aes(ymin=ratioL,
                             ymax=ratioU),width=.25,position = pd)
p4 <- p4 + geom_point(stat="identity",size = 2,fill=ifelse(dfInfTypeYearRatio$ratio == 200,"white","black"),pch=21,col="black")
p4 <- p4 + scale_y_continuous("Subclinical to clinical ratio",
                              breaks = c(0,10,25,50,100,200),labels=expression(0,10,25,50,100, infinity),trans = pseudo_log_trans(10)) 

p4 <- p4 + xlab("Year")
p4 <- p4 + labs(col="")
p4 <- p4 + theme_pubr(base_size = 18,legend="top",x.text.angle = 0)
p4 <- p4 + coord_cartesian(ylim=c(0,210),xlim=c(2014.5,2021.5))


p5 <- ggplot(dfInfTypeAgeRatio,
             aes(x=centralB,
                 y=ratio))
p5 <- p5 + geom_rect(aes(xmin=min(lowerB),xmax=max(upperB),ymin=estlowAllRatio,ymax=estupAllRatio),fill="light grey")
p5 <- p5 + geom_segment(aes(x=min(lowerB),xend=max(upperB),y=estAllRatio,yend=estAllRatio),linetype="dotted")
p5 <- p5 + geom_errorbar(aes(ymin=ratioL,
                             ymax=ratioU),position = pd)
p5 <- p5 + geom_point(stat="identity",size = 2,fill=ifelse(dfInfTypeAgeRatio$ratio == 200,"white","black"),pch=21,col="black")
p5 <- p5 + xlab("Age") 
p5 <- p5 + scale_y_continuous("Subclinical to clinical ratio",
                              breaks = c(0,10,25,50,100,200),labels=expression(0,10,25,50,100, infinity),trans = pseudo_log_trans(10)) 
p5 <- p5 + labs(col="")
p5 <- p5 + theme_pubr(base_size = 18,legend="top",x.text.angle = 0)
p5 <- p5 + coord_cartesian(ylim=c(0,210),xlim=c(min(dfInfTypeAgeRatio$lowerB),max(dfInfTypeAgeRatio$upperB)))

print((p1 + p2+ p3)/( p4 + p5) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'top') )

ggsave(filename="figure3.png",path = figFolder,
       width = 16, height = 8, device='png', dpi=700)

####### figure 4 individual factors #######
## pre titer vs prob
quantLocs <- as.vector(quantile(2^(df$avg_pre),na.rm=T,probs = c(0,seq(.2, 1, by = 0.1))))
p1 <- ggplot(df,aes(x=2^(avg_pre),y=1-exp(-365*(prediction_round/timeDiff)))) + theme_pubr(base_size = 25,legend = "right")
p1 <- p1 + geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = y ~ splines::ns(x, 3),col="black")
p1 <- p1 + xlab("Avg. DENV pre titer (HAI)") + ylab("Annual P(infection)")
p1 <- p1 + coord_cartesian(ylim=c(0, .15),xlim=c(9,2560))
p1 <- p1 + scale_x_continuous(trans=pseudo_log_trans(2),breaks = c(10,80,640,2560),labels = c("<10","80","640","2560"))
p1 <- p1 + stat_summary_bin(geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot")

## pre titer vs sympt
quantLocs <- as.vector(quantile(2^(df$avg_pre),na.rm=T,probs = c(0,seq(.2, 1, by = 0.1))))
p2 <- ggplot(df,aes(x=2^(avg_pre),y=1-exp(-365*(sxInf/timeDiff)))) + theme_pubr(base_size = 25,legend = "right")
p2 <- p2 + geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = y ~ splines::ns(x, 2),col="black")
p2 <- p2 + xlab("Avg. DENV pre titer (HAI)") + ylab("Annual P(symptoms)")
p2 <- p2 + coord_cartesian(ylim=c(0, .025),xlim=c(9,2560))
p2 <- p2 + scale_x_continuous(trans=pseudo_log_trans(2),breaks = c(10,80,640,2560),labels = c("<10","80","640","2560"))
p2 <- p2 + stat_summary_bin(geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot")

## pre titer vs sympt | infection
quantLocs <- as.vector(quantile(2^(df$avg_pre),na.rm=T,probs = c(0,seq(.2, 1, by = 0.1))))
p3 <- ggplot(df[df$prediction > .5,],aes(x=2^(avg_pre),y=sxInf)) + theme_pubr(base_size = 25,legend = "right")
p3 <- p3 + geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = y ~ splines::ns(x, 2),col="black")
p3 <- p3 + xlab("Avg. DENV pre titer (HAI)") + ylab("P(symptoms  |  infection)")
p3 <- p3 + coord_cartesian(ylim=c(0, .2),xlim=c(9,2560))
p3 <- p3 + scale_x_continuous(trans=pseudo_log_trans(2),breaks = c(10,80,640,2560),labels = c("<10","80","640","2560"))
p3 <- p3 + stat_summary_bin(geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot")

p <- (p1 + p2 + p3) + plot_layout(guides = "collect")  + plot_annotation(tag_levels = 'A')
print(p)
ggsave(filename="figure4.png",path = figFolder,
       width = 18, height = 6, device='png', dpi=700)





####### figure S2 figure 2 by serotype #######

## hai pre v post 

x <- c(0:5120)
y <- 4*x
dfline<- data.frame(x = x, y=y)

df$trainInf_forplotting <- ifelse(df$trainInf ==0,ifelse(df$prediction_round,1,0),2)
p1 <- ggplot(df) 
p1 <- p1 + geom_smooth(dfline,mapping=aes(x=x,y=y),col="black",lwd=.65)
p1 <- p1 + geom_jitter(data = df[df$trainInf_forplotting!=2,],aes(x=d1_pre,y=d1_post,fill=as.factor(trainInf_forplotting),col=as.factor(trainInf_forplotting)),pch=21,alpha=.5)
p1 <- p1 + geom_jitter(data=df[df$trainInf_forplotting==2,],aes(x=d1_pre,y=d1_post,fill=as.factor(trainInf_forplotting)),pch=21)
p1 <- p1 + theme_pubr(base_size = 32,legend = "top") + facet_wrap(~ageCatPreFollow,nrow=1)
p1 <- p1 + scale_x_continuous(name = "DENV-1 HAI (pre interval)",trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p1 <- p1 + scale_y_continuous(name = "DENV-1 HAI (post interval)", trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p1 <- p1 + scale_fill_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p1 <- p1 + scale_color_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p1 <- p1 + guides(fill=guide_legend("Predicted Infection"),color="none")

p2 <- ggplot(df) 
p2 <- p2 + geom_smooth(dfline,mapping=aes(x=x,y=y),col="black",lwd=.65)
p2 <- p2 + geom_jitter(data = df[df$trainInf_forplotting!=2,],aes(x=d2_pre,y=d2_post,fill=as.factor(trainInf_forplotting),col=as.factor(trainInf_forplotting)),pch=21,alpha=.5)
p2 <- p2 + geom_jitter(data=df[df$trainInf_forplotting==2,],aes(x=d2_pre,y=d2_post,fill=as.factor(trainInf_forplotting)),pch=21)
p2 <- p2 + theme_pubr(base_size = 32,legend = "top") + facet_wrap(~ageCatPreFollow,nrow=1)
p2 <- p2 + scale_x_continuous(name = "DENV-2 HAI (pre interval)",trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p2 <- p2 + scale_y_continuous(name = "DENV-2 HAI (post interval)", trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p2 <- p2 + scale_fill_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p2 <- p2 + scale_color_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p2 <- p2 + guides(fill=guide_legend("Predicted Infection"),color="none")

p3 <- ggplot(df) 
p3 <- p3 + geom_smooth(dfline,mapping=aes(x=x,y=y),col="black",lwd=.65)
p3 <- p3 + geom_jitter(data = df[df$trainInf_forplotting!=2,],aes(x=d3_pre,y=d3_post,fill=as.factor(trainInf_forplotting),col=as.factor(trainInf_forplotting)),pch=21,alpha=.5)
p3 <- p3 + geom_jitter(data=df[df$trainInf_forplotting==2,],aes(x=d3_pre,y=d3_post,fill=as.factor(trainInf_forplotting)),pch=21)
p3 <- p3 + theme_pubr(base_size = 32,legend = "top") + facet_wrap(~ageCatPreFollow,nrow=1)
p3 <- p3 + scale_x_continuous(name = "DENV-3 HAI (pre interval)",trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p3 <- p3 + scale_y_continuous(name = "DENV-3 HAI (post interval)", trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p3 <- p3 + scale_fill_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p3 <- p3 + scale_color_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p3 <- p3 + guides(fill=guide_legend("Predicted Infection"),color="none")


p4 <- ggplot(df) 
p4 <- p4 + geom_smooth(dfline,mapping=aes(x=x,y=y),col="black",lwd=.65)
p4 <- p4 + geom_jitter(data = df[df$trainInf_forplotting!=2,],aes(x=d4_pre,y=d4_post,fill=as.factor(trainInf_forplotting),col=as.factor(trainInf_forplotting)),pch=21,alpha=.5)
p4 <- p4 + geom_jitter(data=df[df$trainInf_forplotting==2,],aes(x=d4_pre,y=d4_post,fill=as.factor(trainInf_forplotting)),pch=21)
p4 <- p4 + theme_pubr(base_size = 32,legend = "top") + facet_wrap(~ageCatPreFollow,nrow=1)
p4 <- p4 + scale_x_continuous(name = "DENV-4 HAI (pre interval)",trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p4 <- p4 + scale_y_continuous(name = "DENV-4 HAI (post interval)", trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p4 <- p4 + scale_fill_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p4 <- p4 + scale_color_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p4 <- p4 + guides(fill=guide_legend("Predicted Infection"),color="none")


p5 <- ggplot(df) 
p5 <- p5 + geom_smooth(dfline,mapping=aes(x=x,y=y),col="black",lwd=.65)
p5 <- p5 + geom_jitter(data = df[df$trainInf_forplotting!=2,],aes(x=je_pre,y=je_post,fill=as.factor(trainInf_forplotting),col=as.factor(trainInf_forplotting)),pch=21,alpha=.5)
p5 <- p5 + geom_jitter(data=df[df$trainInf_forplotting==2,],aes(x=je_pre,y=je_post,fill=as.factor(trainInf_forplotting)),pch=21)
p5 <- p5 + theme_pubr(base_size = 32,legend = "top") + facet_wrap(~ageCatPreFollow,nrow=1)
p5 <- p5 + scale_x_continuous(name = "JE HAI (pre interval)",trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p5 <- p5 + scale_y_continuous(name = "JE HAI (post interval)", trans=pseudo_log_trans(2),breaks = c(10,80,640,5120),labels = c("<10","80","640","5120"),limits=c(5,5500))
p5 <- p5 + scale_fill_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p5 <- p5 + scale_color_manual(labels = c("1" = "+", "0" = "-","2"="Lab. conf."),values =c("#00AFBB", "#E7B800", "#FC4E07"))
p5 <- p5 + guides(fill=guide_legend("Predicted Infection"),color="none")

print(p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 5,guides = "collect") & theme(legend.position = 'right',
                                                                                legend.direction = 'vertical'))

ggsave(filename="figureS1.png",path = figFolder,
       width = 29, height = 30, device='png', dpi=700)


####### figure S3 individuals with multiple infections #####

df$pos <- 1:length(df$subjectNoAnon)
df$cumsumInf <- NA
for(i in 1:length(unique(df$subjectNoAnon))){
  hold <- df[df$subjectNoAnon == unique(df$subjectNoAnon)[i],]
  hold$cumInf <- cumsum(hold$prediction_round[order(hold$followUp)])[match(hold$followUp,hold$followUp[order(hold$followUp)])]
  df[match(hold$pos,df$pos),]$cumsumInf <- hold$cumInf 
}

# plot how this varies by age
p1 <- ggplot(df[(df$prediction_round),],aes(x=ageAtFollow,fill=(as.factor(cumsumInf)))) 
p1 <- p1 + geom_bar() + scale_x_binned(breaks = seq(5,90,5),limits=c(0,90),name = "Age (years)",
                                       labels = c("","10","","20","","30","","40","","50","","60","","70","","80","","90")) 
p1 <- p1 + theme_pubr(base_size = 24,legend = "top")
p1 <- p1 + labs(fill="Subj. inf. #")
p1 <- p1 + ylab("Number of intervals")
p1 <- p1 + scale_fill_viridis_d(labels = c("TRUE" = "Primary", "FALSE" = "Secondary"))
p1 <- p1 + guides(fill=guide_legend(reverse=TRUE))

# plot how probability varies by age
quantLocs <- as.vector(quantile((df$ageAtFollow),na.rm=T,probs = c(0,seq(.2, 1, by = 0.1))))
p2 <- ggplot(df[df$cumsumInf>0,],aes(x=(ageAtFollow),y=as.numeric(cumsumInf>1))) + theme_pubr(base_size = 24,legend = "top")
p2 <- p2 + geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = y ~ splines::ns(x, 3),col="black")
p2 <- p2 + xlab("Age (years)") + ylab("P(2nd/3rd infection | infection)")#ylab("Annual probability of infection")
#p2 <- p2 + coord_cartesian(ylim=c(0, .2),xlim=c(9,2560))
#p2 <- p2 + scale_x_continuous(breaks=c(0,10,40,160,640,2560),trans = pseudo_log_trans(2))
p2 <- p2 + stat_summary_bin(geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot")


p <- p1 + p2 + plot_annotation(tag_levels = 'A')
print(p)
ggsave(filename="figureS2_mult_infs.png",path = figFolder,
       width = 14, height = 10, device='png', dpi=700)

####### figure S4 sex differences -- foi #######

ageVec <- c(1:60)
dfm <- df[df$gender == "M",]
dff <- df[df$gender == "F",]
pvec17male <- df2prop(dfm,2017,ageVec,40)
pvec17female <- df2prop(dff,2017,ageVec,40)

# run constant FOI model
# 2017
fit=glm(cbind(pvec17male$inf,pvec17male$n-pvec17male$inf)~offset(log(pvec17male$a)),family=binomial(link="cloglog"), data=pvec17male)
phi=exp(coef(fit))
pvec17male$fit = 1-exp(-phi*pvec17male$a)
pvec17male$foi = exp(fit$coef)
pvec17male$gender <- "M"

fit=glm(cbind(pvec17female$inf,pvec17female$n-pvec17female$inf)~offset(log(pvec17female$a)),family=binomial(link="cloglog"), data=pvec17female)
phi=exp(coef(fit))
pvec17female$fit = 1-exp(-phi*pvec17female$a)
pvec17female$foi = exp(fit$coef)
pvec17female$gender <- "F"

pvec17 <- rbind(pvec17female,pvec17male)

# make actual plot
quantLocs <- as.vector(quantile(df[(df$sampDate_pre_year< 2018)&(df$ageAtEnr>1),]$ageAtEnr,na.rm=T,probs = seq(0, 1, by = 0.05)))
pd = position_dodge(1)
p3 <- ggplot()
p3 <- p3 + geom_line(pvec17,mapping=aes(x=a,y=fit,col=gender),size=1.2)
p3 <- p3 + xlab("Age (years)") + ylab("Seroprevalence")
p3 <- p3 + coord_cartesian(xlim = c(0,60),ylim=c(0,1))
p3 <- p3 + theme_pubr(base_size = 18,legend="right") 
p3 <- p3 + stat_summary_bin(df[(df$sampDate_pre_year< 2018)&(df$ageAtEnr>1),],mapping=aes(x=ageAtEnr,y=as.numeric(!((d1_pre<40)&(d2_pre<40)&(d3_pre<40)&(d4_pre<40))),col= gender),geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot",position=pd)
p3 <- p3 + labs("C") +  guides(col=guide_legend("Sex")) 
print(p3)

ggsave(filename="figureS3_age_sex_seroprevalence.png",path = figFolder,
       width = 10, height = 6, device='png', dpi=700)

####### figure S5 sex differences #######
quantLocs <- as.vector(quantile((df$ageAtPreFollow),na.rm=T,probs = seq(0, 1, by = 0.1)))
pd = position_dodge(3)
p1 <- ggplot(df,aes(x=(ageAtPreFollow),y=infection,col=gender,fill=gender)) + theme_pubr(base_size = 25,legend = "right")
p1 <- p1 + geom_smooth(method = "loess")
p1 <- p1 + xlab("Age (years)") + ylab("Annual P(infection)")
p1 <- p1 + coord_cartesian(ylim=c(0, .2))
p1 <- p1 + labs(col = "Sex",fill = "Sex")
p1 <- p1 + stat_summary_bin(geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot",position = pd)
p1 <- p1 + labs("A")

quantLocs <- as.vector(quantile((df$ageAtPreFollow),na.rm=T,probs = seq(0, 1, by = 0.1)))
pd = position_dodge(3)
p2 <- ggplot(df,aes(x=(ageAtPreFollow),y=avg_pre,col=gender,fill=gender)) + theme_pubr(base_size = 25,legend = "right")
p2 <- p2 + geom_smooth(data=df[df$ageAtPreFollow > 0,],mapping=aes(x=(ageAtPreFollow),y=avg_pre,col=gender,fill=gender),
                       method = "loess")
p2 <- p2 + ylab("Avg. DENV pre titer (HAI)") + xlab("Age (years)")
p2 <- p2 + labs(col = "Sex",fill = "Sex")
p2 <- p2 + scale_y_continuous(breaks=c(0,10,20,40,80,160,320,640,1280,2560,5120),trans = pseudo_log_trans(2))
p2 <- p2 + stat_summary_bin(geom = "pointrange",breaks = quantLocs,fun.data = "mean_cl_boot",position = pd)
p2 <- p2  + labs("B")

print(p1 + p2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A'))

ggsave(filename="figureS5.png",path = figFolder,
       width = 16, height = 6, device='png', dpi=700)
# Household Composition Analysis #####


if(max(df$avg_pre,na.rm = T)>100){
  df$avg_pre <- log2(df$avg_pre)
}

if(max(df$houseNoPreTiterMean,na.rm = T)>100){
  df$houseNoPreTiterMean <- log2(df$houseNoPreTiterMean)
}

df$avg_pre_cat <- ifelse(df$avg_pre<log2(20),"<20",
                         ifelse((df$avg_pre>=log2(20))&(df$avg_pre<log2(40)),"20-40",
                                ifelse((df$avg_pre>=log2(40))&(df$avg_pre<log2(80)),"40-80",
                                       ifelse((df$avg_pre>=log2(80))&(df$avg_pre<log2(160)),"80-160",
                                              ifelse((df$avg_pre>=log2(160)),"160+","Odd")))))
df$avg_pre_cat  <- factor(df$avg_pre_cat, levels = c("20-40" ,"<20","40-80", "80-160","160+"))

df$sampDate_post_year <- as.factor(df$sampDate_post_year)
df$sampDate_post_month <- as.factor(df$sampDate_post_month)

df$houseNumSizeupdate <- df$houseNumMalesNewbornupdate + df$houseNumFemaleNewbornupdate + df$houseNumFemales5to18update + df$houseNumMales5to18update +
  df$houseNumFemalesGE18update + df$houseNumMalesGE18update +
  df$houseNumFemalesLT5update + df$houseNumMalesLT5update

#create variable lists for analyses
varlist <- c("houseNumMalesNewbornupdate","houseNumFemaleNewbornupdate","houseNumMalesLT5update",
             "houseNumFemalesLT5update","houseNumMales5to18update","houseNumFemales5to18update",
             "houseNumMalesGE18update","houseNumFemalesGE18update","houseNumNewbornupdate",
             "houseNumLT5update","houseNum5to18update","houseNumGE18update")
varlist1 <- c("houseNumNewbornupdate","houseNumLT5update","houseNum5to18update","houseNumGE18update")
varlist2 <- c("houseNumMalesNewbornupdate","houseNumFemaleNewbornupdate","houseNumMalesLT5update",
              "houseNumFemalesLT5update","houseNumMales5to18update","houseNumFemales5to18update",
              "houseNumMalesGE18update","houseNumFemalesGE18update")

# house size impact on risk
print("House Size")

model <- glmmTMB(prediction_round ~  houseNumSizeupdate + sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                 data=df, 
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = "houseNumSizeupdate"))
print("aOR for house size")
print(cis95)

# all unadjusted analyses ###
modelouts <- c()
{
  model <- glmmTMB(prediction_round ~  houseNumMalesNewbornupdate,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemaleNewbornupdate,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemaleNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesLT5update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesLT5update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMales5to18update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemales5to18update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesGE18update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesGE18update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumNewbornupdate,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumLT5update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNum5to18update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNum5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumGE18update,
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumGE18update"))
  modelouts <- rbind(modelouts,cis95)
}

DV<-varlist # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18",
                "# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A1<-data.frame(DV,IV,ES,LCI,UCI)

A1$IV = factor(A1$IV, levels = c("# Males GE18",
                                 "# Males 5-18",
                                 "# Males 1-5",
                                 "# Males NB",
                                 "# Females GE18",
                                 "# Females 5-18",
                                 "# Females 1-5",
                                 "# Females NB",
                                 "# Indiv. GE18",
                                 "# Indiv. 5-18",
                                 "# Indiv. 1-5",
                                 "# Indiv. NB"))

# all adjusted analyses with household random effects ###

modelouts <- c()
{
  model <- glmmTMB(prediction_round ~  houseNumMalesNewbornupdate + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemaleNewbornupdate + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemaleNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesLT5update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesLT5update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMales5to18update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemales5to18update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesGE18update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesGE18update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumNewbornupdate + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumLT5update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNum5to18update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNum5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumGE18update + (1|houseAnon),
                   data=df, 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumGE18update"))
  modelouts <- rbind(modelouts,cis95)
}

DV<-varlist # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18",
                "# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A1_re<-data.frame(DV,IV,ES,LCI,UCI)

A1_re$IV = factor(A1_re$IV, levels = c("# Males GE18",
                                       "# Males 5-18",
                                       "# Males 1-5",
                                       "# Males NB",
                                       "# Females GE18",
                                       "# Females 5-18",
                                       "# Females 1-5",
                                       "# Females NB",
                                       "# Indiv. GE18",
                                       "# Indiv. 5-18",
                                       "# Indiv. 1-5",
                                       "# Indiv. NB"))

# run for males + females
model <- glmmTMB(prediction_round ~  houseNumNewbornupdate + houseNumLT5update + houseNum5to18update + houseNumGE18update  +  sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                 data=df, 
                 family = binomial(link="logit"))
modelouts <- exp(confint(model,parm = varlist1))

# Data for 1 forest plot ###
DV<-varlist1 # Heading for Facet Wrap
IV<-as.factor(c("# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A2all<-data.frame(DV,IV,ES,LCI,UCI)



# run for males and females
model <- glmmTMB(prediction_round ~ houseNumMalesNewbornupdate + houseNumFemaleNewbornupdate + houseNumMalesLT5update + houseNumFemalesLT5update + houseNumMales5to18update + houseNumFemales5to18update + 
                   houseNumMalesGE18update + houseNumFemalesGE18update +  sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                 data=df, 
                 family = binomial(link="logit"))
modelouts <- exp(confint(model,parm = varlist2))

# Data for 1 forest plot ###
DV<-varlist2 # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A2males_females<-data.frame(DV,IV,ES,LCI,UCI)


## plots 

# all
A1$analysis <- "Univariate"
A1_re$analysis <- "Univariate w/ r.e.s"
A1v1 <- A1
A1v1 <- A1v1[A1v1$DV %in% varlist1,]
A1_re_v1 <- A1_re
A1_re_v1 <- A1_re_v1[A1_re_v1$DV %in% varlist1,]
A2all$analysis <- "Multivariate w/ r.e.s"
A <- rbind(A1v1,A2all)
A <- rbind(A1_re_v1,A)

A$analysis <- factor(A$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))

pd = position_dodge(.5)
p2 <- ggplot(data=A[A$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,3)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("Household Structure") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2,3))+ 
  theme_pubr(base_size = 25,legend = "right")

print("all")
print(A)

print(p2)
ggsave(filename="figureSupp_all.png",path = figFolder,
       width = 16, height = 6, device='png', dpi=700)

# males females
A1$analysis <- "Univariate"
A1_re$analysis <- "Univariate w/ r.e.s"
A1v1 <- A1
A1v1 <- A1v1[A1v1$DV %in% varlist2,]
A1_re_v1 <- A1_re
A1_re_v1 <- A1_re_v1[A1_re_v1$DV %in% varlist2,]
A2males_females$analysis <- "Multivariate w/ r.e.s"
A <- rbind(A1v1,A2males_females)
A <- rbind(A1_re_v1,A)
A$analysis <- factor(A$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))

pd = position_dodge(.5)
p1 <- ggplot(data=A[A$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis)) +
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,3)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("Household Structure") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5,3))+ 
  theme_pubr(base_size = 25,legend = "right")

print("males and females")
print(A)

###### Household HAI and Attack Rate Analysis #####
df$houseFUprePosPropCat <- ifelse(df$houseFUprePosProp==0,0,
                                  ifelse((df$houseFUprePosProp>0)&(df$houseFUprePosProp<.2),1,
                                         ifelse((df$houseFUprePosProp>=.2)&(df$houseFUprePosProp<.4),2,
                                                ifelse((df$houseFUprePosProp>=.4)&(df$houseFUprePosProp<=1),2,NA))))
df$houseFUprePosPropCat <- as.factor(df$houseFUprePosPropCat)

df$houseNoPreTiterMeanCat <- ifelse(df$houseNoPreTiterMean<=log2(40),0,
                                    ifelse((df$houseNoPreTiterMean>log2(40))&(df$houseNoPreTiterMean<=log2(66)),1,
                                           ifelse((df$houseNoPreTiterMean>log2(66))&(df$houseNoPreTiterMean<=log2(100)),2,
                                                  ifelse((df$houseNoPreTiterMean>log2(100))&(df$houseNoPreTiterMean<=log2(152.5)),2,
                                                         ifelse(df$houseNoPreTiterMean>log2(152.5),2,NA)))))
df$houseNoPreTiterMeanCat  <- factor(df$houseNoPreTiterMeanCat, levels = c("0","1","2"))

#  run with full data
model1 <- glmmTMB(prediction_round ~ houseFUprePosPropCat,
                  data=df,
                  family = binomial(link="logit"))

model1_re <- glmmTMB(prediction_round ~ houseFUprePosPropCat + (1|houseAnon),
                     data=df,
                     family = binomial(link="logit"))

model1_mult <- glmmTMB(prediction_round ~ houseFUprePosPropCat + sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df,
                       family = binomial(link="logit"))

model2 <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat,
                  data=df,
                  family = binomial(link="logit"))

model2_re <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat + (1|houseAnon),
                     data=df,
                     family = binomial(link="logit"))

model2_mult <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat + sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df,
                       family = binomial(link="logit"))

model3 <- glmmTMB(prediction_round ~ avg_pre_cat,
                  data=df,
                  family = binomial(link="logit"))

model3_re <- glmmTMB(prediction_round ~ avg_pre_cat + (1|houseAnon),
                     data=df,
                     family = binomial(link="logit"))

model3_mult <- glmmTMB(prediction_round ~ sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df,
                       family = binomial(link="logit"))

# get results
varlist <- c("houseFUprePosPropCat1", "houseFUprePosPropCat2")
modelouts1 <- exp(confint(model1,parm = varlist))

modelouts1_re <- exp(confint(model1_re,parm = varlist))

modelouts1_mult <- exp(confint(model1_mult,parm = varlist))

varlist <- c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2")
modelouts2 <- exp(confint(model2,parm = varlist))

modelouts2_re <- exp(confint(model2_re,parm = varlist))

modelouts2_mult <- exp(confint(model2_mult,parm = varlist))

varlist <- c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+")
modelouts3 <- exp(confint(model3,parm = varlist))

modelouts3_re <- exp(confint(model3_re,parm = varlist))

modelouts3_mult <- exp(confint(model3_mult,parm = varlist))

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2")) # Independent variable names
ES<-as.numeric(modelouts1[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1[,2]) # Upper 95% confidence interval

A1<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1 <- rbind(A1,hold)
A1$analysis <- "Univariate"

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2"))# Independent variable names
ES<-as.numeric(modelouts1_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1_re[,2]) # Upper 95% confidence interval

A1_re<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1_re <- rbind(A1_re,hold)
A1_re$analysis <- "Univariate w/ r.e.s"

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2"))# Independent variable names
ES<-as.numeric(modelouts1_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1_mult[,2]) # Upper 95% confidence interval

A1_multi<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1_multi <- rbind(A1_multi,hold)
A1_multi$analysis <- "Multivariate w/ r.e.s"


A1 <- rbind(A1_re,A1)
A1 <- rbind(A1,A1_multi)
A1$analysis <- factor(A1$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A1$IV <- factor(A1$IV, levels=c("0","0-0.2",">0.2"))

DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2 <- rbind(A2,hold)
A2$analysis <- "Univariate"

DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2_re[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2_re<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2_re <- rbind(A2_re,hold)
A2_re$analysis <- "Univariate w/ r.e.s"


DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2_mult[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2_multi<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2_multi <- rbind(A2_multi,hold)
A2_multi$analysis <- "Multivariate w/ r.e.s"

A2 <- rbind(A2_re,A2)
A2 <- rbind(A2,A2_multi)
A2$analysis <- factor(A2$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A2$IV <- factor(A2$IV, levels=c("<40","40-66",">66"))

DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3 <- rbind(A3,hold)
A3$analysis <- "Univariate"

DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3_re[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3_re<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3_re <- rbind(A3_re,hold)
A3_re$analysis <- "Univariate w/ r.e.s"


DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3_mult[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3_multi<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3_multi <- rbind(A3_multi,hold)
A3_multi$analysis <- "Multivariate w/ r.e.s"

A3 <- rbind(A3_re,A3)
A3 <- rbind(A3,A3_multi)
A3$analysis <- factor(A3$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A3$IV <- factor(A3$IV, levels=c("<20","20-40","40-80","80-160",">160"))


pd = position_dodge(.5)
p3 <- ggplot(data=A1[A1$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,2)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("AR of household (previous year)") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2))+ 
  theme_pubr(base_size = 25,legend = "right")

print(p3)
ggsave(filename="household_attack_analysis.png",path = figFolder,
       width = 16, height = 10, device='png', dpi=700)

print("HH AR")
print(A1)


pd = position_dodge(.5)
p4 <- ggplot(data=A2[A2$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,2)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("Avg. household HAI titers (pre)") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2))+ 
  theme_pubr(base_size = 25,legend = "right")

print(p4)
ggsave(filename="household_hi_analysis.png",path = figFolder,
       width = 16, height = 10, device='png', dpi=700)

print("HH HI")
print(A2)

print(p2 + p1 + p3 + p4 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A'))
ggsave(filename="figure_5_allHH.png",path = figFolder,
       width = 20, height = 14, device='png', dpi=700)

# Household Composition Analysis sensitivity analyses (80% of household sampled) ######
varlist <- c("houseNumMalesNewbornupdate","houseNumFemaleNewbornupdate","houseNumMalesLT5update",
             "houseNumFemalesLT5update","houseNumMales5to18update","houseNumFemales5to18update",
             "houseNumMalesGE18update","houseNumFemalesGE18update","houseNumNewbornupdate",
             "houseNumLT5update","houseNum5to18update","houseNumGE18update")
varlist1 <- c("houseNumNewbornupdate","houseNumLT5update","houseNum5to18update","houseNumGE18update")
varlist2 <- c("houseNumMalesNewbornupdate","houseNumFemaleNewbornupdate","houseNumMalesLT5update",
              "houseNumFemalesLT5update","houseNumMales5to18update","houseNumFemales5to18update",
              "houseNumMalesGE18update","houseNumFemalesGE18update")
# all unadjusted analyses ###
modelouts <- c()
{
  model <- glmmTMB(prediction_round ~  houseNumMalesNewbornupdate,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemaleNewbornupdate,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemaleNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesLT5update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesLT5update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMales5to18update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemales5to18update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesGE18update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesGE18update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumNewbornupdate,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumLT5update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNum5to18update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNum5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumGE18update,
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumGE18update"))
  modelouts <- rbind(modelouts,cis95)
}

DV<-varlist # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18",
                "# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18"))  # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A1<-data.frame(DV,IV,ES,LCI,UCI)

A1$IV = factor(A1$IV, levels = c("# Males GE18",
                                 "# Males 5-18",
                                 "# Males 1-5",
                                 "# Females GE18",
                                 "# Females 5-18",
                                 "# Females 1-5",
                                 "# Indiv. GE18",
                                 "# Indiv. 5-18",
                                 "# Indiv. 1-5"))

# all adjusted analyses with household random effects ###

modelouts <- c()
{
  model <- glmmTMB(prediction_round ~  houseNumMalesNewbornupdate + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemaleNewbornupdate + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemaleNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesLT5update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesLT5update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMales5to18update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemales5to18update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesGE18update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesGE18update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumNewbornupdate + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,],
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumLT5update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNum5to18update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNum5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumGE18update + (1|houseAnon),
                   data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumGE18update"))
  modelouts <- rbind(modelouts,cis95)
}

DV<-varlist # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18",
                "# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A1_re<-data.frame(DV,IV,ES,LCI,UCI)

A1_re$IV = factor(A1_re$IV, levels = c("# Males GE18",
                                       "# Males 5-18",
                                       "# Males 1-5",
                                       "# Males NB",
                                       "# Females GE18",
                                       "# Females 5-18",
                                       "# Females 1-5",
                                       "# Females NB",
                                       "# Indiv. GE18",
                                       "# Indiv. 5-18",
                                       "# Indiv. 1-5",
                                       "# Indiv. NB"))

# run for males + females
model <- glmmTMB(prediction_round ~  houseNumNewbornupdate + houseNumLT5update + houseNum5to18update + houseNumGE18update  +  sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                 data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                 family = binomial(link="logit"))
modelouts <- exp(confint(model,parm = varlist1))

# Data for 1 forest plot ###
DV<-varlist1 # Heading for Facet Wrap
IV<-as.factor(c("# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18"))  # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A2all<-data.frame(DV,IV,ES,LCI,UCI)



# run for males and females
model <- glmmTMB(prediction_round ~  houseNumMalesNewbornupdate + houseNumFemaleNewbornupdate +houseNumMalesLT5update + houseNumFemalesLT5update + houseNumMales5to18update + houseNumFemales5to18update + 
                   houseNumMalesGE18update + houseNumFemalesGE18update +  sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                 data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                 family = binomial(link="logit"))
modelouts <- exp(confint(model,parm = varlist2))

# Data for 1 forest plot ###
DV<-varlist2 # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A2males_females<-data.frame(DV,IV,ES,LCI,UCI)


## plots 

# all
A1$analysis <- "Univariate"
A1_re$analysis <- "Univariate w/ r.e.s"
A1v1 <- A1
A1v1 <- A1v1[A1v1$DV %in% varlist1,]
A1_re_v1 <- A1_re
A1_re_v1 <- A1_re_v1[A1_re_v1$DV %in% varlist1,]
A2all$analysis <- "Multivariate w/ r.e.s"
A <- rbind(A1v1,A2all)
A <- rbind(A1_re_v1,A)

A$analysis <- factor(A$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))

pd = position_dodge(.5)
p2 <- ggplot(data=A[A$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) +  
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,3)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("Household Structure") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2,3))+ 
  theme_pubr(base_size = 25,legend = "right")

print("all subset 80pc")
print(A)

# males females
A1$analysis <- "Univariate"
A1_re$analysis <- "Univariate w/ r.e.s"
A1v1 <- A1
A1v1 <- A1v1[A1v1$DV %in% varlist2,]
A1_re_v1 <- A1_re
A1_re_v1 <- A1_re_v1[A1_re_v1$DV %in% varlist2,]
A2males_females$analysis <- "Multivariate w/ r.e.s"
A <- rbind(A1v1,A2males_females)
A <- rbind(A1_re_v1,A)
A$analysis <- factor(A$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))

pd = position_dodge(.5)
p1 <- ggplot(data=A[A$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis)) +
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+
  coord_flip(ylim=c(0,3)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+
  xlab("Household Structure") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5,3))+ 
  theme_pubr(base_size = 25,legend = "right")

print("males and females subset 80pc")
print(A)


###### Household HAI and Attack Rate Analysis (80% of household sampled)  #####
#  run with full data
model1 <- glmmTMB(prediction_round ~ houseFUprePosPropCat,
                  data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                  family = binomial(link="logit"))

model1_re <- glmmTMB(prediction_round ~ houseFUprePosPropCat + (1|houseAnon),
                     data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                     family = binomial(link="logit"))

model1_mult <- glmmTMB(prediction_round ~ houseFUprePosPropCat + sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                       family = binomial(link="logit"))

model2 <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat,
                  data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                  family = binomial(link="logit"))

model2_re <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat + (1|houseAnon),
                     data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                     family = binomial(link="logit"))

model2_mult <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat + sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                       family = binomial(link="logit"))

model3 <- glmmTMB(prediction_round ~ avg_pre_cat,
                  data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                  family = binomial(link="logit"))

model3_re <- glmmTMB(prediction_round ~ avg_pre_cat + (1|houseAnon),
                     data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                     family = binomial(link="logit"))

model3_mult <- glmmTMB(prediction_round ~ sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df[((df$numSamplesPost)/(df$houseNumSizeupdate+1))>=.8,], 
                       family = binomial(link="logit"))

# get results
varlist <- c("houseFUprePosPropCat1", "houseFUprePosPropCat2")
modelouts1 <- exp(confint(model1,parm = varlist))

modelouts1_re <- exp(confint(model1_re,parm = varlist))

modelouts1_mult <- exp(confint(model1_mult,parm = varlist))

varlist <- c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2")
modelouts2 <- exp(confint(model2,parm = varlist))

modelouts2_re <- exp(confint(model2_re,parm = varlist))

modelouts2_mult <- exp(confint(model2_mult,parm = varlist))

varlist <- c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+")
modelouts3 <- exp(confint(model3,parm = varlist))

modelouts3_re <- exp(confint(model3_re,parm = varlist))

modelouts3_mult <- exp(confint(model3_mult,parm = varlist))

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2")) # Independent variable names
ES<-as.numeric(modelouts1[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1[,2]) # Upper 95% confidence interval

A1<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1 <- rbind(A1,hold)
A1$analysis <- "Univariate"

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2"))# Independent variable names
ES<-as.numeric(modelouts1_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1_re[,2]) # Upper 95% confidence interval

A1_re<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1_re <- rbind(A1_re,hold)
A1_re$analysis <- "Univariate w/ r.e.s"

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2"))# Independent variable names
ES<-as.numeric(modelouts1_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1_mult[,2]) # Upper 95% confidence interval

A1_multi<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1_multi <- rbind(A1_multi,hold)
A1_multi$analysis <- "Multivariate w/ r.e.s"


A1 <- rbind(A1_re,A1)
A1 <- rbind(A1,A1_multi)
A1$analysis <- factor(A1$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A1$IV <- factor(A1$IV, levels=c("0","0-0.2",">0.2"))

DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2 <- rbind(A2,hold)
A2$analysis <- "Univariate"

DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2_re[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2_re<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2_re <- rbind(A2_re,hold)
A2_re$analysis <- "Univariate w/ r.e.s"


DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2_mult[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2_multi<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2_multi <- rbind(A2_multi,hold)
A2_multi$analysis <- "Multivariate w/ r.e.s"

A2 <- rbind(A2_re,A2)
A2 <- rbind(A2,A2_multi)
A2$analysis <- factor(A2$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A2$IV <- factor(A2$IV, levels=c("<40","40-66",">66"))

DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3 <- rbind(A3,hold)
A3$analysis <- "Univariate"

DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3_re[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3_re<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3_re <- rbind(A3_re,hold)
A3_re$analysis <- "Univariate w/ r.e.s"


DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3_mult[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3_multi<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3_multi <- rbind(A3_multi,hold)
A3_multi$analysis <- "Multivariate w/ r.e.s"

A3 <- rbind(A3_re,A3)
A3 <- rbind(A3,A3_multi)
A3$analysis <- factor(A3$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A3$IV <- factor(A3$IV, levels=c("<20","20-40","40-80","80-160",">160"))


pd = position_dodge(.5)
p3 <- ggplot(data=A1[A1$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,2)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("AR of household (previous year)") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2))+ 
  theme_pubr(base_size = 25,legend = "right")

print("HH AR 80% sensitivity")
print(A1)


pd = position_dodge(.5)
p4 <- ggplot(data=A2[A2$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,2)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("Avg. household HAI titers (pre)") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2))+ 
  theme_pubr(base_size = 25,legend = "right")

print(p4)
ggsave(filename="household_hi_analysis_80percHouseData.png",path = figFolder,
       width = 16, height = 10, device='png', dpi=700)

print("HH HI 80% sensitivity")
print(A2)

print(p2 + p1 + p3 + p4 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A'))
ggsave(filename="figure_S9_allHH_80percHouseData.png",path = figFolder,
       width = 20, height = 14, device='png', dpi=700)



# Household Composition Analysis sensitivity analyses (seronaive only) #####
#create variable lists for analyses
varlist <- c("houseNumMalesNewbornupdate","houseNumFemaleNewbornupdate","houseNumMalesLT5update",
             "houseNumFemalesLT5update","houseNumMales5to18update","houseNumFemales5to18update",
             "houseNumMalesGE18update","houseNumFemalesGE18update","houseNumNewbornupdate",
             "houseNumLT5update","houseNum5to18update","houseNumGE18update")
varlist1 <- c("houseNumNewbornupdate","houseNumLT5update","houseNum5to18update","houseNumGE18update")
varlist2 <- c("houseNumMalesNewbornupdate","houseNumFemaleNewbornupdate","houseNumMalesLT5update",
              "houseNumFemalesLT5update","houseNumMales5to18update","houseNumFemales5to18update",
              "houseNumMalesGE18update","houseNumFemalesGE18update")

# all unadjusted analyses ###
modelouts <- c()
{
  model <- glmmTMB(prediction_round ~  houseNumMalesNewbornupdate,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemaleNewbornupdate,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemaleNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesLT5update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesLT5update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMales5to18update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemales5to18update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesGE18update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesGE18update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumNewbornupdate,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumLT5update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNum5to18update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNum5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumGE18update,
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumGE18update"))
  modelouts <- rbind(modelouts,cis95)
}

DV<-varlist # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18",
                "# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A1<-data.frame(DV,IV,ES,LCI,UCI)

A1$IV = factor(A1$IV, levels = c("# Males GE18",
                                 "# Males 5-18",
                                 "# Males 1-5",
                                 "# Males NB",
                                 "# Females GE18",
                                 "# Females 5-18",
                                 "# Females 1-5",
                                 "# Females NB",
                                 "# Indiv. GE18",
                                 "# Indiv. 5-18",
                                 "# Indiv. 1-5",
                                 "# Indiv. NB"))

# all adjusted analyses with household random effects ###

modelouts <- c()
{
  model <- glmmTMB(prediction_round ~  houseNumMalesNewbornupdate + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemaleNewbornupdate + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemaleNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesLT5update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesLT5update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMales5to18update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemales5to18update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemales5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumMalesGE18update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumMalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumFemalesGE18update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumFemalesGE18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumNewbornupdate + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumNewbornupdate"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumLT5update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumLT5update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNum5to18update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNum5to18update"))
  modelouts <- rbind(modelouts,cis95)
  
  model <- glmmTMB(prediction_round ~  houseNumGE18update + (1|houseAnon),
                   data=df[df$max_pre < 40,], 
                   family = binomial(link="logit"))
  cis95 <- exp(confint(model,parm = "houseNumGE18update"))
  modelouts <- rbind(modelouts,cis95)
}

DV<-varlist # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18",
                "# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A1_re<-data.frame(DV,IV,ES,LCI,UCI)

A1_re$IV = factor(A1_re$IV, levels = c("# Males GE18",
                                       "# Males 5-18",
                                       "# Males 1-5",
                                       "# Males NB",
                                       "# Females GE18",
                                       "# Females 5-18",
                                       "# Females 1-5",
                                       "# Females NB",
                                       "# Indiv. GE18",
                                       "# Indiv. 5-18",
                                       "# Indiv. 1-5",
                                       "# Indiv. NB"))

# run for males + females
model <- glmmTMB(prediction_round ~  houseNumNewbornupdate + houseNumLT5update + houseNum5to18update + houseNumGE18update  +  sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                 data=df[df$max_pre < 40,], 
                 family = binomial(link="logit"))
modelouts <- exp(confint(model,parm = varlist1))

# Data for 1 forest plot ###
DV<-varlist1 # Heading for Facet Wrap
IV<-as.factor(c("# Indiv. NB",
                "# Indiv. 1-5",
                "# Indiv. 5-18",
                "# Indiv. GE18"))  # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A2all<-data.frame(DV,IV,ES,LCI,UCI)



# run for males and females
model <- glmmTMB(prediction_round ~  houseNumMalesNewbornupdate + houseNumFemaleNewbornupdate + houseNumMalesLT5update + houseNumFemalesLT5update + houseNumMales5to18update + houseNumFemales5to18update + 
                   houseNumMalesGE18update + houseNumFemalesGE18update +  sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                 data=df[df$max_pre < 40,], 
                 family = binomial(link="logit"))
modelouts <- exp(confint(model,parm = varlist2))

# Data for 1 forest plot ###
DV<-varlist2 # Heading for Facet Wrap
IV<-as.factor(c("# Males NB",
                "# Females NB",
                "# Males 1-5",
                "# Females 1-5",
                "# Males 5-18",
                "# Females 5-18",
                "# Males GE18",
                "# Females GE18")) # Independent variable names
ES<-as.numeric(modelouts[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts[,2]) # Upper 95% confidence interval

A2males_females<-data.frame(DV,IV,ES,LCI,UCI)


## plots 

# all
A1$analysis <- "Univariate"
A1_re$analysis <- "Univariate w/ r.e.s"
A1v1 <- A1
A1v1 <- A1v1[A1v1$DV %in% varlist1,]
A1_re_v1 <- A1_re
A1_re_v1 <- A1_re_v1[A1_re_v1$DV %in% varlist1,]
A2all$analysis <- "Multivariate w/ r.e.s"
A <- rbind(A1v1,A2all)
A <- rbind(A1_re_v1,A)

A$analysis <- factor(A$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))

pd = position_dodge(.5)
p2 <- ggplot(data=A[A$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,3)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("Household Structure") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2,3))+ 
  theme_pubr(base_size = 25,legend = "right")

print("all subset seronaive")
print(A)

# males females
A1$analysis <- "Univariate"
A1_re$analysis <- "Univariate w/ r.e.s"
A1v1 <- A1
A1v1 <- A1v1[A1v1$DV %in% varlist2,]
A1_re_v1 <- A1_re
A1_re_v1 <- A1_re_v1[A1_re_v1$DV %in% varlist2,]
A2males_females$analysis <- "Multivariate w/ r.e.s"
A <- rbind(A1v1,A2males_females)
A <- rbind(A1_re_v1,A)
A$analysis <- factor(A$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))

pd = position_dodge(.5)
p1 <- ggplot(data=A[A$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis)) +
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,3)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("Household Structure") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5,3))+ 
  theme_pubr(base_size = 25,legend = "right")

print("males and females subset seronaive")
print(A)


###### Household HAI and Attack Rate Analysis (seronaive only) #####
#  run with full data
model1 <- glmmTMB(prediction_round ~ houseFUprePosPropCat,
                  data=df[df$max_pre < 40,],
                  family = binomial(link="logit"))

model1_re <- glmmTMB(prediction_round ~ houseFUprePosPropCat + (1|houseAnon),
                     data=df[df$max_pre < 40,],
                     family = binomial(link="logit"))

model1_mult <- glmmTMB(prediction_round ~ houseFUprePosPropCat + sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df[df$max_pre < 40,],
                       family = binomial(link="logit"))

model2 <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat,
                  data=df[df$max_pre < 40,],
                  family = binomial(link="logit"))

model2_re <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat + (1|houseAnon),
                     data=df[df$max_pre < 40,],
                     family = binomial(link="logit"))

model2_mult <- glmmTMB(prediction_round ~ houseNoPreTiterMeanCat + sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df[df$max_pre < 40,],
                       family = binomial(link="logit"))

model3 <- glmmTMB(prediction_round ~ avg_pre_cat,
                  data=df[df$max_pre < 40,],
                  family = binomial(link="logit"))

model3_re <- glmmTMB(prediction_round ~ avg_pre_cat + (1|houseAnon),
                     data=df[df$max_pre < 40,],
                     family = binomial(link="logit"))

model3_mult <- glmmTMB(prediction_round ~ sampDate_post_month + sampDate_post_year + avg_pre_cat + (1|houseAnon),
                       data=df[df$max_pre < 40,],
                       family = binomial(link="logit"))

# get results
varlist <- c("houseFUprePosPropCat1", "houseFUprePosPropCat2")
modelouts1 <- exp(confint(model1,parm = varlist))

modelouts1_re <- exp(confint(model1_re,parm = varlist))

modelouts1_mult <- exp(confint(model1_mult,parm = varlist))

varlist <- c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2")
modelouts2 <- exp(confint(model2,parm = varlist))

modelouts2_re <- exp(confint(model2_re,parm = varlist))

modelouts2_mult <- exp(confint(model2_mult,parm = varlist))

varlist <- c("avg_pre_cat<20")
modelouts3 <- exp(confint(model3,parm = varlist))

modelouts3_re <- exp(confint(model3_re,parm = varlist))

modelouts3_mult <- exp(confint(model3_mult,parm = varlist))

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2")) # Independent variable names
ES<-as.numeric(modelouts1[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1[,2]) # Upper 95% confidence interval

A1<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1 <- rbind(A1,hold)
A1$analysis <- "Univariate"

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2"))# Independent variable names
ES<-as.numeric(modelouts1_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1_re[,2]) # Upper 95% confidence interval

A1_re<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1_re <- rbind(A1_re,hold)
A1_re$analysis <- "Univariate w/ r.e.s"

DV<-c("houseFUprePosPropCat1", "houseFUprePosPropCat2") # Heading for Facet Wrap
IV<-as.factor(c("0-0.2",
                ">0.2"))# Independent variable names
ES<-as.numeric(modelouts1_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts1_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts1_mult[,2]) # Upper 95% confidence interval

A1_multi<-data.frame(DV,IV,ES,LCI,UCI)
hold<-data.frame("houseFUprePosPropCat0",as.factor("0"),1,NA,NA)
names(hold)<-c("DV","IV","ES","LCI","UCI")
A1_multi <- rbind(A1_multi,hold)
A1_multi$analysis <- "Multivariate w/ r.e.s"


A1 <- rbind(A1_re,A1)
A1 <- rbind(A1,A1_multi)
A1$analysis <- factor(A1$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A1$IV <- factor(A1$IV, levels=c("0","0-0.2",">0.2"))

DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2 <- rbind(A2,hold)
A2$analysis <- "Univariate"

DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2_re[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2_re<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2_re <- rbind(A2_re,hold)
A2_re$analysis <- "Univariate w/ r.e.s"


DV<-c("houseNoPreTiterMeanCat1", "houseNoPreTiterMeanCat2") # Heading for Facet Wrap
IV<-as.factor(c("40-66",
                ">66")) # Independent variable names
ES<-as.numeric(modelouts2_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts2_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts2_mult[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A2_multi<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("houseNoPreTiterMeanCat0",as.factor("<40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A2_multi <- rbind(A2_multi,hold)
A2_multi$analysis <- "Multivariate w/ r.e.s"

A2 <- rbind(A2_re,A2)
A2 <- rbind(A2,A2_multi)
A2$analysis <- factor(A2$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A2$IV <- factor(A2$IV, levels=c("<40","40-66",">66"))

DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3 <- rbind(A3,hold)
A3$analysis <- "Univariate"

DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3_re[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3_re[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3_re[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3_re<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3_re <- rbind(A3_re,hold)
A3_re$analysis <- "Univariate w/ r.e.s"


DV<-c("avg_pre_cat<20", "avg_pre_cat40-80","avg_pre_cat80-160","avg_pre_cat160+") # Heading for Facet Wrap
IV<-as.factor(c("<20",
                "40-80",
                "80-160",
                ">160")) # Independent variable names
ES<-as.numeric(modelouts3_mult[,3]) # b Estimate (could be standardized estimate, Odds Ratio, Incident Rate Ratio, etc.)
LCI<-as.numeric(modelouts3_mult[,1]) # Lower 95% confidence interval
UCI<-as.numeric(modelouts3_mult[,2]) # Upper 95% confidence interval
PCHval <- rep(16,length(IV))

A3_multi<-data.frame(DV,IV,ES,LCI,UCI,PCHval)
hold<-data.frame("avg_pre_cat20-40",as.factor("20-40"),1,NA,NA,21)
names(hold)<-c("DV","IV","ES","LCI","UCI","PCHval")
A3_multi <- rbind(A3_multi,hold)
A3_multi$analysis <- "Multivariate w/ r.e.s"

A3 <- rbind(A3_re,A3)
A3 <- rbind(A3,A3_multi)
A3$analysis <- factor(A3$analysis, levels=c('Univariate', 'Univariate w/ r.e.s', 'Multivariate w/ r.e.s'))
A3$IV <- factor(A3$IV, levels=c("<20","20-40","40-80","80-160",">160"))


pd = position_dodge(.5)
p3 <- ggplot(data=A1[A1$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,2)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("AR of household (previous year)") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2))+ 
  theme_pubr(base_size = 25,legend = "right")

print("HH AR seronaive sensitivity")
print(A1)


pd = position_dodge(.5)
p4 <- ggplot(data=A2[A2$analysis == 'Multivariate w/ r.e.s',], aes(x=IV, y=ES, ymin=LCI, ymax=UCI,group=analysis))+
  geom_pointrange(position = pd)+ 
  geom_hline(yintercept=1, lty=2, size =1) + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, cex=1,position = pd)+ 
  coord_flip(ylim=c(0,2)) + 
  geom_point(shape = 15, size = 2,position = pd) + 
  ggtitle("")+ 
  xlab("Avg. household HAI titers (pre)") + 
  ylab("OR (95% CI)") + 
  scale_y_continuous(breaks = c(0,.5,1,1.5,2))+ 
  theme_pubr(base_size = 25,legend = "right")

print(p4)
ggsave(filename="household_hi_analysis_seronaive.png",path = figFolder,
       width = 16, height = 10, device='png', dpi=700)

print("HH HI seronaive sensitivity")
print(A2)

print(p2 + p1 + p3 + p4 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A'))
ggsave(filename="figure_S10_allHH_seronaive.png",path = figFolder,
       width = 20, height = 14, device='png', dpi=700)



###### Individual, household, and temporal factor analysis #####
namevals <- character(0)
estsU <- numeric(0)
estsM <- numeric(0)

modelouts <- c()

# water containers #####
hold <- df

model <- glmmTMB(prediction_round ~ nbWaterCont,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = "nbWaterCont"))
modelouts <- rbind(modelouts,cis95)


# year #####
hold$sampDate_post_year <- as.factor(hold$sampDate_post_year)

model <- glmmTMB(prediction_round ~ sampDate_post_year,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = c("sampDate_post_year2017","sampDate_post_year2018","sampDate_post_year2019","sampDate_post_year2020","sampDate_post_year2021","sampDate_post_year2022")))
modelouts <- rbind(modelouts,cis95)

# month #####
hold$sampDate_post_month <- factor(hold$sampDate_post_month,levels = c("5","1","2","3","4","6","7","8","9","10","11"))
model <- glmmTMB(prediction_round ~ sampDate_post_month,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = c("sampDate_post_month1","sampDate_post_month2","sampDate_post_month3","sampDate_post_month4","sampDate_post_month6",
                                    "sampDate_post_month7","sampDate_post_month8","sampDate_post_month9","sampDate_post_month10","sampDate_post_month11")))
modelouts <- rbind(modelouts,cis95)

# individual pre titer #####
model <- glmmTMB(prediction_round ~ avg_pre,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = "avg_pre"))
modelouts <- rbind(modelouts,cis95)



# individual pre titer Categorical #####
model <- glmmTMB(prediction_round ~ avgHAICatPreFollow,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = c("avgHAICatPreFollow20-40","avgHAICatPreFollow40-80","avgHAICatPreFollow80-160","avgHAICatPreFollow160+")))
modelouts <- rbind(modelouts,cis95)

# water containers plastic #####
model <- glmmTMB(prediction_round ~ nbWaterContPlastic,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = "nbWaterContPlastic"))
modelouts <- rbind(modelouts,cis95)

# water plastic bottles #####
model <- glmmTMB(prediction_round ~ waterContainerSoftDrink,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = "waterContainerSoftDrink"))
modelouts <- rbind(modelouts,cis95)

# house type (removing houses labeled "other") #####
model <- glmmTMB(prediction_round ~ houseType,
                 data=df[df$houseType != "Others",],
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = c("houseTypeSingle","houseTypeTownhouse")))
modelouts <- rbind(modelouts,cis95)


# garbage management #####

model <- glmmTMB(prediction_round ~ as.factor(garbageManagement=="Car"),
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = "as.factor(garbageManagement == \"Car\")TRUE"))
modelouts <- rbind(modelouts,cis95)

# concrete house #####
model <- glmmTMB(prediction_round ~ constructionMaterialsConcrete,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = "constructionMaterialsConcreteNULL"))
modelouts <- rbind(modelouts,cis95)

# zinc roof #####
model <- glmmTMB(prediction_round ~ roofTypeZinc,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm = "roofTypeZincNULL"))
modelouts <- rbind(modelouts,cis95)



# nearby source of water #####
model <- glmmTMB(prediction_round ~ (waterResourcesNearby == 1),
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm="waterResourcesNearby == 1TRUE"))
modelouts <- rbind(modelouts,cis95)

# water supply by pipe #####
model <- glmmTMB(prediction_round ~ (waterTypePipe==1),
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm="waterTypePipe == 1TRUE"))
modelouts <- rbind(modelouts,cis95)



# door screens #####
model <- glmmTMB(prediction_round ~ doorScreens,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm="doorScreens"))
modelouts <- rbind(modelouts,cis95)

# houseNoPreTiterMean #####
model <- glmmTMB(prediction_round ~ houseNoPreTiterMean,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm="houseNoPreTiterMean"))
modelouts <- rbind(modelouts,cis95)


# toilets outside #####
model <- glmmTMB(prediction_round ~ numToiletsOutside,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm="numToiletsOutside"))
modelouts <- rbind(modelouts,cis95)


# occupation  #####
model <- glmmTMB(prediction_round ~ intervalOccupationTop4,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm=c("intervalOccupationTop4Farmer","intervalOccupationTop4Student","intervalOccupationTop4Unemployed")))
modelouts <- rbind(modelouts,cis95)


# gender #####
model <- glmmTMB(prediction_round ~ gender,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm="genderM"))
modelouts <- rbind(modelouts,cis95)

# ageCatPreFollow #####
model <- glmmTMB(prediction_round ~ ageCatPreFollow,
                 data=hold,
                 family = binomial(link="logit"))
cis95 <- exp(confint(model,parm=c("ageCatPreFollow5-18","ageCatPreFollow18-30","ageCatPreFollow30-50","ageCatPreFollow50+")))
modelouts <- rbind(modelouts,cis95)

print("univariate")
print(modelouts)



# multivariate analysis household REs #########

modeloutsMultHouse <- c()

hold$sampDate_post_month <- factor(hold$sampDate_post_month,levels = c("5","1","2","3","4","6","7","8","9","10","11"))
hold$sampDate_post_year <- as.factor(hold$sampDate_post_year)
modelM <- glmmTMB(prediction_round ~  sampDate_post_month + sampDate_post_year + ageCatPreFollow + intervalOccupationTop4  + (1|houseAnon),
                  data=hold,
                  family = binomial(link="logit"))

cis95 <- exp(confint(modelM,parm=c("sampDate_post_year2017","sampDate_post_year2018","sampDate_post_year2019","sampDate_post_year2020","sampDate_post_year2021","sampDate_post_year2022",
                                   "sampDate_post_month1","sampDate_post_month2","sampDate_post_month3","sampDate_post_month4","sampDate_post_month6",
                                   "sampDate_post_month7","sampDate_post_month8","sampDate_post_month9","sampDate_post_month10","sampDate_post_month11",
                                   "ageCatPreFollow5-18","ageCatPreFollow18-30","ageCatPreFollow30-50","ageCatPreFollow50+",
                                   "intervalOccupationTop4Farmer","intervalOccupationTop4Student","intervalOccupationTop4Unemployed")))
modeloutsMultHouse <- rbind(modeloutsMultHouse,cis95)

print("multivariate -- household res")
print(modeloutsMultHouse)


# multivariate analysis household & indiv REs#########
modeloutsMultHouseIndiv <- c()

hold$sampDate_post_year <- as.factor(hold$sampDate_post_year)
hold$sampDate_post_month<- as.factor(hold$sampDate_post_month)
hold$sampDate_post_month <- factor(hold$sampDate_post_month,levels = c("5","1","2","3","4","6","7","8","9","10","11"))
modelM <- glmmTMB(prediction_round ~  sampDate_post_month + sampDate_post_year + ageCatPreFollow + intervalOccupationTop4 +(1|subjectNoAnon) +(1|houseAnon),
                  data=hold,
                  family = binomial(link="logit"))

cis95 <- exp(confint(modelM,parm=c("sampDate_post_year2017","sampDate_post_year2018","sampDate_post_year2019","sampDate_post_year2020","sampDate_post_year2021","sampDate_post_year2022",
                                   "sampDate_post_month1","sampDate_post_month2","sampDate_post_month3","sampDate_post_month4","sampDate_post_month6",
                                   "sampDate_post_month7","sampDate_post_month8","sampDate_post_month9","sampDate_post_month10","sampDate_post_month11",
                                   "ageCatPreFollow5-18","ageCatPreFollow18-30","ageCatPreFollow30-50","ageCatPreFollow50+",
                                   "intervalOccupationTop4Farmer","intervalOccupationTop4Student","intervalOccupationTop4Unemployed")))
modeloutsMultHouseIndiv <- rbind(modeloutsMultHouseIndiv,cis95)

print("multivariate -- indiv and household res")
print(modeloutsMultHouseIndiv)



# multivariate analysis household REs #########

modeloutsMultHouse <- c()
hold$sampDate_post_month <- as.factor(hold$sampDate_post_month)
hold$sampDate_post_month <- factor(hold$sampDate_post_month,levels = c("5","1","2","3","4","6","7","8","9","10","11"))
hold$sampDate_post_year <- as.factor(hold$sampDate_post_year)
modelM <- glmmTMB(prediction_round ~  avgHAICatPreFollow + sampDate_post_month+sampDate_post_year + ageCatPreFollow + intervalOccupationTop4 + (1|houseAnon),
                  data=hold,
                  family = binomial(link="logit"))

cis95 <- exp(confint(modelM,parm=c("sampDate_post_year2017","sampDate_post_year2018","sampDate_post_year2019","sampDate_post_year2020","sampDate_post_year2021","sampDate_post_year2022",
                                   "sampDate_post_month1","sampDate_post_month2","sampDate_post_month3","sampDate_post_month4","sampDate_post_month6",
                                   "sampDate_post_month7","sampDate_post_month8","sampDate_post_month9","sampDate_post_month10","sampDate_post_month11",
                                   "ageCatPreFollow5-18","ageCatPreFollow18-30","ageCatPreFollow30-50","ageCatPreFollow50+",
                                   "intervalOccupationTop4Farmer","intervalOccupationTop4Student","intervalOccupationTop4Unemployed",
                                   "avgHAICatPreFollow20-40","avgHAICatPreFollow40-80","avgHAICatPreFollow80-160","avgHAICatPreFollow160+")))
modeloutsMultHouse <- rbind(modeloutsMultHouse,cis95)

print("multivariate -- household res w/ titers")
print(modeloutsMultHouse)

# create confusion matrix for training data #####

confusion_mat_plot <- function(pred_val,true_val,figure_name,figure_folder,title_str,subtitle_str){
  tableval <- data.frame(confusionMatrix(as.factor(pred_val),true_val)$table)
  
  # create table for plotting with new labels
  plotTable <- tableval %>%
    mutate(goodbad = ifelse(tableval$Prediction == tableval$Reference, "good", "bad")) %>%
    group_by(Reference) %>%
    dplyr::mutate(prop = Freq/sum(Freq))
  
  # fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
  p <- ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
    geom_tile() +
    geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
    scale_fill_manual(values = c(good = "green", bad = "red")) +
    theme_bw() +
    ggtitle(title_str,subtitle = subtitle_str) +
    xlim(rev(levels(tableval$Reference)))
  print(p)
  ggsave(filename=figure_name,path = figure_folder, width = 6, height = 4, device='png', dpi=700)
  
}
confusion_mat_plot(as.factor(round(df[df$training==1,]$prediction)==1),as.factor(df[df$training==1,]$prediction_round==TRUE),paste0("S7.png"),figFolder,"XGboost: Confusion matrix for training data","")


# create distributiion of prediiction plot #####
p <- ggplot(data = df,aes(x=prediction)) + geom_histogram(binwidth = .025)
p <- p + scale_x_continuous("Model prediction")
p <- p + scale_y_continuous("Intervals")
p <- p + theme_pubr()
print(p)
ggsave(filename="S8.png",path = figFolder, width = 6, height = 4, device='png', dpi=700)


