##Spatial and temporal virus prevalence 

#for working on or off the VPN, change the library path to the personal lib 
.libPaths()   # get the path
.libPaths("C:/Program Files/R/R-3.6.2/library") #to change to an off datastore repository

##Packages
library(tidyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(cooccur)
library(scales)
library(dlookr)
library(viridis)
library(pBrackets)
library(extrafont)
if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
  detach("package:plyr", unload=TRUE) 
} 
library(plyr)
library(dplyr)

# Importing data ----------------------------------------------------------

#now updated w species from CO1 barcodes (apart from single mystery one...which I've put as virilis as its probs in the group but not sure what sp)
indiv_PCR_virus_data<-read.table("Indiv_virus_PCR_results.csv",header = TRUE, sep = ",")

indiv_PCR_virus_data$short_code<-as.character(indiv_PCR_virus_data$short_code)
indiv_PCR_virus_data$site<-as.factor(as.character(indiv_PCR_virus_data$site))
indiv_PCR_virus_data$month_year<-as.factor(as.character(indiv_PCR_virus_data$date))
indiv_PCR_virus_data$species<-as.factor(as.character(indiv_PCR_virus_data$species))
indiv_PCR_virus_data$no_flies<-as.numeric(as.character(indiv_PCR_virus_data$no_flies))
indiv_PCR_virus_data$chaq<-as.factor(as.character(indiv_PCR_virus_data$chaq))
indiv_PCR_virus_data$chaq_sh<-as.factor(as.character(indiv_PCR_virus_data$chaq_sh))
indiv_PCR_virus_data$Dimm_SV_L<-as.factor(as.character(indiv_PCR_virus_data$Dimm_SV_L))
indiv_PCR_virus_data$Dimm_SV_N<-as.factor(as.character(indiv_PCR_virus_data$Dimm_SV_N))
indiv_PCR_virus_data$galbut_virus_407<-as.factor(as.character(indiv_PCR_virus_data$galbut_virus_407))
indiv_PCR_virus_data$Dimm_Nora_A<-as.factor(as.character(indiv_PCR_virus_data$Dimm_Nora_A))
indiv_PCR_virus_data$Dimm_Nora_B<-as.factor(as.character(indiv_PCR_virus_data$Dimm_Nora_B))
indiv_PCR_virus_data$muthill_sh<-as.factor(as.character(indiv_PCR_virus_data$muthill_sh))
indiv_PCR_virus_data$muthill_suz<-as.factor(as.character(indiv_PCR_virus_data$muthill_suz))
indiv_PCR_virus_data$mel_nora_6220<-as.factor(as.character(indiv_PCR_virus_data$mel_nora_6220))
indiv_PCR_virus_data$mel_nora_qPCR<-as.factor(as.character(indiv_PCR_virus_data$mel_nora_qPCR))
indiv_PCR_virus_data$prestney_burn_B<-as.factor(as.character(indiv_PCR_virus_data$prestney_burn_B))
indiv_PCR_virus_data$prestney_burn_sh<-as.factor(as.character(indiv_PCR_virus_data$prestney_burn_sh))
indiv_PCR_virus_data$tranent_L<-as.factor(as.character(indiv_PCR_virus_data$tranent_L))
indiv_PCR_virus_data$tranent_M<-as.factor(as.character(indiv_PCR_virus_data$tranent_M))
indiv_PCR_virus_data$tranent_S<-as.factor(as.character(indiv_PCR_virus_data$tranent_S))
indiv_PCR_virus_data$grom_728<-as.factor(as.character(indiv_PCR_virus_data$grom_virus_728))
indiv_PCR_virus_data$motts_mill_221<-as.factor(as.character(indiv_PCR_virus_data$motts_mill_221))
indiv_PCR_virus_data$Dsub_nora_1665<-as.factor(as.character(indiv_PCR_virus_data$Dsub_nora_1665))
indiv_PCR_virus_data$year<-as.factor(as.character(indiv_PCR_virus_data$year))

#ImmSV_L_data<-data.frame(indiv_PCR_virus_data$short_code,
                          #indiv_PCR_virus_data$date,
                          #indiv_PCR_virus_data$year,
                          #indiv_PCR_virus_data$site,
                          #indiv_PCR_virus_data$species,
                          #indiv_PCR_virus_data$no_flies,
                          #indiv_PCR_virus_data$Dimm_SV_L)
#colnames(ImmSV_L_data)<-c("short_code","date","year","site","species","no_flies","Dimm_SV")

###Isolating the mel and imm data
#ImmSV_L_data_imm<-subset(ImmSV_L_data,species=="imm")
#ImmSV_L_data_mel<-subset(ImmSV_L_data,species=="mel")

##Reduced table with results from Chaq, ImmSV, Imm Nora, Prestney burn, Mel Nora, Tranent, Grom, Motts mill, Galbut, Dsub Nora + Muthill (results in 03_20)
#as of 4_20 updated w new sp identifications - now 15 species
Mar_20_data<-data.frame(indiv_PCR_virus_data$short_code,indiv_PCR_virus_data$date, indiv_PCR_virus_data$year,indiv_PCR_virus_data$site,indiv_PCR_virus_data$species,indiv_PCR_virus_data$no_flies,
                        indiv_PCR_virus_data$chaq,
                        indiv_PCR_virus_data$chaq_sh,
                        indiv_PCR_virus_data$Dimm_SV_L,
                        indiv_PCR_virus_data$Dimm_SV_N, 
                        indiv_PCR_virus_data$Dimm_Nora_A, 
                        indiv_PCR_virus_data$Dimm_Nora_B, 
                        indiv_PCR_virus_data$mel_nora_6220, 
                        indiv_PCR_virus_data$mel_nora_qPCR, 
                        indiv_PCR_virus_data$prestney_burn_B, 
                        indiv_PCR_virus_data$prestney_burn_sh, 
                        indiv_PCR_virus_data$muthill_sh,
                        indiv_PCR_virus_data$muthill_suz,
                        indiv_PCR_virus_data$tranent_L, 
                        indiv_PCR_virus_data$tranent_M, 
                        indiv_PCR_virus_data$tranent_S,
                        indiv_PCR_virus_data$grom_728,
                        indiv_PCR_virus_data$motts_mill_221,
                        indiv_PCR_virus_data$galbut_virus_407,
                        indiv_PCR_virus_data$Dsub_nora_1665)
colnames(Mar_20_data)<-c("short_code","date","year","site","species","no_flies","chaq","chaq_sh","Dimm_SV_L","Dimm_SV_N","Dimm_Nora_A","Dimm_Nora_B","mel_nora_6220","mel_nora_qPCR","prest_B_B","prest_B_sh","muthill_sh","muthill_suz","tranent_L","tranent_M","tranent_S","grom_728","motts_mill_221","galbut_407","Dsub_nora_1665")

# Calculating virus prevalence  -------------------------------------------

##DJ0 with help from JJ Welch
##24Sept09, revised 11sep2013

########### Calculating the 2DlogLik bounds on viral prevelance #########

require(compiler)
########### Core Likelihood function #########
calc.log.lik<-cmpfun(function(p,hitFlies,missFlies){
  p.misses<-log((1-p)^missFlies)
  p.hits<-log(1-((1-p)^hitFlies))
  return(sum(c(p.hits,p.misses)))
})

########### Convenience function that is used to identify bounds #########
calc.bound<-cmpfun(function(p,LL,hitFlies,missFlies){
  p.misses<-log((1-p)^missFlies)
  p.hits<-log(1-((1-p)^hitFlies))
  return(
    (LL-sum(c(p.hits,p.misses)))^2
  )
})

######### End-user function that caclulates the estimated prevelance and if requested 
#########(1) interval-log-likelihood bounds, 
#########(2) LRT for for consistency between pooled and non-pooled samples, 
#########(3) Plots the LL surface

InferPrevalence<-function(nFlies,hits,bounds=FALSE,interval=2,test.consistant=FALSE, plot=FALSE){
  #strip out any NAs
  hits<-as.numeric(hits[!is.na(hits)])
  nFlies<-as.numeric(nFlies[!is.na(nFlies)])
  if((length(hits)==0)|(length(nFlies)==0)){return(NA)}
  #seperate out the pots into hits and misses
  hitFlies<-nFlies[which(as.logical(hits))]
  missFlies<-nFlies[which(!hits)]
  
  # maximise the LL (which is negative)
  optimise(f = calc.log.lik, interval = c(0,1), hitFlies, missFlies, maximum = TRUE,tol = .Machine$double.eps)->estimate	
  estimate$maximum->p
  estimate$objective->ML
  
  #Find the bounds (one at a time, seaching above and below the ML estimate
  if(bounds){
    optimise(f = calc.bound, interval = c(0,p), ML-interval, hitFlies, missFlies, maximum = FALSE,tol = .Machine$double.eps)$minimum->lower
    optimise(f = calc.bound, interval = c(p,1), ML-interval, hitFlies, missFlies, maximum = FALSE,tol = .Machine$double.eps)$minimum->upper
  }
  
  #If required to test for consistency between pooled and un-pooled samples
  if(test.consistant){
    #separate the data
    if((sum(nFlies>1)>0)&(sum(nFlies==1)>0)){
      hitFlies[hitFlies>1]->B.hit
      missFlies[missFlies>1]->B.miss
      hitFlies[hitFlies==1]->S.hit
      missFlies[missFlies==1]->S.miss
      #make the two estimates
      optimise(f = calc.log.lik, interval = c(0,1), B.hit, B.miss, maximum = TRUE,tol = .Machine$double.eps)->B.estimate
      optimise(f = calc.log.lik, interval = c(0,1), S.hit, S.miss, maximum = TRUE,tol = .Machine$double.eps)->S.estimate		
      c(B.estimate$maximum,S.estimate$maximum)->two.p
      #calculate their joint LL and do an LRT
      sum(B.estimate$objective,S.estimate$objective)->ML2
      2*(ML2-ML)->TwoDeltaLL
      pchisq(TwoDeltaLL, df =  1,lower.tail=FALSE)->p.value
    }else{
      test.consistant<-FALSE
    }
  }
  
  #Plot, if requested
  if(plot){
    if(!test.consistant){
      surface<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,hitFlies,missFlies))
      plot(seq(0,1,0.0001),surface,type="l",ylim=c((ML-10),ML),xlab="Prevelance",ylab="log Likelihood")
      abline(v=p,col="red",lwd=4)	
      if(bounds){
        abline(h=(ML-interval))
        abline(v=lower,col="red",lty=3)
        abline(v=upper,col="red",lty=3)
      }
    }
    if(test.consistant){
      surface<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,hitFlies,missFlies))
      surface1<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,B.hit,B.miss))
      surface2<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,S.hit,S.miss))
      plot(seq(0,1,0.0001),surface1,type="l",ylim=c(min(c(max(surface),max(surface1),max(surface2)))-10,max(c(surface,surface1,surface2))),xlab="Prevelance",ylab="log Likelihood",col="green")
      points(seq(0,1,0.0001),surface2,type="l",col="blue")
      points(seq(0,1,0.0001),surface,type="l",col="red")
      abline(v=c(p,two.p[1],two.p[2]),col=c("red","green","blue"),lwd=4)
    }
    
  }
  
  #contruct the return list
  result<-list()
  result$prevelance<-p
  result$log.liklihood<-ML
  if(bounds){result$bounds<-c(lower,upper)}
  if(test.consistant){
    result$alt.estimates<-two.p
    result$two.delta.LL<-TwoDeltaLL
    result$p.value<-p.value
    if(p.value<0.05){result$consistent<-FALSE}else{result$consistent<-TRUE}
  }
  return(result)
}


# Colour blind friendly palette (7 values):
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Longer pallette (10 values)
cbPalette_long <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")

# Prevance of DimmSV in Dimm -------------------------------------

#InferPrevalence(nFlies,hits,bounds=TRUE,interval=2,test.consistant=TRUE,plot=TRUE)

##converting P and A to 1 and 0 in DimmSV_imm presence absence col
#ImmSV_L_data_imm$Dimm_SV <- as.character(ImmSV_L_data_imm$Dimm_SV)
#ImmSV_L_data_imm$Dimm_SV[ImmSV_L_data_imm$Dimm_SV == "P"] <- 1
#ImmSV_L_data_imm$Dimm_SV <- as.character(ImmSV_L_data_imm$Dimm_SV)
#ImmSV_L_data_imm$Dimm_SV[ImmSV_L_data_imm$Dimm_SV == "A"] <- 0

#setting hits (vector of 0s and 1s which equate to P:A of DimmSV)
#hits<-ImmSV_L_data_imm$Dimm_SV
#setting nFlies (vector of number of flies in each pool)
#nFlies<-ImmSV_L_data_imm$no_flies

#DimmSV_imm_total<-InferPrevalence(nFlies,hits,bounds=TRUE,interval=2,test.consistant=TRUE,plot=TRUE)


# Prevalence of DimmSV in Imm over time ----------------------------------
#DimmSV_imm_PoT<-ImmSV_L_data_imm %>%
                    #group_by(date) %>%
                    #summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                                      #hits=Dimm_SV,
                                                                      #bounds=TRUE,
                                                                      #interval=2,test.consistant=TRUE,
                                                                      #plot=FALSE)$prevelance)) %>% as.data.frame()

#DimmSV_imm_LLoT<-ImmSV_L_data_imm %>%
                    #group_by(month_year) %>%
                      #summarize(twologlik = as.numeric(InferPrevalence(nFlies=no_flies,
                                                                     #hits=Dimm_SV,
                                                                     #bounds=TRUE,
                                                                     #interval=2,test.consistant=TRUE,
                                                                     #plot=FALSE)$two.delta.LL)) %>% as.data.frame()
############Error: Column `twologlik` must be length 1 (a summary value), not 0

##Plotting prevalence of DimmSV over time in Imm

#glimpse(DimmSV_imm_PoT)
#tbl_df(DimmSV_imm_PoT)

#DimmSV_imm_PoT$date<-as.Date(DimmSV_imm_PoT$date,format = "%d/%m/%Y")
#DimmSV_imm_PoT$prevalence<-as.numeric(as.character(DimmSV_imm_PoT$prevalence))

#creating plot of prevalence over time 
#DimmSV_imm_Potplot <- ggplot(DimmSV_imm_PoT, aes(date, prevalence)) +
  # geom_line(colour="black",lwd=1.5)+
  # ggtitle("Prevalence of D.imm SV in D.imm\n") +
  # labs(x=NULL,y=NULL) +
  # scale_x_date(date_labels ="%b-%y",date_breaks='2 months') + xlab("") +
  # theme(panel.background = element_blank(),
  #       panel.grid.minor = element_blank(), 
  #       panel.grid.major = element_line(color = "gray50", size = 0.5),
  #       panel.grid.major.x = element_blank(), 
  #       text = element_text(size=18),
  #       axis.text.y = element_text(), 
  #       axis.text.x = element_text(),
  #       plot.title = element_text(hjust = -0.05, vjust=2, size = 14))

#subsetting the prevalence over time data into high and low host numbers sections of the year - and calculating prevalence from those subsets.
# DimmSV_imm_PoT_SO16<-DimmSV_imm_PoT %>%
#   filter(date == "2016-09-23" | date == "2016-09-26" | date == "2016-10-27") 
# cat<-c("SO16","SO16","SO16")
# DimmSV_imm_PoT_SO16<-cbind(DimmSV_imm_PoT_SO16,cat)
# 
# DimmSV_imm_PoT_low17<-DimmSV_imm_PoT %>%
#   filter(date == "2017-05-25" | date == "2017-06-13") 
# cat<-c("low17","low17")
# DimmSV_imm_PoT_low17<-cbind(DimmSV_imm_PoT_low17,cat)
# 
# DimmSV_imm_PoT_high17<-DimmSV_imm_PoT %>%
#   filter(date == "2017-07-21" | date == "2017-08-31" | date == "2017-09-28" | date == "2017-10-18") 
# cat<-c("high17","high17","high17","high17")
# DimmSV_imm_PoT_high17<-cbind(DimmSV_imm_PoT_high17,cat)
# 
# DimmSV_imm_PoT_low18<-DimmSV_imm_PoT %>%
#   filter(date == "2018-04-23" | date == "2018-05-24" | date == "2018-06-22") 
# cat<-c("low18","low18","low18")
# DimmSV_imm_PoT_low18<-cbind(DimmSV_imm_PoT_low18,cat)
# 
# DimmSV_imm_PoT_high18<-DimmSV_imm_PoT %>%
#   filter(date == "2018-07-13" | date == "2018-08-21" | date == "2018-09-30" | date == "2018-10-20") 
# cat<-c("high18","high18","high18","high18")
# DimmSV_imm_PoT_high18<-cbind(DimmSV_imm_PoT_high18,cat)
# 
# DimmSV_imm_PoT_subsetted<-rbind(DimmSV_imm_PoT_SO16,DimmSV_imm_PoT_low17,DimmSV_imm_PoT_high17,DimmSV_imm_PoT_low18,DimmSV_imm_PoT_high18)
# 
# DimmSV_imm_PoT_cat<- DimmSV_imm_PoT_subsetted %>%
#   group_by(cat) %>%
#   summarise(mean_prevalence = mean(prevalence),
#             sd_prevalence = sqrt(var(prevalence)))

# Barplot of prevalence of ImmSV by year ---------------------------------
# DimmSV_imm_PbY<-ImmSV_L_data_imm %>%
#   group_by(year) %>%
#   summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
#                                                     hits=Dimm_SV,
#                                                     bounds=TRUE,
#                                                     interval=2,test.consistant=TRUE,
#                                                     plot=FALSE)$prevelance)
#             ) %>% as.data.frame()
# 
# barplot_ImmSV_Imm<-ggplot(DimmSV_imm_PbY, aes(x=year,y=prevalence)) +
#   geom_col(color="black")+theme(text = element_text(size=14),axis.text.y = element_text(), axis.text.x = element_text(size=16),plot.title = element_text(hjust = 0.5),legend.position = "none")+
#   ylim(c(0.0,1.0))+
#   labs(title="Prevalence of ImmSV in D.imm", x="Year", y="Prevalence")


# Barplot of ImmSV prevlance in Species by year ---------------------------

# ##converting P and A to 1 and 0 in DimmSV_imm presence absence col
# ImmSV_L_data$Dimm_SV <- as.character(ImmSV_L_data$Dimm_SV)
# ImmSV_L_data$Dimm_SV[ImmSV_L_data$Dimm_SV == "P"] <- 1
# ImmSV_L_data$Dimm_SV <- as.character(ImmSV_L_data$Dimm_SV)
# ImmSV_L_data$Dimm_SV[ImmSV_L_data$Dimm_SV == "A"] <- 0
# 
# #setting hits (vector of 0s and 1s which equate to P:A of DimmSV)
# hits<-ImmSV_L_data$Dimm_SV
# #setting nFlies (vector of number of flies in each pool)
# nFlies<-ImmSV_L_data$no_flies
# 
# ##calculating prevalence across species and year
# DimmSV_PbY<-ImmSV_L_data %>%
#   group_by(species,year) %>%
#   summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
#                                                     hits=Dimm_SV,
#                                                     bounds=TRUE,
#                                                     interval=2,test.consistant=TRUE,
#                                                     plot=FALSE)$prevelance),
#             lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
#                                                     hits=Dimm_SV,
#                                                     bounds=TRUE,
#                                                     interval=2,test.consistant=TRUE,
#                                                     plot=FALSE)$bounds[1]),
#             upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
#                                                      hits=Dimm_SV,
#                                                      bounds=TRUE,
#                                                      interval=2,test.consistant=TRUE,
#                                                      plot=FALSE)$bounds[2])) %>% 
#   as.data.frame() 
# 
# DimmSV_PbY_imm_mel<-subset(DimmSV_PbY,species=="imm"|species == "mel")

#defining the limits of the errorbars
# limits <- aes(ymax = upper_bound, ymin=lower_bound)
# dodge <- position_dodge(width=0.9)

# barplot_ImmSV_PbY<-ggplot(DimmSV_PbY, aes(x=species,y=prevalence,fill=year)) +
#   geom_bar(position="dodge",stat="identity")+theme(text = element_text(size=14),axis.text.y = element_text(), axis.text.x = element_text(size=16),plot.title = element_text(hjust = 0.5),legend.position = "none")+
#   ylim(c(0.0,1.0))+
#   labs(title="Prevalence of ImmSV", x="Species", y="Prevalence")
# 
# barplot_ImmSV_PbY_imm_mel<-ggplot(DimmSV_PbY_imm_mel, aes(x=year,y=prevalence,fill=species)) +
#   geom_bar(position="dodge",stat="identity")+theme(text = element_text(size=14),axis.text.y = element_text(), axis.text.x = element_text(size=16),plot.title = element_text(hjust = 0.5),legend.position = "none",panel.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_line(color = "gray50", size = 0.5),panel.grid.major.x = element_blank())+
#   ylim(c(0.0,0.8))+
#   labs(title="Prevalence of ImmSV", x="Year", y="Prevalence")+
#   scale_fill_manual(values = c("orange", "skyblue"))+
#   geom_errorbar(limits,position =dodge, width = 0.1)
# 
# #converting to binary response variable in mel dataset for ImmSV
# ImmSV_L_data_mel$Dimm_SV <- as.character(ImmSV_L_data_mel$Dimm_SV)
# ImmSV_L_data_mel$Dimm_SV[ImmSV_L_data_mel$Dimm_SV == "P"] <- 1
# ImmSV_L_data_mel$Dimm_SV <- as.character(ImmSV_L_data_mel$Dimm_SV)
# ImmSV_L_data_mel$Dimm_SV[ImmSV_L_data_mel$Dimm_SV == "A"] <- 0
# 
# DimmSV_mel_total<-ImmSV_L_data_mel %>%
#   summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
#                                                     hits=Dimm_SV,
#                                                     bounds=TRUE,
#                                                     interval=2,test.consistant=TRUE,
#                                                     plot=FALSE)$prevelance),
#             lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
#                                                      hits=Dimm_SV,
#                                                      bounds=TRUE,
#                                                      interval=2,test.consistant=TRUE,
#                                                      plot=FALSE)$bounds[1]),
#             upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
#                                                      hits=Dimm_SV,
#                                                      bounds=TRUE,
#                                                      interval=2,test.consistant=TRUE,
#                                                      plot=FALSE)$bounds[2])) %>% 
#   as.data.frame() 
#   
# DimmSV_mel_total

# Per-species overall prevalence of ImmSV, ImmNora, Prestney Burn, mel nora, muthill, grom and motts mill --------

##converting P and A to 1 and 0 in virus presence absence cols, and combining results from multiple primers

##Chaq virus 

#chaq
Mar_20_data$chaq <- as.character(Mar_20_data$chaq)
Mar_20_data$chaq[Mar_20_data$chaq == "P"] <- 1
Mar_20_data$chaq[Mar_20_data$chaq == "A"] <- 0
Mar_20_data$chaq<-as.numeric(as.character(Mar_20_data$chaq))

#chaq_sh
Mar_20_data$chaq_sh <- as.character(Mar_20_data$chaq_sh)
Mar_20_data$chaq_sh[Mar_20_data$chaq_sh == "P"] <- 1
Mar_20_data$chaq_sh[Mar_20_data$chaq_sh == "A"] <- 0
Mar_20_data$chaq_sh<-as.numeric(as.character(Mar_20_data$chaq_sh))

#Combining
Mar_20_data$chaq_comb<-rowSums(Mar_20_data[,grep("chaq", colnames(Mar_20_data))])
Mar_20_data$chaq_comb[Mar_20_data$chaq_comb > 0] <- 1

##Dimm_SV
#L
Mar_20_data$Dimm_SV_L <- as.character(Mar_20_data$Dimm_SV_L)
Mar_20_data$Dimm_SV_L[Mar_20_data$Dimm_SV_L == "P"] <- 1
Mar_20_data$Dimm_SV_L[Mar_20_data$Dimm_SV_L == "A"] <- 0
Mar_20_data$Dimm_SV_L<-as.numeric(as.character(Mar_20_data$Dimm_SV_L))
#N
Mar_20_data$Dimm_SV_N <- as.character(Mar_20_data$Dimm_SV_N)
Mar_20_data$Dimm_SV_N[Mar_20_data$Dimm_SV_N == "P"] <- 1
Mar_20_data$Dimm_SV_N[Mar_20_data$Dimm_SV_N == "A"] <- 0
Mar_20_data$Dimm_SV_N<-as.numeric(as.character(Mar_20_data$Dimm_SV_N))
#Combining
Mar_20_data$Dimm_SV_comb<-rowSums(Mar_20_data[,grep("Dimm_SV_", colnames(Mar_20_data))])
Mar_20_data$Dimm_SV_comb[Mar_20_data$Dimm_SV_comb > 0] <- 1

##Dimm_Nora
#A
Mar_20_data$Dimm_Nora_A <- as.character(Mar_20_data$Dimm_Nora_A)
Mar_20_data$Dimm_Nora_A[Mar_20_data$Dimm_Nora_A == "P"] <- 1
Mar_20_data$Dimm_Nora_A[Mar_20_data$Dimm_Nora_A == "A"] <- 0
Mar_20_data$Dimm_Nora_A<-as.numeric(as.character(Mar_20_data$Dimm_Nora_A))
#B
Mar_20_data$Dimm_Nora_B <- as.character(Mar_20_data$Dimm_Nora_B)
Mar_20_data$Dimm_Nora_B[Mar_20_data$Dimm_Nora_B == "P"] <- 1
Mar_20_data$Dimm_Nora_B[Mar_20_data$Dimm_Nora_B == "A"] <- 0
Mar_20_data$Dimm_Nora_B<-as.numeric(as.character(Mar_20_data$Dimm_Nora_B))
#Combining
Mar_20_data$Dimm_Nora_comb<-rowSums(Mar_20_data[,grep("Dimm_Nora_", colnames(Mar_20_data))])
Mar_20_data$Dimm_Nora_comb[Mar_20_data$Dimm_Nora_comb > 0] <- 1

##mel_nora
#6220
Mar_20_data$mel_nora_6220 <- as.character(Mar_20_data$mel_nora_6220)
Mar_20_data$mel_nora_6220[Mar_20_data$mel_nora_6220 == "P"] <- 1
Mar_20_data$mel_nora_6220[Mar_20_data$mel_nora_6220 == "A"] <- 0
Mar_20_data$mel_nora_6220<-as.numeric(as.character(Mar_20_data$mel_nora_6220))
#qpcr
Mar_20_data$mel_nora_qPCR <- as.character(Mar_20_data$mel_nora_qPCR)
Mar_20_data$mel_nora_qPCR[Mar_20_data$mel_nora_qPCR == "P"] <- 1
Mar_20_data$mel_nora_qPCR[Mar_20_data$mel_nora_qPCR == "A"] <- 0
Mar_20_data$mel_nora_qPCR<-as.numeric(as.character(Mar_20_data$mel_nora_qPCR))
#Combining
Mar_20_data$mel_nora_comb<-rowSums(Mar_20_data[,grep("mel_nora_", colnames(Mar_20_data))])
Mar_20_data$mel_nora_comb[Mar_20_data$mel_nora_comb > 0] <- 1

##Prestney Burn
#B
Mar_20_data$prest_B_B <- as.character(Mar_20_data$prest_B_B)
Mar_20_data$prest_B_B[Mar_20_data$prest_B_B == "P"] <- 1
Mar_20_data$prest_B_B[Mar_20_data$prest_B_B == "A"] <- 0
Mar_20_data$prest_B_B<-as.numeric(as.character(Mar_20_data$prest_B_B))
#sh
Mar_20_data$prest_B_sh <- as.character(Mar_20_data$prest_B_sh)
Mar_20_data$prest_B_sh[Mar_20_data$prest_B_sh == "P"] <- 1
Mar_20_data$prest_B_sh[Mar_20_data$prest_B_sh == "A"] <- 0
Mar_20_data$prest_B_sh<-as.numeric(as.character(Mar_20_data$prest_B_sh))
#Combining
Mar_20_data$prest_B_comb<-rowSums(Mar_20_data[,grep("prest_B_", colnames(Mar_20_data))])
Mar_20_data$prest_B_comb[Mar_20_data$prest_B_comb > 0] <- 1

##muthill
#6220
Mar_20_data$muthill_sh <- as.character(Mar_20_data$muthill_sh)
Mar_20_data$muthill_sh[Mar_20_data$muthill_sh == "P"] <- 1
Mar_20_data$muthill_sh[Mar_20_data$muthill_sh == "A"] <- 0
Mar_20_data$muthill_sh<-as.numeric(as.character(Mar_20_data$muthill_sh))
#suz 
Mar_20_data$muthill_suz <- as.character(Mar_20_data$muthill_suz)
Mar_20_data$muthill_suz[Mar_20_data$muthill_suz == "P"] <- 1
Mar_20_data$muthill_suz[Mar_20_data$muthill_suz == "A"] <- 0
Mar_20_data$muthill_suz<-as.numeric(as.character(Mar_20_data$muthill_suz))
#Combining
Mar_20_data$muthill_comb<-rowSums(Mar_20_data[,grep("muthill_", colnames(Mar_20_data))])
Mar_20_data$muthill_comb[Mar_20_data$muthill_comb > 0] <- 1

##tranent
#L
Mar_20_data$tranent_L <- as.character(Mar_20_data$tranent_L)
Mar_20_data$tranent_L[Mar_20_data$tranent_L == "P"] <- 1
Mar_20_data$tranent_L[Mar_20_data$tranent_L == "A"] <- 0
Mar_20_data$tranent_L<-as.numeric(as.character(Mar_20_data$tranent_L))
#M
Mar_20_data$tranent_M <- as.character(Mar_20_data$tranent_M)
Mar_20_data$tranent_M[Mar_20_data$tranent_M == "P"] <- 1
Mar_20_data$tranent_M[Mar_20_data$tranent_M == "A"] <- 0
Mar_20_data$tranent_M<-as.numeric(as.character(Mar_20_data$tranent_M))
#S
Mar_20_data$tranent_S <- as.character(Mar_20_data$tranent_S)
Mar_20_data$tranent_S[Mar_20_data$tranent_S == "P"] <- 1
Mar_20_data$tranent_S[Mar_20_data$tranent_S == "A"] <- 0
Mar_20_data$tranent_S<-as.numeric(as.character(Mar_20_data$tranent_S))
#combined
Mar_20_data$tranent_comb<-rowSums(Mar_20_data[,grep("tranent_", colnames(Mar_20_data))])
Mar_20_data$tranent_comb[Mar_20_data$tranent_comb > 0] <- 1

##Grom
#728 primers 
Mar_20_data$grom_728 <- as.character(Mar_20_data$grom_728)
Mar_20_data$grom_728[Mar_20_data$grom_728 == "P"] <- 1
Mar_20_data$grom_728[Mar_20_data$grom_728 == "A"] <- 0
Mar_20_data$grom_728<-as.numeric(as.character(Mar_20_data$grom_728))

#Motts mill
#221 primers 
Mar_20_data$motts_mill_221 <- as.character(Mar_20_data$motts_mill_221)
Mar_20_data$motts_mill_221[Mar_20_data$motts_mill_221 == "P"] <- 1
Mar_20_data$motts_mill_221[Mar_20_data$motts_mill_221 == "A"] <- 0
Mar_20_data$motts_mill_221<-as.numeric(as.character(Mar_20_data$motts_mill_221))

#Galbut
#407 primers
Mar_20_data$galbut_407 <- as.character(Mar_20_data$galbut_407)
Mar_20_data$galbut_407[Mar_20_data$galbut_407 == "P"] <- 1
Mar_20_data$galbut_407[Mar_20_data$galbut_407 == "A"] <- 0
Mar_20_data$galbut_407<-as.numeric(as.character(Mar_20_data$galbut_407))

#Dsub Nora
#1665F primers
Mar_20_data$Dsub_nora_1665 <- as.character(Mar_20_data$Dsub_nora_1665)
Mar_20_data$Dsub_nora_1665[Mar_20_data$Dsub_nora_1665 == "P"] <- 1
Mar_20_data$Dsub_nora_1665[Mar_20_data$Dsub_nora_1665 == "A"] <- 0
Mar_20_data$Dsub_nora_1665<-as.numeric(as.character(Mar_20_data$Dsub_nora_1665))

##Reducing table to DimmSV, Dimm Nora, mel nora, Prestney burn, Muthill and tranent combined cols
Mar_20_data_reduced<-data.frame(Mar_20_data$short_code,Mar_20_data$date, Mar_20_data$year,Mar_20_data$site,Mar_20_data$species,Mar_20_data$no_flies,Mar_20_data$chaq_comb,Mar_20_data$Dimm_SV_comb,Mar_20_data$Dimm_Nora_comb,Mar_20_data$mel_nora_comb,Mar_20_data$prest_B_comb,Mar_20_data$muthill_comb,Mar_20_data$tranent_comb,Mar_20_data$grom_728,Mar_20_data$motts_mill_221,Mar_20_data$galbut_407,Mar_20_data$Dsub_nora_1665) 
                                
colnames(Mar_20_data_reduced)<-c("short_code","date","year","site","species","no_flies","chaq","Dimm_SV","Dimm_Nora","mel_nora","prest_B","muthill","tranent","grom","motts_mill","galbut","sub_nora")

# Calc. prev of viruses by species ----------------------------------------

##Chaq
Chaq_comb_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=chaq,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=chaq,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=chaq,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group
Chaq_comb_obs3_imm_mel<-subset(Chaq_comb_byspecies, 
                                 species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
Chaq_comb_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

##DimmSV
DimmSV_comb_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=Dimm_SV,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=Dimm_SV,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=Dimm_SV,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group
DimmSV_comb_obs3_imm_mel<-subset(DimmSV_comb_byspecies, 
                                 species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
DimmSV_comb_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

##Dimm Nora
Dimm_Nora_comb_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=Dimm_Nora,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=Dimm_Nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=Dimm_Nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group
Dimm_Nora_comb_obs3_imm_mel<-subset(Dimm_Nora_comb_byspecies, 
                                 species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
Dimm_Nora_comb_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

##mel nora
mel_nora_comb_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=mel_nora,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=mel_nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=mel_nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group
mel_nora_comb_obs3_imm_mel<-subset(mel_nora_comb_byspecies, 
                                    species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
mel_nora_comb_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

##Prestney burn
prest_B_comb_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=prest_B,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=prest_B,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=prest_B,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group
prest_B_comb_obs3_imm_mel<-subset(prest_B_comb_byspecies, 
                                   species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
prest_B_comb_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

##muthill
muthill_comb_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=muthill,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=muthill,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=muthill,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group
muthill_comb_obs3_imm_mel<-subset(muthill_comb_byspecies, 
                                  species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
muthill_comb_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

##tranent
tranent_comb_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=tranent,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=tranent,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=tranent,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group
tranent_comb_obs3_imm_mel<-subset(tranent_comb_byspecies, 
                                  species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
tranent_comb_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

#Grom virus - calculating prevalence across species 
grom_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=grom,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=grom,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=grom,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group
grom_obs3_imm_mel<-subset(grom_byspecies,
                          species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
grom_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

#Motts mill virus - calculating prevalence across species 
  motts_mill_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=motts_mill,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=motts_mill,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=motts_mill,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 
  
#separating the prevalence for Imm, mel, and obs group
motts_mill_obs3_imm_mel<-subset(motts_mill_byspecies,
                                species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")  
motts_mill_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

#Galbut virus - calculating prevalence across species 
  galbut_byspecies<-Mar_20_data_reduced %>%
    group_by(species) %>%
    summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                      hits=galbut,
                                                      bounds=TRUE,
                                                      interval=2,test.consistant=TRUE,
                                                      plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                       hits=galbut,
                                                       bounds=TRUE,
                                                       interval=2,test.consistant=TRUE,
                                                       plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                       hits=galbut,
                                                       bounds=TRUE,
                                                       interval=2,test.consistant=TRUE,
                                                       plot=FALSE)$bounds[2])) %>% 
    as.data.frame()  
  
#separating the prevalence for Imm, mel, and obs group
galbut_obs3_imm_mel<-subset(galbut_byspecies,
                                species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")   
galbut_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

#Dsub Nora virus - calculating prevalence across species 
sub_nora_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=sub_nora,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=sub_nora,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=sub_nora,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$bounds[2])) %>% 
    as.data.frame()  
  
#separating the prevalence for Imm, mel, and obs group
sub_nora_obs3_imm_mel<-subset(sub_nora_byspecies,
                          species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
sub_nora_obs3_imm_mel$species<-c("Dimm","Dmel","Dobs","Dsus","Dsub")

##Producing horizontal box plots of prevalence of each virus by species. 
##defining the limits of the error bars
limits <- aes(ymax = upper_bound, ymin=lower_bound)
dodge <- position_dodge(width=0.9)

#Chaq
barplot_chaq_comb_PbS<-ggplot(Chaq_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_blank(), axis.text.x = element_text(size = 19),plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Chaq virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#ImmSV
barplot_ImmSV_comb_PbS<-ggplot(DimmSV_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Immigrans Sigma virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#ImmNora
barplot_Imm_Nora_comb_PbS<-ggplot(Dimm_Nora_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 12), axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Immigrans Nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#mel nora
barplot_mel_nora_comb_PbS<-ggplot(mel_nora_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_blank(), axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Melanogaster Nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#grom
barplot_grom_PbS<-ggplot(grom_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_blank(), axis.text.x = element_text(size = 19) ,plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Grom virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#motts mill
barplot_motts_mill_PbS<-ggplot(motts_mill_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_blank(), axis.text.x = element_text(size = 19), plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Motts Mill virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#galbut
barplot_galbut_PbS<-ggplot(galbut_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_blank(), axis.text.x = element_text(size = 19), plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Galbut virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#Dsub nora 
barplot_sub_nora_PbS<-ggplot(sub_nora_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_blank(), axis.text.x = element_text(size = 19), plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Dsub nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

# FOR BEN - Mel nora prevalence across species ----------------------------
#mel nora
barplot_mel_nora_comb_PbS_all<-ggplot(mel_nora_comb_byspecies, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(), axis.text.x = element_text(),plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Melanogaster Nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)


##Summarising the prevalence and number of infections/number sampled across species 
n_mel_nora_by_species <- Mar_20_data_reduced %>% 
  group_by(species, mel_nora) %>%
  summarise(n_pools=length(mel_nora),
            n_flies=sum(no_flies))

write.csv(mel_nora_comb_byspecies,file = "mel_nora_prev_bysp_all.csv")
write.csv(n_mel_nora_by_species,file = "mel_nora_infection_by_n_pool.csv")

# FOR THIRD YEAR REVIEW - Tranent virus across all species
#tranent
barplot_tranent_comb_PbS_all<-ggplot(tranent_comb_byspecies, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(), axis.text.x = element_text(),plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Tranent virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#for just imm, obs3, mel
##defining the limits of the error bars
limits_log <- aes(ymax = upper_bound*100, ymin=lower_bound*100)
dodge <- position_dodge(width=0.9)

barplot_tranent_PbS<-ggplot(tranent_comb_obs3_imm_mel, aes(x=species,y=log1p(prevalence),fill=species))+
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19), plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0,1))+
  labs(title="Tranent virus in common species", x="", y="")+
  geom_errorbar(limits_log,position =dodge, width = 0.1)

#for just the rarer species 
#excluding the prevalence for Imm, mel, and obs group from the table - selecting the other sp
tranent_comb_rare_sp<-subset(tranent_comb_byspecies, 
                                  species == "bus" | species =="chy" | species == "def" | species == "fun" | species == "hel" | species == "hyd" | species == "pha" | species == "tri" | species == "cam" | species == "vir")
tranent_comb_rare_sp$species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dpha","Dtri","Dvir")

#for rare species
#trying to put a log transformation on this data but keep the scale linear?
barplot_tranent_PbS_rare<-ggplot(tranent_comb_rare_sp, aes(x=species,y=prevalence*100,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette_long)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19), plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  labs(title="Tranent virus in rare(r) species", x="", y="")+
  geom_errorbar(limits_log,position=dodge, width = 0.1)+
  scale_y_continuous(trans = log_trans(), 
                     breaks = trans_breaks("log", function(x) exp(x)),
                     labels = c(-40,-30,-20,-10,0,10))

#have done a log10 transformation on data (as percentages so *100) with the zero values (taken as anything lower than 0.001 proportion) changed to 0.01% so that log10 0.01 = -2
#added to cols in tranent_comb_rare_sp so that can be plotted
limits_log <- aes(ymax =log10(upper_bound100) ,ymin=log10(lower_bound100))
#dodge <- position_dodge(width=0.9)

#prepping the dataframe
tranent_comb_rare_sp$prevalence100<-tranent_comb_rare_sp$prevalence*100
tranent_comb_rare_sp$lower_bound100<-tranent_comb_rare_sp$lower_bound*100
tranent_comb_rare_sp$upper_bound100<-tranent_comb_rare_sp$upper_bound*100
#replacing the low vals w 0.01%
tranent_comb_rare_sp[1:4,5:6]<-0.01
tranent_comb_rare_sp[6:10,5:6]<-0.01
tranent_comb_rare_sp$log10upper_bound100<-log10(tranent_comb_rare_sp$upper_bound100)
tranent_comb_rare_sp$log10prevalence100<-log10(tranent_comb_rare_sp$prevalence100)
tranent_comb_rare_sp$log10lower_bound100<-log10(tranent_comb_rare_sp$lower_bound100)

#plotting
barplot_tranent_PbS_rare<-ggplot(tranent_comb_rare_sp, aes(x=species,y=log10(prevalence100),fill=species)) +
  geom_bar(stat="identity",position = , colour="black")+
  scale_fill_manual(values = cbPalette_long)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=17,colour = "grey25"),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19), plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(0.00005,2))+
  labs(title="Tranent virus in rare(r) species", x="Species", y="log10(prevalence%)")+
  geom_errorbar(limits_log, width = 0.1)

#same method on the common sp
#prepping the dataframe
tranent_comb_obs3_imm_mel$prevalence100<-tranent_comb_obs3_imm_mel$prevalence*100
tranent_comb_obs3_imm_mel$lower_bound100<-tranent_comb_obs3_imm_mel$lower_bound*100
tranent_comb_obs3_imm_mel$upper_bound100<-tranent_comb_obs3_imm_mel$upper_bound*100
#in this case the lowest species w observations is at 0.02 %, I've changed the single v low (-15) % prev to 0.01 (in the second row of the data)
#upper stays the same and is logged
tranent_comb_obs3_imm_mel$log10upper_bound100<-log10(tranent_comb_obs3_imm_mel$upper_bound100)
#lower100 and prev100 have second row changed to 0.01% and then logged
tranent_comb_obs3_imm_mel[2,5:6]<-0.01
tranent_comb_obs3_imm_mel$log10lower_bound100<-log10(tranent_comb_obs3_imm_mel$lower_bound100)
tranent_comb_obs3_imm_mel$log10prevalence100<-log10(tranent_comb_obs3_imm_mel$prevalence100)

barplot_tranent_PbS<-ggplot(tranent_comb_obs3_imm_mel, aes(x=species,y=log10prevalence100,fill=species))+
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=17,colour = "grey25"),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19), plot.title = element_text(hjust = 0.5, size = 18),legend.position = "none")+
  ylim(c(-3,3))+
  labs(title="Tranent virus in common species", x="Species", y="Log10(prevalence%)")+
  geom_errorbar(limits_log,position =dodge, width = 0.1)


##Summarising the prevalence and number of infections/number sampled across species 
n_tranent_by_species <- Mar_20_data_reduced %>% 
  group_by(species, tranent) %>%
  summarise(n_pools=length(tranent),
            n_flies=sum(no_flies))

#prestney burn
barplot_prest_B_comb_PbS<-ggplot(prest_B_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3),text = element_text(size=14),axis.text.y = element_blank(), axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Prestney burn virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#muthill
barplot_muthill_comb_PbS<-ggplot(muthill_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3),text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19),plot.title = element_text(hjust = 0.5,size=18),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Muthill virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

###Creating plot with Imm nora, Mel nora, prestney burn and muthill

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

four_virus_hbars_prev_by_sp<-multiplot(barplot_Imm_Nora_comb_PbS,barplot_mel_nora_comb_PbS, barplot_prest_B_comb_PbS, barplot_muthill_comb_PbS, layout = matrix(c(1,2,3,4), nrow=2, byrow=TRUE))

##also making plot w Imm nora, Mel nora, Prestney burn, Muthill, Grom and Motts Mill - for examining patterns of cross species infections
six_virus_hbars_prev_by_sp<-multiplot(barplot_Imm_Nora_comb_PbS,barplot_mel_nora_comb_PbS, barplot_prest_B_comb_PbS, barplot_muthill_comb_PbS, barplot_grom_PbS, barplot_motts_mill_PbS, layout = matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, byrow=TRUE))

#Plot of the Nora virus prevalences across whole dataset
#ImmNora
barplot_Imm_Nora_comb_PbS_a<-ggplot(Dimm_Nora_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19) ,plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Immigrans Nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

#mel nora
barplot_mel_nora_comb_PbS_a<-ggplot(mel_nora_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_blank(), axis.text.x = element_text(size = 19), plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Melanogaster Nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

nora_viruses_hbars_prev_by_sp<-multiplot(barplot_Imm_Nora_comb_PbS_a,barplot_mel_nora_comb_PbS_a,barplot_sub_nora_PbS, layout = matrix(c(1,2,3), nrow=1, byrow=TRUE))

##Plot of Motts mill, Prestney burn and Grom prevalence 
#prestney burn
barplot_prest_B_comb_PbS<-ggplot(prest_B_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3),text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19),plot.title = element_text(hjust = 0.5,size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  labs(title="Prestney burn virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

gr_mm_pb_hbars_prev_by_sp<-multiplot(barplot_prest_B_comb_PbS, barplot_grom_PbS, barplot_motts_mill_PbS, layout = matrix(c(1,2,3), nrow=1, byrow=TRUE))

#Plot of Chaq and Galbut virus
#Chaq
barplot_chaq_comb_PbS<-ggplot(Chaq_comb_obs3_imm_mel, aes(x=species,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  coord_flip()+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 19),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Chaq virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)

chaq_gal_hbars_prev_by_sp<-multiplot(barplot_chaq_comb_PbS, barplot_galbut_PbS, layout = matrix(c(1,2), nrow=1, byrow=TRUE))

# Quantifying Prestney burn presence/absence over time 

prest_B_over_time <- Mar_20_data_reduced %>%
  group_by(date,species,prest_B) %>%
  summarise(n_pools=length(prest_B),
            n_flies=sum(no_flies))
prest_B_over_time_positive <- subset(prest_B_over_time,prest_B==1)

##same for other species 
mel_nora_over_time <- Mar_20_data_reduced %>%
  group_by(date,species,mel_nora) %>%
  summarise(n_pools=length(mel_nora),
            n_flies=sum(no_flies))
mel_nora_over_time_positive <- subset(mel_nora_over_time,mel_nora==1)

muthill_over_time <- Mar_20_data_reduced %>%
  group_by(date,species,muthill) %>%
  summarise(n_pools=length(muthill),
            n_flies=sum(no_flies))
muthill_over_time_positive <- subset(muthill_over_time,muthill==1)

Dimm_nora_over_time <- Mar_20_data_reduced %>%
  group_by(date,species,Dimm_Nora) %>%
  summarise(n_pools=length(Dimm_Nora),
            n_flies=sum(no_flies))
Dimm_nora_over_time_positive <- subset(Dimm_nora_over_time,Dimm_Nora==1)


Dimm_SV_over_time <- Mar_20_data_reduced %>%
  group_by(date,species,Dimm_SV) %>%
  summarise(n_pools=length(Dimm_SV),
            n_flies=sum(no_flies))
Dimm_SV_over_time_positive <- subset(Dimm_SV_over_time,Dimm_SV==1)

##no of flies collected of each species
no_flies_by_sp <- Mar_20_data_reduced %>%
  group_by(species) %>%
  summarise(n_flies=sum(no_flies))


# Prestney Burn prevalence over time ---------------------------------------

prest_b_PoT<-Mar_20_data_reduced %>%
  group_by(date,species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=prest_B,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance)) %>% as.data.frame()

glimpse(prest_b_PoT)
tbl_df(prest_b_PoT)

prest_b_PoT$date<-as.Date(prest_b_PoT$date,format = "%d/%m/%Y")
prest_b_PoT$prevalence<-as.numeric(as.character(prest_b_PoT$prevalence))

#creating plot of prevalence over time
prest_b_Potplot<- prest_b_PoT %>% 
      select(prevalence,date,species) %>%
      filter(species %in% c("imm","sub","sil", "mel", "obs")) %>% 
      ggplot(aes(date, prevalence, color = species)) +
      geom_line(lwd=1.5)+
      ggtitle("Prevalence of Prestney burn over time\n") +
      labs(x=NULL,y=NULL) +
      scale_x_date(date_labels ="%b-%y",date_breaks='2 months') + xlab("") +
      theme(panel.background = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_line(color = "gray50", size = 0.5),
            panel.grid.major.x = element_blank(), 
            text = element_text(size=18),
            axis.text.y = element_text(), 
            axis.text.x = element_text(),
            plot.title = element_text(hjust = -0.05, vjust=2, size = 14))

# Co-infection? -----------------------------------------------------------

##counting no. of flies collected from each sp.
Mar_20_data_reduced %>%
  group_by(species) %>%
  summarise(no_flies_by_sp = sum(no_flies))

##Looking at the number of co-infections
##creating a number of infections col in the data_frame
Mar_20_data_reduced$n_virus<-rowSums(Mar_20_data_reduced[,7:17])

#looking at only single flies 
Mar_20_data_reduced_onlysingles <- Mar_20_data_reduced %>%
  filter(no_flies==1) 
#changing n_virus to a factor
Mar_20_data_reduced_onlysingles$n_virus<-as.factor(Mar_20_data_reduced_onlysingles$n_virus)

##summarising the single fly data
coinfection_summary <- Mar_20_data_reduced_onlysingles %>%
  group_by(n_virus) %>%
  summarise(no_flies = length(short_code))

##looking at coinfection data distribution 
ggplot(Mar_20_data_reduced_onlysingles, aes(x=n_virus)) +
  geom_bar(colour = "black",fill = "darkorange",width = 0.82)+
  xlab("No. of viral infections") +
  ylab("No. of single flies") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", size = 0.5),
        axis.text.y = element_text(size = 14,colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"),
        text = element_text(size = 15))+
  annotate(geom = "text",x=5.5,y=180,size = 6, label = "Total n = 416 single flies", color = "black")+
  annotate(geom = "text",x=1,y=123,size = 5, label = "113", color = "black")+
  annotate(geom = "text",x=2,y=179,size = 5, label = "169", color = "black")+
  annotate(geom = "text",x=3,y=99,size = 5, label = "89", color = "black")+
  annotate(geom = "text",x=4,y=38,size = 5, label = "28", color = "black")+
  annotate(geom = "text",x=5,y=24,size = 5, label = "14", color = "black")+
  annotate(geom = "text",x=6,y=13,size = 5, label = "3", color = "black")

#looking at it in a different way - how likely are you to be infected with anything...and then how likely to have mutiple infections 
#stacked horizontal plot
#adding col with single variable 
coinfection_summary$n_col<-as.factor("1")
#making n_virus a factor
coinfection_summary$n_virus<-as.factor(coinfection_summary$n_virus)
#adding column with the proportion of the single flies infected with each of the no_virus categories
coinfection_summary$prop<-coinfection_summary$no_flies/416

coinfection_proportion_graph <- coinfection_summary %>%
  arrange(n_virus) %>%
  mutate(n_virus = factor(n_virus,levels = c("5","4","3","2","1","0"))) %>%
ggplot(aes(x=n_col,y=prop,fill=n_virus))+
  geom_bar(position = "stack",stat = "identity",width = 0.7)+
  coord_flip()+
  scale_fill_manual(values = cbPalette)+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 24),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y  = element_blank(),
        text = element_text(size = 14),
        legend.text = element_text(size = 24),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"))

##examining any common coinfection combinations using a contingency table

#first coverting the table of viruses to long form...using tidyr
keycol <- "virus"
valuecol <- "PA"
gathercols <- c("chaq","Dimm_SV","Dimm_Nora","mel_nora","prest_B","muthill","tranent","grom","motts_mill","galbut","sub_nora")
Mar_20_data_reduced_onlysingles_long<-gather_(Mar_20_data_reduced_onlysingles, keycol, valuecol, gathercols)

##counting the number of (single) individuals infected with each virus
singlefly_n_infections<-Mar_20_data_reduced_onlysingles_long %>%
  group_by(virus,PA) %>%
  filter(PA==1) %>%
  summarise(n_infected = length(short_code))
#export this as a table
write.csv(singlefly_n_infections, file = "single_fly_n_infections.csv")

#then removing any rows where PA==0, as that means there is no infection
Mar_20_data_reduced_onlysingles_long_pos<-Mar_20_data_reduced_onlysingles_long %>%
  filter(PA==1)

#now counting the combinations
viral_combinations <-data.frame(t(combn(unique(Mar_20_data_reduced_onlysingles_long_pos$virus), 2)),stringsAsFactors = F) %>%
  mutate(counts = map2_dbl(X1, X2, ~length(intersect(Mar_20_data_reduced_onlysingles_long_pos$short_code[Mar_20_data_reduced_onlysingles_long_pos$virus==.x], 
                                                     Mar_20_data_reduced_onlysingles_long_pos$short_code[Mar_20_data_reduced_onlysingles_long_pos$virus==.y]))))

ordered_viral_combinations <- viral_combinations[with(viral_combinations, order(-counts)),]

#symetric
viral_combinations <- expand.grid(X1=unique(Mar_20_data_reduced_onlysingles_long_pos$virus),
            X2=unique(Mar_20_data_reduced_onlysingles_long_pos$virus), stringsAsFactors = F) %>%
  mutate(counts = map2_dbl(X1, X2, ~length(intersect(Mar_20_data_reduced_onlysingles_long_pos$short_code[Mar_20_data_reduced_onlysingles_long_pos$virus==.x], 
                                                     Mar_20_data_reduced_onlysingles_long_pos$short_code[Mar_20_data_reduced_onlysingles_long_pos$virus==.y])))) %>% 
  filter(X1 != X2) 

#creating a heat plot from the results 
ggplot(viral_combinations, aes(x = X1, y = X2)) + 
  geom_tile(aes(fill = counts),colour = "black") +
  scale_fill_gradient(low = "white", high = "#009E73") +
  theme(axis.ticks = element_blank(), 
        axis.text = element_text(),
        axis.title = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

###A second shot at analysing this using the cooccur package in R 
#supposed to be used for species and sites, but I'll use the individuals as 'sites' and viruses as species 

#trimming off any un-needless columns from the data frame so it only contains the short codes and the viruses
cooccur_input<-Mar_20_data_reduced_onlysingles %>%
  select(chaq,Dimm_SV,Dimm_Nora,mel_nora,prest_B,muthill,tranent,grom,motts_mill,galbut,sub_nora)
rownames(cooccur_input)<-Mar_20_data_reduced_onlysingles$short_code
colnames(cooccur_input)<-c("Chaq virus","Dimm sigma virus","Dimm nora virus","Dmel nora virus","Prestney burn virus","Muthill virus","Tranent virus","Grom virus","Motts mill virus","Galbut virus","Dsub nora virus")
cooccur_input<-as.matrix(cooccur_input)
cooccur_input<-t(cooccur_input)

#create coocurrance object
coocurrance_data<-cooccur(cooccur_input, spp_names = TRUE, thresh = FALSE)

#makes a heat map of the co-occurance matrix highlighting combinations with more or less co-occurances than expected.
plot(coocurrance_data)+theme(legend.key.size = unit(0.1, "cm"),
                             legend.text = element_text(size = 3),
                             text = element_text(size = 20, colour = "black"),
                             legend.background = element_rect(size = 1),
                             axis.text = element_text(size = 36))

#makes a table of the species combinations and the obs and exp probs of combinations, along with whether the observed co-occurances deviate significantly from the expected ones...
prob.table(coocurrance_data)


# Prevalence over time? ----------------------------------------------------

#Can we gain anything meaningful by looking at variation in virus prevalence over time 

#first looking at which viruses have the highest estimated prevalences 

#Prestney burn in Dmel + Dobs species?
#Muthill across all 5 most common species
#Chaq in Dmel
#Galbut in Dmel
#Imm nora in all 5 most common sp?

#First adding a col to Mar_20_data_reduced which denotes the 'season' of the collection - which is actually the groups I did the sequencing in (eg. SO16, DJ1617, JO17, AJ18, JO18)
Mar_20_data_reduced$season <- ifelse(Mar_20_data_reduced$date=='26/09/2016'|Mar_20_data_reduced$date=='23/09/2016'|Mar_20_data_reduced$date=='27/10/2016','SO16',ifelse(Mar_20_data_reduced$date=='19/12/2016'|Mar_20_data_reduced$date=='31/03/2017'|Mar_20_data_reduced$date=='28/04/2017'|Mar_20_data_reduced$date=='25/05/2017'|Mar_20_data_reduced$date=='13/06/2017','DJ1617',
                                            ifelse(Mar_20_data_reduced$date=='21/07/2017'|Mar_20_data_reduced$date=='31/08/2017'|Mar_20_data_reduced$date=='28/09/2017'|Mar_20_data_reduced$date=='18/10/2017','JO17',
                                                   ifelse(Mar_20_data_reduced$date=='23/04/2018'|Mar_20_data_reduced$date=='24/05/2018'|Mar_20_data_reduced$date=='22/06/2018','AJ18',
                                                          ifelse(Mar_20_data_reduced$date=='13/07/2018'|Mar_20_data_reduced$date=='21/08/2018'|Mar_20_data_reduced$date=='30/09/2018'|Mar_20_data_reduced$date=='20/10/2018','JO18','no season')
                                                          )
                                                   )
                                            )
                                     )
#Looking at prestney burn, grom and motts mill together
##Prestney burn
prest_B_comb_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=prest_B,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=prest_B,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=prest_B,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
prest_B_comb_obs3_imm_mel_sp_season<-subset(prest_B_comb_byspecies_season, 
                                  species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
prest_B_comb_obs3_imm_mel_sp_season$species<-ifelse(prest_B_comb_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                    ifelse(prest_B_comb_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                           ifelse(prest_B_comb_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                                  ifelse(prest_B_comb_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                         ifelse(prest_B_comb_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                                         )
                                                                  )
                                                           )
                                                    )

#plotting prestney burn prevalence across seasons 
#ordering the seasons factor levels 
barplot_prest_b_comb_PbS_season<-prest_B_comb_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                           levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Prestney Burn virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

##Grom
grom_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=grom,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=grom,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=grom,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
grom_obs3_imm_mel_sp_season<-subset(grom_byspecies_season, 
                                            species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
grom_obs3_imm_mel_sp_season$species<-ifelse(grom_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                    ifelse(grom_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                           ifelse(grom_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                                  ifelse(grom_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                         ifelse(grom_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                                  )
                                                           )
                                                    )
)

#plotting grom prevalence across seasons 
#ordering the seasons factor levels 
barplot_grom_byspecies_season<-grom_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Grom virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

##Motts Mill
motts_mill_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=motts_mill,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=motts_mill,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=motts_mill,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
motts_mill_obs3_imm_mel_sp_season<-subset(motts_mill_byspecies_season, 
                                    species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
motts_mill_obs3_imm_mel_sp_season$species<-ifelse(motts_mill_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                            ifelse(motts_mill_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                   ifelse(motts_mill_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                          ifelse(motts_mill_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                 ifelse(motts_mill_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                          )
                                                   )
                                            )
)

#plotting motts mill prevalence across seasons 
#ordering the seasons factor levels 
barplot_motts_mill_byspecies_season<-motts_mill_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Motts mill virus", x="", y="")+
  geom_errorbar(limits,position=dodge, width = 0.1)+
  facet_grid(~species)

#combining the plots for PB, Motts mill and Grom 
comb_PB_grom_motts_prev_by_season_sp<-multiplot(barplot_prest_b_comb_PbS_season,barplot_grom_byspecies_season,barplot_motts_mill_byspecies_season,layout = matrix(c(1,2,3),nrow=3,byrow = TRUE))

##trying out a wee model in MCMCglmm for incidence of prestney burn 
library(MCMCglmm)
binary_prior <- list(R=list(V=1, nu=1000), 
                     G=list(G1=list(V=1,nu=1000,alpha.mu=0, alpha.V=1)))
fit_m0_prest_b<-MCMCglmm(prest_B~season + species,
                         random =~site, 
                         data = Mar_20_data_reduced,
                         nitt = 1000100, thin = 100, 
                         prior = binary_prior)
plot(fit_m0_prest_b)
summary(fit_m0_prest_b$Sol)
summary(fit_m0_prest_b$VCV)

###Trying to look at Nora virus prevalences across season - again using the Mar_20_data_reduced dataset
##Imm Nora
Imm_nora_comb_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=Dimm_Nora,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=Dimm_Nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=Dimm_Nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
Imm_nora_comb_obs3_imm_mel_sp_season<-subset(Imm_nora_comb_byspecies_season, 
                                            species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
Imm_nora_comb_obs3_imm_mel_sp_season$species<-ifelse(Imm_nora_comb_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                    ifelse(Imm_nora_comb_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                           ifelse(Imm_nora_comb_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                                  ifelse(Imm_nora_comb_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                         ifelse(Imm_nora_comb_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                                  )
                                                           )
                                                    )
)

#Plotting Imm Nora prevalence across seasons
#ordering the seasons factor levels 
barplot_Imm_nora_comb_PbS_season<-Imm_nora_comb_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Immigrans nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

##Mel Nora
mel_nora_comb_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=mel_nora,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=mel_nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=mel_nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
mel_nora_comb_obs3_imm_mel_sp_season<-subset(mel_nora_comb_byspecies_season, 
                                             species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
mel_nora_comb_obs3_imm_mel_sp_season$species<-ifelse(mel_nora_comb_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                     ifelse(mel_nora_comb_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                            ifelse(mel_nora_comb_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                                   ifelse(mel_nora_comb_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                          ifelse(mel_nora_comb_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                                   )
                                                            )
                                                     )
)

#Plotting mel Nora prevalence across seasons
#ordering the seasons factor levels 
barplot_mel_nora_comb_PbS_season<-mel_nora_comb_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Melanogaster nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

#Dsub nora virus prevalence across seasons 
sub_nora_comb_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=sub_nora,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=sub_nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=sub_nora,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
sub_nora_comb_obs3_imm_mel_sp_season<-subset(sub_nora_comb_byspecies_season, 
                                             species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
sub_nora_comb_obs3_imm_mel_sp_season$species<-ifelse(sub_nora_comb_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                     ifelse(sub_nora_comb_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                            ifelse(sub_nora_comb_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                                   ifelse(sub_nora_comb_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                          ifelse(sub_nora_comb_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                                   )
                                                            )
                                                     )
)

#Plotting mel Nora prevalence across seasons
#ordering the seasons factor levels 
barplot_sub_nora_comb_PbS_season<-sub_nora_comb_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Subobscura nora virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

#combining the three Nora virus prevalence plots across seasons and common species 
comb_3noras_prev_by_season_sp<-multiplot(barplot_Imm_nora_comb_PbS_season,barplot_mel_nora_comb_PbS_season,barplot_sub_nora_comb_PbS_season,layout = matrix(c(1,2,3), nrow=3, byrow=TRUE))

###looking at Chaq and Galbut prevalence over seasons and in comparison to each other
#Chaq virus prevalence across seasons 
chaq_comb_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=chaq,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=chaq,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=chaq,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
chaq_comb_obs3_imm_mel_sp_season<-subset(chaq_comb_byspecies_season, 
                                             species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
chaq_comb_obs3_imm_mel_sp_season$species<-ifelse(chaq_comb_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                     ifelse(chaq_comb_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                            ifelse(chaq_comb_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                                   ifelse(chaq_comb_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                          ifelse(chaq_comb_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                                   )
                                                            )
                                                     )
)

#Plotting chaq prevalence across seasons
#ordering the seasons factor levels 
barplot_chaq_comb_PbS_season<-chaq_comb_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size = 17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Chaq virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

#Galbut virus prevalence across seasons 
galbut_comb_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=galbut,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=galbut,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=galbut,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
galbut_comb_obs3_imm_mel_sp_season<-subset(galbut_comb_byspecies_season, 
                                         species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
galbut_comb_obs3_imm_mel_sp_season$species<-ifelse(galbut_comb_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                 ifelse(galbut_comb_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                        ifelse(galbut_comb_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                               ifelse(galbut_comb_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                      ifelse(galbut_comb_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                               )
                                                        )
                                                 )
)

#Plotting galbut prevalence across seasons
#ordering the seasons factor levels 
barplot_galbut_comb_PbS_season<-galbut_comb_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size=17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Galbut virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

#combining the Chaq and Galbut virus prevalence plots across seasons and common species 
comb_chaq_gal_prev_by_season_sp<-multiplot(barplot_chaq_comb_PbS_season,barplot_galbut_comb_PbS_season,layout = matrix(c(1,2), nrow=2, byrow=TRUE))

#Looking at the prevalence of Muthill virus over time - across seasons
muthill_comb_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=muthill,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=muthill,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=muthill,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
muthill_comb_obs3_imm_mel_sp_season<-subset(muthill_comb_byspecies_season, 
                                           species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
muthill_comb_obs3_imm_mel_sp_season$species<-ifelse(muthill_comb_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                   ifelse(muthill_comb_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                          ifelse(muthill_comb_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                                 ifelse(muthill_comb_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                        ifelse(muthill_comb_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                                 )
                                                          )
                                                   )
)

#Plotting muthill prevalence across seasons
#ordering the seasons factor levels 
barplot_muthill_comb_PbS_season<-muthill_comb_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size=17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Muthill virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

#Imm SV prevalence across seasons
ImmSV_comb_byspecies_season<-Mar_20_data_reduced %>%
  group_by(species,season) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=Dimm_SV,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=Dimm_SV,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=Dimm_SV,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 

#separating the prevalence for Imm, mel, and obs group, and renaming the factor levels
ImmSV_comb_obs3_imm_mel_sp_season<-subset(ImmSV_comb_byspecies_season, 
                                            species == "mel"|species=="imm"|species == "obs"|species == "sub"|species == "sil")
ImmSV_comb_obs3_imm_mel_sp_season$species<-ifelse(ImmSV_comb_obs3_imm_mel_sp_season$species=='imm','Dimm',
                                                    ifelse(ImmSV_comb_obs3_imm_mel_sp_season$species=='mel','Dmel',
                                                           ifelse(ImmSV_comb_obs3_imm_mel_sp_season$species=='obs','Dobs',
                                                                  ifelse(ImmSV_comb_obs3_imm_mel_sp_season$species=='sub','Dsub',
                                                                         ifelse(ImmSV_comb_obs3_imm_mel_sp_season$species=='sil','Dsus','no species')
                                                                  )
                                                           )
                                                    )
)

#Plotting ImmSV prevalence across seasons
#ordering the seasons factor levels 
barplot_ImmSV_comb_PbS_season<-ImmSV_comb_obs3_imm_mel_sp_season %>%
  mutate(season = factor(season, 
                         levels = c('SO16','DJ1617','JO17','AJ18','JO18'))) %>%
  ggplot(aes(x=season,y=prevalence,fill=species)) +
  geom_bar(position="dodge",stat="identity",colour="black")+
  theme(panel.background = element_blank(),panel.grid.major = element_line(color = "lightgrey", size = 0.3), text = element_text(size=14),axis.text.y = element_text(size=17), axis.text.x = element_text(size = 16,angle = 30,hjust = 1),plot.title = element_text(hjust = 0.5, size = 19),legend.position = "none")+
  ylim(c(0.0,1.0))+
  scale_fill_manual(values = cbPalette)+
  labs(title="Dimm sigma virus", x="", y="")+
  geom_errorbar(limits,position =dodge, width = 0.1)+
  facet_grid(~species)

###Host range - viruses with P/A data
##table which counts no of indivs infected w each virus
##Summarising the prevalence and number of infections/number sampled across species 
n_chaq_by_species <- Mar_20_data_reduced %>% 
  group_by(species, chaq) %>%
  summarise(n_pools=length(chaq),
            n_flies=sum(no_flies))
subset(n_chaq_by_species, chaq==1)

n_Dimm_SV_by_species <- Mar_20_data_reduced %>% 
  group_by(species, Dimm_SV) %>%
  summarise(n_pools=length(Dimm_SV),
            n_flies=sum(no_flies))
subset(n_Dimm_SV_by_species, Dimm_SV==1)

n_Dimm_Nora_by_species <- Mar_20_data_reduced %>% 
  group_by(species, Dimm_Nora) %>%
  summarise(n_pools=length(Dimm_Nora),
            n_flies=sum(no_flies))
subset(n_Dimm_Nora_by_species, Dimm_Nora==1)

n_mel_nora_by_species <- Mar_20_data_reduced %>% 
  group_by(species, mel_nora) %>%
  summarise(n_pools=length(mel_nora),
            n_flies=sum(no_flies))
subset(n_mel_nora_by_species, mel_nora==1)

n_prest_B_by_species <- Mar_20_data_reduced %>% 
  group_by(species, prest_B) %>%
  summarise(n_pools=length(prest_B),
            n_flies=sum(no_flies))
subset(n_prest_B_by_species, prest_B==1)

n_muthill_by_species <- Mar_20_data_reduced %>% 
  group_by(species, muthill) %>%
  summarise(n_pools=length(muthill),
            n_flies=sum(no_flies))
subset(n_muthill_by_species, muthill==1)

subset(n_tranent_by_species, tranent==1)

n_grom_by_species <- Mar_20_data_reduced %>% 
  group_by(species, grom) %>%
  summarise(n_pools=length(grom),
            n_flies=sum(no_flies))
subset(n_grom_by_species, grom==1)

n_sub_nora_by_species <- Mar_20_data_reduced %>% 
  group_by(species, sub_nora) %>%
  summarise(n_pools=length(sub_nora),
            n_flies=sum(no_flies))
subset(n_sub_nora_by_species, sub_nora==1)

n_motts_mill_by_species <- Mar_20_data_reduced %>% 
  group_by(species,motts_mill) %>%
  summarise(n_pools=length(motts_mill),
            n_flies=sum(no_flies))
subset(n_motts_mill_by_species, motts_mill==1)

n_galbut_by_species <- Mar_20_data_reduced %>% 
  group_by(species, galbut) %>%
  summarise(n_pools=length(galbut),
            n_flies=sum(no_flies))
subset(n_galbut_by_species, galbut==1)


# Base R prevalence plots -------------------------------------------------
#creating index of which species have the largest total number of indivs collected
n_by_species<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(total_flies = sum(no_flies))
n_by_species$species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir")

#using ML to calculate prevalence from PA
tranent_comb_byspecies<-Mar_20_data_reduced %>%
  group_by(species) %>%
  summarize(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,
                                                    hits=tranent,
                                                    bounds=TRUE,
                                                    interval=2,test.consistant=TRUE,
                                                    plot=FALSE)$prevelance),
            lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=tranent,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[1]),
            upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,
                                                     hits=tranent,
                                                     bounds=TRUE,
                                                     interval=2,test.consistant=TRUE,
                                                     plot=FALSE)$bounds[2])) %>% 
  as.data.frame() 
tranent_comb_byspecies$species<-factor(c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir"), levels=c("Dimm","Dsub","Dobs","Dsus","Dmel","Dfun","Dpha","Dtri","Dhel","Dhyd","Sdef","Ccos","Dbus","Dvir","Hcam"))
#adding the total_species col to the dataframe
tranent_comb_byspecies <- merge(tranent_comb_byspecies,n_by_species,by="species")
#ordering the data by descending total no of flies
tranent_comb_byspecies<-tranent_comb_byspecies[order(-tranent_comb_byspecies$total_flies),]

data<-tranent_comb_byspecies
#using match to reorder the dataframe so that the error bars and indexes are in the correct order
#target <- c("Ccos","Dbus","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsub","Dsus","Dtri","Dvir","Hcam","Sdef")
#data[match(target, data$species),]->data
#converting all the vals to percentages from proportions
data[,2:4] <- data[,2:4]*100
#changing any value < 0.01% to 0.01%
index <-data[,2:4] < 0.01
data[,2:4][index] <- 0.01
prev_index<-index[,1]#for changing colours etc. of sp w prev < 0.01
lb_index<-index[,2]#for changing the lower bounds the error bars when prev is < 0.01 to 0.001
data[,3][lb_index]<-0.001#changing those lower bounds

# a bar or boxplot in base R for prevalence w a log scale
#making some fake data so that the whole range is encompassed
fake_species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir")
fake_prevalence<-as.numeric(seq(0.01,100,length.out = 15))
fake_data<-data.frame(fake_species,fake_prevalence)

#making plot
par(mar = c(5.1,5.1, 2.6, 2.1), # change the margins
    lwd = 1.7,# increase the line thickness
    cex.axis = 1.3, # increase default axis label size
    cex.lab = 1.4)
#initally making the empty plot using some fake data, so that the whole plot space needed is created
barplot(fake_data$fake_prevalence~fake_data$fake_species,axes = FALSE,xaxt='n',yaxt='n', border= NA, col = NA,main="", xlab="", ylab="",log = 'y') # invisible bars - only plot
#using a loop to add custom lines on the log scale
line_list<-c(0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
for(i in 1:length(line_list)){
  list<-line_list
  par(xpd = TRUE)
  lines(x = c(par("usr")[1],par("usr")[2]), y = c(list[i],list[i]),lty = 1,lwd=2,col="grey90")
}
par(xpd=FALSE)#REMEMBER TO RUN THIS AFTER
#adding y axis
axis(2, at = c(0.01,0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100), tick = TRUE, labels = c(0.01,0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100),lwd = 2,lwd.ticks = 2,las = 1,mgp=c(3, 0.75, 0))
#adding y axis label
mtext(side = 2, line = 4, "Virus prevalence (%)", cex = 1.6,padj = 0.6)
#adding x axis
vec<-seq(par("usr")[1]+0.74,par("usr")[2]-0.66,length.out = 15)
axis(1, at = vec,
     tick = FALSE,
     labels = FALSE)
#x coordinates for labels 
xtext_cols<-rep.int("black",15)
xtext_cols[prev_index]<-"grey60"
text(x = vec,
     y = par("usr")[4]-1.994,
     labels = c("Dimm","Dsub","Dobs","Dsus","Dmel","Dfun","Dpha","Dtri","Dhel","Dhyd","Sdef","Ccos","Dbus","Dvir","Hcam"),
     ## Rotate the labels by 35 degrees.
     xpd = NA,
     srt = 0,
     adj = 0.5,
     cex = 1.5,
     col = xtext_cols)
#adding x axis label
mtext(side = 1, line = 3,"Drosophila species", cex = 1.6)

#Setting the amount of space to leave before each bar
bar_spacing<-c(-0.273143,rep.int(0.273143,14)) 
#plotting actual data on to the plot
#making colour vector fr bars
xbar_cols<-rep.int("#009E73",15)
xbar_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bars
xborder_cols<-rep.int("black",15)
xborder_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bar borders
barplot(data$prevalence~data$species,add=TRUE,yaxt ="n",xaxt="n",log = "y",col=xbar_cols,space=bar_spacing,xpd=TRUE,border=xborder_cols)
#adding error bars
barCenters<-barplot(data$prevalence~data$species,plot=FALSE,yaxt ="n",xaxt="n",log = "y",space=bar_spacing)
arrows(barCenters,data$lower_bound,barCenters,data$upper_bound, lwd=1.7, angle=90, code=3,length = 0.1,col = xtext_cols)
#adding text to plot with virus name 
mtext(side=3,line=1,"Tranent Virus",cex=1.8,padj = 0.25)

# All viruses log10 prevalence plots loops ---------------------------------

###Loop which calculates a list of virus prevalence data tables across all species 
##prepping the data
comb_virus_data<-Mar_20_data_reduced

n_by_species<-comb_virus_data %>%
  group_by(species) %>%
  summarize(total_flies = sum(no_flies))
n_by_species$species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir")

viruses<-list(colnames(comb_virus_data[,7:17]))
prev_by_sp_list<-list()

#cycling through the viruses in comb_viruses and calculating a by species data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by species df for",viruses[[1]][i],sep = " "))
  prev_by_sp_list[[i]]<-comb_virus_data %>%
    group_by(species) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2])) %>% 
    as.data.frame() 
  
  prev_by_sp_list[[i]]$species<-factor(c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir"), levels=c("Dimm","Dsub","Dobs","Dsus","Dmel","Dfun","Dpha","Dtri","Dhel","Dhyd","Sdef","Ccos","Dbus","Dvir","Hcam"))
  #adding the total_species col to the dataframe
  prev_by_sp_list[[i]] <- merge(prev_by_sp_list[[i]],n_by_species,by="species")
  #ordering the data by descending total no of flies
  prev_by_sp_list[[i]]<-prev_by_sp_list[[i]][order(-prev_by_sp_list[[i]]$total_flies),]

  names(prev_by_sp_list)[i] <- viruses[[1]][i]
    
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence by species df created",sep = " "))
  
}

#once this loop has been run - check that you have a list of prevalence tables ordered by the number of each species collected
glimpse(prev_by_sp_list)

#Loop for creating a multipage plot series of prevalence by species for each virus

#Make list of viruses in the form you want the title to take...eg. capitalised
virus_titles<-c("Chaq","Dimm Sigma","Dimm Nora","Dmel Nora","Prestney Burn","Muthill","Tranent","Grom","Motts Mill","Galbut","Dsub Nora")
multi_virus_data<-prev_by_sp_list

for (i in 1:length(multi_virus_data)) {
  
  #isolating data for only one virus
  data<-multi_virus_data[[i]]
  
  #preparing the data
  #proportions -> percentages
  data[,2:4] <- data[,2:4]*100
  #changing any value < 0.01% to 0.01%
  index <-data[,2:4] < 0.01
  data[,2:4][index] <- 0.01
  prev_index<-index[,1] #for changing colours etc. of sp w prev < 0.01
  lb_index<-index[,2] #for changing the lower bounds the error bars when prev is < 0.01 to 0.001
  data[,3][lb_index]<-0.001 #changing those lower bounds
  
  ##if none of the non-zero lower bounds are less than 0.1...changing the limit of the graph to 0.1, if none are less than 1, changing to 1
  non_index<-data[,3:4]>0.01
  
  if (min(data[,3:4][non_index])>1) {
    
    data[,2][prev_index]<-1#change prevalence to 1% if <0.01 so that that't the bottom of the graph
    
    line_list<-c(1.5,2,3,5,7.5,10,15,20,30,50,75,100)
    
    #making some fake data so that the whole range is encompassed 
    fake_species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir")
    fake_prevalence<-as.numeric(seq(1,100,length.out = 15))
    fake_data<-data.frame(fake_species,fake_prevalence)
    
    y_axis<-c(1,1.5,2,3,5,7.5,10,15,20,30,50,75,100)
    x_text_y<-par("usr")[4]-1.25
  
  } else if (min(data[,3:4][non_index])>0.1) {
    
    data[,2][prev_index]<-0.1#change prevalence to 0.1% if <0.01 so that that't the bottom of the graph
    
    line_list<-c(0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir")
    fake_prevalence<-as.numeric(seq(0.1,100,length.out = 15))
    fake_data<-data.frame(fake_species,fake_prevalence)
    
    y_axis<-c(0.1,0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    x_text_y<-par("usr")[4]-1.93
    
  } else {
    
    data=data #keep the same, w bottom at 0.01%
    
    line_list<-c(0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir")
    fake_prevalence<-as.numeric(seq(0.01,100,length.out = 15))
    fake_data<-data.frame(fake_species,fake_prevalence)
    
    y_axis<-c(0.01,0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    x_text_y<-par("usr")[4]-1.994
  }
  
  pdf(file = paste("log_prevalence_by_species_plots/log_prevalence_by_species_", sub(" ","_",virus_titles[i]), ".pdf", sep = ""),width = 12,height = 6)
  #setting margins for plotting 
  par(mar = c(5.1,5.1, 2.6, 2.1), # change the margins
      lwd = 1.7,# increase the line thickness
      cex.axis = 1.3, # increase default axis label size
      cex.lab = 1.4)
  
  # starting to prepare the plot
  #initally making the empty plot using some fake data, so that the whole plot space needed is created
  barplot(fake_data$fake_prevalence~fake_data$fake_species,axes = FALSE,xaxt='n',yaxt='n', border= NA, col = NA,main="", xlab="", ylab="",log = 'y') # invisible bars - only plot
  #using a loop to add custom lines on the log scale
  for(j in 1:length(line_list)){
    list<-line_list
    par(xpd = TRUE)
    lines(x = c(par("usr")[1],par("usr")[2]), y = c(list[j],list[j]),lty = 1,lwd=2,col="grey90")
  }
  par(xpd=FALSE)#REMEMBER TO RUN THIS AFTER
  #adding y axis
  axis(2, at = y_axis, tick = TRUE, labels = y_axis, lwd = 2,lwd.ticks = 2,las = 1,mgp=c(3, 0.75, 0))
  #adding y axis label
  mtext(side = 2, line = 4, "Virus prevalence (%)", cex = 1.6,padj = 0.6)
  #adding x axis
  vec<-seq(par("usr")[1]+0.74,par("usr")[2]-0.66,length.out = 15)
  axis(1, at = vec,
       tick = FALSE,
       labels = FALSE)
  #x coordinates for labels 
  xtext_cols<-rep.int("black",15)
  xtext_cols[prev_index]<-"grey60"
  text(x = vec,
       y = x_text_y,
       labels = c("Dimm","Dsub","Dobs","Dsus","Dmel","Dfun","Dpha","Dtri","Dhel","Dhyd","Sdef","Ccos","Dbus","Dvir","Hcam"),
       ## Rotate the labels by 35 degrees.
       xpd = NA,
       srt = 0,
       adj = 0.5,
       cex = 1.5,
       col = xtext_cols)
  #adding x axis label
  mtext(side = 1, line = 3,"Drosophila species", cex = 1.6)
  #Setting the amount of space to leave before each bar
  bar_spacing<-c(-0.273143,rep.int(0.273143,14)) 
  
  #making colour vector fr bars
  xbar_cols<-rep.int("#009E73",15)
  xbar_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bars
  xborder_cols<-rep.int("black",15)
  xborder_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bar borders
  #plotting actual data on to the plot
  barplot(data$prevalence~data$species,add=TRUE,yaxt ="n",xaxt="n",log = "y",col=xbar_cols,space=bar_spacing,xpd=TRUE,border=xborder_cols)
  #adding error bars
  barCenters<-barplot(data$prevalence~data$species,plot=FALSE,yaxt ="n",xaxt="n",log = "y",space=bar_spacing)
  arrows(barCenters,data$lower_bound,barCenters,data$upper_bound, lwd=1.7, angle=90, code=3,length = 0.1,col = xtext_cols)
  #adding text to plot with virus name 
  mtext(side=3,line=1,paste(virus_titles[i],"Virus",sep = " "),cex=1.8,padj = 0.25)
  
  dev.off()
}

# Base R plots separated by season ----------------------------------------

Mar_20_data_reduced -> comb_virus_data
comb_virus_data$species<-revalue(comb_virus_data$species,c("imm"="Dimm","sub"="Dsub","obs"="Dobs","sil"="Dsus","mel"="Dmel","fun"="Dfun","pha"="Dpha","tri"="Dtri","hel"="Dhel","hyd"="Dhyd","def"="Sdef","chy"="Ccos","bus"="Dbus","vir"="Dvir","cam"="Hcam"))
#comb_virus_data$season<-as.factor(comb_virus_data$season)

##this time cycling through the dataframes and creating prevalence estimates separated by both species AND season 
viruses<-list(colnames(comb_virus_data[,7:17]))
prev_by_sp_season_list<-list()

#cycling through the viruses in comb_viruses and calculating a by species data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by species + season df for",viruses[[1]][i],sep = " "))
  prev_by_sp_season_list[[i]]<-comb_virus_data %>%
    dplyr::group_by(species, season) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2])) %>% 
    as.data.frame() 
  
  names(prev_by_sp_season_list)[i] <- viruses[[1]][i]
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence by species + season df created",sep = " "))
}

#once this loop has been run - check that you have a list of prevalence tables ordered by the number of each species collected
glimpse(prev_by_sp_season_list)

##Now producing plots of the log10 prevalence of viruses over the 5 seasons

#making some fake data so that the whole range is encompassed on the (empty) original plot
fake_species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir")
fake_prevalence<-as.numeric(seq(0.01,100,length.out = 15))
fake_data<-data.frame(fake_species,fake_prevalence)

#Make list of viruses in the form you want the title to take...eg. capitalised
virus_titles<-c("Chaq","Dimm Sigma","Dimm Nora","Dmel Nora","Prestney Burn","Muthill","Tranent","Grom","Motts Mill","Galbut","Dsub Nora")
multi_virus_data_sea<-prev_by_sp_season_list

#formatting and subsetting data
##could I separate into years? - eg. three groups
data<-multi_virus_data_sea[[3]]

#same manipulations of data as previous loop eg. percentages

#preparing the data
#proportions -> percentages
data[,3:5] <- data[,3:5]*100
#changing any value < 0.01% to 0.01%
index <-data[,3:5] < 0.01
data[,3:5][index] <- 0.01
prev_index<-index[,1]#for changing colours etc. of sp w prev < 0.01
lb_index<-index[,2]#for changing the lower bounds the error bars when prev is < 0.01 to 0.001
data[,4][lb_index]<-0.001#changing those lower bounds

data[data$season=="SO16",]->SO16_dat
data[data$season=="DJ1617" | data$season=="JO17",]-> y17_dat
data[data$season=="AJ18" | data$season=="JO18",]->y18_dat

barplot(SO16_dat$prevalence~SO16_dat$species)

##TO BE CONTINUED...


# Implementing bivariate model of prevalence on data ----------------------


