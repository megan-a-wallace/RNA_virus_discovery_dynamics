######################################
##RNA virus prevalence calculations ##
######################################

##Megan A. Wallace
##March 2021

##installing packages
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
library(lubridate)
library(vegan)

##Prevalence function 

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

##Colour pallettes

# Colour blind friendly palette (7 values):
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Longer pallette (10 values)
cbPalette_long <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")

# Importing data ----------------------------------------------------------

#now updated w species from CO1 barcodes (apart from single mystery one...which I've put as virilis as its probs in the group but not sure what sp)
indiv_PCR_virus_data<-read.table("Indiv_virus_PCR_results.csv",header = TRUE, sep = ",")

indiv_PCR_virus_data$short_code<-as.character(indiv_PCR_virus_data$short_code)
indiv_PCR_virus_data$trap_date<-as.Date(indiv_PCR_virus_data$trap_date,format = "%d/%m/%y")
indiv_PCR_virus_data$date<-as.Date(indiv_PCR_virus_data$date,format = "%d/%m/%y")
indiv_PCR_virus_data$month_year<-as.factor(as.character(indiv_PCR_virus_data$month_year))
indiv_PCR_virus_data$month<-factor(indiv_PCR_virus_data$month, levels = c("mar","apr","may","jun","jul","aug","sep","oct","dec"))
indiv_PCR_virus_data$year<-factor(as.character(indiv_PCR_virus_data$year),levels = c("2016","2017","2018"))
indiv_PCR_virus_data$site<-as.factor(as.character(gsub(" ", "", indiv_PCR_virus_data$site)))
indiv_PCR_virus_data$collect_time<-as.numeric(hm(indiv_PCR_virus_data$collect_time))/60 #I've converted this to minutes 
indiv_PCR_virus_data$trap_wind<-as.numeric(indiv_PCR_virus_data$trap_wind)
indiv_PCR_virus_data$trap_humidity<-as.numeric(indiv_PCR_virus_data$trap_humidity)
indiv_PCR_virus_data$collect_humidity<-as.numeric(indiv_PCR_virus_data$collect_humidity)
indiv_PCR_virus_data$rain_vol<-as.numeric(indiv_PCR_virus_data$rain_vol)
indiv_PCR_virus_data$trap_duration<-as.numeric(indiv_PCR_virus_data$trap_duration)
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

##Need to deal with the fact that I don't have collection time or rain volume data for the two months of collections - the sep and oct 16

##Reduced table with results from Chaq, ImmSV, Imm Nora, Prestney burn, Mel Nora, Tranent, Grom, Motts mill, Galbut, Dsub Nora + Muthill (results in 03_20)
#as of 4_20 updated w new sp identifications - now 15 species
Mar_20_data<-data.frame(indiv_PCR_virus_data$short_code,
                        indiv_PCR_virus_data$date,
                        indiv_PCR_virus_data$month,
                        indiv_PCR_virus_data$month_year,
                        indiv_PCR_virus_data$year,
                        indiv_PCR_virus_data$site,
                        indiv_PCR_virus_data$lat,
                        indiv_PCR_virus_data$lon,
                        indiv_PCR_virus_data$collect_time,
                        indiv_PCR_virus_data$trap_temp,
                        indiv_PCR_virus_data$trap_wind,
                        indiv_PCR_virus_data$trap_humidity,
                        indiv_PCR_virus_data$collect_temp,
                        indiv_PCR_virus_data$collect_wind,
                        indiv_PCR_virus_data$collect_humidity,
                        indiv_PCR_virus_data$rain_vol,
                        indiv_PCR_virus_data$trap_duration,
                        indiv_PCR_virus_data$species,
                        indiv_PCR_virus_data$no_flies,
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
colnames(Mar_20_data)<-c("short_code","date","month","month_year","year","site","lat","lon","collect_time","trap_temp","trap_wind","trap_humidity","collect_temp","collect_wind","collect_humidity","rain_vol","trap_duration","species","no_flies","chaq","chaq_sh","Dimm_SV_L","Dimm_SV_N","Dimm_Nora_A","Dimm_Nora_B","mel_nora_6220","mel_nora_qPCR","prest_B_B","prest_B_sh","muthill_sh","muthill_suz","tranent_L","tranent_M","tranent_S","grom_728","motts_mill_221","galbut_407","Dsub_nora_1665")

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

##At the moment, just using the average of the temp, wind and humidity at trap laying and collection for these values, later, I'll use the data from the temperature measures in the traps
Mar_20_data$temp<-rowMeans(Mar_20_data[,grep("_temp", colnames(Mar_20_data))],na.rm = TRUE)
Mar_20_data$wind<-rowMeans(Mar_20_data[,grep("_wind", colnames(Mar_20_data))],na.rm = TRUE)
Mar_20_data$humidity<-rowMeans(Mar_20_data[,grep("_humidity", colnames(Mar_20_data))],na.rm = TRUE)

##Reducing table to chaq, DimmSV, Dimm Nora, mel nora, Prestney burn, Muthill, tranent, grom, motts mill, galbut and Dsub nora combined cols
Mar_20_data_reduced<-data.frame(Mar_20_data$short_code,Mar_20_data$date, Mar_20_data$month, Mar_20_data$month_year, Mar_20_data$year,Mar_20_data$site,Mar_20_data$lat,Mar_20_data$lon,Mar_20_data$collect_time,Mar_20_data$temp,Mar_20_data$wind,Mar_20_data$humidity,Mar_20_data$rain_vol,Mar_20_data$trap_duration,Mar_20_data$species,Mar_20_data$no_flies,Mar_20_data$chaq_comb,Mar_20_data$Dimm_SV_comb,Mar_20_data$Dimm_Nora_comb,Mar_20_data$mel_nora_comb,Mar_20_data$prest_B_comb,Mar_20_data$muthill_comb,Mar_20_data$tranent_comb,Mar_20_data$grom_728,Mar_20_data$motts_mill_221,Mar_20_data$galbut_407,Mar_20_data$Dsub_nora_1665) 

colnames(Mar_20_data_reduced)<-c("short_code","date","month","month_year","year","site","lat","lon","collect_time","temp","wind","humidity","rain_vol","trap_duration","species","no_flies","chaq_virus","Dimm_SV_virus","Dimm_Nora_virus","mel_nora_virus","prest_B_virus","muthill_virus","tranent_virus","grom_virus","motts_mill_virus","galbut_virus","sub_nora_virus")

##Adding an extra col specifying the sequencing pool
Mar_20_data_reduced$seq_pool<-ifelse(Mar_20_data_reduced$month_year %in% c('sep_16','oct_16'),'SO16',
                                     ifelse(Mar_20_data_reduced$month_year %in% c('dec_16','mar_17','apr_17','may_17','jun_17'),'DJ1617',
                                            ifelse(Mar_20_data_reduced$month_year %in% c('jul_17','aug_17','sep_17','oct_17'),'JO17',
                                                   ifelse(Mar_20_data_reduced$month_year %in% c('apr_18','may_18','jun_18'),'AJ18',
                                                          ifelse(Mar_20_data_reduced$month_year %in% c('jul_18','aug_18','sep_18','oct_18'),'JO18','not sequenced')
                                                   )
                                            )
                                     )
)

#Creating a table of the number of flies in each sequencing pool, grouped by species - for summary table in sequencing chapter
seq_sp_pools<-Mar_20_data_reduced %>% group_by(seq_pool, species) %>%
  summarise(nflies = sum(no_flies))
##and one just grouped by species
sp_pools<-Mar_20_data_reduced %>% group_by(species) %>%
  summarise(nflies = sum(no_flies))
##Creating the season covariate to model early and late drosophila season virus prevalence 
#In this system I define early as Mar-Jun, and late as Jul-Dec (but not sure whether July should be early or late...)
Mar_20_data_reduced$season<-ifelse(Mar_20_data_reduced$month %in% c('mar','apr','may','jun'),'early',
                                   ifelse(Mar_20_data_reduced$month %in% c('jul','aug','sep','oct','nov','dec'),'late','out_of_season')
)

##For the purposes of this analysis, the pilot site BB is the same as BR - as they are so close to one another...
levels(Mar_20_data_reduced$site) <- c("BR","BG", "BN", "BR", "CE", "CK", "CR", "CS", "DK", "EC", "EL", "GF", "HL", "IB", "IN", "LW", "OX", "SI", "TN", "TR", "VG")

##Adding in altitude of each site to the data 
sampling_site_data<-read.csv(file = "../../../../../thesis/ch2_metagenomic_virus_discovery/sampling_site_coords.csv", header = TRUE)
Mar_20_data_reduced$alt <- sampling_site_data$elevation[match(Mar_20_data_reduced$site, sampling_site_data$code)]

# Global prevalence of each virus ----

###Loop which calculates a list of virus prevalence across the whole dataset
comb_virus_data<-Mar_20_data_reduced

viruses<-list(colnames(comb_virus_data[,17:27]))
prev_list<-list()
##to sum the log likelihoods with no separation by any co-variate
LL_global_vec<-as.numeric(vector(length = length(viruses[[1]])))

#cycling through the viruses in comb_viruses and calculating a by species data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("calculating prevalence for",viruses[[1]][i],sep = " "))
  prev_list[[i]]<-comb_virus_data %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2]),
              log_likelihood = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$log.liklihood)) %>% 
    as.data.frame() 
  
  names(prev_list)[i] <- viruses[[1]][i]
  
  LL_global_vec[i] <- prev_list[[i]]$log_likelihood
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence calculated",sep = " "))
  
}

#total LL for global virus prevalence
sum(LL_global_vec)->total_LL_global_prev

# All viruses by species log10 prevalence plots loops ---------------------------------

###Loop which calculates a list of virus prevalence data tables across all species 
##prepping the data
comb_virus_data<-Mar_20_data_reduced

n_by_species<-comb_virus_data %>%
  group_by(species) %>%
  summarize(total_flies = sum(no_flies))
n_by_species$species<-c("Dbus","Hcam","Ccos","Sdef","Dfun","Dhel","Dhyd","Dimm","Dmel","Dobs","Dpha","Dsus","Dsub","Dtri","Dvir")

viruses<-list(colnames(comb_virus_data[,17:27]))
prev_by_sp_list<-list()

#cycling through the viruses in comb_viruses and calculating a by species data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by species df for",viruses[[1]][i],sep = " "))
  prev_by_sp_list[[i]]<-comb_virus_data %>%
    group_by(species) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2]),
              log_likelihood = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$log.liklihood)) %>%  
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

#And now summing the log likelihoods for each virus divided by species to compare to the overall likelihood without species included. 

LL_table<-as.data.frame(cbind(virus_titles,LL_global_vec,as.numeric(vector(length = 11)),as.numeric(vector(length = 11)),as.numeric(vector(length = 11))))
colnames(LL_table)<-c("virus","global","by_species","by_season","by_site")
LL_table$global<-as.numeric(LL_table$global)

for (i in 1:length(multi_virus_data)) {
  
  #isolating data for only one virus
  data<-multi_virus_data[[i]]
  
  #summing the LLs across separate species
  sum(data$log_likelihood)->virus_sp_LL
  
  #inserting into the LL summary table
  LL_table$by_species[i]<-as.numeric(virus_sp_LL)
}

###Loop which calculates a list of virus prevalence data tables across all seasons and species - for LL calculations
##prepping the data
comb_virus_data<-Mar_20_data_reduced

viruses<-list(colnames(comb_virus_data[,7:17]))
prev_by_sp_season_list<-list()

#cycling through the viruses in comb_viruses and calculating a by species data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by species df for",viruses[[1]][i],sep = " "))
  prev_by_sp_season_list[[i]]<-comb_virus_data %>%
    group_by(species,season) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2]),
              log_likelihood = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$log.liklihood)) %>%  
    as.data.frame()
  
  names(prev_by_sp_season_list)[i] <- viruses[[1]][i]
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence by species and season df created",sep = " "))
  
}

multi_virus_data<-prev_by_sp_season_list

for (i in 1:length(multi_virus_data)) {
  
  #isolating data for only one virus
  data<-multi_virus_data[[i]]
  
  #summing the LLs across separate species
  sum(data$log_likelihood)->virus_sp_season_LL
  
  #inserting into the LL summary table
  LL_table$by_season[i]<-as.numeric(virus_sp_season_LL)
}

###Loop which calculates a list of virus prevalence data tables across all sites, seasons and species - for LL calculations
##prepping the data
comb_virus_data<-Mar_20_data_reduced

viruses<-list(colnames(comb_virus_data[,7:17]))
prev_by_sp_season_site_list<-list()

#cycling through the viruses in comb_viruses and calculating a by species, site, season data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by species, site and season df for",viruses[[1]][i],sep = " "))
  prev_by_sp_season_site_list[[i]]<-comb_virus_data %>%
    group_by(species,season,site) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2]),
              log_likelihood = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$log.liklihood)) %>%  
    as.data.frame()
  
  names(prev_by_sp_season_site_list)[i] <- viruses[[1]][i]
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence by species, site and season df created",sep = " "))
  
}

multi_virus_data<-prev_by_sp_season_site_list

for (i in 1:length(multi_virus_data)) {
  
  #isolating data for only one virus
  data<-multi_virus_data[[i]]
  
  #summing the LLs across separate species
  sum(data$log_likelihood)->virus_sp_season_site_LL
  
  #inserting into the LL summary table
  LL_table$by_site[i]<-as.numeric(virus_sp_season_site_LL)
}

##LRT for each virus global vs. by species, site and season

LRT_results<-data.frame(virus = character(length = 11),twochLL_sp = numeric(length = 11),p_sp = numeric(length = 11),twochLL_sea = numeric(length = 11),p_sea= numeric(length = 11), twochLL_site = numeric(length = 11), p_site = numeric(length = 11))

for (j in 1:length(LL_table$virus)){
  
  virus<-LL_table$virus[j]
  
  as.numeric(LL_table$global[j])->global_LL
  as.numeric(LL_table$by_species[j])->sp_LL
  as.numeric(LL_table$by_season[j])->sea_LL
  as.numeric(LL_table$by_site[j])->site_LL
  
  #comparing global w species
  teststat_sp <- -2 * (global_LL-sp_LL)
  p.val_sp <- pchisq(teststat_sp, df = 1, lower.tail = FALSE)
  
  #comparing species w species + season
  teststat_sea <- -2 * (sp_LL-sea_LL)
  p.val_sea <- pchisq(teststat_sea, df = 1, lower.tail = FALSE)
  
  #comparing species + season w species + season + site
  teststat_site <- -2 * (sea_LL-site_LL)
  p.val_site <- pchisq(teststat_site, df = 1, lower.tail = FALSE)
  
  LRT_results[j,]$virus<-virus
  LRT_results[j,]$twochLL_sp<-teststat_sp
  LRT_results[j,]$p_sp<-p.val_sp
  LRT_results[j,]$twochLL_sea<-teststat_sea
  LRT_results[j,]$p_sea<-p.val_sea
  LRT_results[j,]$twochLL_site<-teststat_site
  LRT_results[j,]$p_site<-p.val_site
}

write.table(LRT_results,file = "../../../../../thesis/ch2_metagenomic_virus_discovery/LRT_prevalence_results.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
# All viruses by month log10 prevalence plots loops ---------------------------------

###Loop which calculates a list of virus prevalence data tables across all species 
##prepping the data
comb_virus_data<-Mar_20_data_reduced

n_by_month<-comb_virus_data %>%
  group_by(month) %>%
  summarize(total_flies = sum(no_flies))
n_by_month$month<-c("March","April","May","June","July","August","September","October","December")

viruses<-list(colnames(comb_virus_data[,17:27]))
prev_by_month_list<-list()

#cycling through the viruses in comb_viruses and calculating a by month data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by month df for",viruses[[1]][i],sep = " "))
  prev_by_month_list[[i]]<-comb_virus_data %>%
    group_by(month) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2])) %>% 
    as.data.frame() 
  
  prev_by_month_list[[i]]$month<-factor(c("March","April","May","June","July","August","September","October","December"), levels=c("March","April","May","June","July","August","September","October","December"))
  #adding the total_species col to the dataframe
  prev_by_month_list[[i]] <- merge(prev_by_month_list[[i]],n_by_month,by="month")
  #ordering the data by chronological months
  prev_by_month_list[[i]]<-prev_by_month_list[[i]][order(prev_by_month_list[[i]]$month),]
  
  names(prev_by_month_list)[i] <- viruses[[1]][i]
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence by species df created",sep = " "))
}

#once this loop has been run - check that you have a list of prevalence tables ordered by the number of each species collected
glimpse(prev_by_month_list)

#Loop for creating a multipage plot series of prevalence by month for each virus

#Make list of viruses in the form you want the title to take...eg. capitalised
virus_titles<-c("Chaq satellite","Dimm Sigma","Dimm Nora","Dmel Nora","Prestney Burn","Muthill","Tranent","Grom","Motts Mill","Galbut","Dsub Nora")
multi_virus_data<-prev_by_month_list

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
    fake_months<-c("March","April","May","June","July","August","September","October","December")
    fake_prevalence<-as.numeric(seq(1,100,length.out = 9))
    fake_data<-data.frame(fake_months,fake_prevalence)
    
    y_axis<-c(1,1.5,2,3,5,7.5,10,15,20,30,50,75,100)
    x_text_y<-par("usr")[4]-1.25
    
  } else if (min(data[,3:4][non_index])>0.1) {
    
    data[,2][prev_index]<-0.1#change prevalence to 0.1% if <0.01 so that that't the bottom of the graph
    
    line_list<-c(0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_months<-c("March","April","May","June","July","August","September","October","December")
    fake_prevalence<-as.numeric(seq(0.1,100,length.out = 9))
    fake_data<-data.frame(fake_months,fake_prevalence)
    
    y_axis<-c(0.1,0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    x_text_y<-par("usr")[4]-1.93
    
  } else {
    
    data=data #keep the same, w bottom at 0.01%
    
    line_list<-c(0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_months<-c("March","April","May","June","July","August","September","October","December")
    fake_prevalence<-as.numeric(seq(0.01,100,length.out = 9))
    fake_data<-data.frame(fake_months,fake_prevalence)
    
    y_axis<-c(0.01,0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    x_text_y<-par("usr")[4]-1.994
  }
  
  pdf(file = paste("log_prevalence_by_month_plots/log_prevalence_by_month_", sub(" ","_",virus_titles[i]), ".pdf", sep = ""),width = 12,height = 6)
  #setting margins for plotting 
  par(mar = c(5.1,5.1, 2.6, 2.1), # change the margins
      lwd = 1.7,# increase the line thickness
      cex.axis = 1.3, # increase default axis label size
      cex.lab = 1.4)
  
  # starting to prepare the plot
  #initally making the empty plot using some fake data, so that the whole plot space needed is created
  barplot(fake_data$fake_prevalence~fake_data$fake_months,axes = FALSE,xaxt='n',yaxt='n', border= NA, col = NA,main="", xlab="", ylab="",log = 'y') # invisible bars - only plot
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
  vec<-seq(par("usr")[1]+0.54,par("usr")[2]-0.748,length.out = 9)
  axis(1, at = vec,
       tick = FALSE,
       labels = FALSE)
  #x coordinates for labels 
  xtext_cols<-rep.int("black",9)
  xtext_cols[prev_index]<-"grey60"
  text(x = vec,
       y = x_text_y,
       labels = c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Dec"),
       ## Rotate the labels by 35 degrees.
       xpd = NA,
       srt = 0,
       adj = 0.5,
       cex = 1.5,
       col = xtext_cols)
  #adding x axis label
  mtext(side = 1, line = 3,"Months", cex = 1.6)
  #Setting the amount of space to leave before each bar
  bar_spacing<-c(-0.184,rep.int(0.27,8)) 
  
  #making colour vector fr bars
  xbar_cols<-rep.int("#CC79A7",9)
  xbar_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bars
  xborder_cols<-rep.int("black",9)
  xborder_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bar borders
  #plotting actual data on to the plot
  barplot(data$prevalence~data$month,add=TRUE,yaxt ="n",xaxt="n",log = "y",col=xbar_cols,xpd=TRUE,border=xborder_cols,space=bar_spacing)
  #
  #adding error bars
  barCenters<-barplot(data$prevalence~data$month,plot=FALSE,yaxt ="n",xaxt="n",log = "y",space=bar_spacing)
  #space=bar_spacing
  arrows(barCenters,data$lower_bound,barCenters,data$upper_bound, lwd=1.7, angle=90, code=3,length = 0.1,col = xtext_cols)
  #adding text to plot with virus name 
  mtext(side=3,line=1,paste(virus_titles[i],"Virus",sep = " "),cex=1.8,padj = 0.25)
  
  dev.off()
}

###Virus by site log 10 prevalence -----

###Loop which calculates a list of virus prevalence data tables across all species 
##prepping the data
comb_virus_data<-Mar_20_data_reduced

n_by_site<-comb_virus_data %>%
  group_by(site) %>%
  summarize(total_flies = sum(no_flies))
n_by_site$site<-c("BR","BG", "BN", "CE", "CK", "CR", "CS", "DK", "EC", "EL", "GF", "HL", "IB", "IN", "LW", "OX", "SI", "TN", "TR", "VG")

viruses<-list(colnames(comb_virus_data[,17:27]))
prev_by_site_list<-list()

#cycling through the viruses in comb_viruses and calculating a by month data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by site df for",viruses[[1]][i],sep = " "))
  prev_by_site_list[[i]]<-comb_virus_data %>%
    group_by(site) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=no_flies,hits=get(viruses[[1]][i]),bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2])) %>% 
    as.data.frame() 
  
  prev_by_site_list[[i]]$site<-factor(c("BR","BG", "BN", "CE", "CK", "CR", "CS", "DK", "EC", "EL", "GF", "HL", "IB", "IN", "LW", "OX", "SI", "TN", "TR", "VG"), levels=c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","DK","IN","VG","OX","CK","TN","EC","GF","EL"))
  #adding the total_species col to the dataframe
  prev_by_site_list[[i]] <- merge(prev_by_site_list[[i]],sampling_site_data,by.x = "site", by.y = "code")
  #ordering the data by longitude
  prev_by_site_list[[i]]<-prev_by_site_list[[i]][order(prev_by_site_list[[i]]$lon),]
  
  names(prev_by_site_list)[i] <- viruses[[1]][i]
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence by site df created",sep = " "))
}

#once this loop has been run - check that you have a list of prevalence tables ordered by the number of each species collected
glimpse(prev_by_site_list)

#Loop for creating a multipage plot series of prevalence by month for each virus

#Make list of viruses in the form you want the title to take...eg. capitalised
virus_titles<-c("Chaq satellite","Dimm Sigma","Dimm Nora","Dmel Nora","Prestney Burn","Muthill","Tranent","Grom","Motts Mill","Galbut","Dsub Nora")
multi_virus_data<-prev_by_site_list

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
    fake_sites<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","DK","IN","VG","OX","CK","TN","EC","GF","EL")
    fake_prevalence<-as.numeric(seq(1,100,length.out = 20))
    fake_data<-data.frame(fake_sites,fake_prevalence)
    
    y_axis<-c(1,1.5,2,3,5,7.5,10,15,20,30,50,75,100)
    x_text_y<-par("usr")[4]-1.24
    
  } else if (min(data[,3:4][non_index])>0.1) {
    
    data[,2][prev_index]<-0.1#change prevalence to 0.1% if <0.01 so that that't the bottom of the graph
    
    line_list<-c(0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_sites<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","DK","IN","VG","OX","CK","TN","EC","GF","EL")
    fake_prevalence<-as.numeric(seq(0.1,100,length.out = 20))
    fake_data<-data.frame(fake_sites,fake_prevalence)
    
    y_axis<-c(0.1,0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    x_text_y<-par("usr")[4]-1.93
    
  } else {
    
    data=data #keep the same, w bottom at 0.01%
    
    line_list<-c(0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_sites<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","DK","IN","VG","OX","CK","TN","EC","GF","EL")
    fake_prevalence<-as.numeric(seq(0.01,100,length.out = 20))
    fake_data<-data.frame(fake_sites,fake_prevalence)
    
    y_axis<-c(0.01,0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    x_text_y<-par("usr")[4]-1.994
  }
  
  pdf(file = paste("log_prevalence_by_site_plots/log_prevalence_by_site_", sub(" ","_",virus_titles[i]), ".pdf", sep = ""),width = 12,height = 6)
  #setting margins for plotting 
  par(mar = c(7.1,5.1, 2.6, 2.1), # change the margins
      lwd = 1.7,# increase the line thickness
      cex.axis = 1.3, # increase default axis label size
      cex.lab = 1.4)
  
  # starting to prepare the plot
  #initally making the empty plot using some fake data, so that the whole plot space needed is created
  barplot(fake_data$fake_prevalence~fake_data$fake_sites,axes = FALSE,xaxt='n',yaxt='n', border= NA, col = NA,main="", xlab="", ylab="",log = 'y') # invisible bars - only plot
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
  vec<-seq(par("usr")[1]+0.54,par("usr")[2]-0.748,length.out = 20)
  axis(1, at = vec,
       tick = FALSE,
       labels = FALSE)
  #x coordinates for labels 
  xtext_cols<-rep.int("black",20)
  xtext_cols[prev_index]<-"grey60"
  text(x = vec,
       y = x_text_y,
       labels = c("Crammond","Sighthill","Corstorphine Hill","Braidburn","Botanic gardens","Bush Estate","Bangholm","Hillend","Burdiehouse burn","Lasswade","Gilmerton","Dalkeith","Inveresk","Vogrie","Oxenfoord","Cockenzie","Tranent","Haddington","Gosford","Elvingston"),
       ## Rotate the labels by 35 degrees.
       xpd = NA,
       srt = 30,
       adj = 1,
       cex = 1.5,
       col = xtext_cols)
  #adding x axis label
  mtext(side = 1, line = 5,"Sites", cex = 1.6)
  #Setting the amount of space to leave before each bar
  bar_spacing<-c(-0.22,rep.int(0.27,19)) 
  
  #making colour vector fr bars
  xbar_cols<-rep.int("#0072B2",20)
  xbar_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bars
  xborder_cols<-rep.int("black",20)
  xborder_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bar borders
  #plotting actual data on to the plot
  barplot(data$prevalence~data$site,add=TRUE,yaxt ="n",xaxt="n",log = "y",col=xbar_cols,xpd=TRUE,border=xborder_cols,space=bar_spacing)
  #
  #adding error bars
  barCenters<-barplot(data$prevalence~data$site,plot=FALSE,yaxt ="n",xaxt="n",log = "y",space=bar_spacing)
  #space=bar_spacing
  arrows(barCenters,data$lower_bound,barCenters,data$upper_bound, lwd=1.7, angle=90, code=3,length = 0.1,col = xtext_cols)
  #adding text to plot with virus name 
  mtext(side=3,line=1,paste(virus_titles[i],"Virus",sep = " "),cex=1.8,padj = 0.25)
  
  dev.off()
}

### Confection ----

#looking at only single flies, 416 observations
Mar_20_data_reduced_onlysingles <- Mar_20_data_reduced %>%
  filter(no_flies==1) 
#Summarizing the number of infections 
Mar_20_data_reduced_onlysingles$n_virus<-as.factor(rowSums(Mar_20_data_reduced_onlysingles[,17:27]))

coinfection_summary <- Mar_20_data_reduced_onlysingles %>%
  group_by(n_virus) %>%
  summarise(no_flies = length(short_code))

##looking at coinfection data distribution 
ggplot(Mar_20_data_reduced_onlysingles, aes(x=n_virus)) +
  geom_bar(colour = "black",fill = "#0072B2",width = 0.82)+
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
  geom_bar(position = "stack",stat = "identity",width = 1)+
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
gathercols<-paste(gathercols,"_virus",sep = "")
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

#trimming off any needless columns from the data frame so it only contains the short codes and the viruses
cooccur_input<-Mar_20_data_reduced_onlysingles %>%
  select(gathercols)
rownames(cooccur_input)<-Mar_20_data_reduced_onlysingles$short_code
colnames(cooccur_input)<-c("Chaq satellite","Dimm sigma virus","Dimm nora virus","Dmel nora virus","Prestney burn virus","Muthill virus","Tranent virus","Grom virus","Motts mill virus","Galbut virus","Dsub nora virus")
cooccur_input<-as.matrix(cooccur_input)
cooccur_input<-t(cooccur_input)

#create coocurrance object
coocurrance_data<-cooccur(cooccur_input, spp_names = TRUE, thresh = FALSE)

#makes a heat map of the co-occurance matrix highlighting combinations with more or less co-occurances than expected.
plot(coocurrance_data)+theme(legend.key.size = unit(0.1, "cm"),
                             legend.text = element_text(size = 3),
                             text = element_text(size = 24, colour = "black"),
                             legend.background = element_rect(size = 1),
                             axis.text = element_text(size = 36))

#makes a table of the species combinations and the obs and exp probs of combinations, along with whether the observed co-occurances deviate significantly from the expected ones...
prob.table(coocurrance_data)->prob_cooccorence_table

##which virus combinations do you see less of than expected? - four combinations 
lt_sig_cooccur<-prob_cooccorence_table %>% filter(p_lt < 0.05)
##which virus combinations do you see more of than expected? - seven combinations (one of which is chaq and galbut)
gt_sig_cooccur<-prob_cooccorence_table %>% filter(p_gt < 0.05)

##combining these tables and exporting them 
rbind(lt_sig_cooccur,gt_sig_cooccur)->sig_coocur_combs

write.table(sig_coocur_combs,file = "sig_cooccur_virus_combs.csv",quote = FALSE,sep = ",",row.names = FALSE,col.names = TRUE)

##Community Diversity -----

flies_by_site<-Mar_20_data_reduced %>% 
  group_by(site) %>% 
  summarise(no_of_flies = sum(no_flies),
            prop_of_flies = (sum(no_flies))/2227)
flies_by_site<-flies_by_site[order(flies_by_site$prop_of_flies),]#increasing order by proportion

##% of collections in 17 and 18 which were from jul-oct, excluded 16 as didn't sample whole year
flies_by_month<-Mar_20_data_reduced %>% 
  group_by(year,month) %>% 
  filter(year == "2017" | year == "2018") %>% 
  summarise(no_of_flies = sum(no_flies),
            prop_of_flies = (sum(no_flies))/1903)
sum(flies_by_month$prop_of_flies[c(5:8,12:15)])

##species composition
flies_by_sp<-Mar_20_data_reduced %>% 
  group_by(species) %>% 
  summarise(no_of_flies = sum(no_flies),
            prop_of_flies = (sum(no_flies))/2227)
flies_by_sp<-flies_by_sp[order(flies_by_sp$prop_of_flies),]

##species composition by season
flies_by_sp_season<-Mar_20_data_reduced %>% 
  group_by(species,season) %>% 
  summarise(no_of_flies = sum(no_flies))

##no of flies in early season
sum(flies_by_sp_season[flies_by_sp_season$season == "early",]$no_of_flies)
(68+62+5+136)/(sum(flies_by_sp_season[flies_by_sp_season$season == "early",]$no_of_flies)) # obs group percentage

##no. of flies in late season
sum(flies_by_sp_season[flies_by_sp_season$season == "late",]$no_of_flies)
(1240)/(sum(flies_by_sp_season[flies_by_sp_season$season == "late",]$no_of_flies)) # imm percentage

###Shannon diversity calculations 

##late season shannon diversity
diversity(flies_by_sp_season[flies_by_sp_season$season == "late",]$no_of_flies, index = "shannon")
##early season shannon diversity
diversity(flies_by_sp_season[flies_by_sp_season$season == "early",]$no_of_flies, index = "shannon")

##effective number of species 
exp(diversity(flies_by_sp_season[flies_by_sp_season$season == "early",]$no_of_flies, index = "shannon"))
exp(diversity(flies_by_sp_season[flies_by_sp_season$season == "late",]$no_of_flies, index = "shannon"))
