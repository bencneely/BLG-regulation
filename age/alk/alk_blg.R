#--------------------------------------------------------------
#Ben Neely
#05/07/2024
#Create ALKs to assign ages to unaged fish in each sample
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("FSA" %in% rownames(installed.packages()) == FALSE) {install.packages("FSA")}
library(FSA)

if("rio" %in% rownames(installed.packages()) == FALSE) {install.packages("rio")}
library(rio)

if("nnet" %in% rownames(installed.packages()) == FALSE) {install.packages("nnet")}
library(nnet)

if("lubridate" %in% rownames(installed.packages()) == FALSE) {install.packages("lubridate")}
library(lubridate)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=22,color="black",face="bold"),
        axis.text=element_text(size=18,color="black"),
        legend.position="none")
options(scipen=999)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/BLG regulation eval/age/")

## Read in pre-regulation data and separate into aged and unaged
pre=import("alk/precatch_raw.csv")%>%
  filter(spp=="Bluegill",
         tl>=50)%>%
  mutate(cmgrp=lencat(tl,10))%>%
  select(impd,year,spp,cmgrp,age)

pre_aged=pre%>%
  filter(age>=0)

pre_unaged=pre%>%
  filter(is.na(age))

## Read in pre-regulation data and separate into aged and unaged
post=import("alk/postcatch_raw.csv")%>%
  filter(spp=="Bluegill",
         tl>=50)%>%
  mutate(cmgrp=lencat(tl,10))%>%
  select(impd,year,spp,cmgrp,age)

post_aged=post%>%
  filter(age>=0)

post_unaged=post%>%
  filter(is.na(age))

## Set seed for reproducible results
set.seed(57)

################################################################################
################################################################################
## Create data for pre-regulation

##############################################################
## GDCL
## Create age-length key from aged fish data using multinomial logistic regression model
GDCL_mlr_pre=multinom(age~cmgrp,data=subset(pre_aged,impd=="GDCL"))

## Predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
GDCL_alk_pre=predict(GDCL_mlr_pre,data.frame(cmgrp=lens),type="probs")
row.names(GDCL_alk_pre)=lens

## Apply ALK to unaged data set and combine with aged fish
GDCL_pre=alkIndivAge(GDCL_alk_pre,age~cmgrp,data=subset(pre_unaged,impd=="GDCL"))%>%
  bind_rows(subset(pre_aged,impd=="GDCL"))%>%
  mutate(trt="con",
         period="pre")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## HOCB
## Create age-length key from aged fish data using multinomial logistic regression model
HOCB_mlr_pre=multinom(age~cmgrp,data=subset(pre_aged,impd=="HOCB"))

## Predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
HOCB_alk_pre=predict(HOCB_mlr_pre,data.frame(cmgrp=lens),type="probs")
row.names(HOCB_alk_pre)=lens

## Apply ALK to unaged data set and combine with aged fish
HOCB_pre=alkIndivAge(HOCB_alk_pre,age~cmgrp,data=subset(pre_unaged,impd=="HOCB"))%>%
  bind_rows(subset(pre_aged,impd=="HOCB"))%>%
  mutate(trt="con",
         period="pre")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## JWSL
## Create age-length key from aged fish data using multinomial logistic regression model
JWSL_mlr_pre=multinom(age~cmgrp,data=subset(pre_aged,impd=="JWSL"))

## Predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
JWSL_alk_pre=predict(JWSL_mlr_pre,data.frame(cmgrp=lens),type="probs")
row.names(JWSL_alk_pre)=lens

## Apply ALK to unaged data set and combine with aged fish
JWSL_pre=alkIndivAge(JWSL_alk_pre,age~cmgrp,data=subset(pre_unaged,impd=="JWSL"))%>%
  bind_rows(subset(pre_aged,impd=="JWSL"))%>%
  mutate(trt="exp",
         period="pre")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## LXLX
## Create age-length key from aged fish data using multinomial logistic regression model
LXLX_mlr_pre=multinom(age~cmgrp,data=subset(pre_aged,impd=="LXLX"))

## Predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
LXLX_alk_pre=predict(LXLX_mlr_pre,data.frame(cmgrp=lens),type="probs")
row.names(LXLX_alk_pre)=lens

## Apply ALK to unaged data set and combine with aged fish
LXLX_pre=alkIndivAge(LXLX_alk_pre,age~cmgrp,data=subset(pre_unaged,impd=="LXLX"))%>%
  bind_rows(subset(pre_aged,impd=="LXLX"))%>%
  mutate(trt="exp",
         period="pre")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## MISL
## Create age-length key from aged fish data using multinomial logistic regression model
MISL_mlr_pre=multinom(age~cmgrp,data=subset(pre_aged,impd=="MISL"))

## Predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
MISL_alk_pre=predict(MISL_mlr_pre,data.frame(cmgrp=lens),type="probs")
row.names(MISL_alk_pre)=lens

## Apply ALK to unaged data set and combine with aged fish
MISL_pre=alkIndivAge(MISL_alk_pre,age~cmgrp,data=subset(pre_unaged,impd=="MISL"))%>%
  bind_rows(subset(pre_aged,impd=="MISL"))%>%
  mutate(trt="exp",
         period="pre")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## POCL
## Create age-length key from aged fish data using multinomial logistic regression model
POCL_mlr_pre=multinom(age~cmgrp,data=subset(pre_aged,impd=="POCL"))

## Predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
POCL_alk_pre=predict(POCL_mlr_pre,data.frame(cmgrp=lens),type="probs")
row.names(POCL_alk_pre)=lens

## Apply ALK to unaged data set and combine with aged fish
POCL_pre=alkIndivAge(POCL_alk_pre,age~cmgrp,data=subset(pre_unaged,impd=="POCL"))%>%
  bind_rows(subset(pre_aged,impd=="POCL"))%>%
  mutate(trt="con",
         period="pre")%>%
  select(impd,spp,trt,period,age,cmgrp)
##############################################################
## PTSB
## Create age-length key from aged fish data using multinomial logistic regression model
PTSB_mlr_pre=multinom(age~cmgrp,data=subset(pre_aged,impd=="PTSB"))

## Predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
PTSB_alk_pre=predict(PTSB_mlr_pre,data.frame(cmgrp=lens),type="probs")
row.names(PTSB_alk_pre)=lens

## Apply ALK to unaged data set and combine with aged fish
PTSB_pre=alkIndivAge(PTSB_alk_pre,age~cmgrp,data=subset(pre_unaged,impd=="PTSB"))%>%
  bind_rows(subset(pre_aged,impd=="PTSB"))%>%
  mutate(trt="exp",
         period="pre")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## ZPCB
## Create age-length key from aged fish data using multinomial logistic regression model
ZPCB_mlr_pre=multinom(age~cmgrp,data=subset(pre_aged,impd=="ZPCB"))

## Predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
ZPCB_alk_pre=predict(ZPCB_mlr_pre,data.frame(cmgrp=lens),type="probs")
row.names(ZPCB_alk_pre)=lens

## Apply ALK to unaged data set and combine with aged fish
ZPCB_pre=alkIndivAge(ZPCB_alk_pre,age~cmgrp,data=subset(pre_unaged,impd=="ZPCB"))%>%
  bind_rows(subset(pre_aged,impd=="ZPCB"))%>%
  mutate(trt="con",
         period="pre")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
##############################################################
## Combine pre-regulation data
allaged_pre=bind_rows(GDCL_pre,HOCB_pre,JWSL_pre,LXLX_pre,
                      MISL_pre,POCL_pre,PTSB_pre,ZPCB_pre)

################################################################################
################################################################################
## Create data for post-regulation

##############################################################
## GDCL
## Create age-length key from aged fish data using multinomial logistic regression model
GDCL_mlr_post=multinom(age~cmgrp,data=subset(post_aged,impd=="GDCL"))

## predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
GDCL_alk_post=predict(GDCL_mlr_post,data.frame(cmgrp=lens),type="probs")
row.names(GDCL_alk_post)=lens

## Apply ALK to unaged data set and combine with aged fish
GDCL_post=alkIndivAge(GDCL_alk_post,age~cmgrp,data=subset(post_unaged,impd=="GDCL"))%>%
  bind_rows(subset(post_aged,impd=="GDCL"))%>%
  mutate(trt="con",
         period="post")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## HOCB
## Create age-length key from aged fish data using multinomial logistic regression model
HOCB_mlr_post=multinom(age~cmgrp,data=subset(post_aged,impd=="HOCB"))

## predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
HOCB_alk_post=predict(HOCB_mlr_post,data.frame(cmgrp=lens),type="probs")
row.names(HOCB_alk_post)=lens

## Apply ALK to unaged data set and combine with aged fish
HOCB_post=alkIndivAge(HOCB_alk_post,age~cmgrp,data=subset(post_unaged,impd=="HOCB"))%>%
  bind_rows(subset(post_aged,impd=="HOCB"))%>%
  mutate(trt="con",
         period="post")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## JWSL
## Create age-length key from aged fish data using multinomial logistic regression model
JWSL_mlr_post=multinom(age~cmgrp,data=subset(post_aged,impd=="JWSL"))

## predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
JWSL_alk_post=predict(JWSL_mlr_post,data.frame(cmgrp=lens),type="probs")
row.names(JWSL_alk_post)=lens

## Apply ALK to unaged data set and combine with aged fish
JWSL_post=alkIndivAge(JWSL_alk_post,age~cmgrp,data=subset(post_unaged,impd=="JWSL"))%>%
  bind_rows(subset(post_aged,impd=="JWSL"))%>%
  mutate(trt="exp",
         period="post")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## LXLX
## Create age-length key from aged fish data using multinomial logistic regression model
LXLX_mlr_post=multinom(age~cmgrp,data=subset(post_aged,impd=="LXLX"))

## predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
LXLX_alk_post=predict(LXLX_mlr_post,data.frame(cmgrp=lens),type="probs")
row.names(LXLX_alk_post)=lens

## Apply ALK to unaged data set and combine with aged fish
LXLX_post=alkIndivAge(LXLX_alk_post,age~cmgrp,data=subset(post_unaged,impd=="LXLX"))%>%
  bind_rows(subset(post_aged,impd=="LXLX"))%>%
  mutate(trt="exp",
         period="post")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## MISL
## Create age-length key from aged fish data using multinomial logistic regression model
MISL_mlr_post=multinom(age~cmgrp,data=subset(post_aged,impd=="MISL"))

## predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
MISL_alk_post=predict(MISL_mlr_post,data.frame(cmgrp=lens),type="probs")
row.names(MISL_alk_post)=lens

## Apply ALK to unaged data set and combine with aged fish
MISL_post=alkIndivAge(MISL_alk_post,age~cmgrp,data=subset(post_unaged,impd=="MISL"))%>%
  bind_rows(subset(post_aged,impd=="MISL"))%>%
  mutate(trt="exp",
         period="post")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## POCL
## Create age-length key from aged fish data using multinomial logistic regression model
POCL_mlr_post=multinom(age~cmgrp,data=subset(post_aged,impd=="POCL"))

## predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
POCL_alk_post=predict(POCL_mlr_post,data.frame(cmgrp=lens),type="probs")
row.names(POCL_alk_post)=lens

## Apply ALK to unaged data set and combine with aged fish
POCL_post=alkIndivAge(POCL_alk_post,age~cmgrp,data=subset(post_unaged,impd=="POCL"))%>%
  bind_rows(subset(post_aged,impd=="POCL"))%>%
  mutate(trt="con",
         period="post")%>%
  select(impd,spp,trt,period,age,cmgrp)
##############################################################
## PTSB
## Create age-length key from aged fish data using multinomial logistic regression model
PTSB_mlr_post=multinom(age~cmgrp,data=subset(post_aged,impd=="PTSB"))

## predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
PTSB_alk_post=predict(PTSB_mlr_post,data.frame(cmgrp=lens),type="probs")
row.names(PTSB_alk_post)=lens

## Apply ALK to unaged data set and combine with aged fish
PTSB_post=alkIndivAge(PTSB_alk_post,age~cmgrp,data=subset(post_unaged,impd=="PTSB"))%>%
  bind_rows(subset(post_aged,impd=="PTSB"))%>%
  mutate(trt="exp",
         period="post")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
## ZPCB
## Create age-length key from aged fish data using multinomial logistic regression model
ZPCB_mlr_post=multinom(age~cmgrp,data=subset(post_aged,impd=="ZPCB"))

## predict probability that a fish in a given cm group is a certain age
lens=seq(50,200,10)
ZPCB_alk_post=predict(ZPCB_mlr_post,data.frame(cmgrp=lens),type="probs")
row.names(ZPCB_alk_post)=lens

## Apply ALK to unaged data set and combine with aged fish
ZPCB_post=alkIndivAge(ZPCB_alk_post,age~cmgrp,data=subset(post_unaged,impd=="ZPCB"))%>%
  bind_rows(subset(post_aged,impd=="ZPCB"))%>%
  mutate(trt="con",
         period="post")%>%
  select(impd,spp,trt,period,age,cmgrp)

##############################################################
##############################################################
## Combine post-regulation data
allaged_post=bind_rows(GDCL_post,HOCB_post,JWSL_post,LXLX_post,
                      MISL_post,POCL_post,PTSB_post,ZPCB_post)

##############################################################
##############################################################
## Combine and export all aging data
allaged=bind_rows(allaged_pre,allaged_post)
export(allaged,"alk/allaged_blg.csv")
