#--------------------------------------------------------------
#Ben Neely
#05/20/2024
#Summary statistics for Bluegill paper
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

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/BLG regulation eval/summary statistics/")

## Read in data with import
regfish=import("regfish.csv")
regsamp=import("regsamp.csv")
preage=import("precatch_raw.csv")
postage=import("postcatch_raw.csv")

## How many fish did we sample
## Clean up data a bit
regfish1=regfish%>%
  filter(gear=="EF",
         spp=="Bluegill")%>%
  expandCounts(~count)%>%
  mutate(period=case_when(year<2019 ~ "pre",
                          year>=2022 ~ "post",
                          TRUE ~ "other"),
         period=factor(period,levels=c("pre","post")))%>%
  select(impd,year,spp,tl,period)


tmp=xtabs(~impd+period,regfish1)%>%
  as_tibble%>%
  mutate(trt=case_when(impd=="LXLX"|impd=="PTSB"|impd=="JWSL"|impd=="MISL" ~ "regulation",
                      impd=="GDCL"|impd=="HOCB"|impd=="ZPCB"|impd=="POCL" ~ "control",
                      TRUE ~ "other"))
sum(tmp$n)
sum(subset(tmp,trt=="control")$n)
sum(subset(tmp,trt=="regulation")$n)
sum(subset(tmp,period=="pre")$n)
sum(subset(tmp,period=="post")$n)

## How much effort did we expend
## Clean up data a bit
regsamp1=regsamp%>%
  filter(gear=="EF")%>%
  mutate(period=case_when(year<2019 ~ "pre",
                          year>=2022 ~ "post",
                          TRUE ~ "other"),
         period=factor(period,levels=c("pre","post")))%>%
  select(impd,year,period,effort)


tmp=xtabs(effort~impd+period,regsamp1)%>%
  as_tibble%>%
  mutate(trt=case_when(impd=="LXLX"|impd=="PTSB"|impd=="JWSL"|impd=="MISL" ~ "regulation",
                       impd=="GDCL"|impd=="HOCB"|impd=="ZPCB"|impd=="POCL" ~ "control",
                       TRUE ~ "other"))
sum(tmp$n)
sum(subset(tmp,trt=="control")$n)
sum(subset(tmp,trt=="regulation")$n)
sum(subset(tmp,period=="pre")$n)
sum(subset(tmp,period=="post")$n)

## How many fish did we age - pre
preage1=preage%>%
  filter(gear=="EF",
         spp=="Bluegill")%>%
  mutate(period=case_when(year<2019 ~ "pre",
                          year>=2022 ~ "post",
                          TRUE ~ "other"),
         period=factor(period,levels=c("pre","post")))%>%
  select(impd,year,spp,tl,aged,period)
  
tmp=xtabs(aged~impd+period,preage1)%>%
  as_tibble%>%
  mutate(trt=case_when(impd=="LXLX"|impd=="PTSB"|impd=="JWSL"|impd=="MISL" ~ "regulation",
                       impd=="GDCL"|impd=="HOCB"|impd=="ZPCB"|impd=="POCL" ~ "control",
                       TRUE ~ "other"))

sum(tmp$n)
sum(subset(tmp,trt=="control")$n)
sum(subset(tmp,trt=="regulation")$n)
sum(subset(tmp,period=="pre")$n)
sum(subset(tmp,period=="post")$n)

## How many fish did we age - post
postage1=postage%>%
  filter(gear=="EF",
         spp=="Bluegill")%>%
  mutate(period=case_when(year<2019 ~ "pre",
                          year>=2022 ~ "post",
                          TRUE ~ "other"),
         period=factor(period,levels=c("pre","post")))%>%
  select(impd,year,spp,tl,aged,period)

tmp=xtabs(aged~impd+period,postage1)%>%
  as_tibble%>%
  mutate(trt=case_when(impd=="LXLX"|impd=="PTSB"|impd=="JWSL"|impd=="MISL" ~ "regulation",
                       impd=="GDCL"|impd=="HOCB"|impd=="ZPCB"|impd=="POCL" ~ "control",
                       TRUE ~ "other"))

sum(tmp$n)
sum(subset(tmp,trt=="control")$n)
sum(subset(tmp,trt=="regulation")$n)
sum(subset(tmp,period=="pre")$n)
sum(subset(tmp,period=="post")$n)

## How many fish did we age with ALK- pre
preage2=preage%>%
  filter(gear=="EF",
         spp=="Bluegill",
         year==2017 | year==2018,
         aged==0)
  
tmp=xtabs(~impd+year,preage2)%>%
  as_tibble%>%
  mutate(trt=case_when(impd=="LXLX"|impd=="PTSB"|impd=="JWSL"|impd=="MISL" ~ "regulation",
                       impd=="GDCL"|impd=="HOCB"|impd=="ZPCB"|impd=="POCL" ~ "control",
                       TRUE ~ "other"))

sum(tmp$n)
sum(subset(tmp,trt=="control")$n)
sum(subset(tmp,trt=="regulation")$n)


## How many fish did we age with ALK- post
postage2=postage%>%
  filter(gear=="EF",
         spp=="Bluegill",
         year==2022,
         aged==0)

tmp=xtabs(~impd+year,postage2)%>%
  as_tibble%>%
  mutate(trt=case_when(impd=="LXLX"|impd=="PTSB"|impd=="JWSL"|impd=="MISL" ~ "regulation",
                       impd=="GDCL"|impd=="HOCB"|impd=="ZPCB"|impd=="POCL" ~ "control",
                       TRUE ~ "other"))

sum(tmp$n)
sum(subset(tmp,trt=="control")$n)
sum(subset(tmp,trt=="regulation")$n)
