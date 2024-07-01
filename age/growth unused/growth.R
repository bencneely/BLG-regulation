#--------------------------------------------------------------
#Ben Neely
#05/20/2024
#Bluegill growth
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

if("patchwork" %in% rownames(installed.packages()) == FALSE) {install.packages("patchwork")}
library(patchwork)

if("nlstools" %in% rownames(installed.packages()) == FALSE) {install.packages("nlstools")}
library(nlstools)

if("car" %in% rownames(installed.packages()) == FALSE) {install.packages("car")}
library(car)

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=18,color="black",face="bold"),
        axis.text=element_text(size=14,color="black"),
        legend.position=c(1,0),
        legend.justification=c("right","bottom"),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.background=element_rect(fill="transparent"),
        legend.margin=margin(c(0,1,1,1)))
options(scipen=999)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/BLG regulation eval/age/growth unused")

## Read in data with import
dat=import("allaged.csv")%>%
  mutate(period=factor(period,
                       levels=c("pre","post")))

blg=filter(dat,spp=="Bluegill")

##################################################################
##################################################################
##################################################################
## Set up vb prediction functions for each population during each time period
## Pre
ZPCB_pre_ages=seq(0,max(subset(blg,impd=="ZPCB")$age),by=0.05)
ZPCB_pre_vbpred=function(x) predict(x,data.frame(age=ZPCB_pre_ages))

GDCL_pre_ages=seq(0,max(subset(blg,impd=="GDCL")$age),by=0.05)
GDCL_pre_vbpred=function(x) predict(x,data.frame(age=GDCL_pre_ages))

HOCB_pre_ages=seq(0,max(subset(blg,impd=="HOCB")$age),by=0.05)
HOCB_pre_vbpred=function(x) predict(x,data.frame(age=HOCB_pre_ages))

POCL_pre_ages=seq(0,max(subset(blg,impd=="POCL")$age),by=0.05)
POCL_pre_vbpred=function(x) predict(x,data.frame(age=POCL_pre_ages))
#
JWSL_pre_ages=seq(0,max(subset(blg,impd=="JWSL")$age),by=0.05)
JWSL_pre_vbpred=function(x) predict(x,data.frame(age=JWSL_pre_ages))

LXLX_pre_ages=seq(0,max(subset(blg,impd=="LXLX")$age),by=0.05)
LXLX_pre_vbpred=function(x) predict(x,data.frame(age=LXLX_pre_ages))

MISL_pre_ages=seq(0,max(subset(blg,impd=="MISL")$age),by=0.05)
MISL_pre_vbpred=function(x) predict(x,data.frame(age=MISL_pre_ages))

PTSB_pre_ages=seq(0,max(subset(blg,impd=="PTSB")$age),by=0.05)
PTSB_pre_vbpred=function(x) predict(x,data.frame(age=PTSB_pre_ages))

## Post
ZPCB_post_ages=seq(0,max(subset(blg,impd=="ZPCB")$age),by=0.05)
ZPCB_post_vbpred=function(x) predict(x,data.frame(age=ZPCB_post_ages))

GDCL_post_ages=seq(0,max(subset(blg,impd=="GDCL")$age),by=0.05)
GDCL_post_vbpred=function(x) predict(x,data.frame(age=GDCL_post_ages))

HOCB_post_ages=seq(0,max(subset(blg,impd=="HOCB")$age),by=0.05)
HOCB_post_vbpred=function(x) predict(x,data.frame(age=HOCB_post_ages))

POCL_post_ages=seq(0,max(subset(blg,impd=="POCL")$age),by=0.05)
POCL_post_vbpred=function(x) predict(x,data.frame(age=POCL_post_ages))
#
JWSL_post_ages=seq(0,max(subset(blg,impd=="JWSL")$age),by=0.05)
JWSL_post_vbpred=function(x) predict(x,data.frame(age=JWSL_post_ages))

LXLX_post_ages=seq(0,max(subset(blg,impd=="LXLX")$age),by=0.05)
LXLX_post_vbpred=function(x) predict(x,data.frame(age=LXLX_post_ages))

MISL_post_ages=seq(0,max(subset(blg,impd=="MISL")$age),by=0.05)
MISL_post_vbpred=function(x) predict(x,data.frame(age=MISL_post_ages))

PTSB_post_ages=seq(0,max(subset(blg,impd=="PTSB")$age),by=0.05)
PTSB_post_vbpred=function(x) predict(x,data.frame(age=PTSB_post_ages))

## Define growth model function and set starting values
vb=vbFuns()
vb_sv=list(Linf=200,k=0.2,t0=-1.5)

##################################################################
## Fit models and calculate confidence intervals
##################################################################
## ZPCB
################
## pre
ZPCB_pre_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="ZPCB" & period=="pre"),start=vb_sv)
ZPCB_pre_vb_boot=nlsBoot(ZPCB_pre_vb_fit,niter=400)
(ZPCB_pre_vb_parms=cbind(ests=coef(ZPCB_pre_vb_fit),confint(ZPCB_pre_vb_boot)))

## Use model to predict length at age for plotting
ZPCB_pre_predboot=Boot(ZPCB_pre_vb_fit,f=ZPCB_pre_vbpred)
ZPCB_pre_preds=data.frame("ZPCB",
                          "pre",
                          ZPCB_pre_ages,
                          ZPCB_pre_vbpred(ZPCB_pre_vb_fit),
                          confint(ZPCB_pre_predboot))
names(ZPCB_pre_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="ZPCB" & period=="pre"))

################
## post
ZPCB_post_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="ZPCB" & period=="post"),start=vb_sv)
ZPCB_post_vb_boot=nlsBoot(ZPCB_post_vb_fit,niter=400)
(ZPCB_post_vb_parms=cbind(ests=coef(ZPCB_post_vb_fit),confint(ZPCB_post_vb_boot)))

## Use model to predict length at age for plotting
ZPCB_post_predboot=Boot(ZPCB_post_vb_fit,f=ZPCB_post_vbpred)
ZPCB_post_preds=data.frame("ZPCB",
                           "post",
                           ZPCB_post_ages,
                           ZPCB_post_vbpred(ZPCB_post_vb_fit),
                           confint(ZPCB_post_predboot))
names(ZPCB_post_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="ZPCB" & period=="post"))

################
## plot
ZPCB_preds=bind_rows(ZPCB_pre_preds,ZPCB_post_preds)%>%
  mutate(period=factor(period,levels=c("pre","post")))
ZPCB_points=subset(blg,impd=="ZPCB")%>%
  mutate(period=factor(period,levels=c("pre","post")))

ZPCB_plot=ggplot()+
  geom_line(ZPCB_preds,mapping=aes(x=age,y=tl,color=period))+
  geom_ribbon(ZPCB_preds,mapping=aes(x=age,ymin=lci95,ymax=uci95,fill=period),alpha=0.5)+
  geom_jitter(ZPCB_points,mapping=aes(x=age,y=cmgrp,color=period,shape=period),size=3,width=0.2,alpha=0.35)+
  scale_color_manual(values=c("#ff474c","#475f94"),
                     labels=c("Pre","Post"))+
  scale_fill_manual(values=c("#ff474c","#475f94"),
                    labels=c("Pre","Post"))+
  scale_shape_manual(values=c(17,19),
                   labels=c("Pre","Post"))+
  scale_x_continuous(breaks=seq(0,8,1),
                   name="Estimated age")+
  scale_y_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  annotate("text",label="Campbell",x=0,y=220,hjust=0,vjust=1,size=6)+
  coord_cartesian(xlim=c(-0.25,8.25),
                  ylim=c(-1,225),
                  expand=F)+
  pubtheme

##################################################################
## GDCL
################
## pre
GDCL_pre_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="GDCL" & period=="pre"),start=vb_sv)
GDCL_pre_vb_boot=nlsBoot(GDCL_pre_vb_fit,niter=400)
(GDCL_pre_vb_parms=cbind(ests=coef(GDCL_pre_vb_fit),confint(GDCL_pre_vb_boot)))

## Use model to predict length at age for plotting
GDCL_pre_predboot=Boot(GDCL_pre_vb_fit,f=GDCL_pre_vbpred)
GDCL_pre_preds=data.frame("GDCL",
                          "pre",
                          GDCL_pre_ages,
                          GDCL_pre_vbpred(GDCL_pre_vb_fit),
                          confint(GDCL_pre_predboot))
names(GDCL_pre_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="GDCL" & period=="pre"))

################
## post
GDCL_post_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="GDCL" & period=="post"),start=vb_sv)
GDCL_post_vb_boot=nlsBoot(GDCL_post_vb_fit,niter=400)
(GDCL_post_vb_parms=cbind(ests=coef(GDCL_post_vb_fit),confint(GDCL_post_vb_boot)))

## Use model to predict length at age for plotting
GDCL_post_predboot=Boot(GDCL_post_vb_fit,f=GDCL_post_vbpred)
GDCL_post_preds=data.frame("GDCL",
                           "post",
                           GDCL_post_ages,
                           GDCL_post_vbpred(GDCL_post_vb_fit),
                           confint(GDCL_post_predboot))
names(GDCL_post_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="GDCL" & period=="post"))

################
## plot
GDCL_preds=bind_rows(GDCL_pre_preds,GDCL_post_preds)%>%
  mutate(period=factor(period,levels=c("pre","post")))
GDCL_points=subset(blg,impd=="GDCL")%>%
  mutate(period=factor(period,levels=c("pre","post")))

GDCL_plot=ggplot()+
  geom_line(GDCL_preds,mapping=aes(x=age,y=tl,color=period))+
  geom_ribbon(GDCL_preds,mapping=aes(x=age,ymin=lci95,ymax=uci95,fill=period),alpha=0.5)+
  geom_jitter(GDCL_points,mapping=aes(x=age,y=cmgrp,color=period,shape=period),size=3,width=0.2,alpha=0.35)+
  scale_color_manual(values=c("#ff474c","#475f94"),
                     labels=c("Pre","Post"))+
  scale_fill_manual(values=c("#ff474c","#475f94"),
                    labels=c("Pre","Post"))+
  scale_shape_manual(values=c(17,19),
                     labels=c("Pre","Post"))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  scale_y_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  annotate("text",label="Gardner",x=0,y=220,hjust=0,vjust=1,size=6)+
  coord_cartesian(xlim=c(-0.25,8.25),
                  ylim=c(-1,225),
                  expand=F)+
  pubtheme

##################################################################
## HOCB
################
## pre
HOCB_pre_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="HOCB" & period=="pre"),start=vb_sv)
HOCB_pre_vb_boot=nlsBoot(HOCB_pre_vb_fit,niter=400)
(HOCB_pre_vb_parms=cbind(ests=coef(HOCB_pre_vb_fit),confint(HOCB_pre_vb_boot)))

## Use model to predict length at age for plotting
HOCB_pre_predboot=Boot(HOCB_pre_vb_fit,f=HOCB_pre_vbpred)
HOCB_pre_preds=data.frame("HOCB",
                          "pre",
                          HOCB_pre_ages,
                          HOCB_pre_vbpred(HOCB_pre_vb_fit),
                          confint(HOCB_pre_predboot))
names(HOCB_pre_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="HOCB" & period=="pre"))

################
## post
HOCB_post_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="HOCB" & period=="post"),start=vb_sv)
HOCB_post_vb_boot=nlsBoot(HOCB_post_vb_fit,niter=400)
(HOCB_post_vb_parms=cbind(ests=coef(HOCB_post_vb_fit),confint(HOCB_post_vb_boot)))

## Use model to predict length at age for plotting
HOCB_post_predboot=Boot(HOCB_post_vb_fit,f=HOCB_post_vbpred)
HOCB_post_preds=data.frame("HOCB",
                           "post",
                           HOCB_post_ages,
                           HOCB_post_vbpred(HOCB_post_vb_fit),
                           confint(HOCB_post_predboot))
names(HOCB_post_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="HOCB" & period=="post"))

################
## plot
HOCB_preds=bind_rows(HOCB_pre_preds,HOCB_post_preds)%>%
  mutate(period=factor(period,levels=c("pre","post")))
HOCB_points=subset(blg,impd=="HOCB")%>%
  mutate(period=factor(period,levels=c("pre","post")))

HOCB_plot=ggplot()+
  geom_line(HOCB_preds,mapping=aes(x=age,y=tl,color=period))+
  geom_ribbon(HOCB_preds,mapping=aes(x=age,ymin=lci95,ymax=uci95,fill=period),alpha=0.5)+
  geom_jitter(HOCB_points,mapping=aes(x=age,y=cmgrp,color=period,shape=period),size=3,width=0.2,alpha=0.35)+
  scale_color_manual(values=c("#ff474c","#475f94"),
                     labels=c("Pre","Post"))+
  scale_fill_manual(values=c("#ff474c","#475f94"),
                    labels=c("Pre","Post"))+
  scale_shape_manual(values=c(17,19),
                     labels=c("Pre","Post"))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  scale_y_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  annotate("text",label="Holton",x=0,y=220,hjust=0,vjust=1,size=6)+
  coord_cartesian(xlim=c(-0.25,8.25),
                  ylim=c(-1,225),
                  expand=F)+
  pubtheme

##################################################################
## POCL
################
## pre
POCL_pre_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="POCL" & period=="pre"),start=vb_sv)
POCL_pre_vb_boot=nlsBoot(POCL_pre_vb_fit,niter=400)
(POCL_pre_vb_parms=cbind(ests=coef(POCL_pre_vb_fit),confint(POCL_pre_vb_boot)))

## Use model to predict length at age for plotting
POCL_pre_predboot=Boot(POCL_pre_vb_fit,f=POCL_pre_vbpred)
POCL_pre_preds=data.frame("POCL",
                          "pre",
                          POCL_pre_ages,
                          POCL_pre_vbpred(POCL_pre_vb_fit),
                          confint(POCL_pre_predboot))
names(POCL_pre_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="POCL" & period=="pre"))

################
## post
POCL_post_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="POCL" & period=="post"),start=vb_sv)
POCL_post_vb_boot=nlsBoot(POCL_post_vb_fit,niter=400)
(POCL_post_vb_parms=cbind(ests=coef(POCL_post_vb_fit),confint(POCL_post_vb_boot)))

## Use model to predict length at age for plotting
POCL_post_predboot=Boot(POCL_post_vb_fit,f=POCL_post_vbpred)
POCL_post_preds=data.frame("POCL",
                           "post",
                           POCL_post_ages,
                           POCL_post_vbpred(POCL_post_vb_fit),
                           confint(POCL_post_predboot))
names(POCL_post_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="POCL" & period=="post"))

################
## plot
POCL_preds=bind_rows(POCL_pre_preds,POCL_post_preds)%>%
  mutate(period=factor(period,levels=c("pre","post")))
POCL_points=subset(blg,impd=="POCL")%>%
  mutate(period=factor(period,levels=c("pre","post")))

POCL_plot=ggplot()+
  geom_line(POCL_preds,mapping=aes(x=age,y=tl,color=period))+
  geom_ribbon(POCL_preds,mapping=aes(x=age,ymin=lci95,ymax=uci95,fill=period),alpha=0.5)+
  geom_jitter(POCL_points,mapping=aes(x=age,y=cmgrp,color=period,shape=period),size=3,width=0.2,alpha=0.35)+
  scale_color_manual(values=c("#ff474c","#475f94"),
                     labels=c("Pre","Post"))+
  scale_fill_manual(values=c("#ff474c","#475f94"),
                    labels=c("Pre","Post"))+
  scale_shape_manual(values=c(17,19),
                     labels=c("Pre","Post"))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  scale_y_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  annotate("text",label="Paola",x=0,y=220,hjust=0,vjust=1,size=6)+
  coord_cartesian(xlim=c(-0.25,8.25),
                  ylim=c(-1,225),
                  expand=F)+
  pubtheme

##################################################################
##################################################################
## JWSL
################
## pre
#JWSL_pre_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="JWSL" & period=="pre"),start=vb_sv)
#JWSL_pre_vb_boot=nlsBoot(JWSL_pre_vb_fit,niter=400)
#(JWSL_pre_vb_parms=cbind(ests=coef(JWSL_pre_vb_fit),confint(JWSL_pre_vb_boot)))

## Use model to predict length at age for plotting
#JWSL_pre_predboot=Boot(JWSL_pre_vb_fit,f=JWSL_pre_vbpred)
#JWSL_pre_preds=data.frame("JWSL",
#                          "pre",
#                          JWSL_pre_ages,
#                          JWSL_pre_vbpred(JWSL_pre_vb_fit),
#                          confint(JWSL_pre_predboot))
#names(JWSL_pre_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
#nrow(subset(blg,impd=="JWSL" & period=="pre"))

################
## post
JWSL_post_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="JWSL" & period=="post"),start=vb_sv)
JWSL_post_vb_boot=nlsBoot(JWSL_post_vb_fit,niter=400)
(JWSL_post_vb_parms=cbind(ests=coef(JWSL_post_vb_fit),confint(JWSL_post_vb_boot)))

## Use model to predict length at age for plotting
JWSL_post_predboot=Boot(JWSL_post_vb_fit,f=JWSL_post_vbpred)
JWSL_post_preds=data.frame("JWSL",
                           "post",
                           JWSL_post_ages,
                           JWSL_post_vbpred(JWSL_post_vb_fit),
                           confint(JWSL_post_predboot))
names(JWSL_post_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="JWSL" & period=="post"))

################
## plot
JWSL_preds=bind_rows(JWSL_post_preds)
JWSL_points=subset(blg,impd=="JWSL")%>%
  mutate(period=factor(period,levels=c("pre","post")))

JWSL_plot=ggplot()+
  geom_line(JWSL_preds,mapping=aes(x=age,y=tl),color="#475f94")+
  geom_ribbon(JWSL_preds,mapping=aes(x=age,ymin=lci95,ymax=uci95),fill="#475f94",alpha=0.5)+
  geom_jitter(JWSL_points,mapping=aes(x=age,y=cmgrp,color=period,shape=period),size=3,width=0.2,alpha=0.35)+
  scale_color_manual(values=c("#ff474c","#475f94"),
                     labels=c("Pre","Post"))+
  scale_fill_manual(values=c("#475f94"),
                    labels=c("Post"))+
  scale_shape_manual(values=c(17,19),
                     labels=c("Pre","Post"))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  scale_y_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  annotate("text",label="Jewell",x=0,y=220,hjust=0,vjust=1,size=6)+
  coord_cartesian(xlim=c(-0.25,8.25),
                  ylim=c(-1,225),
                  expand=F)+
  pubtheme

##################################################################
## LXLX
################
## pre
#LXLX_pre_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="LXLX" & period=="pre"),start=vb_sv)
#LXLX_pre_vb_boot=nlsBoot(LXLX_pre_vb_fit,niter=400)
#(LXLX_pre_vb_parms=cbind(ests=coef(LXLX_pre_vb_fit),confint(LXLX_pre_vb_boot)))

## Use model to predict length at age for plotting
#LXLX_pre_predboot=Boot(LXLX_pre_vb_fit,f=LXLX_pre_vbpred)
#LXLX_pre_preds=data.frame("LXLX",
#                          "pre",
#                          LXLX_pre_ages,
#                          LXLX_pre_vbpred(LXLX_pre_vb_fit),
#                          confint(LXLX_pre_predboot))
#names(LXLX_pre_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
#nrow(subset(blg,impd=="LXLX" & period=="pre"))

################
## post
LXLX_post_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="LXLX" & period=="post"),start=vb_sv)
LXLX_post_vb_boot=nlsBoot(LXLX_post_vb_fit,niter=400)
(LXLX_post_vb_parms=cbind(ests=coef(LXLX_post_vb_fit),confint(LXLX_post_vb_boot)))

## Use model to predict length at age for plotting
LXLX_post_predboot=Boot(LXLX_post_vb_fit,f=LXLX_post_vbpred)
LXLX_post_preds=data.frame("LXLX",
                           "post",
                           LXLX_post_ages,
                           LXLX_post_vbpred(LXLX_post_vb_fit),
                           confint(LXLX_post_predboot))
names(LXLX_post_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="LXLX" & period=="post"))

################
## plot
LXLX_preds=bind_rows(LXLX_post_preds)
LXLX_points=subset(blg,impd=="LXLX")%>%
  mutate(period=factor(period,levels=c("pre","post")))

LXLX_plot=ggplot()+
  geom_line(LXLX_preds,mapping=aes(x=age,y=tl),color="#475f94")+
  geom_ribbon(LXLX_preds,mapping=aes(x=age,ymin=lci95,ymax=uci95),fill="#475f94",alpha=0.5)+
  geom_jitter(LXLX_points,mapping=aes(x=age,y=cmgrp,color=period,shape=period),size=3,width=0.2,alpha=0.35)+
  scale_color_manual(values=c("#ff474c","#475f94"),
                     labels=c("Pre","Post"))+
  scale_fill_manual(values=c("#475f94"),
                    labels=c("Post"))+
  scale_shape_manual(values=c(17,19),
                     labels=c("Pre","Post"))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  scale_y_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  annotate("text",label="Lenexa",x=0,y=220,hjust=0,vjust=1,size=6)+
  coord_cartesian(xlim=c(-0.25,8.25),
                  ylim=c(-1,225),
                  expand=F)+
  pubtheme

##################################################################
## MISL
################
## pre
MISL_pre_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="MISL" & period=="pre"),start=vb_sv)
MISL_pre_vb_boot=nlsBoot(MISL_pre_vb_fit,niter=400)
(MISL_pre_vb_parms=cbind(ests=coef(MISL_pre_vb_fit),confint(MISL_pre_vb_boot)))

## Use model to predict length at age for plotting
MISL_pre_predboot=Boot(MISL_pre_vb_fit,f=MISL_pre_vbpred)
MISL_pre_preds=data.frame("MISL",
                          "pre",
                          MISL_pre_ages,
                          MISL_pre_vbpred(MISL_pre_vb_fit),
                          confint(MISL_pre_predboot))
names(MISL_pre_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="MISL" & period=="pre"))

################
## post
MISL_post_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="MISL" & period=="post"),start=vb_sv)
MISL_post_vb_boot=nlsBoot(MISL_post_vb_fit,niter=400)
(MISL_post_vb_parms=cbind(ests=coef(MISL_post_vb_fit),confint(MISL_post_vb_boot)))

## Use model to predict length at age for plotting
MISL_post_predboot=Boot(MISL_post_vb_fit,f=MISL_post_vbpred)
MISL_post_preds=data.frame("MISL",
                           "post",
                           MISL_post_ages,
                           MISL_post_vbpred(MISL_post_vb_fit),
                           confint(MISL_post_predboot))
names(MISL_post_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
nrow(subset(blg,impd=="MISL" & period=="post"))

################
## plot
MISL_preds=bind_rows(MISL_pre_preds,MISL_post_preds)%>%
  mutate(period=factor(period,levels=c("pre","post")))
MISL_points=subset(blg,impd=="MISL")%>%
  mutate(period=factor(period,levels=c("pre","post")))

MISL_plot=ggplot()+
  geom_line(MISL_preds,mapping=aes(x=age,y=tl,color=period))+
  geom_ribbon(MISL_preds,mapping=aes(x=age,ymin=lci95,ymax=uci95,fill=period),alpha=0.5)+
  geom_jitter(MISL_points,mapping=aes(x=age,y=cmgrp,color=period,shape=period),size=3,width=0.2,alpha=0.35)+
  scale_color_manual(values=c("#ff474c","#475f94"),
                     labels=c("Pre","Post"))+
  scale_fill_manual(values=c("#ff474c","#475f94"),
                    labels=c("Pre","Post"))+
  scale_shape_manual(values=c(17,19),
                     labels=c("Pre","Post"))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  scale_y_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  annotate("text",label="MISL",x=0,y=220,hjust=0,vjust=1,size=6)+
  coord_cartesian(xlim=c(-0.25,8.25),
                  ylim=c(-1,225),
                  expand=F)+
  pubtheme

##################################################################
## PTSB
################
## pre
#PTSB_pre_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="PTSB" & period=="pre"),start=vb_sv)
#PTSB_pre_vb_boot=nlsBoot(PTSB_pre_vb_fit,niter=400)
#(PTSB_pre_vb_parms=cbind(ests=coef(PTSB_pre_vb_fit),confint(PTSB_pre_vb_boot)))

## Use model to predict length at age for plotting
#PTSB_pre_predboot=Boot(PTSB_pre_vb_fit,f=PTSB_pre_vbpred)
#PTSB_pre_preds=data.frame("PTSB",
#                          "pre",
#                          PTSB_pre_ages,
#                          PTSB_pre_vbpred(PTSB_pre_vb_fit),
#                          confint(PTSB_pre_predboot))
#names(PTSB_pre_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
#nrow(subset(blg,impd=="PTSB" & period=="pre"))

################
## post
#PTSB_post_vb_fit=nls(cmgrp~vb(age,Linf,k,t0),data=subset(blg,impd=="PTSB" & period=="post"),start=vb_sv)
#PTSB_post_vb_boot=nlsBoot(PTSB_post_vb_fit,niter=400)
#(PTSB_post_vb_parms=cbind(ests=coef(PTSB_post_vb_fit),confint(PTSB_post_vb_boot)))

## Use model to predict length at age for plotting
#PTSB_post_predboot=Boot(PTSB_post_vb_fit,f=PTSB_post_vbpred)
#PTSB_post_preds=data.frame("PTSB",
#                           "post",
#                           PTSB_post_ages,
#                           PTSB_post_vbpred(PTSB_post_vb_fit),
#                           confint(PTSB_post_predboot))
#names(PTSB_post_preds)=c("impd","period","age","tl","lci95","uci95")

## Number of fish
#nrow(subset(blg,impd=="PTSB" & period=="post"))

################
## plot
#PTSB_preds=bind_rows(PTSB_pre_preds,PTSB_post_preds)%>%
#  mutate(period=factor(period,levels=c("pre","post")))
PTSB_points=subset(blg,impd=="PTSB")%>%
  mutate(period=factor(period,levels=c("pre","post")))

PTSB_plot=ggplot()+
#  geom_line(PTSB_preds,mapping=aes(x=age,y=tl,color=period))+
#  geom_ribbon(PTSB_preds,mapping=aes(x=age,ymin=lci95,ymax=uci95,fill=period),alpha=0.5)+
  geom_jitter(PTSB_points,mapping=aes(x=age,y=cmgrp,color=period,shape=period),size=3,width=0.2,alpha=0.35)+
  scale_color_manual(values=c("#ff474c","#475f94"),
                     labels=c("Pre","Post"))+
#  scale_fill_manual(values=c("#ff474c","#475f94"),
#                    labels=c("Pre","Post"))+
  scale_shape_manual(values=c(17,19),
                     labels=c("Pre","Post"))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  scale_y_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  annotate("text",label="Pott #2",x=0,y=220,hjust=0,vjust=1,size=6)+
  coord_cartesian(xlim=c(-0.25,8.25),
                  ylim=c(-1,225),
                  expand=F)+
  pubtheme

##################################################################
## Combine plots
out=(ZPCB_plot|GDCL_plot|HOCB_plot|POCL_plot)/(JWSL_plot|LXLX_plot|MISL_plot|PTSB_plot)
ggsave(plot=out,"growth.png",width=16,height=8,units="in",bg="white")