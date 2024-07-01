#--------------------------------------------------------------
#Ben Neely
#05/13/2024
#Size structure of Redear Sunfish populations
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

if("car" %in% rownames(installed.packages()) == FALSE) {install.packages("car")}
library(car)

if("patchwork" %in% rownames(installed.packages()) == FALSE) {install.packages("patchwork")}
library(patchwork)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=18,color="black",face="bold"),
        axis.text.y=element_text(size=14,color="black"),
        axis.text.x=element_text(size=14,color="black"),
        legend.position=c(1,0),
        legend.justification=c("right","bottom"),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.background=element_rect(fill="transparent"),
        legend.margin=margin(c(0,1,1,1)))
options(scipen=999)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/BLG regulation eval/size structure/")

## Read in data with import
fish=import("regfish.csv")

## Manipulate data a bit for comparing size structure
fish1=fish%>%
  filter(gear=="EF",
         spp=="Redear Sunfish")%>%
  expandCounts(~count)%>%
  mutate(period=case_when(year<2019 ~ "pre",
                          year>=2019 ~ "post",
                          TRUE ~ "other"),
         period=factor(period,levels=c("pre","post")),
         trt=case_when(impd=="LXLX"|impd=="PTSB"|impd=="JWSL"|impd=="MISL" ~ "regulation",
                       impd=="GDCL"|impd=="HOCB"|impd=="ZPCB"|impd=="POCL" ~ "control",
                       TRUE ~ "other"))%>%
  select(impd,year,trt,spp,tl,period)

################################################################################
################################################################################
## KS tests comparing pre and post and plot
ksdat=fish1%>%
  filter(year==2017|year==2018|year==2022)

## Function to find location of D
dfind=function(x, y) {
  n=length(x); m=length(y)
  w=c(x,y)
  o=order(w)
  z=cumsum(ifelse(o<=n,m,-n))
  i=which.max(abs(z))
  w[o[i]]
}

################################################################################
################################################################################
## ZPCB
ZPCB=ksdat%>%
  filter(impd=="ZPCB")
nrow(subset(ZPCB,period=="pre"))
nrow(subset(ZPCB,period=="post"))

## KS test
ksTest(tl~period,data=ZPCB)
#D = 0.462, P < 0.001

## Location of D on X
ZPCB_dloc=dfind(subset(ZPCB,period=="pre")$tl,subset(ZPCB,period=="post")$tl)

## pre and post cumulative frequency at D
ZPCB_kspre_D=ecdf(subset(ZPCB,period=="pre")$tl) (v=ZPCB_dloc)
ZPCB_kspost_D=ecdf(subset(ZPCB,period=="post")$tl) (v=ZPCB_dloc)

###############################
## ECDF plot
ZPCB_ss=ggplot(ZPCB,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 74",
                              "Post: N = 26"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="Cumulative frequency")+
  scale_x_continuous(breaks=seq(0,260,20),
                     name="")+
  coord_cartesian(xlim=c(0,265),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=ZPCB_dloc,xend=ZPCB_dloc,y=ZPCB_kspre_D,yend=ZPCB_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Campbell",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.462),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme
################################################################################
## GDCL

################################################################################
## HOCB

################################################################################
## POCL

################################################################################
## JWSL
JWSL=ksdat%>%
  filter(impd=="JWSL")
nrow(subset(JWSL,period=="pre"))
nrow(subset(JWSL,period=="post"))

## KS test
ksTest(tl~period,data=JWSL)
#D = 0.396, P < 0.001

## Location of D on X
JWSL_dloc=dfind(subset(JWSL,period=="pre")$tl,subset(JWSL,period=="post")$tl)

## pre and post cumulative frequency at D
JWSL_kspre_D=ecdf(subset(JWSL,period=="pre")$tl) (v=JWSL_dloc)
JWSL_kspost_D=ecdf(subset(JWSL,period=="post")$tl) (v=JWSL_dloc)

###############################
## ECDF plot
JWSL_ss=ggplot(JWSL,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 139",
                              "Post: N = 106"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="Cumulative frequency")+
  scale_x_continuous(breaks=seq(0,260,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,265),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=JWSL_dloc,xend=JWSL_dloc,y=JWSL_kspre_D,yend=JWSL_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Jewell",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.396),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## LXLX

################################################################################
## MISL

################################################################################
## PTSB

#################################################################
#################################################################
## Combine plots and export
plots=ZPCB_ss/JWSL_ss
ggsave(plot=plots,"res size structure.png",width=8,height=8,units="in",bg="white")