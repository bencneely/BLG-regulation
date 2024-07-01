#--------------------------------------------------------------
#Ben Neely
#05/13/2024
#Size structure of Bluegill populations
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
        axis.text.x=element_text(size=11,color="black"),
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
         spp=="Bluegill")%>%
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
#D = 0.684, P < 0.001

## Location of D on X
ZPCB_dloc=dfind(subset(ZPCB,period=="pre")$tl,subset(ZPCB,period=="post")$tl)

## pre and post cumulative frequency at D
ZPCB_kspre_D=ecdf(subset(ZPCB,period=="pre")$tl) (v=ZPCB_dloc)
ZPCB_kspost_D=ecdf(subset(ZPCB,period=="post")$tl) (v=ZPCB_dloc)

###############################
## ECDF plot
ZPCB_ss=ggplot(ZPCB,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 481",
                              "Post: N = 116"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="Cumulative frequency")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=ZPCB_dloc,xend=ZPCB_dloc,y=ZPCB_kspre_D,yend=ZPCB_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Campbell",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.684),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme
################################################################################
## GDCL
GDCL=ksdat%>%
  filter(impd=="GDCL")
nrow(subset(GDCL,period=="pre"))
nrow(subset(GDCL,period=="post"))

## KS test
ksTest(tl~period,data=GDCL)
#D = 0.275, P < 0.001

## Location of D on X
GDCL_dloc=dfind(subset(GDCL,period=="pre")$tl,subset(GDCL,period=="post")$tl)

## pre and post cumulative frequency at D
GDCL_kspre_D=ecdf(subset(GDCL,period=="pre")$tl) (v=GDCL_dloc)
GDCL_kspost_D=ecdf(subset(GDCL,period=="post")$tl) (v=GDCL_dloc)

###############################
## ECDF plot
GDCL_ss=ggplot(GDCL,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 329",
                              "Post: N = 82"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=GDCL_dloc,xend=GDCL_dloc,y=GDCL_kspre_D,yend=GDCL_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Gardner",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.275),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## HOCB
HOCB=ksdat%>%
  filter(impd=="HOCB")
nrow(subset(HOCB,period=="pre"))
nrow(subset(HOCB,period=="post"))

## KS test
ksTest(tl~period,data=HOCB)
#D = 0.614, P < 0.001

## Location of D on X
HOCB_dloc=dfind(subset(HOCB,period=="pre")$tl,subset(HOCB,period=="post")$tl)

## pre and post cumulative frequency at D
HOCB_kspre_D=ecdf(subset(HOCB,period=="pre")$tl) (v=HOCB_dloc)
HOCB_kspost_D=ecdf(subset(HOCB,period=="post")$tl) (v=HOCB_dloc)

###############################
## ECDF plot
HOCB_ss=ggplot(HOCB,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 104",
                              "Post: N = 110"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=HOCB_dloc,xend=HOCB_dloc,y=HOCB_kspre_D,yend=HOCB_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Holton",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.614),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## POCL
POCL=ksdat%>%
  filter(impd=="POCL")
nrow(subset(POCL,period=="pre"))
nrow(subset(POCL,period=="post"))

## KS test
ksTest(tl~period,data=POCL)
#D = 0.222, P < 0.001

## Location of D on X
POCL_dloc=dfind(subset(POCL,period=="pre")$tl,subset(POCL,period=="post")$tl)

## pre and post cumulative frequency at D
POCL_kspre_D=ecdf(subset(POCL,period=="pre")$tl) (v=POCL_dloc)
POCL_kspost_D=ecdf(subset(POCL,period=="post")$tl) (v=POCL_dloc)

###############################
## ECDF plot
POCL_ss=ggplot(POCL,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 200",
                              "Post: N = 183"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=POCL_dloc,xend=POCL_dloc,y=POCL_kspre_D,yend=POCL_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Paola",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.222),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## JWSL
JWSL=ksdat%>%
  filter(impd=="JWSL")
nrow(subset(JWSL,period=="pre"))
nrow(subset(JWSL,period=="post"))

## KS test
ksTest(tl~period,data=JWSL)
#D = 0.136, P = 0.060

## Location of D on X
JWSL_dloc=dfind(subset(JWSL,period=="pre")$tl,subset(JWSL,period=="post")$tl)

## pre and post cumulative frequency at D
JWSL_kspre_D=ecdf(subset(JWSL,period=="pre")$tl) (v=JWSL_dloc)
JWSL_kspost_D=ecdf(subset(JWSL,period=="post")$tl) (v=JWSL_dloc)

###############################
## ECDF plot
JWSL_ss=ggplot(JWSL,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 211",
                              "Post: N = 173"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="Cumulative frequency")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=JWSL_dloc,xend=JWSL_dloc,y=JWSL_kspre_D,yend=JWSL_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Jewell",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.136),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==","'0.060'"),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## LXLX
LXLX=ksdat%>%
  filter(impd=="LXLX")
nrow(subset(LXLX,period=="pre"))
nrow(subset(LXLX,period=="post"))

## KS test
ksTest(tl~period,data=LXLX)
#D = 0.131, P = 0.617

## Location of D on X
LXLX_dloc=dfind(subset(LXLX,period=="pre")$tl,subset(LXLX,period=="post")$tl)

## pre and post cumulative frequency at D
LXLX_kspre_D=ecdf(subset(LXLX,period=="pre")$tl) (v=LXLX_dloc)
LXLX_kspost_D=ecdf(subset(LXLX,period=="post")$tl) (v=LXLX_dloc)

###############################
## ECDF plot
LXLX_ss=ggplot(LXLX,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 433",
                              "Post: N = 36"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=LXLX_dloc,xend=LXLX_dloc,y=LXLX_kspre_D,yend=LXLX_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Lenexa",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.131),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",0.617),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## MISL
MISL=ksdat%>%
  filter(impd=="MISL")
nrow(subset(MISL,period=="pre"))
nrow(subset(MISL,period=="post"))

## KS test
ksTest(tl~period,data=MISL)
#D = 0.230, P < 0.001

## Location of D on X
MISL_dloc=dfind(subset(MISL,period=="pre")$tl,subset(MISL,period=="post")$tl)

## pre and post cumulative frequency at D
MISL_kspre_D=ecdf(subset(MISL,period=="pre")$tl) (v=MISL_dloc)
MISL_kspost_D=ecdf(subset(MISL,period=="post")$tl) (v=MISL_dloc)

###############################
## ECDF plot
MISL_ss=ggplot(MISL,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 747",
                              "Post: N = 272"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=MISL_dloc,xend=MISL_dloc,y=MISL_kspre_D,yend=MISL_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Miami",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==","'0.230'"),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## PTSB
PTSB=ksdat%>%
  filter(impd=="PTSB")
nrow(subset(PTSB,period=="pre"))
nrow(subset(PTSB,period=="post"))

## KS test
ksTest(tl~period,data=PTSB)
#D = 0.215, P = 0.005

## Location of D on X
PTSB_dloc=dfind(subset(PTSB,period=="pre")$tl,subset(PTSB,period=="post")$tl)

## pre and post cumulative frequency at D
PTSB_kspre_D=ecdf(subset(PTSB,period=="pre")$tl) (v=PTSB_dloc)
PTSB_kspost_D=ecdf(subset(PTSB,period=="post")$tl) (v=PTSB_dloc)

###############################
## ECDF plot
PTSB_ss=ggplot(PTSB,aes(x=tl,color=period))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(labels=c("Pre: N = 231",
                              "Post: N = 89"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("segment",x=PTSB_dloc,xend=PTSB_dloc,y=PTSB_kspre_D,yend=PTSB_kspost_D,alpha=1,linewidth=1.5)+
  annotate("text",label="Pott #2",x=1,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(D)","==",0.215),parse=T,x=1,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",0.005),parse=T,x=1,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

#################################################################
#################################################################
## Combine plots and export
plots=(ZPCB_ss|GDCL_ss|HOCB_ss|POCL_ss)/(JWSL_ss|LXLX_ss|MISL_ss|PTSB_ss)
ggsave(plot=plots,"size structure.png",width=16,height=8,units="in",bg="white")