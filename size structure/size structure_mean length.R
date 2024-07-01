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

if("Rmisc" %in% rownames(installed.packages()) == FALSE) {install.packages("Rmisc")}
library(Rmisc)

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

## Manipulate data a bit for comparing mean length
fish1=fish%>%
  filter(gear=="EF",
         spp=="Bluegill",
         tl>=50,
         year==2017|year==2018|year==2022)%>%
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
################################################################################
## ZPCB
ZPCB=filter(fish1,impd=="ZPCB")

## Find mean length and 95% CI
ZPCB%>%
  group_by(period)%>%
  summarize(n=n(),
            mn=mean(tl),
            lci=CI(tl,ci=0.95)['lower'],
            uci=CI(tl,ci=0.95)['upper'])%>%
  ungroup()

## Run t-test
t.test(tl~period,ZPCB)

###############################
## ECDF plot
ZPCB_ss=ggplot()+
  stat_ecdf(ZPCB,mapping=aes(x=tl,color=period),linewidth=1)+
  scale_color_manual(labels=c("Pre: N = 350",
                              "Post: N = 98"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="Cumulative frequency")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("rect",xmin=92,xmax=95,ymin=-Inf,ymax=Inf,fill="#ff474c",alpha=0.5)+
  annotate("rect",xmin=129,xmax=142,ymin=-Inf,ymax=Inf,fill="#475f94",alpha=0.5)+
  annotate("text",label="Campbell",x=5,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(t)","==",13.044),parse=T,x=5,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=5,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme
################################################################################
## GDCL
GDCL=filter(fish1,impd=="GDCL")

## Find mean length and 95% CI
GDCL%>%
  group_by(period)%>%
  summarize(n=n(),
            mn=mean(tl),
            lci=CI(tl,ci=0.95)['lower'],
            uci=CI(tl,ci=0.95)['upper'])%>%
ungroup()

## Run t-test
t.test(tl~period,GDCL)

###############################
## ECDF plot
GDCL_ss=ggplot()+
  stat_ecdf(GDCL,mapping=aes(x=tl,color=period),linewidth=1)+
  scale_color_manual(labels=c("Pre: N = 329",
                              "Post: N = 78"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("rect",xmin=134,xmax=140,ymin=-Inf,ymax=Inf,fill="#ff474c",alpha=0.5)+
  annotate("rect",xmin=120,xmax=139,ymin=-Inf,ymax=Inf,fill="#475f94",alpha=0.5)+
  annotate("text",label="Gardner",x=5,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(t)","==",-1.496),parse=T,x=5,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",0.138),parse=T,x=5,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## HOCB
HOCB=filter(fish1,impd=="HOCB")

## Find mean length and 95% CI
HOCB%>%
  group_by(period)%>%
  summarize(n=n(),
            mn=mean(tl),
            lci=CI(tl,ci=0.95)['lower'],
            uci=CI(tl,ci=0.95)['upper'])%>%
ungroup()

## Run t-test
t.test(tl~period,HOCB)

###############################
## ECDF plot
HOCB_ss=ggplot()+
  stat_ecdf(HOCB,mapping=aes(x=tl,color=period),linewidth=1)+
  scale_color_manual(labels=c("Pre: N = 103",
                              "Post: N = 108"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("rect",xmin=70,xmax=80,ymin=-Inf,ymax=Inf,fill="#ff474c",alpha=0.5)+
  annotate("rect",xmin=110,xmax=120,ymin=-Inf,ymax=Inf,fill="#475f94",alpha=0.5)+
  annotate("text",label="Holton",x=5,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(t)","==",-11.496),parse=T,x=5,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=5,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## POCL
POCL=filter(fish1,impd=="POCL")

## Find mean length and 95% CI
POCL%>%
  group_by(period)%>%
  summarize(n=n(),
            mn=mean(tl),
            lci=CI(tl,ci=0.95)['lower'],
            uci=CI(tl,ci=0.95)['upper'])%>%
  ungroup()

## Run t-test
t.test(tl~period,POCL)

###############################
## ECDF plot
POCL_ss=ggplot()+
  stat_ecdf(POCL,mapping=aes(x=tl,color=period),linewidth=1)+
  scale_color_manual(labels=c("Pre: N = 199",
                              "Post: N = 173"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("rect",xmin=101,xmax=109,ymin=-Inf,ymax=Inf,fill="#ff474c",alpha=0.5)+
  annotate("rect",xmin=109,xmax=117,ymin=-Inf,ymax=Inf,fill="#475f94",alpha=0.5)+
  annotate("text",label="Paola",x=5,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(t)","==",-2.769),parse=T,x=5,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",0.005),parse=T,x=5,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## JWSL
JWSL=filter(fish1,impd=="JWSL")

## Find mean length and 95% CI
JWSL%>%
  group_by(period)%>%
  summarize(n=n(),
            mn=mean(tl),
            lci=CI(tl,ci=0.95)['lower'],
            uci=CI(tl,ci=0.95)['upper'])%>%
  ungroup()

## Run t-test
t.test(tl~period,JWSL)

###############################
## ECDF plot
JWSL_ss=ggplot()+
  stat_ecdf(JWSL,mapping=aes(x=tl,color=period),linewidth=1)+
  scale_color_manual(labels=c("Pre: N = 196",
                              "Post: N = 160"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="Cumulative frequency")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("rect",xmin=93,xmax=100,ymin=-Inf,ymax=Inf,fill="#ff474c",alpha=0.5)+
  annotate("rect",xmin=94,xmax=103,ymin=-Inf,ymax=Inf,fill="#475f94",alpha=0.5)+
  annotate("text",label="Jewell",x=5,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(t)","==",-0.663),parse=T,x=5,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",0.508),parse=T,x=5,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## LXLX
LXLX=filter(fish1,impd=="LXLX")

## Find mean length and 95% CI
LXLX%>%
  group_by(period)%>%
  summarize(n=n(),
            mn=mean(tl),
            lci=CI(tl,ci=0.95)['lower'],
            uci=CI(tl,ci=0.95)['upper'])%>%
  ungroup()

## Run t-test
t.test(tl~period,LXLX)

###############################
## ECDF plot
LXLX_ss=ggplot()+
  stat_ecdf(LXLX,mapping=aes(x=tl,color=period),linewidth=1)+
  scale_color_manual(labels=c("Pre: N = 433",
                              "Post: N = 35"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("rect",xmin=85,xmax=90,ymin=-Inf,ymax=Inf,fill="#ff474c",alpha=0.5)+
  annotate("rect",xmin=78,xmax=94,ymin=-Inf,ymax=Inf,fill="#475f94",alpha=0.5)+
  annotate("text",label="Lenexa",x=5,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(t)","==",0.338),parse=T,x=5,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",0.736),parse=T,x=5,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## MISL
MISL=filter(fish1,impd=="MISL")

## Find mean length and 95% CI
MISL%>%
  group_by(period)%>%
  summarize(n=n(),
            mn=mean(tl),
            lci=CI(tl,ci=0.95)['lower'],
            uci=CI(tl,ci=0.95)['upper'])%>%
  ungroup()

## Run t-test
t.test(tl~period,MISL)

###############################
## ECDF plot
MISL_ss=ggplot()+
  stat_ecdf(MISL,mapping=aes(x=tl,color=period),linewidth=1)+
  scale_color_manual(labels=c("Pre: N = 747",
                              "Post: N = 269"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("rect",xmin=119,xmax=123,ymin=-Inf,ymax=Inf,fill="#ff474c",alpha=0.5)+
  annotate("rect",xmin=126,xmax=134,ymin=-Inf,ymax=Inf,fill="#475f94",alpha=0.5)+
  annotate("text",label="Miami",x=5,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(t)","==",-3.889),parse=T,x=5,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","<",0.001),parse=T,x=5,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

################################################################################
## PTSB
PTSB=filter(fish1,impd=="PTSB")

## Find mean length and 95% CI
PTSB%>%
  group_by(period)%>%
  summarize(n=n(),
            mn=mean(tl),
            lci=CI(tl,ci=0.95)['lower'],
            uci=CI(tl,ci=0.95)['upper'])%>%
  ungroup()

## Run t-test
t.test(tl~period,PTSB)

###############################
## ECDF plot
PTSB_ss=ggplot()+
  stat_ecdf(PTSB,mapping=aes(x=tl,color=period),linewidth=1)+
  scale_color_manual(labels=c("Pre: N = 230",
                              "Post: N = 87"),
                     values=c("#ff474c","#475f94"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     name="")+
  scale_x_continuous(breaks=seq(0,220,20),
                     name="Total length (mm)")+
  coord_cartesian(xlim=c(0,225),
                  ylim=c(-0.01,1.01),
                  expand=F)+
  annotate("rect",xmin=70,xmax=76,ymin=-Inf,ymax=Inf,fill="#ff474c",alpha=0.5)+
  annotate("rect",xmin=75,xmax=85,ymin=-Inf,ymax=Inf,fill="#475f94",alpha=0.5)+
  annotate("text",label="Pott #2",x=5,y=1,hjust=0,vjust=1,size=6)+
  annotate("text",label=paste("italic(t)","==",-2.211),parse=T,x=5,y=0.92,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",0.029),parse=T,x=5,y=0.88,hjust=0,vjust=1,size=4)+
  pubtheme

#################################################################
#################################################################
## Combine plots and export
plots=(ZPCB_ss|GDCL_ss|HOCB_ss|POCL_ss)/(JWSL_ss|LXLX_ss|MISL_ss|PTSB_ss)
ggsave(plot=plots,"size structure_mean length.png",width=16,height=8,units="in",bg="white")