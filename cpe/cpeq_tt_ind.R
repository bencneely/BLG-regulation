#--------------------------------------------------------------
#Ben Neely
#05/15/2024
#Relative abundance of Bluegill populations
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

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=26,color="black",face="bold"),
        axis.text.y=element_text(size=20,color="black"),
        axis.text.x=element_text(size=20,color="black"),
        legend.position=c(1,1),
        legend.justification=c("right","top"),
        legend.key.width=unit(2,"cm"),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        legend.background=element_rect(fill="transparent"))
options(scipen=999)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/BLG regulation eval/cpe/")

## Read in data with import
fish=import("regfish.csv")
samp=import("regsamp.csv")

## Manipulate data a bit for comparing relative abundance
## Fish data
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
                       TRUE ~ "other"),
         impd=factor(impd,levels=c("JWSL","LXLX","MISL","PTSB",
                                   "ZPCB","GDCL","HOCB","POCL")))%>%
  select(uid,impd,spp,tl,period)

## Sample data
samp1=samp%>%
  filter(gear=="EF")%>%
  mutate(period=case_when(year<2019 ~ "pre",
                          year>=2019 ~ "post",
                          TRUE ~ "other"),
         period=factor(period,levels=c("pre","post")),
         trt=case_when(impd=="LXLX"|impd=="PTSB"|impd=="JWSL"|impd=="MISL" ~ "regulation",
                       impd=="GDCL"|impd=="HOCB"|impd=="ZPCB"|impd=="POCL" ~ "control",
                       TRUE ~ "other"),
         impd=factor(impd,levels=c("JWSL","LXLX","MISL","PTSB",
                                   "ZPCB","GDCL","HOCB","POCL")))%>%
  select(uid,impd,gear,effort,year,trt,period)

################################################################################
################################################################################
## Group fish to calculate CPEq
fish2=fish1%>%
  filter(tl>=150)%>%
  group_by(uid,impd,spp,period)%>%
  summarize(tot=n())%>%
  ungroup()

## Join with sample data
tmp=samp1%>%
  left_join(fish2,by="uid")%>%
  complete(spp,tot,fill=list(spp="Bluegill",tot=0))%>%
  select(uid,impd=impd.x,trt,year,period=period.x,gear,effort,spp,tot)%>%
  drop_na()

## Calculate CPEq per sample, add 0.1, and log10 transform
## Only keep one sample per impd/period
cpedat=tmp%>%
  mutate(cpeq=(tot/effort),
         log_cpeq=log10(cpeq))%>%
  filter(year==2017|year==2018|year==2022)

################################################################################
## Two-sample t-tests
## One pre (2017 or 2018) and one post (2022) estimate for each impoundment

################
## Control
con=cpedat%>%
  filter(trt=="control")%>%
  group_by(impd,trt,period)%>%
  summarize(cpeq=mean(cpeq))%>%
  ungroup()%>%
  mutate(log_cpeq=log10(cpeq+1))

## Two-sample t-test (n=4)
con_tt=t.test(log_cpeq~period,paired=F,data=con)

################
## Regulation
reg=cpedat%>%
  filter(trt=="regulation")%>%
  group_by(impd,trt,period)%>%
  summarize(cpeq=mean(cpeq))%>%
  ungroup()%>%
  mutate(log_cpeq=log10(cpeq+1))

## Paired t-test (n=4)
reg_tt=t.test(log_cpeq~period,paired=F,data=reg)

################################################################################
## Create plots
plotdat=bind_rows(con,reg)%>%
  mutate(trt=factor(trt,levels=c("regulation","control")),
         period=factor(period,levels=c("post","pre")),
         lab=c("C","C","G","G","H","H","O","O",
               "J","J","L","L","M","M","P","P"))
## Create plot
ggplot(plotdat)+
  geom_point(aes(x=log_cpeq,y=trt,color=period),position=position_dodge(width=0.3),size=10)+
  geom_path(aes(x=log_cpeq,y=trt,color=period),position=position_dodge(width=0.3),linewidth=2)+
  geom_text(aes(x=log_cpeq,y=trt,label=lab,group=period),color="white",position=position_dodge(width=0.3),size=8)+
  scale_color_manual(values=c("#475f94","#ff474c"),
                    labels=c("Post","Pre"),
                    guide=guide_legend(reverse=T))+
  scale_x_continuous(breaks=seq(0,2.5,0.5),
                     labels=c(1,3,10,30,100,300),
                     name=expression("CPE"*italic(q)))+
  scale_y_discrete(labels=c("Regulation","Control"),
                   name="")+
  coord_cartesian(xlim=c(-0.07,2.6),
                  ylim=c(0.5,2.5),
                  expand=F)+
  pubtheme

## Export plot
ggsave(plot=last_plot(),"CPEq_tt_ind.png",height=6,width=10,units="in",bg="white")
