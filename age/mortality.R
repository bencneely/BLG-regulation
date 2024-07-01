#--------------------------------------------------------------
#Ben Neely
#05/13/2024
#Annual mortality of Bluegill populations
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
        axis.title=element_text(size=22,color="black",face="bold"),
        axis.text=element_text(size=18,color="black"),
        legend.position=c(1,1),
        legend.justification=c("right","top"),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        legend.background=element_rect(fill="transparent"),
        legend.margin=margin(c(0,1,1,1)))
options(scipen=999)

## Set working directory
setwd("C:/Users/Ben.Neely/OneDrive - State of Kansas, OITS/Desktop/BLG regulation eval/age/")

## Read in data with import
dat=import("allaged.csv")%>%
  mutate(period=factor(period,
                       levels=c("pre","post")))

blg=filter(dat,spp=="Bluegill")

##################################################################
##################################################################
##################################################################
## Annual mortality
## Summarize data to show number caught per impd/period/age class
blg_cc=blg%>%
  filter(age>0)%>%
  group_by(impd,trt,period,age)%>%
  summarize(ct=n()+1,
            lnct=log(ct))%>%
  ungroup()

#################################################################
#################################################################
## GDCL
gdcl_cc=blg_cc%>%
  filter(impd=="GDCL")

###################
## Set up pre data
gdcl_pre=blg_cc%>%
  filter(impd=="GDCL" & period=="pre")

## Fit weighted catch curve to pre data
gdcl_pre_cc=catchCurve(ct~age,gdcl_pre,weighted=T)

## Extract variables for plotting pre data
gdcl_pre_ccvars=bind_cols(age=gdcl_pre_cc$age.e,
                          ln_catch=gdcl_pre_cc$log.catch.e,
                          weight=gdcl_pre_cc$weights.e,
                          impd="gdcl",
                          period="pre")

###################
## Set up post data
gdcl_post=blg_cc%>%
  filter(impd=="GDCL" & period=="post")

## Fit weighted catch curve to post data
gdcl_post_cc=catchCurve(ct~age,gdcl_post,weighted=T)

## Extract variables for plotting post data
gdcl_post_ccvars=bind_cols(age=gdcl_post_cc$age.e,
                           ln_catch=gdcl_post_cc$log.catch.e,
                           weight=gdcl_post_cc$weights.e,
                           impd="gdcl",
                           period="post")

###################
## Combine pre and post data
gdcl_comb=bind_rows(gdcl_pre_ccvars,gdcl_post_ccvars)%>%
  mutate(catch=exp(ln_catch),
         period=factor(period,levels=c("pre","post")),
         trt="control")

###################
## Weighted catch curve using DVR
## Page 211 IFAR
gdcl_lm=lm(ln_catch~age*period,weights=weight,data=gdcl_comb)
Anova(gdcl_lm)
## Interaction term P=0.224 suggesting no difference in slope (Z)

## Extract mortality estimates
gdcl_coef=coef(gdcl_lm)
gdcl_pre_Z=-1*gdcl_coef[["age"]]
gdcl_post_Z=-1*(gdcl_coef[["age"]]+gdcl_coef[["age:periodpost"]])
gdcl_pre_A=1-exp(-gdcl_pre_Z)
gdcl_post_A=1-exp(-gdcl_post_Z)

## Summarize into table
gdcl_out=bind_cols(impd="gdcl",
                   trt="control",
                   Fval=1.681,
                   Pval=0.224,
                   pre_Z=gdcl_pre_Z,
                   pre_A=gdcl_pre_A,
                   post_Z=gdcl_post_Z,
                   post_A=gdcl_post_A)

#################################################################
#################################################################
## HOCB
hocb_cc=blg_cc%>%
  filter(impd=="HOCB")

###################
## Set up pre data
hocb_pre=blg_cc%>%
  filter(impd=="HOCB" & period=="pre")

## Fit weighted catch curve to pre data
hocb_pre_cc=catchCurve(ct~age,hocb_pre,weighted=T)

## Extract variables for plotting pre data
hocb_pre_ccvars=bind_cols(age=hocb_pre_cc$age.e,
                          ln_catch=hocb_pre_cc$log.catch.e,
                          weight=hocb_pre_cc$weights.e,
                          impd="hocb",
                          period="pre")

###################
## Set up post data
hocb_post=blg_cc%>%
  filter(impd=="HOCB" & period=="post")

## Fit weighted catch curve to post data
hocb_post_cc=catchCurve(ct~age,hocb_post,weighted=T)

## Extract variables for plotting post data
hocb_post_ccvars=bind_cols(age=hocb_post_cc$age.e,
                           ln_catch=hocb_post_cc$log.catch.e,
                           weight=hocb_post_cc$weights.e,
                           impd="hocb",
                           period="post")

###################
## Combine pre and post data
hocb_comb=bind_rows(hocb_pre_ccvars,hocb_post_ccvars)%>%
  mutate(catch=exp(ln_catch),
         period=factor(period,levels=c("pre","post")),
         trt="control")

###################
## Weighted catch curve using DVR
## Page 211 IFAR
hocb_lm=lm(ln_catch~age*period,weights=weight,data=hocb_comb)
Anova(hocb_lm)
## Interaction term P=0.265 suggesting no difference in slope (Z)

## Extract mortality estimates
hocb_coef=coef(hocb_lm)
hocb_pre_Z=-1*hocb_coef[["age"]]
hocb_post_Z=-1*(hocb_coef[["age"]]+hocb_coef[["age:periodpost"]])
hocb_pre_A=1-exp(-hocb_pre_Z)
hocb_post_A=1-exp(-hocb_post_Z)

## Summarize into table
hocb_out=bind_cols(impd="hocb",
                   trt="control",
                   Fval=1.572,
                   Pval=0.265,
                   pre_Z=hocb_pre_Z,
                   pre_A=hocb_pre_A,
                   post_Z=hocb_post_Z,
                   post_A=hocb_post_A)

#################################################################
#################################################################
## ZPCB
zpcb_cc=blg_cc%>%
  filter(impd=="ZPCB")

###################
## Set up pre data
zpcb_pre=blg_cc%>%
  filter(impd=="ZPCB" & period=="pre")

## Fit weighted catch curve to pre data
zpcb_pre_cc=catchCurve(ct~age,zpcb_pre,weighted=T)

## Extract variables for plotting pre data
zpcb_pre_ccvars=bind_cols(age=zpcb_pre_cc$age.e,
                          ln_catch=zpcb_pre_cc$log.catch.e,
                          weight=zpcb_pre_cc$weights.e,
                          impd="zpcb",
                          period="pre")

###################
## Set up post data
zpcb_post=blg_cc%>%
  filter(impd=="ZPCB" & period=="post")

## Fit weighted catch curve to post data
zpcb_post_cc=catchCurve(ct~age,zpcb_post,weighted=T)

## Extract variables for plotting post data
zpcb_post_ccvars=bind_cols(age=zpcb_post_cc$age.e,
                           ln_catch=zpcb_post_cc$log.catch.e,
                           weight=zpcb_post_cc$weights.e,
                           impd="zpcb",
                           period="post")

###################
## Combine pre and post data
zpcb_comb=bind_rows(zpcb_pre_ccvars,zpcb_post_ccvars)%>%
  mutate(catch=exp(ln_catch),
         period=factor(period,levels=c("pre","post")),
         trt="control")

###################
## Weighted catch curve using DVR
## Page 211 IFAR
zpcb_lm=lm(ln_catch~age*period,weights=weight,data=zpcb_comb)
Anova(zpcb_lm)
## Interaction term P=0.046 suggesting difference in slope (Z)

## Extract mortality estimates
zpcb_coef=coef(zpcb_lm)
zpcb_pre_Z=-1*zpcb_coef[["age"]]
zpcb_post_Z=-1*(zpcb_coef[["age"]]+zpcb_coef[["age:periodpost"]])
zpcb_pre_A=1-exp(-zpcb_pre_Z)
zpcb_post_A=1-exp(-zpcb_post_Z)

## Summarize into table
zpcb_out=bind_cols(impd="zpcb",
                   trt="control",
                   Fval=8.211,
                   Pval=0.046,
                   pre_Z=zpcb_pre_Z,
                   pre_A=zpcb_pre_A,
                   post_Z=zpcb_post_Z,
                   post_A=zpcb_post_A)

#################################################################
#################################################################
## POCL
pocl_cc=blg_cc%>%
  filter(impd=="POCL")

###################
## Set up pre data
pocl_pre=blg_cc%>%
  filter(impd=="POCL" & period=="pre")

## Fit weighted catch curve to pre data
pocl_pre_cc=catchCurve(ct~age,pocl_pre,weighted=T)

## Extract variables for plotting pre data
pocl_pre_ccvars=bind_cols(age=pocl_pre_cc$age.e,
                          ln_catch=pocl_pre_cc$log.catch.e,
                          weight=pocl_pre_cc$weights.e,
                          impd="pocl",
                          period="pre")

###################
## Set up post data
pocl_post=blg_cc%>%
  filter(impd=="POCL" & period=="post")

## Fit weighted catch curve to post data
pocl_post_cc=catchCurve(ct~age,pocl_post,weighted=T)

## Extract variables for plotting post data
pocl_post_ccvars=bind_cols(age=pocl_post_cc$age.e,
                           ln_catch=pocl_post_cc$log.catch.e,
                           weight=pocl_post_cc$weights.e,
                           impd="pocl",
                           period="post")

###################
## Combine pre and post data
pocl_comb=bind_rows(pocl_pre_ccvars,pocl_post_ccvars)%>%
  mutate(catch=exp(ln_catch),
         period=factor(period,levels=c("pre","post")),
         trt="control")

###################
## Weighted catch curve using DVR
## Page 211 IFAR
pocl_lm=lm(ln_catch~age*period,weights=weight,data=pocl_comb)
Anova(pocl_lm)
## Interaction term P=0.050 suggesting difference in slope (Z)

## Extract mortality estimates
pocl_coef=coef(pocl_lm)
pocl_pre_Z=-1*pocl_coef[["age"]]
pocl_post_Z=-1*(pocl_coef[["age"]]+pocl_coef[["age:periodpost"]])
pocl_pre_A=1-exp(-pocl_pre_Z)
pocl_post_A=1-exp(-pocl_post_Z)

## Summarize into table
pocl_out=bind_cols(impd="pocl",
                   trt="control",
                   Fval=7.719,
                   Pval=0.050,
                   pre_Z=pocl_pre_Z,
                   pre_A=pocl_pre_A,
                   post_Z=pocl_post_Z,
                   post_A=pocl_post_A)

#################################################################
#################################################################
## LXLX
lxlx_cc=blg_cc%>%
  filter(impd=="LXLX")

###################
## Set up pre data
lxlx_pre=blg_cc%>%
  filter(impd=="LXLX" & period=="pre")

## Fit weighted catch curve to pre data
lxlx_pre_cc=catchCurve(ct~age,lxlx_pre,weighted=T)

## Extract variables for plotting pre data
lxlx_pre_ccvars=bind_cols(age=lxlx_pre_cc$age.e,
                          ln_catch=lxlx_pre_cc$log.catch.e,
                          weight=lxlx_pre_cc$weights.e,
                          impd="lxlx",
                          period="pre")

###################
## Set up post data
lxlx_post=blg_cc%>%
  filter(impd=="LXLX" & period=="post")

## Fit weighted catch curve to post data
lxlx_post_cc=catchCurve(ct~age,lxlx_post,weighted=T)

## Extract variables for plotting post data
lxlx_post_ccvars=bind_cols(age=lxlx_post_cc$age.e,
                           ln_catch=lxlx_post_cc$log.catch.e,
                           weight=lxlx_post_cc$weights.e,
                           impd="lxlx",
                           period="post")

###################
## Combine pre and post data
lxlx_comb=bind_rows(lxlx_pre_ccvars,lxlx_post_ccvars)%>%
  mutate(catch=exp(ln_catch),
         period=factor(period,levels=c("pre","post")),
         trt="regulation")

###################
## Weighted catch curve using DVR
## Page 211 IFAR
lxlx_lm=lm(ln_catch~age*period,weights=weight,data=lxlx_comb)
Anova(lxlx_lm)
## Interaction term P=0.241 suggesting no difference in slope (Z)

## Extract mortality estimates
lxlx_coef=coef(lxlx_lm)
lxlx_pre_Z=-1*lxlx_coef[["age"]]
lxlx_post_Z=-1*(lxlx_coef[["age"]]+lxlx_coef[["age:periodpost"]])
lxlx_pre_A=1-exp(-lxlx_pre_Z)
lxlx_post_A=1-exp(-lxlx_post_Z)

## Summarize into table
lxlx_out=bind_cols(impd="lxlx",
                   trt="regulation",
                   Fval=2.119,
                   Pval=0.241,
                   pre_Z=lxlx_pre_Z,
                   pre_A=lxlx_pre_A,
                   post_Z=lxlx_post_Z,
                   post_A=lxlx_post_A)

#################################################################
#################################################################
## PTSB
ptsb_cc=blg_cc%>%
  filter(impd=="PTSB")

###################
## Set up pre data
ptsb_pre=blg_cc%>%
  filter(impd=="PTSB" & period=="pre")

## Fit weighted catch curve to pre data
ptsb_pre_cc=catchCurve(ct~age,ptsb_pre,weighted=T)

## Extract variables for plotting pre data
ptsb_pre_ccvars=bind_cols(age=ptsb_pre_cc$age.e,
                          ln_catch=ptsb_pre_cc$log.catch.e,
                          weight=ptsb_pre_cc$weights.e,
                          impd="ptsb",
                          period="pre")

###################
## Set up post data
ptsb_post=blg_cc%>%
  filter(impd=="PTSB" & period=="post")

## Fit weighted catch curve to post data
ptsb_post_cc=catchCurve(ct~age,ptsb_post,weighted=T)

## Extract variables for plotting post data
ptsb_post_ccvars=bind_cols(age=ptsb_post_cc$age.e,
                           ln_catch=ptsb_post_cc$log.catch.e,
                           weight=ptsb_post_cc$weights.e,
                           impd="ptsb",
                           period="post")

###################
## Combine pre and post data
ptsb_comb=bind_rows(ptsb_pre_ccvars,ptsb_post_ccvars)%>%
  mutate(catch=exp(ln_catch),
         period=factor(period,levels=c("pre","post")),
         trt="regulation")

###################
## Weighted catch curve using DVR
## Page 211 IFAR
ptsb_lm=lm(ln_catch~age*period,weights=weight,data=ptsb_comb)
Anova(ptsb_lm)
## Interaction term P=0.665 suggesting no difference in slope (Z)

## Extract mortality estimates
ptsb_coef=coef(ptsb_lm)
ptsb_pre_Z=-1*ptsb_coef[["age"]]
ptsb_post_Z=-1*(ptsb_coef[["age"]]+ptsb_coef[["age:periodpost"]])
ptsb_pre_A=1-exp(-ptsb_pre_Z)
ptsb_post_A=1-exp(-ptsb_post_Z)

## Summarize into table
ptsb_out=bind_cols(impd="ptsb",
                   trt="regulation",
                   Fval=0.208,
                   Pval=0.665,
                   pre_Z=ptsb_pre_Z,
                   pre_A=ptsb_pre_A,
                   post_Z=ptsb_post_Z,
                   post_A=ptsb_post_A)

#################################################################
#################################################################
## JWSL
jwsl_cc=blg_cc%>%
  filter(impd=="JWSL")

###################
## Set up pre data
jwsl_pre=blg_cc%>%
  filter(impd=="JWSL" & period=="pre")

## Fit weighted catch curve to pre data
jwsl_pre_cc=catchCurve(ct~age,jwsl_pre,weighted=T)

## Extract variables for plotting pre data
jwsl_pre_ccvars=bind_cols(age=jwsl_pre_cc$age.e,
                          ln_catch=jwsl_pre_cc$log.catch.e,
                          weight=jwsl_pre_cc$weights.e,
                          impd="jwsl",
                          period="pre")

###################
## Set up post data
jwsl_post=blg_cc%>%
  filter(impd=="JWSL" & period=="post")

## Fit weighted catch curve to post data
jwsl_post_cc=catchCurve(ct~age,jwsl_post,weighted=T)

## Extract variables for plotting post data
jwsl_post_ccvars=bind_cols(age=jwsl_post_cc$age.e,
                           ln_catch=jwsl_post_cc$log.catch.e,
                           weight=jwsl_post_cc$weights.e,
                           impd="jwsl",
                           period="post")

###################
## Combine pre and post data
jwsl_comb=bind_rows(jwsl_pre_ccvars,jwsl_post_ccvars)%>%
  mutate(catch=exp(ln_catch),
         period=factor(period,levels=c("pre","post")),
         trt="regulation")

###################
## Weighted catch curve using DVR
## Page 211 IFAR
jwsl_lm=lm(ln_catch~age*period,weights=weight,data=jwsl_comb)
Anova(jwsl_lm)
## Interaction term P=0.032 suggesting difference in slope (Z)

## Extract mortality estimates
jwsl_coef=coef(jwsl_lm)
jwsl_pre_Z=-1*jwsl_coef[["age"]]
jwsl_post_Z=-1*(jwsl_coef[["age"]]+jwsl_coef[["age:periodpost"]])
jwsl_pre_A=1-exp(-jwsl_pre_Z)
jwsl_post_A=1-exp(-jwsl_post_Z)

## Summarize into table
jwsl_out=bind_cols(impd="jwsl",
                   trt="regulation",
                   Fval=6.469,
                   Pval=0.032,
                   pre_Z=jwsl_pre_Z,
                   pre_A=jwsl_pre_A,
                   post_Z=jwsl_post_Z,
                   post_A=jwsl_post_A)

#################################################################
#################################################################
## MISL
misl_cc=blg_cc%>%
  filter(impd=="MISL")

###################
## Set up pre data
misl_pre=blg_cc%>%
  filter(impd=="MISL" & period=="pre")

## Fit weighted catch curve to pre data
misl_pre_cc=catchCurve(ct~age,misl_pre,weighted=T)

## Extract variables for plotting pre data
misl_pre_ccvars=bind_cols(age=misl_pre_cc$age.e,
                          ln_catch=misl_pre_cc$log.catch.e,
                          weight=misl_pre_cc$weights.e,
                          impd="misl",
                          period="pre")

###################
## Set up post data
misl_post=blg_cc%>%
  filter(impd=="MISL" & period=="post")

## Fit weighted catch curve to post data
misl_post_cc=catchCurve(ct~age,misl_post,weighted=T)

## Extract variables for plotting post data
misl_post_ccvars=bind_cols(age=misl_post_cc$age.e,
                           ln_catch=misl_post_cc$log.catch.e,
                           weight=misl_post_cc$weights.e,
                           impd="misl",
                           period="post")

###################
## Combine pre and post data
misl_comb=bind_rows(misl_pre_ccvars,misl_post_ccvars)%>%
  mutate(catch=exp(ln_catch),
         period=factor(period,levels=c("pre","post")),
         trt="regulation")

###################
## Weighted catch curve using DVR
## Page 211 IFAR
misl_lm=lm(ln_catch~age*period,weights=weight,data=misl_comb)
Anova(misl_lm)
## Interaction term P=0.032 suggesting difference in slope (Z)

## Extract mortality estimates
misl_coef=coef(misl_lm)
misl_pre_Z=-1*misl_coef[["age"]]
misl_post_Z=-1*(misl_coef[["age"]]+misl_coef[["age:periodpost"]])
misl_pre_A=1-exp(-misl_pre_Z)
misl_post_A=1-exp(-misl_post_Z)

## Summarize into table
misl_out=bind_cols(impd="misl",
                   trt="regulation",
                   Fval=9.018,
                   Pval=0.024,
                   pre_Z=misl_pre_Z,
                   pre_A=misl_pre_A,
                   post_Z=misl_post_Z,
                   post_A=misl_post_A)

#################################################################
## Combined annual mortality estimates with DVR test statistics
out=bind_rows(gdcl_out,hocb_out,zpcb_out,pocl_out,
              lxlx_out,ptsb_out,jwsl_out,misl_out)%>%
  mutate(across(where(is.numeric),~num(.,digits = 3)))

#################################################################
#################################################################
## Create plots
##GDCL
gdcl_plot=ggplot(gdcl_comb,aes(x=age,y=ln_catch))+
  geom_smooth(aes(color=period,fill=period),method="lm")+
  geom_point(aes(shape=period,color=period),size=4)+
  scale_shape_manual(labels=c("Pre: A = 46.2%",
                              "Post: A = 33.3%"),
                     values=c(15,19))+
  scale_color_manual(labels=c("Pre: A = 46.2%",
                              "Post: A = 33.3%"),
                     values=c("#FF474C","#475F94"))+
  scale_fill_manual(labels=c("Pre: A = 46.2%",
                             "Post: A = 33.3%"),
                    values=c("#FF474C","#475F94"))+
  scale_y_continuous(breaks=seq(0,8,1),
                     name="")+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="")+
  coord_cartesian(xlim=c(0.9,8.1),
                  ylim=c(0,8.01),
                  expand=F)+
  annotate("text",label="Gardner",x=1,y=8,hjust=0,vjust=1,size=8)+
  annotate("text",label=paste("italic(F)","==",gdcl_out$Fval),parse=T,x=1,y=7.4,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",gdcl_out$Pval),parse=T,x=1,y=7.05,hjust=0,vjust=1,size=4)+
  pubtheme

## HOCB
hocb_plot=ggplot(hocb_comb,aes(x=age,y=ln_catch))+
  geom_smooth(aes(color=period,fill=period),method="lm")+
  geom_point(aes(shape=period,color=period),size=4)+
  scale_shape_manual(labels=c("Pre: A = 65.6%",
                              "Post: A = 37.8%"),
                     values=c(15,19))+
  scale_color_manual(labels=c("Pre: A = 65.6%",
                              "Post: A = 37.8%"),
                     values=c("#FF474C","#475F94"))+
  scale_fill_manual(labels=c("Pre: A = 65.6%",
                             "Post: A = 37.8%"),
                    values=c("#FF474C","#475F94"))+
  scale_y_continuous(breaks=seq(0,8,1),
                     name="")+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="")+
  coord_cartesian(xlim=c(0.9,8.1),
                  ylim=c(0,8.01),
                  expand=F)+
  annotate("text",label="Holton",x=1,y=8,hjust=0,vjust=1,size=8)+
  annotate("text",label=paste("italic(F)","==",hocb_out$Fval),parse=T,x=1,y=7.4,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",hocb_out$Pval),parse=T,x=1,y=7.05,hjust=0,vjust=1,size=4)+
  pubtheme

## ZPCB
zpcb_plot=ggplot(zpcb_comb,aes(x=age,y=ln_catch))+
  geom_smooth(aes(color=period,fill=period),method="lm")+
  geom_point(aes(shape=period,color=period),size=4)+
  scale_shape_manual(labels=c("Pre: A = 92.9%",
                              "Post: A = 12.7%"),
                     values=c(15,19))+
  scale_color_manual(labels=c("Pre: A = 92.9%",
                              "Post: A = 12.7%"),
                     values=c("#FF474C","#475F94"))+
  scale_fill_manual(labels=c("Pre: A = 92.9%",
                             "Post: A = 12.7%"),
                    values=c("#FF474C","#475F94"))+
  scale_y_continuous(breaks=seq(0,8,1),
                     name=expression(log[10]*(catch+1)))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="")+
  coord_cartesian(xlim=c(0.9,8.1),
                  ylim=c(0,8.01),
                  expand=F)+
  annotate("text",label="Campbell",x=1,y=8,hjust=0,vjust=1,size=8)+
  annotate("text",label=paste("italic(F)","==",zpcb_out$Fval),parse=T,x=1,y=7.4,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",zpcb_out$Pval),parse=T,x=1,y=7.05,hjust=0,vjust=1,size=4)+
  pubtheme

## POCL
pocl_plot=ggplot(pocl_comb,aes(x=age,y=ln_catch))+
  geom_smooth(aes(color=period,fill=period),method="lm")+
  geom_point(aes(shape=period,color=period),size=4)+
  scale_shape_manual(labels=c("Pre: A = 60.4%",
                              "Post: A = 29.8%"),
                     values=c(15,19))+
  scale_color_manual(labels=c("Pre: A = 60.4%",
                              "Post: A = 29.8%"),
                     values=c("#FF474C","#475F94"))+
  scale_fill_manual(labels=c("Pre: A = 60.4%",
                             "Post: A = 29.8%"),
                    values=c("#FF474C","#475F94"))+
  scale_y_continuous(breaks=seq(0,8,1),
                     name="")+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="")+
  coord_cartesian(xlim=c(0.9,8.1),
                  ylim=c(0,8.01),
                  expand=F)+
  annotate("text",label="Paola",x=1,y=8,hjust=0,vjust=1,size=8)+
  annotate("text",label=paste("italic(F)","==",pocl_out$Fval),parse=T,x=1,y=7.4,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","=='0.050'"),parse=T,x=1,y=7.05,hjust=0,vjust=1,size=4)+
  pubtheme

## LXLX
lxlx_plot=ggplot(lxlx_comb,aes(x=age,y=ln_catch))+
  geom_smooth(aes(color=period,fill=period),method="lm")+
  geom_point(aes(shape=period,color=period),size=4)+
  scale_shape_manual(labels=c("Pre: A = 72.0%",
                              "Post: A = 51.1%"),
                     values=c(15,19))+
  scale_color_manual(labels=c("Pre: A = 72.0%",
                              "Post: A = 51.1%"),
                     values=c("#FF474C","#475F94"))+
  scale_fill_manual(labels=c("Pre: A = 72.0%",
                             "Post: A = 51.1%"),
                    values=c("#FF474C","#475F94"))+
  scale_y_continuous(breaks=seq(0,8,1),
                     name="")+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  coord_cartesian(xlim=c(0.9,8.1),
                  ylim=c(0,8.01),
                  expand=F)+
  annotate("text",label="Lenexa",x=1,y=8,hjust=0,vjust=1,size=8)+
  annotate("text",label=paste("italic(F)","==",lxlx_out$Fval),parse=T,x=1,y=7.4,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",lxlx_out$Pval),parse=T,x=1,y=7.05,hjust=0,vjust=1,size=4)+
  pubtheme

## PTSB
ptsb_plot=ggplot(ptsb_comb,aes(x=age,y=ln_catch))+
  geom_smooth(aes(color=period,fill=period),method="lm")+
  geom_point(aes(shape=period,color=period),size=4)+
  scale_shape_manual(labels=c("Pre: A = 60.5%",
                              "Post: A = 52.5%"),
                     values=c(15,19))+
  scale_color_manual(labels=c("Pre: A = 60.5%",
                              "Post: A = 52.5%"),
                     values=c("#FF474C","#475F94"))+
  scale_fill_manual(labels=c("Pre: A = 60.5%",
                             "Post: A = 52.5%"),
                    values=c("#FF474C","#475F94"))+
  scale_y_continuous(breaks=seq(0,8,1),
                     name="")+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  coord_cartesian(xlim=c(0.9,8.1),
                  ylim=c(0,8.01),
                  expand=F)+
  annotate("text",label="Pott #2",x=1,y=8,hjust=0,vjust=1,size=8)+
  annotate("text",label=paste("italic(F)","==",ptsb_out$Fval),parse=T,x=1,y=7.4,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",ptsb_out$Pval),parse=T,x=1,y=7.05,hjust=0,vjust=1,size=4)+
  pubtheme

## JWSL
jwsl_plot=ggplot(jwsl_comb,aes(x=age,y=ln_catch))+
  geom_smooth(aes(color=period,fill=period),method="lm")+
  geom_point(aes(shape=period,color=period),size=4)+
  scale_shape_manual(labels=c("Pre: A = 6.1%",
                              "Post: A = 51.0%"),
                     values=c(15,19))+
  scale_color_manual(labels=c("Pre: A = 6.1%",
                              "Post: A = 51.0%"),
                     values=c("#FF474C","#475F94"))+
  scale_fill_manual(labels=c("Pre: A = 6.1%",
                             "Post: A = 51.0%"),
                    values=c("#FF474C","#475F94"))+
  scale_y_continuous(breaks=seq(0,8,1),
                     name=expression(log[10]*(catch+1)))+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  coord_cartesian(xlim=c(0.9,8.1),
                  ylim=c(0,8.01),
                  expand=F)+
  annotate("text",label="Jewell",x=1,y=8,hjust=0,vjust=1,size=8)+
  annotate("text",label=paste("italic(F)","==",jwsl_out$Fval),parse=T,x=1,y=7.4,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",jwsl_out$Pval),parse=T,x=1,y=7.05,hjust=0,vjust=1,size=4)+
  pubtheme

## MISL
misl_plot=ggplot(misl_comb,aes(x=age,y=ln_catch))+
  geom_smooth(aes(color=period,fill=period),method="lm")+
  geom_point(aes(shape=period,color=period),size=4)+
  scale_shape_manual(labels=c("Pre: A = 61.1%",
                              "Post: A = 38.3%"),
                     values=c(15,19))+
  scale_color_manual(labels=c("Pre: A = 61.1%",
                              "Post: A = 38.3%"),
                     values=c("#FF474C","#475F94"))+
  scale_fill_manual(labels=c("Pre: A = 61.1%",
                             "Post: A = 38.3%"),
                    values=c("#FF474C","#475F94"))+
  scale_y_continuous(breaks=seq(0,8,1),
                     name="")+
  scale_x_continuous(breaks=seq(0,8,1),
                     name="Estimated age")+
  coord_cartesian(xlim=c(0.9,8.1),
                  ylim=c(0,8.01),
                  expand=F)+
  annotate("text",label="Miami",x=1,y=8,hjust=0,vjust=1,size=8)+
  annotate("text",label=paste("italic(F)","==",misl_out$Fval),parse=T,x=1,y=7.4,hjust=0,vjust=1,size=4)+
  annotate("text",label=paste("italic(P)","==",misl_out$Pval),parse=T,x=1,y=7.05,hjust=0,vjust=1,size=4)+
  pubtheme

#################################################################
#################################################################
## Combine plots and export
plots=(zpcb_plot|gdcl_plot|hocb_plot|pocl_plot)/(jwsl_plot|lxlx_plot|misl_plot|ptsb_plot)
ggsave(plot=plots,"mortality.png",width=14,height=8,units="in",bg="white")