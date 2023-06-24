library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)
library(tidyverse)
library(lubridate)
library(forestplot)
library(extrafont)


mfr2=fread('/Applications/plab/data/swed/ver/mfr_150202_recored_main_filtered_F.csv')


mfr2=subset(mfr2, AR>1989)

grav1=subset(mfr2, parity_clean=="1")
grav1=subset(grav1, GRDBS>="259")	

grav2=subset(mfr2, parity_clean=="2")

temp=subset(grav2, OFRIABEF=="1"|OFRISTIM=="1"|OFRIKIRU=="1"|OFRIICSI=="1"|OFRIANN=="1")
grav2=subset(grav2, !lpnr_mor %in% temp$lpnr_mor)
rm(temp)

grav1=subset(grav1, lpnr_mor!="NA")			
grav2=subset(grav2, lpnr_mor!="NA")

grav1=grav1[!duplicated(grav1$lpnr_mor),]
grav2=grav2[!duplicated(grav2$lpnr_mor),]

grav1=semi_join(grav1,grav2,by="lpnr_mor")
grav2=semi_join(grav2,grav1,by="lpnr_mor")

grav1=grav1[order(grav1$lpnr_mor),]
grav2=grav2[order(grav2$lpnr_mor),]

mydata=data.frame(lpnr_mor=grav1$lpnr_mor)

mydata=mydata %>% add_column(GRDBSpreg1=grav1$GRDBS, GRDBSpreg2=grav2$GRDBS)

temp1=subset(grav1,SECMARK=="1")
temp2=subset(grav2,TSECTIO=="1")
temp3=rbind(temp1,temp2)
temp4=temp3[!duplicated(temp3$lpnr_mor),]

mydata$CSpreg1 = mydata$lpnr_mor %in% temp4$lpnr_mor
rm(temp1,temp2,temp3,temp4)

temp=subset(grav2, SECMARK=="1")
mydata$CSpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav1, FLSPONT=="1")
mydata$spontpreg1= mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, FLSPONT=="1")
mydata$spontpreg2= mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)


mydata$PTL = mydata$GRDBSpreg2<259
mydata$PTL1 = mydata$GRDBSpreg2<196
mydata$PTL2 = mydata$GRDBSpreg2>=196 & mydata$GRDBSpreg2<224
mydata$PTL3 = mydata$GRDBSpreg2>=224 & mydata$GRDBSpreg2<238
mydata$PTL4 = mydata$GRDBSpreg2>=238 & mydata$GRDBSpreg2<259
mydata$w37 = mydata$GRDBSpreg2>=259 & mydata$GRDBSpreg2<266
mydata$w38 = mydata$GRDBSpreg2>=266 & mydata$GRDBSpreg2<273
mydata$opt = mydata$GRDBSpreg2>=273 & mydata$GRDBSpreg2<287
mydata$w41 = mydata$GRDBSpreg2>=287 & mydata$GRDBSpreg2<294
mydata$w42 = mydata$GRDBSpreg2>=294


#below is for adjusting for confounders
#year of birth, maternal age, dates of birth, birth weight

mydata=mydata %>% add_column(ARpreg1=grav1$AR, ARpreg2=grav2$AR,
			     MALDERpreg1=grav1$MALDER, MALDERpreg2=grav2$MALDER,
			     BFODDATpreg1=grav1$BFODDAT, BFODDATpreg2=grav2$BFODDAT,
			     BVIKTpreg1=grav1$BVIKT, BVIKTpreg2=grav2$BVIKT)

#interpregnancy interval
mydata$BFODDATpreg1=ymd(mydata$BFODDATpreg1)
mydata$BFODDATpreg2=ymd(mydata$BFODDATpreg2)
mydata=mydata %>% add_column(intpregint=difftime(mydata$BFODDATpreg2, mydata$BFODDATpreg1, days)-mydata$GRDBSpreg2)
mydata$intpregint=as.numeric(mydata$intpregint, units="days")
mydata=mydata %>% add_column(intpregintmonths=mydata$intpregint/30.43)


#BMI

mydata=mydata %>% add_column(BMIpreg1=grav1$MVIKT/((grav1$MLANGD/100)^2), BMIpreg2=grav2$MVIKT/((grav2$MLANGD/100)^2))

mydata$BMIpreg1.2 = ifelse(is.na(mydata$BMIpreg1), grav1$MVIKT/((grav2$MLANGD/100)^2), mydata$BMIpreg1)
mydata$BMIpreg2.2 = ifelse(is.na(mydata$BMIpreg2), grav2$MVIKT/((grav1$MLANGD/100)^2), mydata$BMIpreg2)
mydata=subset(mydata, select=-c(BMIpreg1,BMIpreg2))
mydata=mydata %>% rename(BMIpreg1=BMIpreg1.2)
mydata=mydata %>% rename(BMIpreg2=BMIpreg2.2)
mydata$BMIpreg1=formatC(mydata$BMIpreg1, digits=1, format="f")
mydata$BMIpreg2=formatC(mydata$BMIpreg2, digits=1, format="f")

mydata$BMIpreg1[mydata$BMIpreg1=="NA"]=23.47
mydata$BMIpreg2[mydata$BMIpreg2=="NA"]=24.12
mydata$BMIpreg1=as.numeric(mydata$BMIpreg1)
mydata$BMIpreg2=as.numeric(mydata$BMIpreg2)


#BMI groups

mydata$BMIpreg1.1 = mydata$BMIpreg1<18.5
mydata$BMIpreg1.2 = mydata$BMIpreg1>=18.5 & mydata$BMIpreg1<25.0
mydata$BMIpreg1.3 = mydata$BMIpreg1>=25.0 & mydata$BMIpreg1<30.0
mydata$BMIpreg1.4 = mydata$BMIpreg1>=30.0 & mydata$BMIpreg1<35.0
mydata$BMIpreg1.5 = mydata$BMIpreg1>=35.0 & mydata$BMIpreg1<40.0
mydata$BMIpreg1.6 = mydata$BMIpreg1>=40.0 & mydata$BMIpreg1<99

mydata$BMIpreg2.1 = mydata$BMIpreg2<18.5
mydata$BMIpreg2.2 = mydata$BMIpreg2>=18.5 & mydata$BMIpreg2<25.0
mydata$BMIpreg2.3 = mydata$BMIpreg2>=25.0 & mydata$BMIpreg2<30.0
mydata$BMIpreg2.4 = mydata$BMIpreg2>=30.0 & mydata$BMIpreg2<35.0
mydata$BMIpreg2.5 = mydata$BMIpreg2>=35.0 & mydata$BMIpreg2<40.0
mydata$BMIpreg2.6 = mydata$BMIpreg2>=40.0 & mydata$BMIpreg2<99


#BMI increase and BMI decrease

temp1=subset(mydata, BMIpreg1.1|BMIpreg1.2|BMIpreg1.3|BMIpreg1.4|BMIpreg1.5|BMIpreg1.6)
temp2=subset(mydata, BMIpreg2.1|BMIpreg2.2|BMIpreg2.3|BMIpreg2.4|BMIpreg2.5|BMIpreg2.6)

temp1=semi_join(temp1,temp2,by="lpnr_mor")
temp2=semi_join(temp2,temp1,by="lpnr_mor")

# TODO check about these 99 placeholders here!!
temp3=data.frame(lpnr_mor=temp1$lpnr_mor,
                                  BMIscorepreg1=ifelse(temp1$BMIpreg1.1,1,
						ifelse(temp1$BMIpreg1.2,2,
						ifelse(temp1$BMIpreg1.3,3,
						ifelse(temp1$BMIpreg1.4,4,
						ifelse(temp1$BMIpreg1.5,5,
						ifelse(temp1$BMIpreg1.6,6,99)))))),
                                  BMIscorepreg2=ifelse(temp2$BMIpreg2.1,1,
						ifelse(temp2$BMIpreg2.2,2,
						ifelse(temp2$BMIpreg2.3,3,
						ifelse(temp2$BMIpreg2.4,4,
						ifelse(temp2$BMIpreg2.5,5,
						ifelse(temp2$BMIpreg2.6,6,99))))))
)

temp3=temp3 %>% add_column(BMIdiff=temp3$BMIscorepreg2-temp3$BMIscorepreg1)
temp3=temp3 %>% add_column(BMIinc=temp3$BMIdiff>0, BMIdec=temp3$BMIdiff<0)

temp4=subset(temp3, BMIinc)
temp5=subset(temp3, BMIdec)

mydata$BMIinc = mydata$lpnr_mor %in% temp4$lpnr_mor
mydata$BMIdec = mydata$lpnr_mor %in% temp5$lpnr_mor

rm(temp1,temp2,temp3,temp4,temp5)


#not born in sweden, unemployed, smoking
temp=subset(grav1, MFODLAND!="SVERIGE")
mydata$notborninswepreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, MFODLAND!="SVERIGE")
mydata$notborninswepreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav1, ARBETE=="3")
mydata$unemployedpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, ARBETE=="3")
mydata$unemployedpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav1, ROK1!="1")
mydata$smokepreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, ROK1!="1")
mydata$smokepreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)


#diabetes, hypertension, gestational hypertension, preeclampsia, allHTN

# TODO fix this filter to not use any_vars!
temp1=subset(grav1, DIABETES=="1"|DIABETES=="2")
temp2=grav1 %>% filter_all(any_vars(. %in% c('E107','O240','O241','O243','O244','O249','O240B','O240C','O240D','O240E','O240F','O240X','O244A','O244B','O244X','25000','25009','76110','250A','250B','250C','250D','250E','250F','250X','648A')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$diabetespreg1 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)

temp1=subset(grav2, DIABETES=="1"|DIABETES=="2")
temp2=grav2 %>% filter_all(any_vars(. %in% c('E107','O240','O241','O243','O244','O249','O240B','O240C','O240D','O240E','O240F','O240X','O244A','O244B','O244X','25000','25009','76110','250A','250B','250C','250D','250E','250F','250X','648A')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$diabetespreg2 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)

temp1=subset(grav1, HYPERTON=="1"|HYPERTON=="2")
temp2=grav1 %>% filter_all(any_vars(. %in% c('I10','I109','401','401X','40199')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$HTNpreg1 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)

temp1=subset(grav2, HYPERTON=="1"|HYPERTON=="2")
temp2=grav2 %>% filter_all(any_vars(. %in% c('I10','I109','401','401X','40199')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$HTNpreg2 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)


temp=grav1 %>% filter_all(any_vars(. %in% c('O139','642','642A','642B','642C','642D','642X','63701')))
mydata$GHTNpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=grav2 %>% filter_all(any_vars(. %in% c('O139','642','642A','642B','642C','642D','642X','63701')))
mydata$GHTNpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=grav1 %>% filter_all(any_vars(. %in% c('O14','O140','O141','O142','O149','O141A','O141B','O141X','642E','642F','642G','642H','63703','63704','63709','63710','63799','76210','76220','76230')))
mydata$PEpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=grav2 %>% filter_all(any_vars(. %in% c('O14','O140','O141','O142','O149','O141A','O141B','O141X','642E','642F','642G','642H','63703','63704','63709','63710','63799','76210','76220','76230')))
mydata$PEpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)


mydata$allHTNpreg1 = mydata$HTNpreg1 | mydata$GHTNpreg1 | mydata$PEpreg1

mydata$allHTNpreg2 = mydata$HTNpreg2 | mydata$GHTNpreg2 | mydata$PEpreg2


#boy child, fetal malformation

temp=subset(grav1, KON=="1")
mydata$boychildpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, KON=="1")
mydata$boychildpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav1, MISSB=="1")
mydata$malformpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, MISSB=="1")
mydata$malformpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)



#SGA and LGA according to Marsal

marsalpreg1 = function(GRDBSpreg1,BVIKTpreg1,boychildpreg1){
    if(boychildpreg1){
        mw=-1.907345e-6*GRDBSpreg1^4 + 1.140644e-3*GRDBSpreg1^3 - 1.336265e-1*GRDBSpreg1^2 + 1.976961*GRDBSpreg1 + 2.410053e+2
    } else {
        mw=-2.761948e-6*GRDBSpreg1^4 + 1.744841e-3*GRDBSpreg1^3 - 2.893626e-1*GRDBSpreg1^2 + 18.91197*GRDBSpreg1 - 4.135122e+2
    }
    pnorm(BVIKTpreg1, mw, abs(0.12*mw))*100
}

mydata$PCTmarsalpreg1 = sapply(1:nrow(mydata), function(x) marsalpreg1(mydata$GRDBSpreg1[x],mydata$BVIKTpreg1[x],mydata$boychildpreg1[x]))

mydata$SGApreg1marsal= mydata$PCTmarsalpreg1 < pnorm(-2) *100
mydata$LGApreg1marsal= mydata$PCTmarsalpreg1 > pnorm(2) *100


mydata$PCTmarsalpreg2 = sapply(1:nrow(mydata), function(x) marsalpreg1(mydata$GRDBSpreg2[x],mydata$BVIKTpreg2[x],mydata$boychildpreg2[x]))

mydata$SGApreg2marsal= mydata$PCTmarsalpreg2 < pnorm(-2) *100
mydata$LGApreg2marsal= mydata$PCTmarsalpreg2 > pnorm(2) *100

rm(marsalpreg1,marsalpreg2)


#GA according to ultrasound
temp=subset(grav2, GRMETOD==1|GRMETOD==5|GRMETOD==6|GRMETOD==7)
mydata$GAUSpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)


#df2 descriptive

temp1=subset(mydata,spontpreg2)
mydatavag=subset(temp1,!CSpreg1)
mydatasec=subset(temp1,CSpreg1)


#continuous variables


#maternal age preg 1

mean(mydatavag$MALDERpreg1)
sd(mydatavag$MALDERpreg1)
mean(mydatasec$MALDERpreg1)
sd(mydatasec$MALDERpreg1)


#maternal age preg 2

mean(mydatavag$MALDERpreg2)
sd(mydatavag$MALDERpreg2)
mean(mydatasec$MALDERpreg2)
sd(mydatasec$MALDERpreg2)


#gestational age preg 1

mean(mydatavag$GRDBSpreg1)
sd(mydatavag$GRDBSpreg1)
mean(mydatasec$GRDBSpreg1)
sd(mydatasec$GRDBSpreg1)


#gestational age preg 2

mean(mydatavag$GRDBSpreg2)
sd(mydatavag$GRDBSpreg2)
mean(mydatasec$GRDBSpreg2)
sd(mydatasec$GRDBSpreg2)


#interpregnancy interval (months)

mean(mydatavag$intpregintmonths)
sd(mydatavag$intpregintmonths)
mean(mydatasec$intpregintmonths)
sd(mydatasec$intpregintmonths)


#BMI preg 1

mean(mydatavag$BMIpreg1)
sd(mydatavag$BMIpreg1)
mean(mydatasec$BMIpreg1)
sd(mydatasec$BMIpreg1)

#BMI preg 2

mean(mydatavag$BMIpreg2)
sd(mydatavag$BMIpreg2)
mean(mydatasec$BMIpreg2)
sd(mydatasec$BMIpreg2)


#categorical variables

#Preterm birth

sum(mydatavag$PTL)/nrow(mydatavag)
sum(mydatasec$PTL)/nrow(mydatasec)


#PTL <28+0

sum(mydatavag$PTL1)/nrow(mydatavag)
sum(mydatasec$PTL1)/nrow(mydatasec)


#PTL 28+0 - 31+6

sum(mydatavag$PTL2)/nrow(mydatavag)
sum(mydatasec$PTL2)/nrow(mydatasec)


#PTL 32+0 - 33+6

sum(mydatavag$PTL3)/nrow(mydatavag)
sum(mydatasec$PTL3)/nrow(mydatasec)


#PTL 34+0 - 36+6

sum(mydatavag$PTL4)/nrow(mydatavag)
sum(mydatasec$PTL4)/nrow(mydatasec)


#w37

sum(mydatavag$w37)/nrow(mydatavag)
sum(mydatasec$w37)/nrow(mydatasec)


#w38

sum(mydatavag$w38)/nrow(mydatavag)
sum(mydatasec$w38)/nrow(mydatasec)


#opt

sum(mydatavag$opt)/nrow(mydatavag)
sum(mydatasec$opt)/nrow(mydatasec)


#w41

sum(mydatavag$w41)/nrow(mydatavag)
sum(mydatasec$w41)/nrow(mydatasec)


#w42

sum(mydatavag$w42)/nrow(mydatavag)
sum(mydatasec$w42)/nrow(mydatasec)



#not born in Sweden preg 1

sum(mydatavag$notborninswepreg1)/nrow(mydatavag)
sum(mydatasec$notborninswepreg1)/nrow(mydatasec)




#smoking preg 1

sum(mydatavag$smokepreg1)/nrow(mydatavag)
sum(mydatasec$smokepreg1)/nrow(mydatasec)



#smoking preg 2

sum(mydatavag$smokepreg2)/nrow(mydatavag)
sum(mydatasec$smokepreg2)/nrow(mydatasec)



#diabetes preg 1

sum(mydatavag$diabetespreg1)/nrow(mydatavag)
sum(mydatasec$diabetespreg1)/nrow(mydatasec)



#diabetes preg 2

sum(mydatavag$diabetespreg2)/nrow(mydatavag)
sum(mydatasec$diabetespreg2)/nrow(mydatasec)



#allHTN preg 1

sum(mydatavag$allHTNpreg1)/nrow(mydatavag)
sum(mydatasec$allHTNpreg1)/nrow(mydatasec)



#allHTN preg 2

sum(mydatavag$allHTNpreg2)/nrow(mydatavag)
sum(mydatasec$allHTNpreg2)/nrow(mydatasec)



#boy child preg 1

sum(mydatavag$boychildpreg1)/nrow(mydatavag)
sum(mydatasec$boychildpreg1)/nrow(mydatasec)



#boy child preg 2

sum(mydatavag$boychildpreg2)/nrow(mydatavag)
sum(mydatasec$boychildpreg2)/nrow(mydatasec)



#fetal malformation preg 1

sum(mydatavag$malformpreg1)/nrow(mydatavag)
sum(mydatasec$malformpreg1)/nrow(mydatasec)



#fetal malformation preg 2

sum(mydatavag$malformpreg2)/nrow(mydatavag)
sum(mydatasec$malformpreg2)/nrow(mydatasec)



#SGA preg 1

sum(mydatavag$SGApreg1marsal)/nrow(mydatavag)
sum(mydatasec$SGApreg1marsal)/nrow(mydatasec)



#SGA preg 2

sum(mydatavag$SGApreg2marsal)/nrow(mydatavag)
sum(mydatasec$SGApreg2marsal)/nrow(mydatasec)



#LGA preg 1

sum(mydatavag$LGApreg1marsal)/nrow(mydatavag)
sum(mydatasec$LGApreg1marsal)/nrow(mydatasec)



#LGA preg 2

sum(mydatavag$LGApreg2marsal)/nrow(mydatavag)
sum(mydatasec$LGApreg2marsal)/nrow(mydatasec)



#CS preg 2

sum(mydatavag$CSpreg2)/nrow(mydatavag)
sum(mydatasec$CSpreg2)/nrow(mydatasec)



#GA acc US preg 2

sum(mydatavag$GAUSpreg2)/nrow(mydatavag)
sum(mydatasec$GAUSpreg2)/nrow(mydatasec)



#cleaning

rm(temp1,mydatavag,mydatasec)




#univariable analysis

mydata3=subset(mydata,spontpreg2)



glm1=glm(GRDBSpreg2 ~ MALDERpreg2, family=gaussian(), data=mydata3)
summary(glm1)
confint(glm1)

glm2=glm(GRDBSpreg2 ~ BMIpreg2, family=gaussian(), data=mydata3)
summary(glm2)
confint(glm2)

glm3=glm(GRDBSpreg2 ~ GRDBSpreg1, family=gaussian(), data=mydata3)
summary(glm3)
confint(glm3)

glm4=glm(GRDBSpreg2 ~ notborninswepreg2, family=gaussian(), data=mydata3)
summary(glm4)
confint(glm4)

glm5=glm(GRDBSpreg2 ~ smokepreg2, family=gaussian(), data=mydata3)
summary(glm5)
confint(glm5)

glm6=glm(GRDBSpreg2 ~ diabetespreg2, family=gaussian(), data=mydata3)
summary(glm6)
confint(glm6)

glm7=glm(GRDBSpreg2 ~ allHTNpreg2, family=gaussian(), data=mydata3)
summary(glm7)
confint(glm7)

glm8=glm(GRDBSpreg2 ~ boychildpreg2, family=gaussian(), data=mydata3)
summary(glm8)
confint(glm8)

glm9=glm(GRDBSpreg2 ~ malformpreg2, family=gaussian(), data=mydata3)
summary(glm9)
confint(glm9)

glm10=glm(GRDBSpreg2 ~ SGApreg2marsal, family=gaussian(), data=mydata3)
summary(glm10)
confint(glm10)

glm11=glm(GRDBSpreg2 ~ LGApreg2marsal, family=gaussian(), data=mydata3)
summary(glm11)
confint(glm11)

glm12=glm(GRDBSpreg2 ~ CSpreg1, family=gaussian(), data=mydata3)
summary(glm12)
confint(glm12)


#cleaning

rm(glm1,glm2,glm3,glm4,glm5,glm6,glm7,glm8,glm9,glm10,glm11,glm12)



#df2 multivar and results

mydata3=subset(mydata,spontpreg2)

tempPTL=subset(mydata3,opt|PTL)
tempPTL1=subset(mydata3,opt|PTL1)
tempPTL2=subset(mydata3,opt|PTL2)
tempPTL3=subset(mydata3,opt|PTL3)
tempPTL4=subset(mydata3,opt|PTL4)
temp37=subset(mydata3,opt|w37)
temp38=subset(mydata3,opt|w38)
temp41=subset(mydata3,opt|w41)
temp42=subset(mydata3,opt|w42)

temp1=subset(mydata3,PTL)
temp2=subset(mydata3,PTL1)
temp3=subset(mydata3,PTL2)
temp4=subset(mydata3,PTL3)
temp5=subset(mydata3,PTL4)
temp6=subset(mydata3,w37)
temp7=subset(mydata3,w38)
temp8=subset(mydata3,w41)
temp9=subset(mydata3,w42)



glm1=glm(PTL ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))

glm2=glm(PTL1 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal, family=binomial(), data=tempPTL1)
summary(glm2)
exp(cbind(OR=coef(glm2),confint(glm2)))

glm3=glm(PTL2 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=tempPTL2)
summary(glm3)
exp(cbind(OR=coef(glm3),confint(glm3)))

glm4=glm(PTL3 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=tempPTL3)
summary(glm4)
exp(cbind(OR=coef(glm4),confint(glm4)))

glm5=glm(PTL4 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=tempPTL4)
summary(glm5)
exp(cbind(OR=coef(glm5),confint(glm5)))

glm6=glm(w37 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=temp37)
summary(glm6)
exp(cbind(OR=coef(glm6),confint(glm6)))

glm7=glm(w38 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=temp38)
summary(glm7)
exp(cbind(OR=coef(glm7),confint(glm7)))

glm8=glm(w41 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=temp41)
summary(glm8)
exp(cbind(OR=coef(glm8),confint(glm8)))

glm9=glm(w42 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=temp42)
summary(glm9)
exp(cbind(OR=coef(glm9),confint(glm9)))




fonts()

#forest plot of above results (table 3a)


tabletext <- cbind(
    c("Gestational age","","<37+0 w","","<28+0 w","28+0 - 31+6 w","32+0 - 33+6 w","34+0 - 36+6 w","37+0 - 37+6 w","38+0 - 38+6 w","","41+0 - 41+6 w",">41+6 w"),
    c("n","","10 912","","408","631","959","8 914","16 348","53 534","","101 617","22 504"),
    c("% of cohort","","2.13%","","0.080%","0.12%","0.19%","1.74%","3.19%","10.45%","","19.83%","4.39%"),
    c("aOR (95% CI)","","1.67 (1.58, 1.77)","","1.77 (1.32, 2.34)","1.64 (1.28, 2.07)","2.16 (1.80, 2.57)","1.62 (1.52, 1.73)","1.26 (1.19, 1.32)","1.10 (1.06, 1.14)","","1.24 (1.21, 1.27)","1.55 (1.49, 1.62)"),
    c("P-value","","<1E-15","","7.2E-5","5.1E-5","<1E-15","<1E-15","<1E-15","2.2E-8","","<1E-15","<1E-15"))

df_c <- data.frame(mean = c(NA,NA,1.67,NA,1.77,1.64,2.16,1.62,1.26,1.10,NA,1.24,1.55),
                                      lower = c(NA,NA,1.58,NA,1.32,1.28,1.80,1.52,1.19,1.06,NA,1.21,1.49),
                                      upper = c(NA,NA,1.77,NA,2.34,2.07,2.57,1.73,1.32,1.14,NA,1.27,1.62))

forestplot(tabletext,
                      txt_gp = fpTxtGp(label = list(gpar(cex=1),
                                                    gpar(cex=1),
                                                    gpar(cex=1),
                                                    gpar(cex=1),
                                                    gpar(cex=1))),
                      df_c,new_page = TRUE,
                      boxsize = 0.2,
                      is.summary = c(TRUE, rep(FALSE, 12)),
                      clip = c(0.7,3.0),
                      xlab = "aOR",
                      xlog = FALSE,
                      zero=1,
                      colgap=unit(5, "mm"),
                      col = fpColors(box = "royalblue",
                                                    line = "darkblue"),
                      vertices = TRUE,
                      lineheight='lines')

rm(tabletext,df_c)


rm(glm1,glm2,glm3,glm4,glm5,glm6,glm7,glm8,glm9,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,tempPTL,tempPTL1,tempPTL2,tempPTL3,tempPTL4,temp37,temp38,temp41,temp42)




#df2 multivar and results only GA according to US

mydata3=subset(mydata,spontpreg2)
mydata3=subset(mydata3,GAUSpreg2)

tempPTL=subset(mydata3,opt|PTL)
tempPTL1=subset(mydata3,opt|PTL1)
tempPTL2=subset(mydata3,opt|PTL2)
tempPTL3=subset(mydata3,opt|PTL3)
tempPTL4=subset(mydata3,opt|PTL4)
temp37=subset(mydata3,opt|w37)
temp38=subset(mydata3,opt|w38)
temp41=subset(mydata3,opt|w41)
temp42=subset(mydata3,opt|w42)

temp1=subset(mydata3,PTL)
temp2=subset(mydata3,PTL1)
temp3=subset(mydata3,PTL2)
temp4=subset(mydata3,PTL3)
temp5=subset(mydata3,PTL4)
temp6=subset(mydata3,w37)
temp7=subset(mydata3,w38)
temp8=subset(mydata3,w41)
temp9=subset(mydata3,w42)



glm1=glm(PTL ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))

glm2=glm(PTL1 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal, family=binomial(), data=tempPTL1)
summary(glm2)
exp(cbind(OR=coef(glm2),confint(glm2)))

glm3=glm(PTL2 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=tempPTL2)
summary(glm3)
exp(cbind(OR=coef(glm3),confint(glm3)))

glm4=glm(PTL3 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=tempPTL3)
summary(glm4)
exp(cbind(OR=coef(glm4),confint(glm4)))

glm5=glm(PTL4 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=tempPTL4)
summary(glm5)
exp(cbind(OR=coef(glm5),confint(glm5)))

glm6=glm(w37 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=temp37)
summary(glm6)
exp(cbind(OR=coef(glm6),confint(glm6)))

glm7=glm(w38 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=temp38)
summary(glm7)
exp(cbind(OR=coef(glm7),confint(glm7)))

glm8=glm(w41 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=temp41)
summary(glm8)
exp(cbind(OR=coef(glm8),confint(glm8)))

glm9=glm(w42 ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, family=binomial(), data=temp42)
summary(glm9)
exp(cbind(OR=coef(glm9),confint(glm9)))




fonts()

#forest plot of above results (table 4a)


tabletext <- cbind(
    c("Gestational age","","<37+0 w","","<28+0 w","28+0 - 31+6 w","32+0 - 33+6 w","34+0 - 36+6 w","37+0 - 37+6 w","38+0 - 38+6 w","","41+0 - 41+6 w",">41+6 w"),
    c("n","","9 086","","295","490","778","7 523","14 188","47 031","","87 476","17 779"),
    c("% of cohort","","2.04%","","0.066%","0.11%","0.17%","1.69%","3.19%","10.56%","","19.64%","3.99%"),
    c("aOR (95% CI)","","1.67 (1.57, 1.78)","","1.57 (1.09, 2.19)","1.78 (1.35, 2.29)","2.24 (1.84, 2.71)","1.61 (1.50, 1.73)","1.26 (1.19, 1.33)","1.10 (1.06, 1.14)","","1.25 (1.22, 1.28)","1.61 (1.53, 1.68)"),
    c("P-value","","<1E-15","","0.011","1.9E-5","<1E-15","<1E-15","7.6E-15","6.7E-8","","<1E-15","<1E-15"))

df_c <- data.frame(mean = c(NA,NA,1.67,NA,1.57,1.78,2.24,1.61,1.26,1.10,NA,1.25,1.61),
                                      lower = c(NA,NA,1.57,NA,1.09,1.35,1.84,1.50,1.19,1.06,NA,1.22,1.53),
                                      upper = c(NA,NA,1.78,NA,2.19,2.29,2.71,1.73,1.33,1.14,NA,1.28,1.68))

forestplot(tabletext,
                      txt_gp = fpTxtGp(label = list(gpar(cex=1),
                                                    gpar(cex=1),
                                                    gpar(cex=1),
                                                    gpar(cex=1),
                                                    gpar(cex=1))),
                      df_c,new_page = TRUE,
                      boxsize = 0.2,
                      is.summary = c(TRUE, rep(FALSE, 12)),
                      clip = c(1,3.5),
                      xticks = c(1,1.5,2,2.5,3),
                      xlab = "aOR",
                      xlog = FALSE,
                      zero=1,
                      colgap=unit(5, "mm"),
                      col = fpColors(box = "royalblue",
                                                    line = "darkblue"),
                      vertices = TRUE,
                      lineheight='lines')

rm(tabletext,df_c)


rm(glm1,glm2,glm3,glm4,glm5,glm6,glm7,glm8,glm9,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,tempPTL,tempPTL1,tempPTL2,tempPTL3,tempPTL4,temp37,temp38,temp41,temp42)



#density plot

theme_set(theme_light())

mydata3=subset(mydata,spontpreg2)

g <- ggplot(mydata3, aes(GRDBSpreg2))

g + geom_density(aes(fill=factor(CSpreg1)), alpha=0.4, size=0.2) + 
    labs(x="Gestational duration of the second pregnancy (days)",
              y="Proportion",
              yaxt="n",
              fill="Mode of delivery the first pregnancy") +
    scale_fill_discrete(labels=c("Vaginal","Cesarean section")) +
    geom_vline(xintercept = 259, linetype=3, color="red") +
    annotate(geom="text", x=235, y=0.05, label="37+0 weeks", color="red") +
    theme(axis.text.y = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black"))

rm(g)


#qq-plot

mydata3=subset(mydata,spontpreg2)

temp1=subset(mydata3,CSpreg1)
qqnorm(temp1$GRDBSpreg2, pch=1, frame=FALSE)
qqline(temp1$GRDBSpreg2, lwd=1, col="red")

temp2=subset(mydata3,!CSpreg1)
qqnorm(temp2$GRDBSpreg2, pch=1, frame=FALSE, size=I(0.2))
qqline(temp2$GRDBSpreg2, lwd=1, col="red")

rm(temp1,temp2,mydata3)




#survival analyses


#Model 2


df2temp1=subset(mydata, spontpreg2)
df2surv1=data.frame(GRDBSpreg2=df2temp1$GRDBSpreg2, delivery=TRUE, CSpreg1=df2temp1$CSpreg1)

df2temp2=subset(mydata, !spontpreg2)
df2surv2=data.frame(GRDBSpreg2=df2temp2$GRDBSpreg2, delivery=FALSE, CSpreg1=df2temp2$CSpreg1)

df2surv=rbind(df2surv1,df2surv2)


survobject=Surv(time=df2surv$GRDBSpreg2, event=df2surv$delivery)
fit=survfit(survobject~CSpreg1, data=df2surv)

ggsurvplot(fit, data=df2surv)

mysurvplot=ggsurvplot(fit, size=0.5, xlim=c(154,308), break.x.by=28, ylab="", xlab="Gestational duration of the second pregnancy (days)", pval=TRUE, risk.table=TRUE, risk.table.title="", risk.table.height=0.30, legend.labs=c("Previous vaginal delivery","Previous CS"), legend.title="")

mysurvplot$plot <- mysurvplot$plot + 
    geom_segment(aes(x=259,y=0,xend=259,yend=1),linetype="dashed",color="red") +
    annotate(geom="text", x=245, y=0.75, label="37+0 weeks", color="red")

mysurvplot


rm(df2temp1,df2temp2,df2surv1,df2surv2,df2surv,survobject,fit,mysurvplot,ggsurvplot)




#Cox model 2

df2temp1=subset(mydata, spontpreg2)
df2surv1=data.frame(GRDBSpreg2=df2temp1$GRDBSpreg2, delivery=TRUE, CSpreg1=df2temp1$CSpreg1, MALDERpreg2=df2temp1$MALDERpreg2, BMIpreg2=df2temp1$BMIpreg2, notborninswepreg2=df2temp1$notborninswepreg2, smokepreg2=df2temp1$smokepreg2, diabetespreg2=df2temp1$diabetespreg2, allHTNpreg2=df2temp1$allHTNpreg2, boychildpreg2=df2temp1$boychildpreg2, malformpreg2=df2temp1$malformpreg2, SGApreg2marsal=df2temp1$SGApreg2marsal, LGApreg2marsal=df2temp1$LGApreg2marsal)

df2temp2=subset(mydata, !spontpreg2)
df2surv2=data.frame(GRDBSpreg2=df2temp2$GRDBSpreg2, delivery=FALSE, CSpreg1=df2temp2$CSpreg1, MALDERpreg2=df2temp2$MALDERpreg2, BMIpreg2=df2temp2$BMIpreg2, notborninswepreg2=df2temp2$notborninswepreg2, smokepreg2=df2temp2$smokepreg2, diabetespreg2=df2temp2$diabetespreg2, allHTNpreg2=df2temp2$allHTNpreg2, boychildpreg2=df2temp2$boychildpreg2, malformpreg2=df2temp2$malformpreg2, SGApreg2marsal=df2temp2$SGApreg2marsal, LGApreg2marsal=df2temp2$LGApreg2marsal)

df2surv=rbind(df2surv1,df2surv2)

cox = coxph(Surv(GRDBSpreg2, delivery) ~ CSpreg1 + MALDERpreg2 + BMIpreg2 + notborninswepreg2 + smokepreg2 + diabetespreg2 + allHTNpreg2 + boychildpreg2 + malformpreg2 + SGApreg2marsal + LGApreg2marsal, data = df2surv)

summary(cox)




rm(fit,df2temp1,df2temp2,df2surv1,df2surv2,df2surv,survobject,cox,mysurvplot)

rm(grav1,grav2,mydata,mydata3)





#df3


grav1=subset(mfr2, parity_clean=="1")
grav2=subset(mfr2, parity_clean=="2")
grav3=subset(mfr2, parity_clean=="3")

grav1=subset(grav1, GRDBS>="259")	
grav2=subset(grav2, GRDBS>="259")	

temp=subset(grav3, OFRIABEF=="1"|OFRISTIM=="1"|OFRIKIRU=="1"|OFRIICSI=="1"|OFRIANN=="1")
grav3=subset(grav3, !lpnr_mor %in% temp$lpnr_mor)
rm(temp)

grav1=subset(grav1, lpnr_mor!="NA")			
grav2=subset(grav2, lpnr_mor!="NA")
grav3=subset(grav3, lpnr_mor!="NA")

grav1=grav1[!duplicated(grav1$lpnr_mor),]
grav2=grav2[!duplicated(grav2$lpnr_mor),]
grav3=grav3[!duplicated(grav3$lpnr_mor),]

grav1=semi_join(grav1,grav2,by="lpnr_mor")
grav1=semi_join(grav1,grav3,by="lpnr_mor")
grav2=semi_join(grav2,grav1,by="lpnr_mor")
grav2=semi_join(grav2,grav3,by="lpnr_mor")
grav3=semi_join(grav3,grav1,by="lpnr_mor")
grav3=semi_join(grav3,grav2,by="lpnr_mor")

grav1=grav1[order(grav1$lpnr_mor),]
grav2=grav2[order(grav2$lpnr_mor),]
grav3=grav3[order(grav3$lpnr_mor),]

mydata=data.frame(lpnr_mor=grav1$lpnr_mor)

mydata=mydata %>% add_column(GRDBSpreg1=grav1$GRDBS, GRDBSpreg2=grav2$GRDBS, GRDBSpreg3=grav3$GRDBS)

temp1=subset(grav1,SECMARK=="1")
temp2=subset(grav2,TSECTIO=="1")
temp3=rbind(temp1,temp2)
temp4=temp3[!duplicated(temp3$lpnr_mor),]
mydata$CSpreg1 = mydata$lpnr_mor %in% temp4$lpnr_mor
rm(temp1,temp2,temp3,temp4)

temp=subset(grav2, SECMARK=="1")
mydata$CSpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav3, SECMARK=="1")
mydata$CSpreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav1, FLSPONT=="1")
mydata$spontpreg1= mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, FLSPONT=="1")
mydata$spontpreg2= mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav3, FLSPONT=="1")
mydata$spontpreg3= mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)


mydata$PTL = mydata$GRDBSpreg3<259
mydata$PTL1 = mydata$GRDBSpreg3<196
mydata$PTL2 = mydata$GRDBSpreg3>=196 & mydata$GRDBSpreg3<224
mydata$PTL3 = mydata$GRDBSpreg3>=224 & mydata$GRDBSpreg3<238
mydata$PTL4 = mydata$GRDBSpreg3>=238 & mydata$GRDBSpreg3<259
mydata$w37 = mydata$GRDBSpreg3>=259 & mydata$GRDBSpreg3<266
mydata$w38 = mydata$GRDBSpreg3>=266 & mydata$GRDBSpreg3<273
mydata$opt = mydata$GRDBSpreg3>=273 & mydata$GRDBSpreg3<287
mydata$w41 = mydata$GRDBSpreg3>=287 & mydata$GRDBSpreg3<294
mydata$w42 = mydata$GRDBSpreg3>=294




#categorization based on number of previous CSs

temp0=subset(mydata, !CSpreg1 & !CSpreg2)
temp2=subset(mydata, CSpreg1 & CSpreg2)

temp1.1=subset(mydata, CSpreg1 & !CSpreg2)
temp1.2=subset(mydata, !CSpreg1 & CSpreg2)
temp1.3=rbind(temp1.1,temp1.2)

mydata$zeroCS = mydata$lpnr_mor %in% temp0$lpnr_mor
mydata$twoCS = mydata$lpnr_mor %in% temp2$lpnr_mor
mydata$oneCS = mydata$lpnr_mor %in% temp1.3$lpnr_mor

rm(temp0,temp2,temp1.1,temp1.2,temp1.3)



#categorization based on sequence of previous delivery modes

vagvag=subset(mydata, !CSpreg1 & !CSpreg2)
vagCS=subset(mydata, !CSpreg1 & CSpreg2)
CSvag=subset(mydata, CSpreg1 & !CSpreg2)
CSCS=subset(mydata, CSpreg1 & CSpreg2)

mydata$vagvag = mydata$lpnr_mor %in% vagvag$lpnr_mor
mydata$vagCS = mydata$lpnr_mor %in% vagCS$lpnr_mor
mydata$CSvag = mydata$lpnr_mor %in% CSvag$lpnr_mor
mydata$CSCS = mydata$lpnr_mor %in% CSCS$lpnr_mor

rm(vagvag,vagCS,CSvag,CSCS)



#below is for adjusting for confounders

#year of birth, maternal age, dates of birth, birth weight

mydata=mydata %>% add_column(ARpreg1=grav1$AR, ARpreg2=grav2$AR, ARpreg3=grav3$AR, MALDERpreg1=grav1$MALDER, MALDERpreg2=grav2$MALDER, MALDERpreg3=grav3$MALDER, BFODDATpreg1=grav1$BFODDAT, BFODDATpreg2=grav2$BFODDAT, BFODDATpreg3=grav3$BFODDAT, BVIKTpreg1=grav1$BVIKT, BVIKTpreg2=grav2$BVIKT, BVIKTpreg3=grav3$BVIKT)



#BMI

mydata=mydata %>% add_column(BMIpreg1=(grav1$MVIKT/((grav1$MLANGD/100)^2)), BMIpreg2=(grav2$MVIKT/((grav2$MLANGD/100)^2)), BMIpreg3=(grav3$MVIKT/((grav3$MLANGD/100)^2)))

mydata$BMIpreg1.2 = ifelse(is.na(mydata$BMIpreg1), grav1$MVIKT/((grav2$MLANGD/100)^2), mydata$BMIpreg1)
mydata$BMIpreg2.2 = ifelse(is.na(mydata$BMIpreg2), grav2$MVIKT/((grav1$MLANGD/100)^2), mydata$BMIpreg2)
mydata$BMIpreg3.2 = ifelse(is.na(mydata$BMIpreg3), grav3$MVIKT/((grav2$MLANGD/100)^2), mydata$BMIpreg3)
mydata=subset(mydata, select=-c(BMIpreg1,BMIpreg2,BMIpreg3))
mydata=mydata %>% rename(BMIpreg1=BMIpreg1.2)
mydata=mydata %>% rename(BMIpreg2=BMIpreg2.2)
mydata=mydata %>% rename(BMIpreg3=BMIpreg3.2)
mydata$BMIpreg1=formatC(mydata$BMIpreg1, digits=1, format="f")
mydata$BMIpreg2=formatC(mydata$BMIpreg2, digits=1, format="f")
mydata$BMIpreg3=formatC(mydata$BMIpreg3, digits=1, format="f")

mydata$BMIpreg1[mydata$BMIpreg1=="NA"]=23.41
mydata$BMIpreg2[mydata$BMIpreg2=="NA"]=24.10
mydata$BMIpreg3[mydata$BMIpreg3=="NA"]=24.98
mydata$BMIpreg1=as.numeric(mydata$BMIpreg1)
mydata$BMIpreg2=as.numeric(mydata$BMIpreg2)
mydata$BMIpreg3=as.numeric(mydata$BMIpreg3)



#BMI groups

mydata$BMIpreg1.1 = mydata$BMIpreg1<18.5
mydata$BMIpreg1.2 = mydata$BMIpreg1>=18.5 & mydata$BMIpreg1<25.0
mydata$BMIpreg1.3 = mydata$BMIpreg1>=25.0 & mydata$BMIpreg1<30.0
mydata$BMIpreg1.4 = mydata$BMIpreg1>=30.0 & mydata$BMIpreg1<35.0
mydata$BMIpreg1.5 = mydata$BMIpreg1>=35.0 & mydata$BMIpreg1<40.0
mydata$BMIpreg1.6 = mydata$BMIpreg1>=40.0 & mydata$BMIpreg1<99

mydata$BMIpreg2.1 = mydata$BMIpreg2<18.5
mydata$BMIpreg2.2 = mydata$BMIpreg2>=18.5 & mydata$BMIpreg2<25.0
mydata$BMIpreg2.3 = mydata$BMIpreg2>=25.0 & mydata$BMIpreg2<30.0
mydata$BMIpreg2.4 = mydata$BMIpreg2>=30.0 & mydata$BMIpreg2<35.0
mydata$BMIpreg2.5 = mydata$BMIpreg2>=35.0 & mydata$BMIpreg2<40.0
mydata$BMIpreg2.6 = mydata$BMIpreg2>=40.0 & mydata$BMIpreg2<99

mydata$BMIpreg3.1 = mydata$BMIpreg3<18.5
mydata$BMIpreg3.2 = mydata$BMIpreg3>=18.5 & mydata$BMIpreg3<25.0
mydata$BMIpreg3.3 = mydata$BMIpreg3>=25.0 & mydata$BMIpreg3<30.0
mydata$BMIpreg3.4 = mydata$BMIpreg3>=30.0 & mydata$BMIpreg3<35.0
mydata$BMIpreg3.5 = mydata$BMIpreg3>=35.0 & mydata$BMIpreg3<40.0
mydata$BMIpreg3.6 = mydata$BMIpreg3>=40.0 & mydata$BMIpreg3<99


#not born in sweden, unemployed, smoking

temp=subset(grav1, MFODLAND!="SVERIGE")
mydata$notborninswepreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, MFODLAND!="SVERIGE")
mydata$notborninswepreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav3, MFODLAND!="SVERIGE")
mydata$notborninswepreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav1, ARBETE=="3")
mydata$unemployedpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, ARBETE=="3")
mydata$unemployedpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav3, ARBETE=="3")
mydata$unemployedpreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav1, ROK1!="1")
mydata$smokepreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, ROK1!="1")
mydata$smokepreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav3, ROK1!="1")
mydata$smokepreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)


#diabetes, hypertension, gestational hypertension, preeclampsia

temp1=subset(grav1, DIABETES=="1"|DIABETES=="2")
temp2=grav1 %>% filter_all(any_vars(. %in% c('E107','O240','O241','O243','O244','O249','O240B','O240C','O240D','O240E','O240F','O240X','O244A','O244B','O244X','25000','25009','76110','250A','250B','250C','250D','250E','250F','250X','648A')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$diabetespreg1 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)

temp1=subset(grav2, DIABETES=="1"|DIABETES=="2")
temp2=grav2 %>% filter_all(any_vars(. %in% c('E107','O240','O241','O243','O244','O249','O240B','O240C','O240D','O240E','O240F','O240X','O244A','O244B','O244X','25000','25009','76110','250A','250B','250C','250D','250E','250F','250X','648A')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$diabetespreg2 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)

temp1=subset(grav3, DIABETES=="1"|DIABETES=="2")
temp2=grav3 %>% filter_all(any_vars(. %in% c('E107','O240','O241','O243','O244','O249','O240B','O240C','O240D','O240E','O240F','O240X','O244A','O244B','O244X','25000','25009','76110','250A','250B','250C','250D','250E','250F','250X','648A')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$diabetespreg3 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)


temp1=subset(grav1, HYPERTON=="1"|HYPERTON=="2")
temp2=grav1 %>% filter_all(any_vars(. %in% c('I10','I109','401','401X','40199')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$HTNpreg1 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)

temp1=subset(grav2, HYPERTON=="1"|HYPERTON=="2")
temp2=grav2 %>% filter_all(any_vars(. %in% c('I10','I109','401','401X','40199')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$HTNpreg2 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)

temp1=subset(grav3, HYPERTON=="1"|HYPERTON=="2")
temp2=grav3 %>% filter_all(any_vars(. %in% c('I10','I109','401','401X','40199')))
temp3=rbind(temp1,temp2)
temp3=temp3[!duplicated(temp3$lpnr_mor),]
mydata$HTNpreg3 = mydata$lpnr_mor %in% temp3$lpnr_mor
rm(temp1,temp2,temp3)



temp=grav1 %>% filter_all(any_vars(. %in% c('O139','642','642A','642B','642C','642D','642X','63701')))
mydata$GHTNpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=grav2 %>% filter_all(any_vars(. %in% c('O139','642','642A','642B','642C','642D','642X','63701')))
mydata$GHTNpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=grav3 %>% filter_all(any_vars(. %in% c('O139','642','642A','642B','642C','642D','642X','63701')))
mydata$GHTNpreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)


temp=grav1 %>% filter_all(any_vars(. %in% c('O14','O140','O141','O142','O149','O141A','O141B','O141X','642E','642F','642G','642H','63703','63704','63709','63710','63799','76210','76220','76230')))
mydata$PEpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=grav2 %>% filter_all(any_vars(. %in% c('O14','O140','O141','O142','O149','O141A','O141B','O141X','642E','642F','642G','642H','63703','63704','63709','63710','63799','76210','76220','76230')))
mydata$PEpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=grav3 %>% filter_all(any_vars(. %in% c('O14','O140','O141','O142','O149','O141A','O141B','O141X','642E','642F','642G','642H','63703','63704','63709','63710','63799','76210','76220','76230')))
mydata$PEpreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)


mydata$allHTNpreg1 = mydata$HTNpreg1 | mydata$GHTNpreg1 | mydata$PEpreg1

mydata$allHTNpreg2 = mydata$HTNpreg2 | mydata$GHTNpreg2 | mydata$PEpreg2

mydata$allHTNpreg3 = mydata$HTNpreg3 | mydata$GHTNpreg3 | mydata$PEpreg3



#boy child, fetal malformation

temp=subset(grav1, KON=="1")
mydata$boychildpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, KON=="1")
mydata$boychildpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav3, KON=="1")
mydata$boychildpreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav1, MISSB=="1")
mydata$malformpreg1 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav2, MISSB=="1")
mydata$malformpreg2 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)

temp=subset(grav3, MISSB=="1")
mydata$malformpreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)



#SGA and LGA according to Marsal


marsalpreg1 = function(GRDBSpreg1,BVIKTpreg1,boychildpreg1){
    if(boychildpreg1){
        mw=-1.907345e-6*GRDBSpreg1^4 + 1.140644e-3*GRDBSpreg1^3 - 1.336265e-1*GRDBSpreg1^2 + 1.976961*GRDBSpreg1 + 2.410053e+2
    } else {
        mw=-2.761948e-6*GRDBSpreg1^4 + 1.744841e-3*GRDBSpreg1^3 - 2.893626e-1*GRDBSpreg1^2 + 18.91197*GRDBSpreg1 - 4.135122e+2
    }
    pnorm(BVIKTpreg1, mw, abs(0.12*mw))*100
}

mydata$PCTmarsalpreg1 = sapply(1:nrow(mydata), function(x) marsalpreg1(mydata$GRDBSpreg1[x],mydata$BVIKTpreg1[x],mydata$boychildpreg1[x]))

mydata$SGApreg1marsal=mydata$PCTmarsalpreg1 < pnorm(-2) *100
mydata$LGApreg1marsal=mydata$PCTmarsalpreg1 > pnorm(2) *100


mydata$PCTmarsalpreg2 = sapply(1:nrow(mydata), function(x) marsalpreg1(mydata$GRDBSpreg2[x],mydata$BVIKTpreg2[x],mydata$boychildpreg2[x]))

mydata$SGApreg2marsal=mydata$PCTmarsalpreg2 < pnorm(-2) *100
mydata$LGApreg2marsal=mydata$PCTmarsalpreg2 > pnorm(2) *100

mydata$PCTmarsalpreg3 = sapply(1:nrow(mydata), function(x) marsalpreg1(mydata$GRDBSpreg3[x],mydata$BVIKTpreg3[x],mydata$boychildpreg3[x]))

mydata$SGApreg3marsal=mydata$PCTmarsalpreg3 < pnorm(-2) *100
mydata$LGApreg3marsal=mydata$PCTmarsalpreg3 > pnorm(2) *100


rm(marsalpreg1,marsalpreg2,marsalpreg3)

mydata$seq = ifelse(mydata$vagvag,1,ifelse(mydata$vagCS,2,ifelse(mydata$CSvag,3,4)))


#GA according to ultrasound

temp=subset(grav3, GRMETOD==1|GRMETOD==5|GRMETOD==6|GRMETOD==7)
mydata$GAUSpreg3 = mydata$lpnr_mor %in% temp$lpnr_mor
rm(temp)



#df3 descriptive

#table 2


vagvag=subset(mydata,vagvag)
vagCS=subset(mydata,vagCS)
CSvag=subset(mydata,CSvag)
CSCS=subset(mydata,CSCS)


#continuous variables



#gestational age preg 1

mean(vagvag$GRDBSpreg1)
sd(vagvag$GRDBSpreg1)
mean(vagCS$GRDBSpreg1)
sd(vagCS$GRDBSpreg1)
mean(CSvag$GRDBSpreg1)
sd(CSvag$GRDBSpreg1)
mean(CSCS$GRDBSpreg1)
sd(CSCS$GRDBSpreg1)

#gestational age preg 2

mean(vagvag$GRDBSpreg2)
sd(vagvag$GRDBSpreg2)
mean(vagCS$GRDBSpreg2)
sd(vagCS$GRDBSpreg2)
mean(CSvag$GRDBSpreg2)
sd(CSvag$GRDBSpreg2)
mean(CSCS$GRDBSpreg2)
sd(CSCS$GRDBSpreg2)

#gestational age preg 3

mean(vagvag$GRDBSpreg3)
sd(vagvag$GRDBSpreg3)
mean(vagCS$GRDBSpreg3)
sd(vagCS$GRDBSpreg3)
mean(CSvag$GRDBSpreg3)
sd(CSvag$GRDBSpreg3)
mean(CSCS$GRDBSpreg3)
sd(CSCS$GRDBSpreg3)



#categorical variables


#Preterm birth

sum(vagvag$PTL)/nrow(vagvag)
sum(vagCS$PTL)/nrow(vagCS)
sum(CSvag$PTL)/nrow(CSvag)
sum(CSCS$PTL)/nrow(CSCS)


#CSpreg3

sum(vagvag$CSpreg3)/nrow(vagvag)
sum(vagCS$CSpreg3)/nrow(vagCS)
sum(CSvag$CSpreg3)/nrow(CSvag)
sum(CSCS$CSpreg3)/nrow(CSCS)


#spontpreg3

sum(vagvag$spontpreg3)/nrow(vagvag)
sum(vagCS$spontpreg3)/nrow(vagCS)
sum(CSvag$spontpreg3)/nrow(CSvag)
sum(CSCS$spontpreg3)/nrow(CSCS)

#GA according to US preg3

sum(vagvag$GAUSpreg3)/nrow(vagvag)
sum(vagCS$GAUSpreg3)/nrow(vagCS)
sum(CSvag$GAUSpreg3)/nrow(CSvag)
sum(CSCS$GAUSpreg3)/nrow(CSCS)


#spontaneous onset preg3 only

temp1=subset(mydata,spontpreg3)
vagvag=subset(temp1,vagvag)
vagCS=subset(temp1,vagCS)
CSvag=subset(temp1,CSvag)
CSCS=subset(temp1,CSCS)


#gestational age preg 3

mean(vagvag$GRDBSpreg3)
sd(vagvag$GRDBSpreg3)
mean(vagCS$GRDBSpreg3)
sd(vagCS$GRDBSpreg3)
mean(CSvag$GRDBSpreg3)
sd(CSvag$GRDBSpreg3)
mean(CSCS$GRDBSpreg3)
sd(CSCS$GRDBSpreg3)

#Preterm birth

sum(vagvag$PTL)/nrow(vagvag)
sum(vagCS$PTL)/nrow(vagCS)
sum(CSvag$PTL)/nrow(CSvag)
sum(CSCS$PTL)/nrow(CSCS)

#CSpreg3

sum(vagvag$CSpreg3)/nrow(vagvag)
sum(vagCS$CSpreg3)/nrow(vagCS)
sum(CSvag$CSpreg3)/nrow(CSvag)
sum(CSCS$CSpreg3)/nrow(CSCS)

#GA acc US preg3

sum(vagvag$GAUSpreg3)/nrow(vagvag)
sum(vagCS$GAUSpreg3)/nrow(vagCS)
sum(CSvag$GAUSpreg3)/nrow(CSvag)
sum(CSCS$GAUSpreg3)/nrow(CSCS)



#cleaning


rm(temp1,vagvag,vagCS,CSvag,CSCS)




#multivariable models


mydata2=subset(mydata,spontpreg3)


#vag-CS vs vag-vag

mydata3=subset(mydata2,!CSCS)
mydata3=subset(mydata3, vagCS | vagvag)
tempPTL=subset(mydata3,opt | PTL)

sum(tempPTL$vagCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ vagCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#vag-CS vs CS-vag

mydata3=subset(mydata2,!CSCS)
mydata3=subset(mydata3, vagCS | CSvag)
tempPTL=subset(mydata3,opt | PTL)

sum(tempPTL$vagCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ vagCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#CS-vag vs vag-vag

mydata3=subset(mydata2, !CSCS)
mydata3=subset(mydata3, CSvag | vagvag)
tempPTL=subset(mydata3, opt | PTL)

sum(tempPTL$CSvag)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSvag + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)




#CS-vag vs vag-CS

mydata3=subset(mydata2, !CSCS)
mydata3=subset(mydata3, vagCS | CSvag)
tempPTL=subset(mydata3,opt|PTL)

sum(tempPTL$CSvag)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSvag + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)


#CS-CS vs vag-vag

mydata3=subset(mydata2, CSCS | vagvag)
tempPTL=subset(mydata3,opt | PTL)

sum(tempPTL$CSCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#CS-CS vs CS-vag

mydata3=subset(mydata2, CSCS | CSvag)
tempPTL=subset(mydata3,opt | PTL)

sum(tempPTL$CSCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#CS-CS vs vag-CS

mydata3=subset(mydata2, CSCS | vagCS)
tempPTL=subset(mydata3, opt | PTL)

sum(tempPTL$CSCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#forest plot of above results (table 3b)


tabletext <- cbind(
    c("Sequence","","CS - vaginal","vaginal - CS","vaginal - CS"),
    c("","","vs","vs","vs"),
    c("","","vaginal - vaginal","vaginal - vaginal","CS - vaginal"),
    c("aOR (95% CI)","","1.15 (0.96, 1.36)","2.51 (2.01, 3.09)","2.25 (1.71, 2.95)"),
    c("P-value","","0.13","<1E-15","5.5E-9"))

df_c <- data.frame(mean = c(NA,NA,1.15,2.51,2.25),
                                      lower = c(NA,NA,0.96,2.01,1.71),
                                      upper = c(NA,NA,1.36,3.09,2.95))

forestplot(tabletext,
                      txt_gp = fpTxtGp(label = list(gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1))),
                      df_c,new_page = TRUE,
                      boxsize = 0.3,
                      is.summary = c(TRUE, rep(FALSE, 5)),
                      clip = c(0.5,3.5),
                      zero=0.5,
                      xlab = "aOR",
                      xticks = c(0.5,1,1.5,2,2.5,3, 3.5),
                      colgap=unit(5, "mm"),
                      col = fpColors(box = "royalblue",
                                                    line = "darkblue"),
                      vertices = TRUE,
                      lineheight='lines')

rm(tabletext,df_c)



tabletext <- cbind(
    c("Sequence","","CS - CS","CS - CS","CS - CS"),
    c("","","vs","vs","vs"),
    c("","","vaginal - vaginal","CS - vaginal","vaginal - CS"),
    c("aOR (95% CI)","","28.50 (21.73, 37.34)","25.50 (18.46, 35.33)","11.09 (7.85, 15.73)"),
    c("P-value","","<1E-15","<1E-15","<1E-15"))

df_c <- data.frame(mean = c(NA,NA,28.50,25.50,11.09),
                                      lower = c(NA,NA,21.73,18.46,7.85),
                                      upper = c(NA,NA,37.34,35.33,15.73))

forestplot(tabletext,
                      txt_gp = fpTxtGp(label = list(gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1))),
                      df_c,new_page = TRUE,
                      boxsize = 0.3,
                      is.summary = c(TRUE, rep(FALSE, 5)),
                      clip = c(1,40),
                      xlab = "aOR",
                      xlog = FALSE,
                      xticks = c(1,10,20,30,40),
                      zero=1,
                      colgap=unit(5, "mm"),
                      col = fpColors(box = "royalblue",
                                                    line = "darkblue"),
                      vertices = TRUE,
                      lineheight='lines')

rm(tabletext,df_c)



#multivariable models of only GA according to ultrasound


mydata2=subset(mydata,spontpreg3)
mydata2=subset(mydata2,GAUSpreg3)


#vag-CS vs vag-vag

mydata3=subset(mydata2,!CSCS)
mydata3=subset(mydata3, vagCS | vagvag)
tempPTL=subset(mydata3,opt|PTL)

sum(tempPTL$vagCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ vagCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#vag-CS vs CS-vag

mydata3=subset(mydata2,!CSCS)
mydata3=subset(mydata3, vagCS | CSvag)
tempPTL=subset(mydata3,opt|PTL)

sum(tempPTL$vagCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ vagCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#CS-vag vs vag-vag

mydata3=subset(mydata2,!CSCS)
mydata3=subset(mydata3, CSvag | vagvag)
tempPTL=subset(mydata3,opt|PTL)

sum(tempPTL$CSvag)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSvag + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#CS-vag vs vag-CS

mydata3=subset(mydata2,!CSCS)
mydata3=subset(mydata3, vagCS | CSvag)
tempPTL=subset(mydata3,opt|PTL)

sum(tempPTL$CSvag)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSvag + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)


#CS-CS vs vag-vag

mydata3=subset(mydata2, CSCS | vagvag)
tempPTL=subset(mydata3,opt | PTL)

sum(tempPTL$CSCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#CS-CS vs CS-vag

mydata3=subset(mydata2, CSCS | CSvag)
tempPTL=subset(mydata3,opt|PTL)

sum(tempPTL$CSCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#CS-CS vs vag-CS

mydata3=subset(mydata2, CSCS | vagCS)
tempPTL=subset(mydata3,opt|PTL)

sum(tempPTL$CSCS)
sum(tempPTL$notborninswepreg3)
sum(tempPTL$unemployedpreg3)
sum(tempPTL$smokepreg3)
sum(tempPTL$diabetespreg3)
sum(tempPTL$allHTNpreg3)
sum(tempPTL$boychildpreg3)
sum(tempPTL$malformpreg3)
sum(tempPTL$SGApreg3marsal)
sum(tempPTL$LGApreg3marsal)

#all >10

glm1=glm(PTL ~ CSCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, family=binomial(), data=tempPTL)
summary(glm1)
exp(cbind(OR=coef(glm1),confint(glm1)))
rm(mydata3,tempPTL,glm1)



#forest plot of above results (table 4b)


tabletext <- cbind(
    c("Sequence","","CS - vaginal","vaginal - CS","vaginal - CS"),
    c("","","vs","vs","vs"),
    c("","","vaginal - vaginal","vaginal - vaginal","CS - vaginal"),
    c("aOR (95% CI)","","1.12 (0.91, 1.35)","2.27 (1.77, 2.87)","2.10 (1.54, 2.84)"),
    c("P-value","","0.28","3.1E-11","2.1E-6"))

df_c <- data.frame(mean = c(NA,NA,1.12,2.27,2.10),
                                      lower = c(NA,NA,0.91,1.77,1.54),
                                      upper = c(NA,NA,1.35,2.87,2.84))

forestplot(tabletext,
                      txt_gp = fpTxtGp(label = list(gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1))),
                      df_c,new_page = TRUE,
                      boxsize = 0.3,
                      is.summary = c(TRUE, rep(FALSE, 5)),
                      clip = c(0.5,3.0),
                      zero=0.5,
                      xlab = "aOR",
                      colgap=unit(5, "mm"),
                      col = fpColors(box = "royalblue",
                                                    line = "darkblue"),
                      vertices = TRUE,
                      lineheight='lines')

rm(tabletext,df_c)


tabletext <- cbind(
    c("Sequence","","CS - CS","CS - CS","CS - CS"),
    c("","","vs","vs","vs"),
    c("","","vaginal - vaginal","CS - vaginal","vaginal - CS"),
    c("aOR (95% CI)","","29.03 (21.63, 38.90)","26.45 (18.60, 37.73)","12.61 (8.61, 18.59)"),
    c("P-value","","<1E-15","<1E-15","<1E-15"))

df_c <- data.frame(mean = c(NA,NA,29.03,26.45,12.61),
                                      lower = c(NA,NA,21.63,18.60,8.61),
                                      upper = c(NA,NA,38.90,37.73,18.59))

forestplot(tabletext,
                      txt_gp = fpTxtGp(label = list(gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1),
                                                                                  gpar(cex=1))),
                      df_c,new_page = TRUE,
                      boxsize = 0.3,
                      is.summary = c(TRUE, rep(FALSE, 5)),
                      clip = c(1,40),
                      zero=0.5,
                      xticks = c(1,10,20,30,40),
                      xlab = "aOR",
                      colgap=unit(5, "mm"),
                      col = fpColors(box = "royalblue",
                                                    line = "darkblue"),
                      vertices = TRUE,
                      lineheight='lines')

rm(tabletext,df_c)




#density plots


#only spontaneous onsets

mydata3=subset(mydata, spontpreg3)

theme_set(theme_light())

g <- ggplot(mydata3, aes(GRDBSpreg3))

g + geom_density(aes(fill=factor(seq)), alpha=0.4, size=0.2) + 
    labs(x="Gestational duration of the third pregnancy (days)",
              y="Proportion",
              yaxt="n",
              fill="Mode of delivery previous pregnancies") +
    scale_fill_discrete(labels=c("Vaginal - vaginal","Vaginal - CS","CS - vaginal", "CS - CS")) +
    geom_vline(xintercept = 259, linetype=3, color="red") +
    annotate(geom="text", x=235, y=0.05, label="37+0 weeks", color="red") +
    theme(axis.text.y = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black"))

rm(g)




#survival analyses


#Model 2


df3temp1=subset(mydata, spontpreg3)
df3surv1=data.frame(GRDBSpreg3=df3temp1$GRDBSpreg3, delivery=TRUE, seq=df3temp1$seq)

df3temp2=subset(mydata, !spontpreg3)
df3surv2=data.frame(GRDBSpreg3=df3temp2$GRDBSpreg3, delivery=FALSE, seq=df3temp2$seq)

df3surv=rbind(df3surv1,df3surv2)


survobject=Surv(time=df3surv$GRDBSpreg3, event=df3surv$delivery)
fit=survfit(survobject~seq, data=df3surv)

ggsurvplot(fit, data=df3surv)

mysurvplot=ggsurvplot(fit, size=0.5, xlim=c(154,308), break.x.by=28, ylab="", xlab="Gestational duration of the third pregnancy (days)", pval=TRUE, risk.table=TRUE, risk.table.title="", risk.table.height=0.30, legend.labs=c("Vaginal - vaginal","Vaginal - CS","CS - vaginal", "CS - CS"), legend.title="Previous modes of delivery")

mysurvplot$plot <- mysurvplot$plot + 
    geom_segment(aes(x=259,y=0,xend=259,yend=1),linetype="dashed",color="red") +
    annotate(geom="text", x=245, y=0.75, label="37+0 weeks", color="red")

mysurvplot



#Cox model 2 on variable seq

df3temp1=subset(mydata, spontpreg3)
df3surv1=data.frame(GRDBSpreg3=df3temp1$GRDBSpreg3, delivery=TRUE, seq=df3temp1$seq, MALDERpreg3=df3temp1$MALDERpreg3, BMIpreg3=df3temp1$BMIpreg3, notborninswepreg3=df3temp1$notborninswepreg3, smokepreg3=df3temp1$smokepreg3, diabetespreg3=df3temp1$diabetespreg3, allHTNpreg3=df3temp1$allHTNpreg3, boychildpreg3=df3temp1$boychildpreg3, malformpreg3=df3temp1$malformpreg3, SGApreg3marsal=df3temp1$SGApreg3marsal, LGApreg3marsal=df3temp1$LGApreg3marsal)

df3temp2=subset(mydata, !spontpreg3)
df3surv2=data.frame(GRDBSpreg3=df3temp2$GRDBSpreg3, delivery=FALSE, seq=df3temp2$seq, MALDERpreg3=df3temp2$MALDERpreg3, BMIpreg3=df3temp2$BMIpreg3, notborninswepreg3=df3temp2$notborninswepreg3, smokepreg3=df3temp2$smokepreg3, diabetespreg3=df3temp2$diabetespreg3, allHTNpreg3=df3temp2$allHTNpreg3, boychildpreg3=df3temp2$boychildpreg3, malformpreg3=df3temp2$malformpreg3, SGApreg3marsal=df3temp2$SGApreg3marsal, LGApreg3marsal=df3temp2$LGApreg3marsal)

df3surv=rbind(df3surv1,df3surv2)

cox = coxph(Surv(GRDBSpreg3, delivery) ~ seq + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, data = df3surv)

summary(cox)



#Cox model 2 on variable vagvag

df3temp1=subset(mydata, spontpreg3)
df3surv1=data.frame(GRDBSpreg3=df3temp1$GRDBSpreg3, delivery=TRUE, vagvag=df3temp1$vagvag, MALDERpreg3=df3temp1$MALDERpreg3, BMIpreg3=df3temp1$BMIpreg3, notborninswepreg3=df3temp1$notborninswepreg3, smokepreg3=df3temp1$smokepreg3, diabetespreg3=df3temp1$diabetespreg3, allHTNpreg3=df3temp1$allHTNpreg3, boychildpreg3=df3temp1$boychildpreg3, malformpreg3=df3temp1$malformpreg3, SGApreg3marsal=df3temp1$SGApreg3marsal, LGApreg3marsal=df3temp1$LGApreg3marsal)

df3temp2=subset(mydata, !spontpreg3)
df3surv2=data.frame(GRDBSpreg3=df3temp2$GRDBSpreg3, delivery=FALSE, vagvag=df3temp2$vagvag, MALDERpreg3=df3temp2$MALDERpreg3, BMIpreg3=df3temp2$BMIpreg3, notborninswepreg3=df3temp2$notborninswepreg3, smokepreg3=df3temp2$smokepreg3, diabetespreg3=df3temp2$diabetespreg3, allHTNpreg3=df3temp2$allHTNpreg3, boychildpreg3=df3temp2$boychildpreg3, malformpreg3=df3temp2$malformpreg3, SGApreg3marsal=df3temp2$SGApreg3marsal, LGApreg3marsal=df3temp2$LGApreg3marsal)

df3surv=rbind(df3surv1,df3surv2)

cox = coxph(Surv(GRDBSpreg3, delivery) ~ vagvag + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, data = df3surv)

summary(cox)



#Cox model 2 on variable vagCS

df3temp1=subset(mydata, spontpreg3)
df3surv1=data.frame(GRDBSpreg3=df3temp1$GRDBSpreg3, delivery=TRUE, vagCS=df3temp1$vagCS, MALDERpreg3=df3temp1$MALDERpreg3, BMIpreg3=df3temp1$BMIpreg3, notborninswepreg3=df3temp1$notborninswepreg3, smokepreg3=df3temp1$smokepreg3, diabetespreg3=df3temp1$diabetespreg3, allHTNpreg3=df3temp1$allHTNpreg3, boychildpreg3=df3temp1$boychildpreg3, malformpreg3=df3temp1$malformpreg3, SGApreg3marsal=df3temp1$SGApreg3marsal, LGApreg3marsal=df3temp1$LGApreg3marsal)

df3temp2=subset(mydata, !spontpreg3)
df3surv2=data.frame(GRDBSpreg3=df3temp2$GRDBSpreg3, delivery=FALSE, vagCS=df3temp2$vagCS, MALDERpreg3=df3temp2$MALDERpreg3, BMIpreg3=df3temp2$BMIpreg3, notborninswepreg3=df3temp2$notborninswepreg3, smokepreg3=df3temp2$smokepreg3, diabetespreg3=df3temp2$diabetespreg3, allHTNpreg3=df3temp2$allHTNpreg3, boychildpreg3=df3temp2$boychildpreg3, malformpreg3=df3temp2$malformpreg3, SGApreg3marsal=df3temp2$SGApreg3marsal, LGApreg3marsal=df3temp2$LGApreg3marsal)

df3surv=rbind(df3surv1,df3surv2)

cox = coxph(Surv(GRDBSpreg3, delivery) ~ vagCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, data = df3surv)

summary(cox)



#Cox model 2 on variable CSvag

df3temp1=subset(mydata, spontpreg3)
df3surv1=data.frame(GRDBSpreg3=df3temp1$GRDBSpreg3, delivery=TRUE, CSvag=df3temp1$CSvag, MALDERpreg3=df3temp1$MALDERpreg3, BMIpreg3=df3temp1$BMIpreg3, notborninswepreg3=df3temp1$notborninswepreg3, smokepreg3=df3temp1$smokepreg3, diabetespreg3=df3temp1$diabetespreg3, allHTNpreg3=df3temp1$allHTNpreg3, boychildpreg3=df3temp1$boychildpreg3, malformpreg3=df3temp1$malformpreg3, SGApreg3marsal=df3temp1$SGApreg3marsal, LGApreg3marsal=df3temp1$LGApreg3marsal)

df3temp2=subset(mydata, !spontpreg3)
df3surv2=data.frame(GRDBSpreg3=df3temp2$GRDBSpreg3, delivery=FALSE, CSvag=df3temp2$CSvag, MALDERpreg3=df3temp2$MALDERpreg3, BMIpreg3=df3temp2$BMIpreg3, notborninswepreg3=df3temp2$notborninswepreg3, smokepreg3=df3temp2$smokepreg3, diabetespreg3=df3temp2$diabetespreg3, allHTNpreg3=df3temp2$allHTNpreg3, boychildpreg3=df3temp2$boychildpreg3, malformpreg3=df3temp2$malformpreg3, SGApreg3marsal=df3temp2$SGApreg3marsal, LGApreg3marsal=df3temp2$LGApreg3marsal)

df3surv=rbind(df3surv1,df3surv2)

cox = coxph(Surv(GRDBSpreg3, delivery) ~ CSvag + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, data = df3surv)

summary(cox)



#Cox model 2 on variable CSCS

df3temp1=subset(mydata, spontpreg3)
df3surv1=data.frame(GRDBSpreg3=df3temp1$GRDBSpreg3, delivery=TRUE, CSCS=df3temp1$CSCS, MALDERpreg3=df3temp1$MALDERpreg3, BMIpreg3=df3temp1$BMIpreg3, notborninswepreg3=df3temp1$notborninswepreg3, smokepreg3=df3temp1$smokepreg3, diabetespreg3=df3temp1$diabetespreg3, allHTNpreg3=df3temp1$allHTNpreg3, boychildpreg3=df3temp1$boychildpreg3, malformpreg3=df3temp1$malformpreg3, SGApreg3marsal=df3temp1$SGApreg3marsal, LGApreg3marsal=df3temp1$LGApreg3marsal)

df3temp2=subset(mydata, !spontpreg3)
df3surv2=data.frame(GRDBSpreg3=df3temp2$GRDBSpreg3, delivery=FALSE, CSCS=df3temp2$CSCS, MALDERpreg3=df3temp2$MALDERpreg3, BMIpreg3=df3temp2$BMIpreg3, notborninswepreg3=df3temp2$notborninswepreg3, smokepreg3=df3temp2$smokepreg3, diabetespreg3=df3temp2$diabetespreg3, allHTNpreg3=df3temp2$allHTNpreg3, boychildpreg3=df3temp2$boychildpreg3, malformpreg3=df3temp2$malformpreg3, SGApreg3marsal=df3temp2$SGApreg3marsal, LGApreg3marsal=df3temp2$LGApreg3marsal)

df3surv=rbind(df3surv1,df3surv2)

cox = coxph(Surv(GRDBSpreg3, delivery) ~ CSCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, data = df3surv)

summary(cox)



#Cox model 2 on variables vagvag, vagCS, CSvag, CSCS

df3temp1=subset(mydata, spontpreg3)
df3surv1=data.frame(GRDBSpreg3=df3temp1$GRDBSpreg3, delivery=TRUE, vagvag=df3temp1$vagvag, vagCS=df3temp1$vagCS, CSvag=df3temp1$CSvag, CSCS=df3temp1$CSCS, MALDERpreg3=df3temp1$MALDERpreg3, BMIpreg3=df3temp1$BMIpreg3, notborninswepreg3=df3temp1$notborninswepreg3, smokepreg3=df3temp1$smokepreg3, diabetespreg3=df3temp1$diabetespreg3, allHTNpreg3=df3temp1$allHTNpreg3, boychildpreg3=df3temp1$boychildpreg3, malformpreg3=df3temp1$malformpreg3, SGApreg3marsal=df3temp1$SGApreg3marsal, LGApreg3marsal=df3temp1$LGApreg3marsal)

df3temp2=subset(mydata, !spontpreg3)
df3surv2=data.frame(GRDBSpreg3=df3temp2$GRDBSpreg3, delivery=FALSE, vagvag=df3temp2$vagvag, vagCS=df3temp2$vagCS, CSvag=df3temp2$CSvag, CSCS=df3temp2$CSCS, MALDERpreg3=df3temp2$MALDERpreg3, BMIpreg3=df3temp2$BMIpreg3, notborninswepreg3=df3temp2$notborninswepreg3, smokepreg3=df3temp2$smokepreg3, diabetespreg3=df3temp2$diabetespreg3, allHTNpreg3=df3temp2$allHTNpreg3, boychildpreg3=df3temp2$boychildpreg3, malformpreg3=df3temp2$malformpreg3, SGApreg3marsal=df3temp2$SGApreg3marsal, LGApreg3marsal=df3temp2$LGApreg3marsal)

df3surv=rbind(df3surv1,df3surv2)

cox = coxph(Surv(GRDBSpreg3, delivery) ~ vagvag + vagCS + CSvag + CSCS + MALDERpreg3 + BMIpreg3 + notborninswepreg3 + smokepreg3 + diabetespreg3 + allHTNpreg3 + boychildpreg3 + malformpreg3 + SGApreg3marsal + LGApreg3marsal, data = df3surv)

summary(cox)




rm(df3temp1,df3temp2,df3surv1,df3surv2,df3surv,survobject,fit,mysurvplot,ggsurvplot,cox)

rm(grav1,grav2,grav3,mydata,mydata2,mydata3)



