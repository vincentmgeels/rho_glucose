library(geepack)
library(gee)
library(ggpubr)
library(ggplot2)
library(splines2)
library(scales)
library(dplyr)

dir <- "/Users/shen/Documents/research/samsi/data/"
setwd(dir)

data<-readRDS(file = "data_clean_days_0719.rds")
data.one<-subset(data,data$group==0)
data.two<-subset(data,data$group==1)
##############
# clean data #
##############
# demo<- readRDS(file = "ds_subject.rds")
# cgm<- readRDS(file = "ds_cgm.rds")
# data<-merge(cgm, demo, by="pt_id", all.x=TRUE)
# data$group<-ifelse(data$trt_group=="CGM Only",0, 1)
# data$gender.ind<-ifelse(data$gender=="F",0, 1)
# 
# 
# ### add days in each visit###
# total.id<-unique(cgm$pt_id)
# data.test<-rep()
# for (i in 1:length(total.id)) {
#  for (j in c(0:5)) {
#    test<- subset(data, data$pt_id==total.id[i] & data$visitnum==j)
#    start.day<-min(unique(test$device_dt_tm_days_from_enroll))
#    test$day.in.visit<-test$device_dt_tm_days_from_enroll-start.day+1
#    data.test<-rbind(data.test,test)
#  } 
#   print(i)
# }
# 
# ### add recording order in every day ###
# data.test$device_tm_order_in_day<-round(as.numeric(data.test$device_tm)/300)
# 
# ### add recording order in every visit ###
# data.test$device_tm_order_in_visit<-round(as.numeric(data.test$device_tm)*data.test$day.in.visit/300)
# 
# saveRDS(data.test, file = "data_clean_days_0719.rds")
# 


#######
# GEE #
#######

data.subset<-data
#summary(gee(glucose_value ~ group, family=gaussian,data=data.subset, id=pt_id, corstr = "independence"))$coefficients

result.gee1<-matrix(NA, nrow = 6,ncol = 4)
for (i in 1:6) {
  data.used<-subset(data.subset,data.subset$visitnum==i-1)
  gee1<-summary(gee(glucose_value ~ group, family=gaussian,data=data.used, id=pt_id, corstr = "independence"))$coefficients
  result.gee1[i,1]<-i-1
  result.gee1[i,2]<-gee1[2,1]
  result.gee1[i,3]<-gee1[2,1]-1.96*gee1[2,4]
  result.gee1[i,4]<-gee1[2,1]+1.96*gee1[2,4]
  print(i)
}

colnames(result.gee1)<-c("visit","effect","low.ci","high.ci")
result.gee1<-as.data.frame(result.gee1)
ggplot(result.gee1, aes(x=visit, y=effect, group=1)) + geom_hline(yintercept=0, linetype="dashed", color = "blue")+
  geom_point(size=10,alpha=0.52) +
  geom_errorbar(width=.1, aes(ymin=low.ci, ymax=high.ci), colour="darkred") + 
  labs(x="visit",y= "group effect difference", title="Group difference estimation and CI with GEE", caption="source: A Randomized Trial Comparing CGM With and Without Routine BGM in Adults with T1D.") 

