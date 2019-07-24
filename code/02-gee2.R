# Current version of GEE code

library(geepack)
library(gee)
library(ggpubr)
library(ggplot2)
library(splines2)
library(scales)
library(dplyr)

dir <- "/Users/shen/Documents/research/samsi/data/"
setwd(dir)

data<-readRDS(file = "ds_cgm_complete_3days.rds")

used.id<-c(5, 7, 11, 18, 21, 24, 25, 26, 29, 32, 34, 35, 41, 46, 47, 49, 50, 53, 55, 69, 74, 82, 83, 88, 89, 91, 94, 96, 101, 118, 120, 124, 127, 129, 133, 135, 150, 153, 156, 158, 159, 164, 167, 168, 171, 176, 184, 190, 191, 192, 194, 203, 204, 205, 207, 212, 213, 219, 220, 225)
data.subset<-subset(data, data$pt_id %in% used.id)
data.subset<-subset(data.subset, data.subset$visitnum %in% c(1,2,3))
data.subset$day.in.visit<-data.subset$day_in_subrange
data.subset$group<-ifelse(data.subset$trt_group=="CGM Only",0, 1)
### in specific visit####

choose.visit<-3
data.used<-as.data.frame(subset(data.subset,data.subset$visitnum==choose.visit))
data.used$device_tm_order_in_day<-data.used$device_tm_bin


t <-data.used$device_tm_order_in_day
knots.choose<- round(288*(c(1:6)/7))
bsMat <- bSpline(t, knots = knots.choose, degree = 3, intercept=T)
data.used$group_spline_day<-data.used$group*bsMat
data.used$bsMat<-bsMat
gee2<-gee(glucose_value ~ bsMat+group_spline_day-1, family=gaussian,data=data.used, id=pt_id,corstr="exchangeable")

spline_beta0<-matrix((summary(gee2)$coefficients)[c(1:10),1],nrow=10,ncol = 1)
spline_beta<-matrix((summary(gee2)$coefficients)[c(11:20),1],nrow=10,ncol = 1)
se_beta0<-gee2$naive.variance[c(1:10),c(1:10)]
se_beta1<-gee2$naive.variance[c(1:10),c(1:10)]+gee2$naive.variance[c(11:20),c(11:20)]+2*gee2$naive.variance[c(1:10),c(11:20)]

estimated_matrix0<-bsMat%*%spline_beta0
se_fitted0<-sqrt(diag(bsMat%*%se_beta0%*%t(bsMat)))

estimated_matrix<-bsMat%*%spline_beta

estimated_matrix1<-estimated_matrix0+estimated_matrix
se_fitted1<-sqrt(diag(bsMat%*%se_beta1%*%t(bsMat)))

uplimit0<-estimated_matrix0+se_fitted0*1.96
lowlimit0<-estimated_matrix0-se_fitted0*1.96

uplimit1<-estimated_matrix1+se_fitted1*1.96
lowlimit1<-estimated_matrix1-se_fitted1*1.96

plot.data.1<-as.data.frame(cbind(data.used$device_tm,estimated_matrix0,uplimit0,lowlimit0))
plot.data.2<-as.data.frame(cbind(data.used$device_tm,estimated_matrix1,uplimit1,lowlimit1))
plot.data.1$group.name<-"CGM"
plot.data.2$group.name<-"CGM+BGM"
plot.data<-rbind(plot.data.1,plot.data.2)

colnames(plot.data)<-c("t","fit","upper","lower","group.name")

plot.new<-plot.data %>%
  mutate(tclock = hms::as_hms(t)) 

plot.new$tclock<-as.POSIXct(plot.new$tclock)

ggplot(plot.new, aes(x=tclock, y=fit,group = group.name,color=group.name))+ geom_line(aes(color=group.name))+
  geom_smooth(aes(ymin = lower, 
                  ymax = upper),
              stat = "identity") +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours") +
  xlab("time in mentoring(in visit 3)") +
  ylab("treatment effect")

### every day in each visit ####
choose.visit<-3
data.used<-subset(data.subset,data.subset$visitnum==choose.visit)

t <-data.used$device_tm_order_in_day

knots.choose<- round(288*(c(1:6)/7))
bsMat <- bSpline(t, knots = knots.choose, degree = 3, intercept=T)
data.used$day1.ind<-ifelse(data.used$day.in.visit==1,1,0)
data.used$day2.ind<-ifelse(data.used$day.in.visit==2,1,0)
data.used$day3.ind<-ifelse(data.used$day.in.visit==3,1,0)
#data.used$day4.ind<-ifelse(data.used$day.in.visit==4,1,0)
#data.used$day5.ind<-ifelse(data.used$day.in.visit==5,1,0)
#data.used$day6.ind<-ifelse(data.used$day.in.visit==6,1,0)
#data.used$day7.ind<-ifelse(data.used$day.in.visit==7,1,0)
#data.used$day8.ind<-ifelse(data.used$day.in.visit==8,1,0)

data.used$group_spline_day1<-data.used$group*bsMat*data.used$day1.ind
data.used$group_spline_day2<-data.used$group*bsMat*data.used$day2.ind
data.used$group_spline_day3<-data.used$group*bsMat*data.used$day3.ind
#data.used$group_spline_day4<-data.used$group*bsMat*data.used$day4.ind
#data.used$group_spline_day5<-data.used$group*bsMat*data.used$day5.ind
#data.used$group_spline_day6<-data.used$group*bsMat*data.used$day6.ind
#data.used$group_spline_day7<-data.used$group*bsMat*data.used$day7.ind
#data.used$group_spline_day8<-data.used$group*bsMat*data.used$day8.ind

data.used$bsMat<-bsMat
gee2<-gee(glucose_value ~ bsMat+group_spline_day1+group_spline_day2+group_spline_day3-1, family=gaussian,data=data.used, id=pt_id, corstr = "exchangeable")
#gee2<-gee(glucose_value ~ bsMat+group_spline_day1+group_spline_day2+group_spline_day3+group_spline_day4+group_spline_day5+group_spline_day6+group_spline_day7+group_spline_day8-1, family=gaussian,data=data.used, id=pt_id, corstr = "independence")


spline_beta_1<-matrix((summary(gee2)$coefficients)[c(11:20),1],nrow=10,ncol = 1)
spline_beta_2<-matrix((summary(gee2)$coefficients)[c(21:30),1],nrow=10,ncol = 1)
spline_beta_3<-matrix((summary(gee2)$coefficients)[c(31:40),1],nrow=10,ncol = 1)
#spline_beta_4<-matrix((summary(gee2)$coefficients)[c(41:50),1],nrow=10,ncol = 1)
#spline_beta_5<-matrix((summary(gee2)$coefficients)[c(51:60),1],nrow=10,ncol = 1)
#spline_beta_6<-matrix((summary(gee2)$coefficients)[c(61:70),1],nrow=10,ncol = 1)
#spline_beta_7<-matrix((summary(gee2)$coefficients)[c(71:80),1],nrow=10,ncol = 1)
#spline_beta_8<-matrix((summary(gee2)$coefficients)[c(81:90),1],nrow=10,ncol = 1)


se_beta_1<-gee2$naive.variance[c(11:20),c(11:20)]
se_beta_2<-gee2$naive.variance[c(21:30),c(21:30)]
se_beta_3<-gee2$naive.variance[c(31:40),c(31:40)]
#se_beta_4<-gee2$naive.variance[c(41:50),c(41:50)]
#se_beta_5<-gee2$naive.variance[c(51:60),c(51:60)]
#se_beta_6<-gee2$naive.variance[c(61:70),c(61:70)]
#se_beta_7<-gee2$naive.variance[c(71:80),c(71:80)]
#se_beta_8<-gee2$naive.variance[c(81:90),c(81:90)]

estimated_matrix_1<-bsMat%*%spline_beta_1
estimated_matrix_2<-bsMat%*%spline_beta_2
estimated_matrix_3<-bsMat%*%spline_beta_3
#estimated_matrix_4<-bsMat%*%spline_beta_4
#estimated_matrix_5<-bsMat%*%spline_beta_5
#estimated_matrix_6<-bsMat%*%spline_beta_6
#estimated_matrix_7<-bsMat%*%spline_beta_7
#estimated_matrix_8<-bsMat%*%spline_beta_8

se_fitted_1<-sqrt(diag(bsMat%*%se_beta_1%*%t(bsMat)))
se_fitted_2<-sqrt(diag(bsMat%*%se_beta_2%*%t(bsMat)))
se_fitted_3<-sqrt(diag(bsMat%*%se_beta_3%*%t(bsMat)))
#se_fitted_4<-sqrt(diag(bsMat%*%se_beta_4%*%t(bsMat)))
#se_fitted_5<-sqrt(diag(bsMat%*%se_beta_5%*%t(bsMat)))
#se_fitted_6<-sqrt(diag(bsMat%*%se_beta_6%*%t(bsMat)))
#se_fitted_7<-sqrt(diag(bsMat%*%se_beta_7%*%t(bsMat)))
#se_fitted_8<-sqrt(diag(bsMat%*%se_beta_8%*%t(bsMat)))


uplimit_1<-estimated_matrix_1+se_fitted_1*1.96
lowlimit_1<-estimated_matrix_1-se_fitted_1*1.96
uplimit_2<-estimated_matrix_2+se_fitted_2*1.96
lowlimit_2<-estimated_matrix_2-se_fitted_2*1.96
uplimit_3<-estimated_matrix_3+se_fitted_2*1.96
lowlimit_3<-estimated_matrix_3-se_fitted_2*1.96
#uplimit_4<-estimated_matrix_4+se_fitted_2*1.96
#lowlimit_4<-estimated_matrix_4-se_fitted_2*1.96
#uplimit_5<-estimated_matrix_5+se_fitted_2*1.96
#lowlimit_5<-estimated_matrix_5-se_fitted_2*1.96
#uplimit_6<-estimated_matrix_6+se_fitted_2*1.96
#lowlimit_6<-estimated_matrix_6-se_fitted_2*1.96
#uplimit_7<-estimated_matrix_7+se_fitted_2*1.96
#lowlimit_7<-estimated_matrix_7-se_fitted_2*1.96
#uplimit_8<-estimated_matrix_8+se_fitted_2*1.96
#lowlimit_8<-estimated_matrix_8-se_fitted_2*1.96


plot.data<-as.data.frame(cbind(data.used$device_tm,estimated_matrix_1,uplimit_1,lowlimit_1,estimated_matrix_2,uplimit_2,lowlimit_2,estimated_matrix_3,uplimit_3,lowlimit_3))

colnames(plot.data)<-c("t","y1","uplimit1","lowlimit1","y2","uplimit2","lowlimit2","y3","uplimit3","lowlimit3")
plot.new<-plot.data %>%
  mutate(tclock = hms::as_hms(t)) 

plot.new$tclock<-as.POSIXct(plot.new$tclock)

ggplot(data = plot.new, aes(x = tclock)) + 
  geom_hline(yintercept = 0, linetype="dashed", color = "red") + 
  geom_smooth(aes(y = y1, ymin = lowlimit_1, ymax = uplimit_1), stat = "identity", color = "blue") + 
  geom_smooth(aes(y = y2, ymin = lowlimit_2, ymax = uplimit_2), stat = "identity", color = "orange") +
  geom_smooth(aes(y = y3, ymin = lowlimit_3, ymax = uplimit_3), stat = "identity", color = "red") +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours") +
  ylim(-55,55) +
  xlab("time in mentoring (in one visit)") +
  ylab("group difference")






