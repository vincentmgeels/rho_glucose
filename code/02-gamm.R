library(nlme)
library(mgcv)
library(ggplot2)
library(splines2)
library(scales)
library(dplyr)
#dir <- "/Users/shen/Documents/research/samsi/data/"
#setwd(dir)

data<-readRDS(file = "./data/ds_cgm_complete_3days.rds")

used.id<-c(5, 7, 11, 18, 21, 24, 25, 26, 29, 32, 34, 35, 41, 46, 47, 49, 50, 53, 55, 69, 74, 82, 83, 88, 89, 91, 94, 96, 101, 118, 120, 124, 127, 129, 133, 135, 150, 153, 156, 158, 159, 164, 167, 168, 171, 176, 184, 190, 191, 192, 194, 203, 204, 205, 207, 212, 213, 219, 220, 225)
data.subset<-as.data.frame(subset(data, data$pt_id %in% used.id))
data.subset<-subset(data.subset, data.subset$visitnum %in% c(1,2,3))
data.subset$day.in.visit<-data.subset$day_in_subrange
data.subset$group<-ifelse(data.subset$trt_group=="CGM Only",0, 1)
data.subset$device_tm_order_in_day<-data.subset$device_tm_bin

data.used<-data.subset
data.used$group<-as.factor(data.used$group)
gammm1<-gamm(glucose_value ~ group+s(device_tm_order_in_day, bs= "cc", by=group,k=7) ,family=gaussian, data=data.used, correlation=corAR1(form=~1|visitnum), random = list(pt_id=~1,visitnum = ~1))

mm3 <- expand.grid(device_tm_order_in_day=seq(0,288,length=1000),group=levels(data.used$group))
pp3 <- predict.gam(gammm1$gam,mm3,se.fit=T,type='response')

rr3 <- data.frame(dv="2 exac2",mm3,pp3[1],pp3[2])
rr3$lower <- with(rr3,fit- qnorm(0.95)*se.fit)
rr3$upper <- with(rr3,fit+ qnorm(0.95)*se.fit)
rr3$group.name<-ifelse(rr3$group==0,"CGM","CGM+BGM") 

rr3$t<-rr3$device_tm_order_in_day*5*60
plot.new<-rr3 %>%
  mutate(tclock = hms::as_hms(t)) 

plot.new$tclock<-as.POSIXct(plot.new$tclock)

ggplot(plot.new, aes(x=tclock, y=fit,group = group.name,color=group.name))+ geom_line(aes(color=group.name))+
  geom_smooth(aes(ymin = lower, 
                  ymax = upper),
              stat = "identity") +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours") +
  xlab("time in mentoring(day)") +
  ylab("treatment effect")

