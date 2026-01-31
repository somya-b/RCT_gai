
{ # load libraries
  library(dplyr)
  library(readxl)
  library(fixest)
}

start <- read.csv('Data/StartDates_Sectors.csv')
data_control <- read_xlsx('Data/RCT_FAlbo_FAeg.xlsx',sheet=3)
data_control$GAI <- data_control$F.Ae.aegypti/data_control$FuncTrap
# data_control$GAI <- data_control$F.Ae.albopictus/data_control$FuncTrap

data <- data_control %>% subset((Eyear<2024)|(Eyear==2024 & Eweek<38))
data$EYEW <- paste0(data$Eyear, " EW",data$Eweek)

ctrl_sec <- unique(data[c('StudyArea','Sector_ID')])
ctrl_area <- unique(data$StudyArea)

times <- data.frame(EYEW=unique(data$EYEW), time=1000:1290)
data <- left_join(data,times,by='EYEW')

int_start <- 1179
nreps <- 1000

tp <- as.character(seq(0,111,by=13))

att_point <- data.frame(matrix(NA, ncol=24, nrow=nreps))
colnames(att_point) <- c("-179",rev(seq(-13,-179,by=-13)),"-1",tp)
att_lb <- att_ub <- att_point
gai_binned <- att_point



for (i in 1:nreps){
  
  int.sites <- sample(ctrl_sec$Sector_ID, size=30) #sample control sectors to be treated
  df <- data %>% mutate(start_time=ifelse(Sector_ID %in% int.sites,int_start,Inf))
  df$event_time <- df$time - df$start_time
  
  model = feglm(
    GAI ~ sunab(start_time,event_time, bin.rel="bin::13") | Sector_ID + time, 
    data = df, vcov = ~Sector_ID)
  
  p <- iplot(model)$prms
  att_point[i,] <- p[,"estimate"]
  att_lb[i,] <- p[,"ci_low"]
  att_ub[i,] <- p[,"ci_high"]
  
  # extracting GAI to calulate IEs
  binned <-bin(sort(unique(df$event_time)),"bin::13")
  binned <- binned[-1] #removing -inf
  agg_data = df %>% filter(Sector_ID %in% int.sites) %>% group_by(event_time) %>%
    select(GAI, event_time) %>% summarise(gai_point=mean(GAI))
  agg_data$bin <- binned
  agg_data$bin[agg_data$event_time == -1 ] = -1
  binned_gai <- agg_data %>% group_by(bin) %>% summarise(gai=mean(gai_point))
  
  gai_binned[i,] <- binned_gai$gai
  
}

den <- gai_binned - att_point
ie <- -100*att_point/den

cm=1/2.54
pdf(paste0("event-time study/",folder,"/control_histograms.pdf"), width = 18.5*cm, height = 18.5*cm)
bins <- seq(0,111,by=13)
par(mfrow=c(3,3), mar=c(4,4,3,3),cex.lab=1.2,las=1)
k=1;
for(i in 16:ncol(ie)){
  if(i<24){
    hist(ie[,i], main=paste0("Bin: ", bins[k]," to ", bins[k+1]-1, " weeks"), 
         xlab='Intervention Efficacy (%)')
    avg = mean(ie[,i])
    lines(x=c(avg,avg),y=c(-1,600),col='red',lwd=1.5)
  } else {
    hist(ie[,i], main="Bin: 104+ weeks", xlab='Intervention Efficacy (%)')
    avg = mean(ie[,i])
    lines(x=c(avg,avg),y=c(-1,600),col='red',lwd=1.5)
  }
  k=k+1
}
dev.off()




