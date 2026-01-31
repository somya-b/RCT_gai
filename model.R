library(fixest)
library(readxl)
library(dplyr)

#Load data
start <- read.csv('Data/StartDates_Sectors.csv')

data_treated <- read_xlsx('Data/RCT_FAlbo_FAeg.xlsx',sheet=1)
data_treated$GAI <- data_treated$F.Ae.aegypti/data_treated$FuncTrap
# data_treated$GAI <- data_treated$F.Ae.albopictus/data_treated$FuncTrap

data_control <- read_xlsx('Data/RCT_FAlbo_FAeg.xlsx',sheet=3)
data_control$GAI <- data_control$F.Ae.aegypti/data_control$FuncTrap
# data_control$GAI <- data_control$F.Ae.albopictus/data_control$FuncTrap

data <- rbind(data_treated[,c("Eyear","Eweek","StudyArea","Sector_ID","GAI")], 
              data_control[,c("Eyear","Eweek","StudyArea","Sector_ID","GAI")])

#RCT only till 13Sept2024 EW37 2024
data <- data %>% subset((Eyear<2024)|(Eyear==2024 & Eweek<38))
data$EYEW <- paste0(data$Eyear, " EW",data$Eweek)

int.sites <- start$Sector_ID 
data <- data %>% mutate(treated=ifelse(Sector_ID %in% int.sites,1,0))

times <- data.frame(EYEW=unique(data$EYEW), time=1000:1290)
start <- left_join(start,times,by='EYEW')
data <- left_join(data,times,by='EYEW')

start_30 <- start$Sector_ID[start$EYEW=="2022 EW30"]
start_35 <- start$Sector_ID[start$EYEW=="2022 EW35"]
start_37 <- start$Sector_ID[start$EYEW=="2022 EW37"]

data <- data %>% mutate(
  release = case_when(
    Sector_ID %in% start_30 & ((Eyear == 2022 & Eweek  >29) | Eyear > 2022) ~ 1,
    Sector_ID %in% start_35 & ((Eyear == 2022 & Eweek  >34) | Eyear > 2022) ~ 1,
    Sector_ID %in% start_37  & ((Eyear == 2022 & Eweek  >36) | Eyear > 2022) ~ 1,
    .default = 0
  )
)

data <- data %>% left_join(start %>% select(Sector_ID,time), by='Sector_ID') %>% 
  rename(start_time='time.y', time='time.x')
data$start_time[is.na(data$start_time)] <- Inf
data$event_time <- data$time - data$start_time

# Join Covariates Data
vars <- read.csv('Data/SG_EnvironmentVar.csv') 
vars <- subset(vars, select = c(Sector_ID, Eweek, Eyear, EYEW, MeanT, MeanWS,NDVI_A,X300m.P,
                                A_HDB_A_New,A_HDB_H_New,HDB_3R_P, HDB_4R_P,HDB_5R_P,
                                length_D, Forest_P, Grass_P,
                                MVege_P, tp, r,MeanT_L1,MeanT_L2,MeanT_L3,MeanT_L4,
                                tp_L1,tp_L2,tp_L3,tp_L4,r_L1,r_L2,r_L3,r_L4))
vars$a_HDB_P <- rowMeans(vars[,c('HDB_3R_P', 'HDB_4R_P','HDB_5R_P')])
vars[,c('HDB_3R_P', 'HDB_4R_P','HDB_5R_P')] <- NULL
vars <- vars %>% rename(A_HDB_A = A_HDB_A_New, A_HDB_H=A_HDB_H_New)

data <- data %>% left_join(vars[,-c(2,3)],by=c('Sector_ID','EYEW')) 


## ------ TWFE without Covariates -----
m_glm_bin = feglm(
  GAI ~ sunab(start_time,event_time, bin.rel="bin::13") | Sector_ID + time, 
  data = data, vcov = ~Sector_ID)
cm=1/2.54


### ----- Calculate IE -----

#Binning GAI to match ATE bins
binned <- bin(sort(unique(data$event_time)),"bin::13")
binned <- binned[-1] #removing -inf
agg_data = data %>% filter(Sector_ID %in% int.sites) %>% group_by(event_time) %>%
  select(GAI, event_time) %>% summarise(gai_point=mean(GAI))
agg_data$bin <- binned
binned_gai <- agg_data %>% group_by(bin) %>% summarise(gai=mean(gai_point))

#IE
ate=m_glm_bin$coeftable[,1]
se_ate=m_glm_bin$coeftable[,2]
g_der<-1/(binned_gai$gai-ate)+ ate/(binned_gai$gai-ate)^2
se=-100*g_der*se_ate
ie <- data.frame(estimate=-100*ate/(binned_gai$gai-ate))
ie$ci_low <- ie$estimate - 1.96*se
ie$ci_high <- ie$estimate + 1.96*se


## ------ TWFE with Covariates -----
# Only time-variant covariates can be used
m_glm_covs_bin = feglm(
  GAI ~ MeanT+ MeanWS+MeanT_L1+MeanT_L2+MeanT_L3+MeanT_L4+tp_L1+tp_L2+tp_L3+tp_L4+r_L1+r_L2+r_L3+r_L4+tp+ r+ sunab(start_time, event_time, bin.rel="bin::13") | Sector_ID + time, 
  data = data, vcov = ~Sector_ID)

pdf('event-time study/plots/aegypti/ate_aegypti_withCov.pdf', width = 12*cm, height =12*cm)
iplot(m_glm_covs_bin, xlab = "Time to treatment", main="Effect on GAI (GLM with covs)" , ref.line = 0)
dev.off()
etable(m_glm_covs_bin)


### ----- Calculate IE -----

# Binning GAI data
binned <-bin(sort(unique(data$event_time)),"bin::13")
binned <- binned[-1] #removing -inf
agg_data = data %>% filter(Sector_ID %in% int.sites) %>% group_by(event_time) %>%
  select(GAI, event_time) %>% summarise(gai_point=mean(GAI))
agg_data$bin <- binned
binned_gai <- agg_data %>% group_by(bin) %>% summarise(gai=mean(gai_point))

#IE
ate_cov=m_glm_covs_bin$coeftable[-c(1:16),1]
se_ate=m_glm_covs_bin$coeftable[-c(1:16),2]
g_der<-1/(binned_gai$gai-ate_cov)+ ate_cov/(binned_gai$gai-ate_cov)^2
se_ie=-100*g_der*se_ate
ie_cov <- data.frame(estimate=-100*ate_cov/(binned_gai$gai-ate_cov))
ie_cov$ci_low <- ie_cov$estimate - 1.96*se_ie
ie_cov$ci_high <- ie_cov$estimate + 1.96*se_ie



## ---- In time check ------
df = data %>% filter(time < 1179) # removing all treatment period
pseudo_int = 52 #pseudo intervention weeks used: 26, 52
df <- df %>% mutate(start_time = case_when(treated == 1 ~ 1178-pseudo_int,
                                           treated == 0 ~ Inf))
df$event_time = df$time - df$start_time

### ----- without Covariates -----
m_glm = feglm(
  GAI ~ sunab(start_time,event_time, bin.rel="bin::13") | Sector_ID + time, 
  data = df, vcov = ~Sector_ID)

cm=1/2.54
folder="aegypti"

pdf(paste0('plots/',folder,'/ae_intime_', pseudo_int,'.pdf'), width = 12*cm, height =12*cm)
iplot(m_glm,xlab = "Time to treatment", 
      main=paste0("Pseudo intervention: ",pseudo_int," weeks"),
      ref.line = 0)
dev.off()

### ----- with Covariates -----
m_glm = feglm(
  GAI ~ MeanT+ MeanWS+MeanT_L1+MeanT_L2+MeanT_L3+MeanT_L4+tp_L1+tp_L2+tp_L3+tp_L4+r_L1+r_L2+r_L3+r_L4+tp+ r+ sunab(start_time,event_time, bin.rel="bin::13") | Sector_ID + time, 
  data = df, vcov = ~Sector_ID)

cm=1/2.54
pdf(paste0('plots/',folder,'/ae_intim_withCov_', pseudo_int,'.pdf'), width = 12*cm, height =12*cm)
iplot(m_glm,xlab = "Time to treatment", 
      main=paste0("Pseudo intervention: ",pseudo_int," weeks"),
      ref.line = 0)
dev.off()



