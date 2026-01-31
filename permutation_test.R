#Load data
start <- read.csv('Data/StartDates_Sectors.csv')

data_treated <- read_xlsx('Data/RCT_FAlbo_FAeg.xlsx',sheet=1)
data_treated$GAI <- data_treated$F.Ae.aegypti/data_treated$FuncTrap


data_control <- read_xlsx('Data/RCT_FAlbo_FAeg.xlsx',sheet=3)
data_control$GAI <- data_control$F.Ae.aegypti/data_control$FuncTrap


data <- rbind(data_treated[,c("Eyear","Eweek","StudyArea","Sector_ID","GAI")], 
              data_control[,c("Eyear","Eweek","StudyArea","Sector_ID","GAI")])

data$EYEW <- paste0(data$Eyear, " EW",data$Eweek)

int.sites <- start$Sector_ID 
data <- data %>% mutate(treated=ifelse(Sector_ID %in% int.sites,1,0))

times <- data.frame(EYEW=unique(data$EYEW), time=1000:1290)
start <- left_join(start,times,by='EYEW')
data <- left_join(data,times,by='EYEW')

# Join Covariates Data
vars <- read.csv('Data/SG_EnvironmentVar_13Sep2024_wLags.csv') 
vars <- subset(vars, select = c(Sector_ID, Eweek, Eyear, EYEW, MeanT, MeanWS,NDVI_A,X300m.P,
                                A_HDB_A_New,A_HDB_H_New,HDB_3R_P, HDB_4R_P,HDB_5R_P,
                                length_D, Forest_P, Grass_P,
                                MVege_P, tp, r,MeanT_L1,MeanT_L2,MeanT_L3,MeanT_L4,
                                tp_L1,tp_L2,tp_L3,tp_L4,r_L1,r_L2,r_L3,r_L4))
vars$a_HDB_P <- rowMeans(vars[,c('HDB_3R_P', 'HDB_4R_P','HDB_5R_P')])
vars[,c('HDB_3R_P', 'HDB_4R_P','HDB_5R_P')] <- NULL
vars <- vars %>% rename(A_HDB_A = A_HDB_A_New, A_HDB_H=A_HDB_H_New)

data <- data %>% left_join(vars[,-c(2,3)],by=c('Sector_ID','EYEW')) 

#Sampling clusters
sites <- unique(data$StudyArea)
int_samples <- vector('list',10000)
for (n in 1:length(int_samples)) {
  int_samples[[n]] <- sample(sites,8)
}

#Checking balance
md<-vector(length=10000)
pre_int_time <- 1179 #earliest release
for (n in 1:length(int_samples)) {
  int <- int_samples[[n]]
  con <- sites[!sites %in% int]
  int_df <- data %>% filter(StudyArea %in% int & time<pre_int_time) %>% summarise(avg=mean(GAI), var=var(GAI))
  con_df <- data %>% filter(StudyArea %in% con & time<pre_int_time)%>% summarise(avg=mean(GAI), var=var(GAI))
  
  md[n] <- int_df$avg-con_df$avg
}

sum(abs(md)<.01)

i<-which(abs(md)<.01)
per_set <- int_samples[i]

atts <- vector('list',length(per_set))
for (s in 1:length(per_set)) {
  int_clusters <- per_set[[s]]
  tmp_data <- data %>% mutate(
    release = case_when(
      StudyArea %in% int_clusters ~ 1,
      .default = 0
    ),
    start_time = case_when(
      StudyArea %in% int_clusters ~pre_int_time,
      .default = Inf
    )
  )
  
  tmp_data$event_time <- tmp_data$time - tmp_data$start_time
  
  mod_out = feglm(
    GAI ~ MeanT+ MeanWS+MeanT_L1+MeanT_L2+MeanT_L3+MeanT_L4+tp_L1+tp_L2+tp_L3+tp_L4+r_L1+r_L2+r_L3+r_L4+tp+ r+ sunab(start_time, event_time, bin.rel="bin::13") | Sector_ID + time, 
    data = tmp_data, vcov = ~Sector_ID)
  ates[[s]] <- mod_out$coefficients
  atts[[s]] <- aggregate(mod_out, c("ATT" = "event_time::[^-]"))
  
  binned <- bin(sort(unique(tmp_data$event_time)),"bin::13")
  binned <- binned[-1] #removing -inf
  agg_data = tmp_data %>% filter(StudyArea %in% int_clusters) %>% group_by(event_time) %>%
    select(GAI, event_time) %>% summarise(gai_point=mean(GAI))
  agg_data$bin <- binned
  binned_gai <- agg_data %>% group_by(bin) %>% summarise(gai=mean(gai_point))
  IEs[[s]] <- -100*tail(mod_out$coefficients,23)/(binned_gai$gai-tail(mod_out$coefficients,23))
}

ates_main <- tail(m_glm_covs_bin$coeftable[,1],9) 
sim_ates <- lapply(ates, function(x){tail(x,9)})

for (n in 1:9) {
  t <- unlist(lapply(sim_ates,function(x)x[n]))
  p <- sum(t<ates_main[n])/length(t)
  print(paste0(names(ates_main[n])," p-value= ",p))
}

cm=1/2.54

library(stringi)

pdf('event-time study/plots/perm_test_ates_ae.pdf',width=18*cm,height = 20*cm)
par(mfrow=c(3,3))
for (n in 1:9) {
  t1  <- as.integer(stri_split_fixed(names(ates_main),"::")[[n]][2])
  main <- ifelse(n<length(ates_main),paste0('Exposure weeks ',t1,' to ',t1+12),"Exposure weeks 104+")
  hist(unlist(lapply(sim_ates,function(x)x[[n]])),xlim=c(-.3,.3), main=main,xlab='ATE' )
  lines(x=rep(ates_main[n],2),y=c(-1,3000),col='red',lwd=2)
}
dev.off()


