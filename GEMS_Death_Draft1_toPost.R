#######################################################
# Code for death prediction
#######################################################
rm(list=ls()) #this clears the workspace
graphics.off(); #close all graphics windows


library(tidyverse)
library(lubridate)
library(gridExtra)
library(pROC)
library(data.table)
library(fields)
library(zoo)
library(ks)
library(KernSmooth)
library(ranger)
library(viridis)
library(purrr)
library(broom)
library(profvis)
library(furrr)
library(mice)
library(glmnet)
library(glmnetUtils)
library(cvAUC) #this is a wrapper of AUC inside the ROCR package
library(table1)
library(data.table)
library(bit64)
library(pdp)

setwd('')

################### FUNCTion for var screening, regression fitting, AUC calc ####
CPR.funct <- function(data,outcome,iter,nvars_opts){
  out=ranger(as.formula(paste(outcome,'~',paste(names,collapse="+"),sep="")),data=data,num.trees=1000,importance="impurity")
  imps=importance(out)
  df_imps_full=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
  
  
  result=data.frame(iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)
  test_record <- NA
  train_record <- NA
  
  for (each in 1:iter){
    print(each)
    train=data %>% sample_frac(.80,replace=F)
    
    test=data[-which(data$index %in% train$index),]
    
    train_record <- c(train_record,table(train[,outcome])[["1"]])
    test_record <- c(test_record,table(test[,outcome])[["1"]])
    
    out=ranger(as.formula(paste(outcome,'~',paste(names,collapse="+"),sep="")),data=train,num.trees=1000,importance="impurity")
    df_imps=data.frame(names=names(ranger::importance(out)),imps=ranger::importance(out)) %>% arrange(desc(imps))
    for (nvars in nvars_opts){
      
      print(nvars)
      out1=glm(as.formula(paste(outcome,'~',paste(df_imps$names[1:nvars],collapse="+"),sep="")),data=train,family="binomial",control=glm.control(maxit=50))
      out2=ranger(as.formula(paste(outcome,'~',paste(df_imps$names[1:nvars],collapse="+"),sep="")),data=train,num.trees=1000)
      
      df=data.frame(iter=each,nvar=nvars,true=test[[outcome]],pred_glm=as.numeric(predict(out1,newdata=test,type="response")),pred_RF=as.numeric(predict(out2,data=test,type="response")$predictions))
      result=rbind(result,df)
    }
  }
  result<-result[-1,]
  train_record<-train_record[-1]
  test_record<-test_record[-1]
  
  AUCs<-result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_glm,.$true,folds=.$iter))
  AUCs2<-result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_RF,.$true,folds=.$iter))
  
  AUC_df<-rbind(bind_rows(AUCs %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="LR"))),
                bind_rows(AUCs2 %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="RF"))))
  AUC_df$nvar<-rep(nvars_opts,2)
  AUC_df
  
  
  calib_fits=data.frame(nvar=NA,iter=NA,intc=NA,intc_LCI=NA,intc_UCI=NA,slope=NA,slope_LCI=NA,slope_UCI=NA)
  for (nvars in nvars_opts){
    for (each in 1:iter){
      data.temp <- result %>% filter(nvar==nvars & iter==each)
      
      intercept <- glm(true~1,offset=log(pred_glm/(1-pred_glm)),family="binomial",data=data.temp)
      #summary(intercept) #Should have intercept of 0. intercept is calibration intercept
      slope <- glm(true~log(pred_glm/(1-pred_glm)),family="binomial",data=data.temp)
      #summary(slope) #Should have slope of 1. beta coefficient = slope=calibration slope
      
      df=data.frame(nvar=nvars,iter=each,
                    intc=coef(intercept),intc_LCI=confint(intercept)[1],intc_UCI=confint(intercept)[2],
                    slope=coef(slope)[2],slope_LCI=confint(slope)[2,1],slope_UCI=confint(slope)[2,2])
      calib_fits=rbind(calib_fits,df)
      
    }
  }
  calib_fits<-calib_fits[-1,]
  calib <- calib_fits %>% group_by(nvar) %>% summarize(mean(intc),mean(intc_LCI),mean(intc_UCI),
                                                       mean(slope),mean(slope_LCI),mean(slope_UCI))
  names(calib) <- c("nvar","intc","intc_LCI","intc_UCI","slope","slope_LCI","slope_UCI")  #renaming
  
  decilesCC <- result %>% split(.,list(.$iter,.$nvar),drop=TRUE) %>% purrr::map(. %>% arrange(pred_glm)) %>% #now have a df for each iteration of each nvar
    purrr::map(~mutate(.x, decile_glm=ntile(pred_glm,10))) %>% #create predicted glm decile groups; equivalent to: purrr::map(list_resB, ~mutate(.x, decile_glm=ntile(pred_glm,10))); str(temp3)
    bind_rows(.) %>% split(.,f=.$nvar) %>% #a list of df for each nvar which contains all iter for that nvar. "nest" might be better for this
    purrr::map(., . %>% group_by(decile_glm) %>% summarize(mean(true),mean(pred_glm))) #for each decile in each nvar, have an avg true and avg predicted
  
  output<-list(df_imps=df_imps_full,result=result,train_record=train_record,test_record=test_record,AUC_df=AUC_df,decilesCC=decilesCC,calib=calib,iter=iter,nvars_opts=nvars_opts)
  
}


####################
#GEMS cases
################### import GEMS data ####
#starting with full dataset
gems1_orig <- read.csv("/gems1.csv", header=T)
#22567 observations, has type=Case/Control

gems1=gems1_orig %>% select(site,	type, f3_gender,	f3_drh_turgor	,	f3_drh_iv	,	f3_drh_hosp,
                            f4a_relationship,	f4a_dad_live	,					
                            f4a_prim_schl,	f4a_ppl_house	,	f4a_yng_children	,			
                            f4a_slp_rooms,	f4a_floor	,	f4a_house_elec	,  			
                            f4a_house_bike,	f4a_house_phone	,	f4a_house_tele	,    			
                            f4a_house_car	,	f4a_house_cart	,	f4a_house_scoot	,    			
                            f4a_house_fridge	,	f4a_house_agland	,	f4a_house_radio	,   			
                            f4a_house_boat	,	f4a_house_none	,	f4a_fuel_elec	,   			
                            f4a_fuel_biogas	,	f4a_fuel_grass	,	f4a_fuel_propane	,   			
                            f4a_fuel_coal	,	f4a_fuel_dung	,	f4a_fuel_natgas	, 			
                            f4a_fuel_charcoal	,	f4a_fuel_crop	,	f4a_fuel_kero	,   			
                            f4a_fuel_wood	,	f4a_fuel_other	,	f4a_ani_goat	,     			
                            f4a_ani_sheep	,	f4a_ani_dog	,	f4a_ani_cat	,      			
                            f4a_ani_cow	,	f4a_ani_rodents	,	f4a_ani_fowl	,       			
                            f4a_ani_other	,	f4a_ani_no	,	f4a_water_house	,    			
                            f4a_water_covwell	,	f4a_water_yard	,	f4a_water_covpwell	, 			
                            f4a_water_pubtap	,	f4a_water_prospring	,	f4a_water_well	,			
                            f4a_water_unspring	,	f4a_water_pubwell	,	f4a_water_river	,    			
                            f4a_water_pond	,	f4a_water_deepwell	,	f4a_water_rain	,   			
                            f4a_water_shallwell	,	f4a_water_bought	,	f4a_water_othr	,    			
                            f4a_water_bore	,	f4a_ms_water	,	f4a_fetch_water	,    			
                            f4a_trip_day	,	f4a_trip_week	,	f4a_water_avail	,   			
                            f4a_store_water	,	f4a_trt_water	,	f4a_trt_method	,   			
                            f4a_notrt_water	,	f4a_disp_feces	,    					
                            f4a_fac_waste	,	f4a_share_fac	,	f4a_wash_eat	,    			
                            f4a_wash_cook	,	f4a_wash_nurse	,	f4a_wash_def	,      			
                            f4a_wash_animal	,	f4a_wash_child	,	f4a_wash_othr	,      			
                            f4a_wash_use	,	f4a_breastfed	,	f4a_drh_days	,     			
                            f4a_max_stools	,	f4a_drh_blood	,	f4a_drh_vomit	,      			
                            f4a_drh_thirst	,	f4a_drh_lessdrink	,    					
                            f4a_drh_bellypain	,	f4a_drh_restless	,   					
                            f4a_drh_lethrgy	,	f4a_drh_consc	,	f4a_drh_strain	,  			
                            f4a_drh_prolapse	,	f4a_drh_cough	,    					
                            f4a_drh_conv	,	f4a_cur_thirsty	,	f4a_cur_skin	,    			
                            f4a_cur_restless	,	f4a_cur_drymouth	,   					
                            f4a_cur_fastbreath	,	f4a_hometrt_ors	,	f4a_hometrt_maize	,  			
                            f4a_hometrt_milk	,	f4a_hometrt_herb	,	f4a_hometrt_zinc	, 			
                            f4a_hometrt_none	,	f4a_hometrt_othrliq	,	f4a_hometrt_ab	,  			
                            f4a_hometrt_othr1	,	f4a_hometrt_othr2	,	f4a_offr_drink	,    			
                            f4a_seek_outside	,	f4a_seek_pharm	,	f4a_seek_friend	,    			
                            f4a_seek_healer	,	f4a_seek_doc	,	f4a_seek_privdoc	,   			
                            f4a_seek_remdy	,	f4a_seek_other	,	f4b_haz	,  			
                            f4b_muac	,	f4b_temp	,	f4b_resp	,           			
                            f4b_chest_indrw	,	f4b_eyes	,	f4b_mouth	,          			
                            f4b_skin	,	f4b_mental	,	f4b_rectal	,         			
                            f4b_bipedal	,	f4b_abn_hair	,	f4b_under_nutr	,     			
                            f4b_skin_flaky	,	f4b_observe_stool	,	f4b_nature_stool	,   			
                            f4b_recommend	,	f4b_volume	,  					
                            f4b_admit	,	f9_memory_aid	,	wealth_index	,       			
                            wiq	,	base_age	,	
                            f5_exp_drh	, 	f5_exp_dys	, 	f5_exp_cou	, 	f5_exp_fever	, 	
                            f5_diag_typ	, 	f5_diag_mal	, 	f5_diag_pne	,			
                            f5_exp_rectal	, 	f5_exp_convul	, 	f5_exp_arthritis	,			
                            f5_rectal	, 	f5_bipedal	, 	f5_abn_hair	, 	f5_under_nutr	, 	f5_skin_flaky,
                            f5_ms_water	, 	f5_main_cont	, 	f5_treat_water	, 	f5_trt_meth	, 	
                            f5_wash_where	, 	f5_wash_piped	, 	f5_wash_noptap	, 	f5_wash_tap	, 	f5_wash_basin,
                            f5_wash_soap	, 	f5_wash_ash	, 					
                            f5_child_feces	, 	f5_feces_visible	, 	f5_feces_else	, 	f5_house_feces,
                            f5_haz, f4b_outcome, f5_status, agegroup, f4b_height, f5_height,
                            f4b_date, f5_date, f4b_haz_f, f7_date, f7_haz, f7_haz_f, f7_height,
                            f5_diag_othr,f5_child_health,
                            f7_relation	,	f7_dad_live	,		
                            f7_prim_schl	,	f7_ppl_house	,	f7_yng_childrn	,
                            f7_slp_rooms	,	f7_floor	,	f7_house_elec	,  
                            f7_house_bike	,	f7_house_phone	,	f7_house_tele	,    
                            f7_house_car	,	f7_house_cart	,	f7_house_scoot	,    
                            f7_house_fridge	,	f7_house_agland	,	f7_house_radio	,   
                            f7_house_boat	,	f7_house_none	,	f7_fuel_elec	,   
                            f7_fuel_biogas	,	f7_fuel_grass	,	f7_fuel_propane	,   
                            f7_fuel_coal	,	f7_fuel_dung	,	f7_fuel_natgas	, 
                            f7_fuel_charcoal	,	f7_fuel_crop	,	f7_fuel_kero	,   
                            f7_fuel_wood	,	f7_fuel_other	,	f7_ani_goat	,     
                            f7_ani_sheep	,	f7_ani_dog	,	f7_ani_cat	,      
                            f7_ani_cow	,	f7_ani_rodents	,	f7_ani_fowl	,       
                            f7_ani_other	,	f7_ani_no	,	f7_water_house	,    
                            f7_water_covwell	,	f7_water_yard	,	f7_water_covpwell	, 
                            f7_water_pubtap	,	f7_water_prospring	,	f7_water_well	,
                            f7_water_unspring	,	f7_water_pubwell	,	f7_water_river	,    
                            f7_water_pond	,	f7_water_deepwell	,	f7_water_rain	,   
                            f7_water_shallwell	,	f7_water_bought	,	f7_water_othr	,    
                            f7_water_bore	,	f7_ms_water	,	f7_water_avail	,   
                            f7_store_water	,	f7_trt_water	,	f7_trt_method	,   
                            f7_disp_feces	,    				
                            f7_fac_waste	,	f7_share_fac	,	f7_wash_eat	,    
                            f7_wash_cook	,	f7_wash_nurse	,	f7_wash_def	,      
                            f7_wash_animal	,	f7_wash_child	,	f7_wash_othr	,      
                            f7_wash_use	,	f7_breastfed	,   		
                            f7_seekcare	,				
                            f7_height	,	f7_muac	,	f7_haz	,  
                            f7_temp	,	f7_resp	,   		
                            f7_bipedal	,	f7_abn_hair	,	f7_under_nutr	,     
                            f7_skin_flaky, caseid,
                            
                            f7_med_cotr,f7_med_gent,f7_med_chlor,
                            f7_med_eryth,f7_med_azith,f7_med_omacr,
                            f7_med_peni,f7_med_amoxy,f7_med_ampi,
                            f7_med_nalid,f7_med_cipro,f7_med_sele,
                            f7_med_otherant,
                            f11_antibiotic,
                            f11_anti_ampi,f11_anti_nali,f11_anti_cotr,
                            f11_anti_cipr,f11_anti_sele,f11_anti_gent,
                            f11_anti_chlo,f11_anti_eryt,f11_anti_azit,
                            f11_anti_macr,f11_anti_peni,f11_anti_amox,
                            f11_anti_other,f4a_hometrt_ab,f4a_seek_remdy,f4a_seek_remdy_spec,
                            f4b_trt_give_cxl,f4b_trt_give_gent,f4b_trt_give_chlor,
                            f4b_trt_give_ery,f4b_trt_give_azi,f4b_trt_give_macr,
                            f4b_trt_give_pen,f4b_trt_give_amox,f4b_trt_give_ampi,
                            f4b_trt_give_nalid,f4b_trt_give_cpnr,f4b_trt_give_slpy,
                            f4b_trt_give_othr,f4b_trt_pres_cxl,f4b_trt_pres_gent,
                            f4b_trt_pres_chlor,
                            f4b_trt_pres_ery,f4b_trt_pres_azi,f4b_trt_pres_macr,
                            f4b_trt_pres_pen,f4b_trt_pres_amox,f4b_trt_pres_ampi,
                            f4b_trt_pres_nalid,f4b_trt_pres_cpnr,f4b_trt_pres_slpy,
                            f4b_trt_pres_othr
)


#variables w/ small cell sizes, dropping for now: "f4a_chlorine","f4a_primcare","f4a_mom_live",
# "f4a_drh_undrink","f4a_drh_fever","f4a_drh_breath","f4b_skin_pinch",
# "f5_diag_meng"

gems1=gems1 %>% mutate(any_breast_fed=factor(case_when((f4a_breastfed==1|f4a_breastfed==2)~1,TRUE~0))) #SMA dichotomizing breastfeeding
gems1=gems1 %>% mutate(any_breast_fed2=factor(case_when((f4a_breastfed==0|f4a_breastfed==1)~0,TRUE~1)))
gems1=gems1 %>% mutate(cont=case_when(site %in% c(1,2,3,4) ~ 1,
                                      TRUE ~ 2)) #SMA creating "cont"inent variable
gems1$site=as.factor(gems1$site)
gems1$index=1:dim(gems1)[1] #SMA creating an ID for each observation

cases_orig <- gems1 %>% filter(type=="Case")
#9439

################### define variables ####
cases_orig=cases_orig %>% mutate(haz_dif = f4b_haz - f5_haz) #dat_joined$haz_dif <- dat_joined$f4b_haz - dat_joined$f5_haz
# summary(cases$haz_dif)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -9.1100 -0.0800  0.1800  0.2027  0.4800  8.3500    1051 
cases_orig=cases_orig %>% mutate(haz_1.0=(case_when(haz_dif>=1.0 ~ 1, TRUE~0)))
# cases=cases %>% mutate(haz_1.0=(case_when(haz_dif<=-1.0 ~ 1, TRUE~0))) #dat_joined$haz_0.1 <- ifelse(dat_joined$haz_dif<=-1.0,1,0)
# # table(cases$haz_1.0)
# # 0    1 
# # 9270  169 
cases_orig=cases_orig %>% mutate(haz_0.5=(case_when(haz_dif>=0.5 ~ 1, TRUE~0)))
# cases=cases %>% mutate(haz_0.5=(case_when(haz_dif<=-0.5 ~ 1, TRUE~0))) #dat_joined$haz_0.5 <- ifelse(dat_joined$haz_dif<=-0.5,1,0)
# # table(cases$haz_0.5)
# # 0    1 
# # 8907  532 
#see pg 51 of notebook for manual checks

cases_orig=cases_orig %>% mutate(change_ht = f5_height - f4b_height)
# summary(cases$change_ht)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -25.033   0.800   1.633   1.772   2.500  41.833    1051 

# #f4b_outcome: 1=resolved; 2=improved; 3=no better; 4=worse; 5=died at hosp; 6=unknown/LTF
# #f5_child_health: 1=appears healthy; 2=improved but not back to normal; 3=no better; 4=worse; 5=died
#death_all: any death after hospital admit
cases_orig=cases_orig %>% mutate(death_all=(case_when(
                                    f4b_outcome==5 ~ 1,
                                    f5_child_health==5 ~ 1, 
                                    TRUE~0)),
                                 death_hosp=(case_when(
                                    f4b_outcome==5 ~ 1,
                                    TRUE~0)),
                                 death_home=(case_when(
                                    (f4b_outcome!=5 & f5_child_health==5) ~1,
                                    TRUE~0)))
# table(cases_orig$death_all)
# 0    1 
# 9249  190
# table(cases_orig$death_hosp)
# 0    1 
# 9390   49 
# table(cases_orig$death_home)
# 0    1 
# 9298  141 

cases_orig$f4b_date_date <- as.Date(as.character(cases_orig$f4b_date))
cases_orig$f5_date_date <- as.Date(as.character(cases_orig$f5_date))
cases_orig$fup_days <- as.numeric(cases_orig$f5_date_date - cases_orig$f4b_date_date)

#create variable where NA missing are set to value 9
cases_orig=cases_orig %>% mutate(f4b_outcome_miss=replace_na(f4b_outcome,9),
                                 f5_child_health_miss=replace_na(f5_child_health,9))
#only a couple missing and one 9, set those to no(0)
cases_orig$f4a_drh_lethrgy <- as.numeric(cases_orig$f4a_drh_lethrgy)
cases_orig <- cases_orig %>% mutate(f4a_drh_lethrgy_miss=case_when((f4a_drh_lethrgy==9) ~ as.numeric(NA), TRUE~f4a_drh_lethrgy))
cases_orig$f4a_drh_restless <- as.numeric(cases_orig$f4a_drh_restless)
cases_orig <- cases_orig %>% mutate(f4a_drh_restless=case_when((f4a_drh_restless==9) ~ as.numeric(NA), TRUE~f4a_drh_restless))

#too few in 14(dam or earth pan), set to other(18)
#combine 9 and 10 (both are covered well)
cases_orig$f4a_ms_water <- as.numeric(cases_orig$f4a_ms_water)
cases_orig <- cases_orig %>% mutate(f4a_ms_water=case_when((f4a_ms_water==14)~18,(f4a_ms_water==9)~10,TRUE~f4a_ms_water))

#assuming f5_child_health==0 are also missing since aren't in codebook, so set to 9
#this doesn't work since updating R, assuming is a tidyverse version issue. go base
# test=cases_orig %>% mutate(f5_child_health_miss=case_when(f5_child_health==0 ~ 9, 
#                                                           TRUE~f5_child_health_miss))
# test<-cases_orig
# test$f5_child_health_miss <- ifelse(test$f5_child_health_miss==0,9,test$f5_child_health_miss)
cases_orig$f5_child_health_miss <- ifelse(cases_orig$f5_child_health_miss==0,9,cases_orig$f5_child_health_miss)
# table(cases$f4b_outcome_miss)
# # 1    2    3    4    5    6    9 
# # 614 4179 4187   11   49  398    1
# table(is.na(cases$f4b_outcome))
# # FALSE  TRUE 
# # 9438     1
# table(cases$f5_child_health_miss)
# # 1    2    3    4    5    9
# # 6266  834  839  470  186  844
# table(is.na(cases$f5_child_health))
# # FALSE  TRUE 
# # 8604   835 

#very few 5 so put in other category
cases_orig$f4a_disp_feces <- as.numeric(cases_orig$f4a_disp_feces)
cases_orig=cases_orig %>% mutate(f4a_disp_feces=(case_when(f4a_disp_feces==5 ~ 6, TRUE~f4a_disp_feces)))

#combining categories. not sure what 0 is, isn't in dictionary
table(cases_orig$f4a_trt_method)
#    0    1    2    3    4    5    6    7 
# 6740    3  656  751 1197   73   14    5
cases_orig$f4a_trt_method <- as.numeric(cases_orig$f4a_trt_method)
cases_orig=cases_orig %>% mutate(f4a_trt_method=(case_when(f4a_trt_method==1 | f4a_trt_method==6 | f4a_trt_method==7 ~ 7, TRUE~f4a_trt_method)))
table(cases_orig$test)
# 0    2    3    4    5    7 
# 6740  656  751 1197   73   22 

#creating antibiotic variables - want one var for each of: any antibiotic at any time, b/f seek care, when present care
#stool collection antibiotic variables (can be from cases or controls)
abx_stool <- c("f11_antibiotic",
               "f11_anti_ampi","f11_anti_nali","f11_anti_cotr",
               "f11_anti_cipr","f11_anti_sele","f11_anti_gent",
               "f11_anti_chlo","f11_anti_eryt","f11_anti_azit",
               "f11_anti_macr","f11_anti_peni","f11_anti_amox",
               "f11_anti_other")
#cases b/f coming health center (non-medical questionnaire)
abx_cases_pre <- c("f4a_hometrt_ab","f4a_seek_remdy","f4a_seek_remdy_spec")
#cases given at health center
abx_cases_at <- c("f4b_trt_give_cxl","f4b_trt_give_gent","f4b_trt_give_chlor",
                  "f4b_trt_give_ery","f4b_trt_give_azi","f4b_trt_give_macr",
                  "f4b_trt_give_pen","f4b_trt_give_amox","f4b_trt_give_ampi",
                  "f4b_trt_give_nalid","f4b_trt_give_cpnr","f4b_trt_give_slpy",
                  "f4b_trt_give_othr")
#cases given at health center for tx at home
abx_cases_home <- c("f4b_trt_pres_cxl","f4b_trt_pres_gent","f4b_trt_pres_chlor",
                    "f4b_trt_pres_ery","f4b_trt_pres_azi","f4b_trt_pres_macr",
                    "f4b_trt_pres_pen","f4b_trt_pres_amox","f4b_trt_pres_ampi",
                    "f4b_trt_pres_nalid","f4b_trt_pres_cpnr","f4b_trt_pres_slpy",
                    "f4b_trt_pres_othr")


# test<-cases[abx_stool]
# test$sum<-rowSums(test[,abx_stool])
# table(test$sum,test$f11_antibiotic)
# 0    1    9
# 0 9179    0    0
# 2    0  131    0
# 3    0  119    0
# 4    0    4    0
# 9    0    0    6
#>>> f11_antibiotic==1 is an accurate indicator variable of if kid received antibiotics after arrive health center but before stool sample collected

#f4a_hometrt_ab==1 b/f coming to health  clinic, child was treated w/ antitiobitcs for diarrhea

# # number of antibiotics each child given at health center
# table(rowSums(cases[,abx_cases_at]==1))
# # 0    1    2    3    4    5 
# # 7388 1322  597  103   21    8
# 
# #number of antibiotics each child prescribed for tx at home
# table(rowSums(cases[,abx_cases_home]==1))
# # 0    1    2    3    4 
# # 2765 5040 1613   19    1
# #one child is missing, but won't change results. come back to if important

cases_orig$abx_bf <- cases_orig$f4a_hometrt_ab #if received abx before present for care
cases_orig$abx_at <- ifelse(rowSums(cases_orig[,abx_cases_at]==1)>0,1,0) #if received abx during care
cases_orig$abx_home <- ifelse(rowSums(cases_orig[,abx_cases_home]==1)>0,1,0) #if received abx rx to take home after care
#f4b_trt_pres_ampi missing for one person, assuming just a data entry error, set to 0
cases_orig$abx_home <- ifelse(is.na(cases_orig$abx_home),0,cases_orig$abx_home)
cases_orig$abx_ever <- ifelse(cases_orig$abx_bf==1 | cases_orig$abx_at==1 | cases_orig$abx_home==1,1,0)  

cases_orig <- cases_orig %>% mutate(month=month(f4b_date_date))

table(cases_orig$f4a_floor)
#    1    2    3    4    6    7    8    9   10 
# 2490  746    4    4   18  262 5813   89   13 
cases_orig=cases_orig %>% mutate(f4a_floor=(case_when((f4a_floor==1 | f4a_floor==2 |
                                                         f4a_floor==3 | f4a_floor==4 | f4a_floor==10) ~ 0, #natural, rudimentary, other floor
                                                      TRUE~1))) #finished floor
# table(cases_orig$test)
# # 0    1 
# # 3257 6182 

#combine to fewer categories: f4a_ms_water
table(cases_orig$f4a_ms_water)
#   1    2    3    4    5    6    7    8   10   11   12   13   15   16   17   18 
# 668  933 3481   72  272  109  972  703  127   58   56  374  529  658  336   91 
cases_orig=cases_orig %>% mutate(f4a_ms_water=(case_when((f4a_ms_water==6 | f4a_ms_water==13 | f4a_ms_water==14) ~ 0, #surface 
                                                 (f4a_ms_water==4 | f4a_ms_water==5 | f4a_ms_water==12 | f4a_ms_water==16) ~ 1, #unimproved
                                                 (f4a_ms_water==3 | f4a_ms_water==7 | f4a_ms_water==8 | f4a_ms_water==9 |
                                                    f4a_ms_water==10 | f4a_ms_water==11 | f4a_ms_water==15 | f4a_ms_water==17) ~ 2, #other improved
                                                 (f4a_ms_water==1 | f4a_ms_water==2) ~ 3, #piped
                                                 TRUE~4))) #other
# table(cases_orig$test)
# 0    1    2    3    4 
# 483 1058 6206 1601   91 
# #use JMP drinking water services ladder
# #surface (subset of unimproved)
# 6-pond/lake
# 13-river/stream
# 14-dam/earth pan
# #unimproved
# 4-open well in house/yard
# 5-open public well
# 12-unprotected spring
# 16-bought
# #other improved
# 3-public tap
# 7-deep tube well
# 8-shallow tube well
# 9-covered well in house/yard
# 10-covered public well
# 11-protected spring
# 15-rainwater
# 17-bore hole
# #safely managed/ piped into dwelling/plot/yard
# 1-piped into house
# 2-piped into yard
# #other
# 18-other

summary(cases_orig$f4a_relationship)
#    1    2    3    4    5    6    7    8    9   10 
# 8941   81   49   12  206    6  125    6    2   11 
#creating a combined category for non-father male relation OR non relation (4,6,8,9)
cases_orig$f4a_relationship <- as.numeric(cases_orig$f4a_relationship)
cases_orig=cases_orig %>% mutate(f4a_relationship=(case_when((f4a_relationship==4 | f4a_relationship==6 | f4a_relationship==8 | f4a_relationship==9) ~ 9, #non-father male relation OR 
                                                         TRUE~f4a_relationship))) #what was originally


str(cases_orig$f4a_drh_cough)
table(cases_orig$f4a_drh_cough)
cases_orig$f4a_drh_cough <- ifelse(cases_orig$f4a_drh_cough==9,0,cases_orig$f4a_drh_cough)

#convert these to factors
vars <- c("f4a_ms_water","f4a_fac_waste","f4a_dad_live","f4b_recommend",
          "f4a_relationship","f4a_prim_schl","f4a_floor","f4a_disp_feces",
          "f4a_wash_use","f4a_water_avail","f4a_trt_method","f4a_drh_blood",
          "f4a_drh_vomit","f4a_drh_thirst","f4a_drh_lessdrink","f4a_drh_bellypain",
          "f4a_drh_restless","f4a_drh_lethrgy_miss","f4a_drh_consc","f4a_drh_strain",
          "f4a_drh_prolapse","f4a_drh_cough","f4a_drh_conv","f4a_cur_thirsty",
          "f4a_cur_skin","f4a_cur_restless","f4a_cur_drymouth","f4a_cur_fastbreath",
          "f4b_nature_stool","month"
)
cases_orig[vars] <- lapply(cases_orig[vars], factor)
#this works too: test <- cases_orig %>% mutate_at(vars, list(~factor(.)))



#collected as categorical, but ordinal so leaving as numeric for now: 
#f4a_prim_schl, f4a_offr_drink, f4a_max_stools, f4a_breastfed, f4b_mouth, f4b_skin, f4b_mental

################### inclusion/exclusion ####
#cases_orig n=9439

#limit to fup time period criteria
# #table(cases_orig$fup_days,cases_orig$death_all)
# #do any of these outlier fup dates include kids who died in hospital?
# temp <- cases_orig %>% filter(f4b_outcome_miss==5)
# table(temp$fup_days,temp$death_all)
# #>>>yes. kids who died in hospital per f4b_outcome_miss had weird fup dates so don't exclude deaths based on that. 
# #but makes sense to exclude survivors based on those cutoffs since want consistent time at risk
# table(is.na(cases_orig$f5_status),cases_orig$death_all)
# # 0    1
# # FALSE 9241  190
# # TRUE     8    0
# table(is.na(cases_orig$f5_status),cases_orig$f5_child_health_miss)
# # 1    2    3    4    5    9
# # FALSE 6266  834  839  470  186  836
# # TRUE     0    0    0    0    0    8
# #>>> are 8 kids who not reported dead, but missing fup info. should drop these too
# #>>> also, f5_status redundant for inclusion/exclusion since consider f5_child_health_miss, and they're in full agreement on who didn't have fup
#if not dead and fup days outside bounds, don't want
#if dead, any fup days acceptable
#f5_status ==1 60d f-up conducted, 0==not conducted
  

# #f4b_outcome: 1=resolved; 2=improved; 3=no better; 4=worse; 5=died at hosp; 6=unknown/LTF
#f5_status: 1=60d f-up conducted; 0=not conducted)
# #f5_child_health_miss: 1=appears healthy; 2=improved but not back to normal; 3=no better; 4=worse; 5=died; 9=missing
#death_all: 0=didn't die; 1=did die; any death after hospital admit
      # cases_orig=cases_orig %>% mutate(death_all=(case_when(
      #   f4b_outcome==5 ~ 1,
      #   f5_child_health==5 ~ 1, 
      #   TRUE~0)),
      #   death_hosp=(case_when(
      #     f4b_outcome==5 ~ 1,
      #     TRUE~0)),
      #   death_home=(case_when(
      #     (f4b_outcome!=5 & f5_child_health==5) ~1,
      #     TRUE~0)))

#in co-author review, decide this fup_days criteria doesn't work for all observations, don't want to apply to all so incorporate fup_days inclusion criteria to each individual scenario below
#cases_orig$drop <- ifelse(is.na(cases_orig$f5_status) | (cases_orig$death_all==0 & (cases_orig$fup_days<49 | cases_orig$fup_days>91)),1,0)

# #cases at hospital discharge:
# table(cases_orig$f4b_outcome_miss,cases_orig$f5_child_health_miss)
# #      1    2    3    4    5    9
# # 1  420   16  107   44    8   19
# # 2 2869  526  229  171   84  300
# # 3 2684  247  503  254   38  461
# # 4    4    2    0    0    2    3
# # 5    0    0    0    0   45    4
# # 6  289   43    0    1    9   56
# # 9    0    0    0    0    0    1

#@ hosp discharge and @ 60d f-up

# #drop from all analyses if:
# #missing missing
# #dead alive
# #alive missing 
# miss.miss <- cases_orig %>% filter((f4b_outcome_miss==6 | f4b_outcome_miss==9) & f5_child_health_miss==9) #n=1+56=57, what expect
# dead.alive <- cases_orig %>% filter(f4b_outcome_miss==5 & (f5_child_health_miss!=5 & f5_child_health_miss!=9)) #n=0+0+0+0=0, what expect
# alive.miss <- cases_orig %>% filter((f4b_outcome_miss!=5 & f4b_outcome_miss!=6 & f4b_outcome_miss!=9) & f5_child_health_miss==9) #n=19+300+461+3=783, what expect

#for dead_all (died after hosp admit)
#missing dead
#missing alive
#dead missing
#dead dead
#alive dead
#alive alive
#need to be documented alive within timeframe if not dead at hospital:
miss.dead <- cases_orig %>% filter((f4b_outcome_miss==6 | f4b_outcome_miss==9) & f5_child_health_miss==5 &
                                     (cases_orig$fup_days>=49 & cases_orig$fup_days<=91)) #n=0+9-1-1=7, what expect
  # miss.dead <- cases_orig %>% filter((f4b_outcome_miss==6 | f4b_outcome_miss==9) & f5_child_health_miss==5) #n=0+9=9, what expect
  # table(miss.dead$fup_days)
  # 56  61  62  65  66  69  92 208 
  # 1   1   2   1   1   1   1   1 
miss.alive <- cases_orig %>% filter((f4b_outcome_miss==6 | f4b_outcome_miss==9) & (f5_child_health_miss!=5 & f5_child_health_miss!=9) & 
                                      (cases_orig$fup_days>=49 & cases_orig$fup_days<=91)) #n=289+43+1-1-3-1=328, what expect
  # miss.alive <- cases_orig %>% filter((f4b_outcome_miss==6 | f4b_outcome_miss==9) & (f5_child_health_miss!=5 & f5_child_health_miss!=9)) #n=289+43+1=333, what expect
  # table(miss.alive$fup_days)
  # 47 51 52 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 86 87 88 89 90 91 92 95 
  # 1  1  1  2  5  4  7 19  8  7 13  7 18 24 12 10 10 18 17 10  9 12  9 12  6  7 15 10  7  6  2  4  5  5  1  4  5  3  8  5  3  1
dead.miss <- cases_orig %>% filter(f4b_outcome_miss==5 & f5_child_health_miss==9 &
                                     (cases_orig$fup_days>=49 & cases_orig$fup_days<=91)) #n=4, what expect
  # dead.miss <- cases_orig %>% filter(f4b_outcome_miss==5 & f5_child_health_miss==9) #n=4, what expect
  # table(dead.miss$fup_days)
  # 64 83 87 90 
  # 1  1  1  1
dead.dead <- cases_orig %>% filter(f4b_outcome_miss==5 & f5_child_health_miss==5) #n=45, what expect
alive.dead <- cases_orig %>% filter((f4b_outcome_miss!=5 & f4b_outcome_miss!=6 & f4b_outcome_miss!=9) & f5_child_health_miss==5 &
                                      (cases_orig$fup_days>=49 & cases_orig$fup_days<=91)) #n=8+84+38+2-9=123, what expect
  # alive.dead <- cases_orig %>% filter((f4b_outcome_miss!=5 & f4b_outcome_miss!=6 & f4b_outcome_miss!=9) & f5_child_health_miss==5) #n=8+84+38+2=132, what expect
  # table(alive.dead$fup_days)
  # 49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69 
  # 5   4   5   1   3   6   2   6   7   6  12   8  10   4   8   5   2   1   3   3   3 
  # 70  72  73  74  75  76  77  78  80  83  84  86  88  90  91  92  93 134 136 138 144 
  # 1   1   1   1   2   2   1   2   1   1   1   1   1   2   1   1   1   1   1   1   1 
  # 164 187 193 
  # 1   1   1 
alive.alive <- cases_orig %>% filter((f4b_outcome_miss!=5 & f4b_outcome_miss!=6 & f4b_outcome_miss!=9) & (f5_child_health_miss!=5 & f5_child_health_miss!=9) &
                                       (cases_orig$fup_days>=49 & cases_orig$fup_days<=91)) #n=420+2869+2684+4+16+526+247+2+107+229+503+0+44+171+254+0-1-1-1-3-6-14-3-2-1-7-3-4-1-1-3-3-1-1-1-2-1-1-1-1=8013, what expect
  # alive.alive <- cases_orig %>% filter((f4b_outcome_miss!=5 & f4b_outcome_miss!=6 & f4b_outcome_miss!=9) & (f5_child_health_miss!=5 & f5_child_health_miss!=9)) #n=420+2869+2684+4+16+526+247+2+107+229+503+0+44+171+254+0=8076, what expect
  # table(alive.alive$fup_days)
  # 34  40  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66 
  # 1   1   1   3   6 219 207  98  80  75  88 102 106 219 181 748 711 683 537 449 315 256 277 
  # 67  68  69  70  71  72  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89 
  # 216 189 183 162 147 122 105  91  88  77  84  85  83  81  87  95 114 123 127 105 121  61  58 
  # 90  91  92  93  94  95  96  97  98  99 100 101 102 105 107 109 110 111 113 118 128 
  # 50   8  14   3   2   1   7   3   4   1   1   3   3   1   1   1   2   1   1   1   1 
cases_anydeath <- rbind(miss.dead,miss.alive,dead.miss,dead.dead,alive.dead,alive.alive) #n=7+328+4+45+123+8013=8520, what expect
#n=8520
# table(cases_anydeath$death_all)
# # 0    1 
# # 8341  179
# table(cases_anydeath$f5_status,cases_anydeath$death_all) # f5_status ==1 (60d f-up conducted, 0==not conducted)
# # 0    1
# # 0    2    5
# # 1 8339  174
# table(cases_anydeath$f4b_outcome_miss,cases_anydeath$death_all)
# 0    1
# 1  583    7
# 2 3744   77
# 3 3680   38
# 4    6    1
# 5    0   49
# 6  328    7
# table(cases_anydeath$f5_child_health_miss,cases_anydeath$death_all)
# 0    1
# 1 6219    0
# 2  816    0
# 3  837    0
# 4  469    0
# 5    0  175
# 9    0    4

#choosing not to drop observations for HAZ plausibility purposes since those who already dead not eligible to be excluded for such reasons


#for dead_hosp (died after hosp admit, before discharge)
cases_hospdeath <- rbind(miss.dead,miss.alive,dead.miss,dead.dead,alive.dead,alive.alive) #n=7+328+4+45+123+8013=8520, what expect
#n=8520
# cases_hospdeath <- rbind(miss.alive,dead.miss,dead.dead,alive.dead,alive.alive) #n=328+4+45+123+8013=8513, what expect
# #>>> anydeath includes miss.dead, hospdeath does not include them (7 ppl)


# #for dead_home (died after hosp discharge, ie @60d f-up)
# #missing alive
# #alive dead
# #alive alive
# miss.alive <- cases %>% filter((f4b_outcome_miss==6 | f4b_outcome_miss==9) & (f5_child_health_miss!=5 & f5_child_health_miss!=9)) #n=289+43+1=333, what expect
# alive.dead <- cases %>% filter((f4b_outcome_miss!=5 & f4b_outcome_miss!=6 & f4b_outcome_miss!=9) & f5_child_health_miss==5) #n=8+84+38+2=132, what expect
# alive.alive <- cases %>% filter((f4b_outcome_miss!=5 & f4b_outcome_miss!=6 & f4b_outcome_miss!=9) & (f5_child_health_miss!=5 & f5_child_health_miss!=9)) #n=420+2869+2684+4+16+526+247+2+107+229+503+44+171+254+0+0=8076, what expect
cases_homedeath <- rbind(miss.dead,miss.alive,alive.dead,alive.alive) #n=7+328+123+8013=8471, what expect
#n=8471
# cases_homedeath <- rbind(miss.alive,alive.dead,alive.alive) #n=328+123+8013=8464, what expect
# #>>> anydeath includes miss.dead, dead.miss, dead.dead, homedeath does not include them (7+4+45=56 ppl)
# #>>> hospdeath includes dead.miss, dead.dead, homedeath does not include them (4+45=49 ppl)


################### drop missing since can't have missing in RF, define "names" variables interested in ####
#DEATH
complete_anydeath <- cases_anydeath %>% filter(!is.na(f4a_dad_live)&!is.na(f4a_prim_schl)&!is.na(f4a_slp_rooms)&!is.na(f4a_water_avail)&
                                           !is.na(f4a_disp_feces)&!is.na(f4a_share_fac)&!is.na(f4a_wash_use)&!is.na(f4a_drh_days)&
                                           !is.na(f4a_drh_thirst)&!is.na(f4a_drh_restless)&!is.na(f4a_drh_lethrgy_miss)&!is.na(f4a_drh_conv)&
                                           !is.na(f4a_cur_skin)&!is.na(f4a_cur_restless)&!is.na(f4a_cur_drymouth)&!is.na(f4a_cur_fastbreath)&
                                           !is.na(f4a_offr_drink)&!is.na(f4b_temp)&!is.na(f4b_chest_indrw)&!is.na(f4b_mouth)&
                                           !is.na(f4b_skin)&!is.na(f4b_under_nutr)&!is.na(f4b_nature_stool)&!is.na(f4a_ppl_house)&
                                           !is.na(f4b_resp)&!is.na(f4b_abn_hair)&!is.na(f3_drh_iv)&!is.na(f4a_store_water)&
                                           !is.na(f4b_mental)&(f4a_drh_vomit!="9")&
                                           !is.na(f4b_haz)&!is.na(f4a_breastfed)&!is.na(f4a_drh_consc)&(f4a_drh_consc!="9")
)
#8520 to 8060 observations

complete_hospdeath <- cases_hospdeath %>% filter(!is.na(f4a_dad_live)&!is.na(f4a_prim_schl)&!is.na(f4a_slp_rooms)&!is.na(f4a_water_avail)&
                                                 !is.na(f4a_disp_feces)&!is.na(f4a_share_fac)&!is.na(f4a_wash_use)&!is.na(f4a_drh_days)&
                                                 !is.na(f4a_drh_thirst)&!is.na(f4a_drh_restless)&!is.na(f4a_drh_lethrgy_miss)&!is.na(f4a_drh_conv)&
                                                 !is.na(f4a_cur_skin)&!is.na(f4a_cur_restless)&!is.na(f4a_cur_drymouth)&!is.na(f4a_cur_fastbreath)&
                                                 !is.na(f4a_offr_drink)&!is.na(f4b_temp)&!is.na(f4b_chest_indrw)&!is.na(f4b_mouth)&
                                                 !is.na(f4b_skin)&!is.na(f4b_under_nutr)&!is.na(f4b_nature_stool)&!is.na(f4a_ppl_house)&
                                                 !is.na(f4b_resp)&!is.na(f4b_abn_hair)&!is.na(f3_drh_iv)&!is.na(f4a_store_water)&
                                                 !is.na(f4b_mental)&(f4a_drh_vomit!="9")&
                                                 !is.na(f4b_haz)&!is.na(f4a_breastfed)&!is.na(f4a_drh_consc)&(f4a_drh_consc!="9")
)
#8520 to 8060 observations

complete_homedeath <- cases_homedeath %>% filter(!is.na(f4a_dad_live)&!is.na(f4a_prim_schl)&!is.na(f4a_slp_rooms)&!is.na(f4a_water_avail)&
                                                 !is.na(f4a_disp_feces)&!is.na(f4a_share_fac)&!is.na(f4a_wash_use)&!is.na(f4a_drh_days)&
                                                 !is.na(f4a_drh_thirst)&!is.na(f4a_drh_restless)&!is.na(f4a_drh_lethrgy_miss)&!is.na(f4a_drh_conv)&
                                                 !is.na(f4a_cur_skin)&!is.na(f4a_cur_restless)&!is.na(f4a_cur_drymouth)&!is.na(f4a_cur_fastbreath)&
                                                 !is.na(f4a_offr_drink)&!is.na(f4b_temp)&!is.na(f4b_chest_indrw)&!is.na(f4b_mouth)&
                                                 !is.na(f4b_skin)&!is.na(f4b_under_nutr)&!is.na(f4b_nature_stool)&!is.na(f4a_ppl_house)&
                                                 !is.na(f4b_resp)&!is.na(f4b_abn_hair)&!is.na(f3_drh_iv)&!is.na(f4a_store_water)&
                                                 !is.na(f4b_mental)&(f4a_drh_vomit!="9")&
                                                 !is.na(f4b_haz)&!is.na(f4a_breastfed)&!is.na(f4a_drh_consc)&(f4a_drh_consc!="9")
)
#8471 to 8017 observations
# table(complete_anydeath$death_all)
# 0    1 
# 7895  165 
# table(complete_anydeath$death_hosp)
# 0    1 
# 8017   43 
# table(complete_anydeath$death_home)
# 0    1 
# 7938  122 
# 43+122
# [1] 165
# table(complete_hospdeath$death_hosp)
# 0    1 
# 8017   43 
# table(complete_hospdeath$death_home)
# 0    1 
# 7938  122 
# table(complete_homedeath$death_home)
# 0    1 
# 7895  122 
# 8017+43
# [1] 8060
#>>>165 deaths total, 43 in hosp, 122 after discharge


# select variables we're interested in. have checked all appropriate ones are factor
names <- c("site","f3_gender","f3_drh_turgor","f3_drh_iv","f3_drh_hosp",
           "f4a_relationship","f4a_dad_live",
           "f4a_prim_schl","f4a_ppl_house","f4a_yng_children",
           "f4a_slp_rooms","f4a_floor","f4a_house_elec",  
           "f4a_house_bike","f4a_house_phone","f4a_house_tele",    
           "f4a_house_car","f4a_house_cart","f4a_house_scoot",    
           "f4a_house_fridge","f4a_house_agland","f4a_house_radio",   
           "f4a_house_boat","f4a_house_none","f4a_fuel_elec",   
           "f4a_fuel_biogas","f4a_fuel_grass","f4a_fuel_propane",   
           "f4a_fuel_coal","f4a_fuel_dung","f4a_fuel_natgas", 
           "f4a_fuel_charcoal","f4a_fuel_crop","f4a_fuel_kero",   
           "f4a_fuel_wood","f4a_fuel_other","f4a_ani_goat",     
           "f4a_ani_sheep","f4a_ani_dog","f4a_ani_cat",      
           "f4a_ani_cow","f4a_ani_rodents","f4a_ani_fowl",       
           "f4a_ani_other","f4a_ani_no","f4a_water_house",    
           "f4a_water_covwell","f4a_water_yard","f4a_water_covpwell", 
           "f4a_water_pubtap","f4a_water_prospring","f4a_water_well",
           "f4a_water_unspring","f4a_water_pubwell","f4a_water_river",    
           "f4a_water_pond","f4a_water_deepwell","f4a_water_rain",   
           "f4a_water_shallwell","f4a_water_bought","f4a_water_othr",    
           "f4a_water_bore","f4a_ms_water","f4a_water_avail",   
           "f4a_store_water","f4a_trt_water","f4a_trt_method",   
           "f4a_disp_feces",    
           "f4a_fac_waste","f4a_share_fac","f4a_wash_eat",    
           "f4a_wash_cook","f4a_wash_nurse","f4a_wash_def",      
           "f4a_wash_animal","f4a_wash_child","f4a_wash_othr",      
           "f4a_wash_use","f4a_breastfed","f4a_drh_days",     
           "f4a_max_stools","f4a_drh_blood","f4a_drh_vomit",      
           "f4a_drh_thirst","f4a_drh_lessdrink",    
           "f4a_drh_bellypain","f4a_drh_restless",   
           "f4a_drh_lethrgy_miss","f4a_drh_consc","f4a_drh_strain",  
           "f4a_drh_prolapse","f4a_drh_cough",    
           "f4a_drh_conv","f4a_cur_thirsty","f4a_cur_skin",    
           "f4a_cur_restless","f4a_cur_drymouth",   
           "f4a_cur_fastbreath","f4a_hometrt_ors","f4a_hometrt_maize",  
           "f4a_hometrt_milk","f4a_hometrt_herb","f4a_hometrt_zinc", 
           "f4a_hometrt_none","f4a_hometrt_othrliq","f4a_hometrt_ab",  
           "f4a_hometrt_othr1","f4a_hometrt_othr2","f4a_offr_drink",    
           "f4a_seek_outside","f4a_seek_pharm","f4a_seek_friend",    
           "f4a_seek_healer","f4a_seek_doc","f4a_seek_privdoc",   
           "f4a_seek_remdy","f4a_seek_other","f4b_haz",  
           #"f4b_muac",
           "f4b_temp","f4b_resp",           
           "f4b_chest_indrw","f4b_eyes","f4b_mouth",          
           "f4b_skin","f4b_mental","f4b_rectal",         
           "f4b_bipedal","f4b_abn_hair","f4b_under_nutr",     
           "f4b_skin_flaky",
           # "f4b_observe_stool","f4b_nature_stool",   
           "f4b_recommend",  
           "f4b_admit", 
           #"wealth_index", 
           "base_age"
           #f5_ variables are at f-up, don't have to inform prediction at hosp admit; haven't checked if any of these f5_ should be factors
           # 'f5_exp_drh', 'f5_exp_dys', 'f5_exp_cou', 'f5_exp_fever', 
           # 'f5_diag_typ', 'f5_diag_mal', 'f5_diag_pne',
           # 'f5_exp_rectal', 'f5_exp_convul', 'f5_exp_arthritis',
           # 'f5_rectal', 'f5_bipedal', 'f5_abn_hair', 'f5_under_nutr', 'f5_skin_flaky',
           # 'f5_ms_water_factor', 'f5_main_cont', 'f5_treat_water', 'f5_trt_meth', 
           # 'f5_wash_where', 'f5_wash_noptap', 'f5_wash_tap', 
           # 'f5_wash_basin', 'f5_wash_soap', 'f5_wash_ash', 
           # 'f5_child_feces_factor', 'f5_feces_visible', 'f5_house_feces_factor',
)

#table(cases$f5_child_health_miss,is.na(cases$f5_skin_flaky))
#f4a_share_fac leads to some being dropped, keep considering for now
#variables w/ small cell sizes, dropping for now: "f4a_chlorine","f4a_primcare","f4a_mom_live",
# "f4a_drh_undrink","f4a_drh_fever","f4a_drh_breath","f4b_skin_pinch",
# "f5_diag_meng","f5_wash_piped","f4a_fetch_water","f4a_trip_day","f4a_trip_week",
# f4a_notrt_water, f4b_volume, 
#these lead to lots being dropped >> are physical exam at f-up, not useful for predicting what happens at enrollment, so not analyze anymore
#f5_rectal
#f5_bipedal
#f5_abn_hair
#f5_under_nutr
#f5_skin_flaky

#checking all appropriates ones are factorized
# temp<-temp_gf[names]
# str(temp)
# summary(temp)

################### cases into age and other sub groups ####
complete_anydeath_age1 <- complete_anydeath %>% filter(agegroup == 1)
complete_anydeath_age2 <- complete_anydeath %>% filter(agegroup == 2)
complete_anydeath_age3 <- complete_anydeath %>% filter(agegroup == 3)
complete_anydeath_age4 <- complete_anydeath %>% filter(agegroup == 1 | agegroup==2)
complete_anydeath_age5 <- complete_anydeath %>% filter(agegroup == 2 | agegroup==3)

complete_hospdeath_age1 <- complete_hospdeath %>% filter(agegroup == 1)
complete_hospdeath_age2 <- complete_hospdeath %>% filter(agegroup == 2)
complete_hospdeath_age3 <- complete_hospdeath %>% filter(agegroup == 3)
complete_hospdeath_age4 <- complete_hospdeath %>% filter(agegroup == 1 | agegroup==2)

complete_homedeath_age1 <- complete_homedeath %>% filter(agegroup == 1)
complete_homedeath_age2 <- complete_homedeath %>% filter(agegroup == 2)
complete_homedeath_age3 <- complete_homedeath %>% filter(agegroup == 3)
complete_homedeath_age4 <- complete_homedeath %>% filter(agegroup == 1 | agegroup==2)

# summary(complete_anydeath$site)
# 1    2    3    4    5    6    7 
# 861 1783  525 1124 1473 1348  959 
#Gambia,Mali,Mozambique,Kenya,India,Bdesh,Pakistan
cases_anydeath_Afr <- complete_anydeath %>% filter(site==1 | site==2 | site==3 | site==4)
cases_anydeath_SEAsia <- complete_anydeath %>% filter(site==5 | site==6 | site==7)

complete_anydeath_Gambia <- complete_anydeath %>% filter(site == 1)
complete_anydeath_Mali <- complete_anydeath %>% filter(site == 2)
complete_anydeath_Mozam <- complete_anydeath %>% filter(site == 3)
complete_anydeath_Kenya <- complete_anydeath %>% filter(site == 4)
complete_anydeath_India <- complete_anydeath %>% filter(site == 5)
complete_anydeath_Bang <- complete_anydeath %>% filter(site == 6)
complete_anydeath_Pak <- complete_anydeath %>% filter(site == 7)

################### descriptive for pub ####
dim(complete_anydeath)
# [1] 8060  349
table(complete_anydeath$death_all)
# 0    1 
# 7895  165 
dim(complete_hospdeath)
# [1] 8060  349
table(complete_anydeath$death_hosp)
# 0    1 
# 8017   43
table(complete_hospdeath$death_hosp)
# 0    1 
# 8017   43 
dim(complete_homedeath)
# [1] 8017  349
table(complete_anydeath$death_home)
# 0    1 
# 7938  122 
table(complete_homedeath$death_home)
# 0    1 
# 7895  122 

table(complete_anydeath$agegroup)
table(complete_anydeath$death_hosp,complete_anydeath$agegroup)
table(complete_anydeath$death_home,complete_anydeath$agegroup)
table(complete_anydeath$death_all,complete_anydeath$agegroup)

table(complete_anydeath$site)
table(complete_anydeath$death_hosp,complete_anydeath$site)
table(complete_anydeath$death_home,complete_anydeath$site)
table(complete_anydeath$death_all,complete_anydeath$site)


#overlapping histograms comparing distributions of top predictive variables
p1 <- hist(complete_anydeath[which(complete_anydeath$death_all==0),]$f4b_muac,freq=F)
p2 <- hist(complete_anydeath[which(complete_anydeath$death_all==1),]$f4b_muac,freq=F)
p3 <- hist(kilifi[which(kilifi$death_all==0),]$f4b_muac,freq=F)
p4 <- hist(kilifi[which(kilifi$death_all==1),]$f4b_muac,freq=F)
#jpeg("/overlap_hist_muac.jpg")
par(mfrow=c(2,1))
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,25),ylim=c(0,0.3),freq=F,xlab="baseline MUAC",main="Death any time by baseline MUAC in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,25),ylim=c(0,0.3),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
plot( p3, col=rgb(1,0,0,1/4), xlim=c(0,25),ylim=c(0,0.3),freq=F,xlab="baseline MUAC",main="Death any time by baseline MUAC in Kilifi")
plot( p4, col=rgb(0,0,1,1/4),  xlim=c(0,25),ylim=c(0,0.3),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()

p1 <- hist(complete_anydeath[which(complete_anydeath$death_all==0),]$f4b_resp,freq=F)
p2 <- hist(complete_anydeath[which(complete_anydeath$death_all==1),]$f4b_resp,freq=F)
p3 <- hist(kilifi[which(kilifi$death_all==0),]$f4b_resp,freq=F)
p4 <- hist(kilifi[which(kilifi$death_all==1),]$f4b_resp,freq=F)
#jpeg("/overlap_hist_resp.jpg")
par(mfrow=c(2,1))
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,120),ylim=c(0,0.08),freq=F,xlab="Respiratory Rate",main="Death any time by Respiratory Rate in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,120),ylim=c(0,0.08),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
plot( p3, col=rgb(1,0,0,1/4), xlim=c(0,120),ylim=c(0,0.08),freq=F,xlab="Respiratory Rate",main="Death any time by Respiratory Rate in Kilifi")
plot( p4, col=rgb(0,0,1,1/4),  xlim=c(0,120),ylim=c(0,0.08),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()

p1 <- hist(complete_anydeath[which(complete_anydeath$death_all==0),]$f4b_temp,freq=F)
p2 <- hist(complete_anydeath[which(complete_anydeath$death_all==1),]$f4b_temp,freq=F)
p3 <- hist(kilifi[which(kilifi$death_all==0),]$f4b_temp,freq=F)
p4 <- hist(kilifi[which(kilifi$death_all==1),]$f4b_temp,freq=F)
#jpeg("/overlap_hist_temp.jpg")
par(mfrow=c(2,1))
plot( p1, col=rgb(1,0,0,1/4), xlim=c(30,45),ylim=c(0,0.6),freq=F,xlab="Temperature (?C)",main="Death any time by Temperature in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(30,45),ylim=c(0,0.6),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
plot( p3, col=rgb(1,0,0,1/4), xlim=c(30,45),ylim=c(0,0.6),freq=F,xlab="Temperature (?C)",main="Death any time by Temperature in Kilifi")
plot( p4, col=rgb(0,0,1,1/4),  xlim=c(30,45),ylim=c(0,0.6),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()

p1 <- hist(complete_anydeath[which(complete_anydeath$death_all==0),]$base_age,freq=F)
p2 <- hist(complete_anydeath[which(complete_anydeath$death_all==1),]$base_age,freq=F)
p3 <- hist(kilifi[which(kilifi$death_all==0),]$base_age,freq=F)
p4 <- hist(kilifi[which(kilifi$death_all==1),]$base_age,freq=F)
#jpeg("/overlap_hist_age.jpg")
par(mfrow=c(2,1))
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,0.08),freq=F,xlab="Age (months)",main="Death any time by baseline Age in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,0.08),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
plot( p3, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,0.08),freq=F,xlab="Age (months)",main="Death any time by baseline Age in Kilifi")
plot( p4, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,008),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()


#dates 
options(max.print=1500)
table(gems1_orig$enrolldate) #Dec 2007 - Mar 2011


names<-append(x=names, values=c("death_all"))
temp <- complete_anydeath[names]


#table comparing included and excluded data
# cases_anydeath (n=8520)
# complete_anydeath (n=8060)
# 8520-8060 = 460

#creating dataset of what dropped for missing predictors
missing <- anti_join(cases_anydeath,complete_anydeath)
dim(missing)

missing <- cbind(missing,included=rep(0,dim(missing)[1]))
data.to.table<-complete_anydeath
data.to.table <- cbind(data.to.table,included=rep(1,dim(data.to.table)[1]))
data.to.table<-rbind(missing,data.to.table)
table1(~ as.factor(death_all) + base_age + site + f4b_muac + f4b_resp + 
         f4b_temp + f4a_ppl_house | as.factor(included), 
       data=data.to.table)


################### descriptive Se/Sp for different regimens ####
#proportion who die among each subset below
#0-6mo
# table(complete_anydeath[which(complete_anydeath$base_age<=6),]$death_all)[2]
# dim(complete_anydeath[which(complete_anydeath$base_age<=6),])[1]
# round(((table(complete_anydeath[which(complete_anydeath$base_age<=6),]$death_all)[2])/(dim(complete_anydeath[which(complete_anydeath$base_age<=6),])[1]))*100,2)

complete_anydeath$test_mo0_6 <- ifelse(complete_anydeath$base_age<=6,1,0)
TP_mo0_6 <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_6)[4]
TN_mo0_6 <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_6)[1]
FP_mo0_6 <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_6)[3]
FN_mo0_6 <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_6)[2]
table(complete_anydeath$death_all,complete_anydeath$test_mo0_6)
# 0    1
# 0 6515 1380
# 1  116   49
#rows are true, columns are test

Se_mo0_6 = round(TP_mo0_6/(TP_mo0_6+FN_mo0_6),2) #test positive | truly positive
Sp_mo0_6 = round(TN_mo0_6/(TN_mo0_6+FP_mo0_6),2) #test negative | truly negative
PPV_mo0_6 = round(TP_mo0_6/(TP_mo0_6+FP_mo0_6),2) #truly positive | test positive
NPV_mo0_6 = round(TN_mo0_6/(TN_mo0_6+FN_mo0_6),2) #truly negative | test negative
FPR_mo0_6 = round(FP_mo0_6/(FP_mo0_6+TN_mo0_6),2) #test positive | truly negative
FNR_mo0_6 = round(FN_mo0_6/(TP_mo0_6+FN_mo0_6),2) #test negative | truly positive

#proportion who screen positive
round((TP_mo0_6+FP_mo0_6)/(TP_mo0_6+TN_mo0_6+FP_mo0_6+FN_mo0_6)*100,1)

#7-59mo
# table(complete_anydeath[which(complete_anydeath$base_age>6),]$death_all)[2]
# dim(complete_anydeath[which(complete_anydeath$base_age>6),])[1]
# round(((table(complete_anydeath[which(complete_anydeath$base_age>6),]$death_all)[2])/(dim(complete_anydeath[which(complete_anydeath$base_age>6),])[1]))*100,2)

complete_anydeath$test_mo7_59 <- ifelse(complete_anydeath$base_age>6,1,0)
TP_mo7_59 <- table(complete_anydeath$death_all,complete_anydeath$test_mo7_59)[4]
TN_mo7_59 <- table(complete_anydeath$death_all,complete_anydeath$test_mo7_59)[1]
FP_mo7_59 <- table(complete_anydeath$death_all,complete_anydeath$test_mo7_59)[3]
FN_mo7_59 <- table(complete_anydeath$death_all,complete_anydeath$test_mo7_59)[2]

Se_mo7_59 = round(TP_mo7_59/(TP_mo7_59+FN_mo7_59),2) #test positive | truly positive
Sp_mo7_59 = round(TN_mo7_59/(TN_mo7_59+FP_mo7_59),2) #test negative | truly negative
PPV_mo7_59 = round(TP_mo7_59/(TP_mo7_59+FP_mo7_59),2) #truly positive | test positive
NPV_mo7_59 = round(TN_mo7_59/(TN_mo7_59+FN_mo7_59),2) #truly negative | test negative
FPR_mo7_59 = round(FP_mo7_59/(FP_mo7_59+TN_mo7_59),2) #test positive | truly negative
FNR_mo7_59 = round(FN_mo7_59/(TP_mo7_59+FN_mo7_59),2) #test negative | truly positive

#proportion who screen positive
round((TP_mo7_59+FP_mo7_59)/(TP_mo7_59+TN_mo7_59+FP_mo7_59+FN_mo7_59)*100,1)

#0-59mo
# table(complete_anydeath$death_all)
# dim(complete_anydeath)
# round(((table(complete_anydeath$death_all)[2])/(dim(complete_anydeath)[1]))*100,2)

# Se #test positive | truly positive
# Sp #test negative | truly negative
# PPV #truly positive | test positive
# NPV #truly negative | test negative
# FPR #test positive | truly negative
# FNR #test negative | truly positive

#0-6mo AND MUAC <12.5
# table(complete_anydeath[which(complete_anydeath$base_age<=6&complete_anydeath$f4b_muac<12.5),]$death_all)[2]
# dim(complete_anydeath[which(complete_anydeath$base_age<=6&complete_anydeath$f4b_muac<12.5),])[1]
# round(((table(complete_anydeath[which(complete_anydeath$base_age<=6&complete_anydeath$f4b_muac<12.5),]$death_all)[2])/(dim(complete_anydeath[which(complete_anydeath$base_age<=6&complete_anydeath$f4b_muac<12.5),])[1]))*100,2)

complete_anydeath$test_mo0_6MUAC <- ifelse((complete_anydeath$base_age<=6&complete_anydeath$f4b_muac<12.5),1,0)
TP_mo0_6MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_6MUAC)[4]
TN_mo0_6MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_6MUAC)[1]
FP_mo0_6MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_6MUAC)[3]
FN_mo0_6MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_6MUAC)[2]

Se_mo0_6MUAC = round(TP_mo0_6MUAC/(TP_mo0_6MUAC+FN_mo0_6MUAC),2) #test positive | truly positive
Sp_mo0_6MUAC = round(TN_mo0_6MUAC/(TN_mo0_6MUAC+FP_mo0_6MUAC),2) #test negative | truly negative
PPV_mo0_6MUAC = round(TP_mo0_6MUAC/(TP_mo0_6MUAC+FP_mo0_6MUAC),2) #truly positive | test positive
NPV_mo0_6MUAC = round(TN_mo0_6MUAC/(TN_mo0_6MUAC+FN_mo0_6MUAC),2) #truly negative | test negative
FPR_mo0_6MUAC = round(FP_mo0_6MUAC/(FP_mo0_6MUAC+TN_mo0_6MUAC),2) #test positive | truly negative
FNR_mo0_6MUAC = round(FN_mo0_6MUAC/(TP_mo0_6MUAC+FN_mo0_6MUAC),2) #test negative | truly positive

#proportion who screen positive
round((TP_mo0_6MUAC+FP_mo0_6MUAC)/(TP_mo0_6MUAC+TN_mo0_6MUAC+FP_mo0_6MUAC+FN_mo0_6MUAC)*100,1)

#7-mo AND MUAC <12.5
# table(complete_anydeath[which(complete_anydeath$base_age>6&complete_anydeath$f4b_muac<12.5),]$death_all)[2]
# dim(complete_anydeath[which(complete_anydeath$base_age>6&complete_anydeath$f4b_muac<12.5),])[1]
# round(((table(complete_anydeath[which(complete_anydeath$base_age>6&complete_anydeath$f4b_muac<12.5),]$death_all)[2])/(dim(complete_anydeath[which(complete_anydeath$base_age>6&complete_anydeath$f4b_muac<12.5),])[1]))*100,2)

complete_anydeath$test_mo7_59MUAC <- ifelse((complete_anydeath$base_age>6&complete_anydeath$f4b_muac<12.5),1,0)
TP_mo7_59MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo7_59MUAC)[4]
TN_mo7_59MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo7_59MUAC)[1]
FP_mo7_59MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo7_59MUAC)[3]
FN_mo7_59MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo7_59MUAC)[2]

Se_mo7_59MUAC = round(TP_mo7_59MUAC/(TP_mo7_59MUAC+FN_mo7_59MUAC),2) #test positive | truly positive
Sp_mo7_59MUAC = round(TN_mo7_59MUAC/(TN_mo7_59MUAC+FP_mo7_59MUAC),2) #test negative | truly negative
PPV_mo7_59MUAC = round(TP_mo7_59MUAC/(TP_mo7_59MUAC+FP_mo7_59MUAC),2) #truly positive | test positive
NPV_mo7_59MUAC = round(TN_mo7_59MUAC/(TN_mo7_59MUAC+FN_mo7_59MUAC),2) #truly negative | test negative
FPR_mo7_59MUAC = round(FP_mo7_59MUAC/(FP_mo7_59MUAC+TN_mo7_59MUAC),2) #test positive | truly negative
FNR_mo7_59MUAC = round(FN_mo7_59MUAC/(TP_mo7_59MUAC+FN_mo7_59MUAC),2) #test negative | truly positive

#proportion who screen positive
round((TP_mo7_59MUAC+FP_mo7_59MUAC)/(TP_mo7_59MUAC+TN_mo7_59MUAC+FP_mo7_59MUAC+FN_mo7_59MUAC)*100,1)

#0-59mo AND MUAC <12.5
# table(complete_anydeath[which(complete_anydeath$f4b_muac<12.5),]$death_all)[2]
# dim(complete_anydeath[which(complete_anydeath$f4b_muac<12.5),])[1]
# round(((table(complete_anydeath[which(complete_anydeath$f4b_muac<12.5),]$death_all)[2])/(dim(complete_anydeath[which(complete_anydeath$f4b_muac<12.5),])[1]))*100,2)

complete_anydeath$test_mo0_59MUAC <- ifelse(complete_anydeath$f4b_muac<12.5,1,0)
TP_mo0_59MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_59MUAC)[4]
TN_mo0_59MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_59MUAC)[1]
FP_mo0_59MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_59MUAC)[3]
FN_mo0_59MUAC <- table(complete_anydeath$death_all,complete_anydeath$test_mo0_59MUAC)[2]

Se_mo0_59MUAC = round(TP_mo0_59MUAC/(TP_mo0_59MUAC+FN_mo0_59MUAC),2) #test positive | truly positive
Sp_mo0_59MUAC = round(TN_mo0_59MUAC/(TN_mo0_59MUAC+FP_mo0_59MUAC),2) #test negative | truly negative
PPV_mo0_59MUAC = round(TP_mo0_59MUAC/(TP_mo0_59MUAC+FP_mo0_59MUAC),2) #truly positive | test positive
NPV_mo0_59MUAC = round(TN_mo0_59MUAC/(TN_mo0_59MUAC+FN_mo0_59MUAC),2) #truly negative | test negative
FPR_mo0_59MUAC = round(FP_mo0_59MUAC/(FP_mo0_59MUAC+TN_mo0_59MUAC),2) #test positive | truly negative
FNR_mo0_59MUAC = round(FN_mo0_59MUAC/(TP_mo0_59MUAC+FN_mo0_59MUAC),2) #test negative | truly positive

#proportion who screen positive
round((TP_mo0_59MUAC+FP_mo0_59MUAC)/(TP_mo0_59MUAC+TN_mo0_59MUAC+FP_mo0_59MUAC+FN_mo0_59MUAC)*100,1)

#single CPR model from 0-59mo, 2var (MUAC, respiratory rate)
#predicted prob of death_all
GEMS_glm_death <- glm(death_all~f4b_muac+f4b_resp,
                               data=cases_anydeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.767397   0.690270   6.907 4.97e-12 ***
#   f4b_muac    -0.751158   0.047493 -15.816  < 2e-16 ***
#   f4b_resp     0.025933   0.006945   3.734 0.000189 ***

round(exp(coef(GEMS_glm_death)),2)
# (Intercept)    f4b_muac    f4b_resp 
# 117.61        0.47        1.03 
round(exp(confint(GEMS_glm_death)),2)
# 2.5 % 97.5 %
#   (Intercept) 30.54 457.97
# f4b_muac     0.43   0.52
# f4b_resp     1.01   1.04

GEMS_2var_fit<-cases_anydeath %>% select(death_all,f4b_muac,f4b_resp)
GEMS_2var_fit$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_2var_fit,type="response"))
GEMS_2var_fit <- na.omit(GEMS_2var_fit)
#8516
GEMS_2var_fit$test_CPR0.05 <- ifelse(GEMS_2var_fit$pred_glm>=0.05,1,0)
#table(signif(GEMS_2var_fit$pred_glm,2),GEMS_2var_fit$test_0.05)
GEMS_2var_fit$test_CPR0.1 <- ifelse(GEMS_2var_fit$pred_glm>=0.1,1,0)
#table(signif(GEMS_2var_fit$pred_glm,2),GEMS_2var_fit$test_0.1)
GEMS_2var_fit$test_CPR0.15 <- ifelse(GEMS_2var_fit$pred_glm>=0.15,1,0)
#table(signif(GEMS_2var_fit$pred_glm,2),GEMS_2var_fit$test_CPR0.15)
GEMS_2var_fit$test_CPR0.2 <- ifelse(GEMS_2var_fit$pred_glm>=0.2,1,0)
GEMS_2var_fit$test_CPR0.25 <- ifelse(GEMS_2var_fit$pred_glm>=0.25,1,0)

# table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.05)
# 0    1
# 0 7714  623
# 1   94   85
#rows are true, columns are test

TP_CPR0.05 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.05)[4]
TN_CPR0.05 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.05)[1]
FP_CPR0.05 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.05)[3]
FN_CPR0.05 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.05)[2]

Se_CPR0.05 = round(TP_CPR0.05/(TP_CPR0.05+FN_CPR0.05),2) #test positive | truly positive
Sp_CPR0.05 = round(TN_CPR0.05/(TN_CPR0.05+FP_CPR0.05),2) #test negative | truly negative
PPV_CPR0.05 = round(TP_CPR0.05/(TP_CPR0.05+FP_CPR0.05),2) #truly positive | test positive
NPV_CPR0.05 = round(TN_CPR0.05/(TN_CPR0.05+FN_CPR0.05),2) #truly negative | test negative
FPR_CPR0.05 = round(FP_CPR0.05/(FP_CPR0.05+TN_CPR0.05),2) #test positive | truly negative
FNR_CPR0.05 = round(FN_CPR0.05/(TP_CPR0.05+FN_CPR0.05),2) #test negative | truly positive

#proportion who screen positive
round((TP_CPR0.05+FP_CPR0.05)/(TP_CPR0.05+TN_CPR0.05+FP_CPR0.05+FN_CPR0.05)*100,1)

TP_CPR0.1 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.1)[4]
TN_CPR0.1 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.1)[1]
FP_CPR0.1 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.1)[3]
FN_CPR0.1 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.1)[2]

Se_CPR0.1 = round(TP_CPR0.1/(TP_CPR0.1+FN_CPR0.1),2) #test positive | truly positive
Sp_CPR0.1 = round(TN_CPR0.1/(TN_CPR0.1+FP_CPR0.1),2) #test negative | truly negative
PPV_CPR0.1 = round(TP_CPR0.1/(TP_CPR0.1+FP_CPR0.1),2) #truly positive | test positive
NPV_CPR0.1 = round(TN_CPR0.1/(TN_CPR0.1+FN_CPR0.1),2) #truly negative | test negative
FPR_CPR0.1 = round(FP_CPR0.1/(FP_CPR0.1+TN_CPR0.1),2) #test positive | truly negative
FNR_CPR0.1 = round(FN_CPR0.1/(TP_CPR0.1+FN_CPR0.1),2) #test negative | truly positive

#proportion who screen positive
round((TP_CPR0.1+FP_CPR0.1)/(TP_CPR0.1+TN_CPR0.1+FP_CPR0.1+FN_CPR0.1)*100,1)

TP_CPR0.15 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.15)[4]
TN_CPR0.15 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.15)[1]
FP_CPR0.15 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.15)[3]
FN_CPR0.15 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.15)[2]

Se_CPR0.15 = round(TP_CPR0.15/(TP_CPR0.15+FN_CPR0.15),2) #test positive | truly positive
Sp_CPR0.15 = round(TN_CPR0.15/(TN_CPR0.15+FP_CPR0.15),2) #test negative | truly negative
PPV_CPR0.15 = round(TP_CPR0.15/(TP_CPR0.15+FP_CPR0.15),2) #truly positive | test positive
NPV_CPR0.15 = round(TN_CPR0.15/(TN_CPR0.15+FN_CPR0.15),2) #truly negative | test negative
FPR_CPR0.15 = round(FP_CPR0.15/(FP_CPR0.15+TN_CPR0.15),2) #test positive | truly negative
FNR_CPR0.15 = round(FN_CPR0.15/(TP_CPR0.15+FN_CPR0.15),2) #test negative | truly positive

#proportion who screen positive
round((TP_CPR0.15+FP_CPR0.15)/(TP_CPR0.15+TN_CPR0.15+FP_CPR0.15+FN_CPR0.15)*100,1)

TP_CPR0.2 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.2)[4]
TN_CPR0.2 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.2)[1]
FP_CPR0.2 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.2)[3]
FN_CPR0.2 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.2)[2]

Se_CPR0.2 = round(TP_CPR0.2/(TP_CPR0.2+FN_CPR0.2),2) #test positive | truly positive
Sp_CPR0.2 = round(TN_CPR0.2/(TN_CPR0.2+FP_CPR0.2),2) #test negative | truly negative
PPV_CPR0.2 = round(TP_CPR0.2/(TP_CPR0.2+FP_CPR0.2),2) #truly positive | test positive
NPV_CPR0.2 = round(TN_CPR0.2/(TN_CPR0.2+FN_CPR0.2),2) #truly negative | test negative
FPR_CPR0.2 = round(FP_CPR0.2/(FP_CPR0.2+TN_CPR0.2),2) #test positive | truly negative
FNR_CPR0.2 = round(FN_CPR0.2/(TP_CPR0.2+FN_CPR0.2),2) #test negative | truly positive

#proportion who screen positive
round((TP_CPR0.2+FP_CPR0.2)/(TP_CPR0.2+TN_CPR0.2+FP_CPR0.2+FN_CPR0.2)*100,1)

TP_CPR0.25 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.25)[4]
TN_CPR0.25 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.25)[1]
FP_CPR0.25 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.25)[3]
FN_CPR0.25 <- table(GEMS_2var_fit$death_all,GEMS_2var_fit$test_CPR0.25)[2]

Se_CPR0.25 = round(TP_CPR0.25/(TP_CPR0.25+FN_CPR0.25),2) #test positive | truly positive
Sp_CPR0.25 = round(TN_CPR0.25/(TN_CPR0.25+FP_CPR0.25),2) #test negative | truly negative
PPV_CPR0.25 = round(TP_CPR0.25/(TP_CPR0.25+FP_CPR0.25),2) #truly positive | test positive
NPV_CPR0.25 = round(TN_CPR0.25/(TN_CPR0.25+FN_CPR0.25),2) #truly negative | test negative
FPR_CPR0.25 = round(FP_CPR0.25/(FP_CPR0.25+TN_CPR0.25),2) #test positive | truly negative
FNR_CPR0.25 = round(FN_CPR0.25/(TP_CPR0.25+FN_CPR0.25),2) #test negative | truly positive

#proportion who screen positive
round((TP_CPR0.25+FP_CPR0.25)/(TP_CPR0.25+TN_CPR0.25+FP_CPR0.25+FN_CPR0.25)*100,1)

#make two color hist of predicted prob based on if do or do not die
p1 <- hist(GEMS_2var_fit[which(GEMS_2var_fit$death_all==0),]$pred_glm,freq=F)
p2 <- hist(GEMS_2var_fit[which(GEMS_2var_fit$death_all==1),]$pred_glm,freq=F)
#jpeg("/overlap_hist_PredProb.jpg")
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,1),ylim=c(0,20),freq=F,xlab="predicted probability",main="Death any time by predicted probability of death")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,1),ylim=c(0,20),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()


#how many ppl don't actually die above given cutoff
table(GEMS_2var_fit[which(GEMS_2var_fit$pred_glm>=0.25),]$death_all)
table(GEMS_2var_fit[which(GEMS_2var_fit$pred_glm>=0.75),]$death_all)
#>>> are still some that survive even at high predicted prob of death
#>>> decision making is how many want to have nutritional intervention on


################### MAIN (0-59mo) any death all cases -HAZ +MUAC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

main <- CPR.funct(data=complete_anydeath,outcome="death_all",iter=100,nvars_opts=c(1:10,15,20,30,40,50))
main[["df_imps"]]
main[["AUC_df"]]
main[["calib"]]

# names      var_red
# 1               f4b_muac 12.967363484
# 2               f4b_resp  5.859098155
# 3               f4b_temp  5.442182853
# 4               base_age  5.279787333
# 5          f4a_ppl_house  3.769922710
# 6           f4a_drh_days  2.938228500
# 7           f4b_abn_hair  2.543875823
# 8       f4a_yng_children  2.526648853
# 9         f4a_offr_drink  2.460260605
# 10         f4a_slp_rooms  2.404415055
# 11        f4b_under_nutr  2.207884416
# 12          f4a_dad_live  2.103694807
# 13        f4a_disp_feces  1.977644900
# 14          f4a_ms_water  1.740106672
# 15             f3_drh_iv  1.618439359
# 16       f4b_chest_indrw  1.580620985
# 17              f4b_skin  1.566358115
# 18                  site  1.552852656
# 19       f4a_water_avail  1.548925549
# 20      f4a_relationship  1.477645553
# 21         f4a_prim_schl  1.473594444
# 22            f4b_mental  1.472513303
# 23         f4b_recommend  1.468801742
# 24         f4a_breastfed  1.461976024
# 25         f4a_share_fac  1.427480032
# 26           f3_drh_hosp  1.417129740
# 27        f4b_skin_flaky  1.348903766
# 28         f3_drh_turgor  1.270248902
# 29             f4b_admit  1.184790359
# 30     f4a_drh_bellypain  1.109433291
# 31    f4a_cur_fastbreath  1.104088902
# 32        f4a_trt_method  1.087808104
# 33          f4a_cur_skin  1.045628317
# 34      f4a_drh_prolapse  1.014995156
# 35      f4a_hometrt_herb  1.013604010
# 36          f4a_wash_def  0.986443005
# 37        f4a_house_elec  0.984153547
# 38        f4a_max_stools  0.972485154
# 39         f4a_drh_vomit  0.950335431
# 40         f4a_drh_cough  0.948494063
# 41  f4a_drh_lethrgy_miss  0.897771472
# 42          f4a_ani_goat  0.842484887
# 43         f4a_wash_cook  0.836300819
# 44     f4a_water_pubwell  0.831211938
# 45        f4a_water_bore  0.822438999
# 46        f4a_wash_nurse  0.792682375
# 47             f4b_mouth  0.790127738
# 48           f4a_ani_cow  0.766835715
# 49       f4a_cur_thirsty  0.756445474
# 50           f4b_bipedal  0.751450073
# 51       f4a_ani_rodents  0.751342378
# 52      f4a_hometrt_none  0.746108398
# 53          f4a_wash_use  0.727154542
# 54             f3_gender  0.708671919
# 55     f4a_drh_lessdrink  0.707029928
# 56           f4a_ani_cat  0.707004732
# 57      f4a_drh_restless  0.704240118
# 58      f4a_seek_outside  0.702924316
# 59      f4a_house_agland  0.699184589
# 60       f4a_hometrt_ors  0.697000704
# 61      f4a_cur_restless  0.695302055
# 62        f4a_drh_thirst  0.687487088
# 63       f4a_seek_healer  0.674943068
# 64         f4a_ani_sheep  0.673808739
# 65        f4a_house_bike  0.671487165
# 66       f4a_house_radio  0.670568576
# 67         f4a_drh_consc  0.666112011
# 68          f4a_drh_conv  0.662763726
# 69     f4a_water_covwell  0.655322548
# 70         f4a_drh_blood  0.654328813
# 71        f4a_wash_child  0.642083574
# 72         f4a_fac_waste  0.633957390
# 73         f4a_trt_water  0.632226053
# 74     f4a_fuel_charcoal  0.630350417
# 75          f4a_ani_fowl  0.626565534
# 76        f4a_house_tele  0.623427425
# 77         f4a_ani_other  0.604246977
# 78             f4a_floor  0.604237322
# 79           f4a_ani_dog  0.594021205
# 80      f4a_seek_privdoc  0.593990314
# 81          f4a_wash_eat  0.569292740
# 82     f4a_hometrt_othr1  0.565145218
# 83       f4a_house_phone  0.564905675
# 84      f4a_water_pubtap  0.560800945
# 85       f4a_store_water  0.558311889
# 86      f4a_cur_drymouth  0.546597783
# 87       f4a_house_scoot  0.535603410
# 88        f4a_house_cart  0.522101706
# 89       f4a_fuel_biogas  0.521120147
# 90      f4a_hometrt_milk  0.516245849
# 91        f4a_hometrt_ab  0.507886481
# 92    f4a_water_deepwell  0.482295311
# 93        f4a_drh_strain  0.479440074
# 94         f4a_fuel_wood  0.478298483
# 95              f4b_eyes  0.460732841
# 96       f4a_water_river  0.458186986
# 97     f4a_hometrt_maize  0.425575293
# 98   f4a_water_prospring  0.417849402
# 99        f4a_water_pond  0.415310523
# 100       f4a_water_yard  0.403275987
# 101        f4a_wash_othr  0.378840785
# 102   f4a_water_covpwell  0.376964358
# 103       f4a_water_rain  0.369609895
# 104     f4a_water_bought  0.360860624
# 105       f4a_water_well  0.356693223
# 106     f4a_house_fridge  0.342458672
# 107           f4a_ani_no  0.341452379
# 108         f4a_seek_doc  0.338144805
# 109       f4a_seek_pharm  0.335474649
# 110  f4a_hometrt_othrliq  0.325464981
# 111        f4a_fuel_kero  0.284948924
# 112       f4a_seek_other  0.278235595
# 113    f4a_hometrt_othr2  0.266613441
# 114       f4a_fuel_grass  0.260864233
# 115        f4a_house_car  0.259911656
# 116           f4b_rectal  0.249900300
# 117      f4a_fuel_natgas  0.243547096
# 118      f4a_wash_animal  0.236173278
# 119       f4a_fuel_other  0.218483853
# 120        f4a_fuel_elec  0.197190054
# 121        f4a_fuel_crop  0.172633789
# 122       f4a_house_none  0.142680074
# 123  f4a_water_shallwell  0.140708857
# 124      f4a_water_house  0.132768201
# 125     f4a_hometrt_zinc  0.121402425
# 126       f4a_water_othr  0.102867407
# 127        f4a_fuel_dung  0.099084738
# 128      f4a_seek_friend  0.065405733
# 129       f4a_seek_remdy  0.064526568
# 130   f4a_water_unspring  0.062048688
# 131        f4a_fuel_coal  0.006515613
# 132     f4a_fuel_propane  0.004823056
# 133       f4a_house_boat  0.002883586

# AUC          SE     lower     upper level Model nvar
# 1  0.8328280 0.011322805 0.8106357 0.8550203  0.95    LR    1
# 2  0.8496393 0.010156834 0.8297323 0.8695464  0.95    LR    2
# 3  0.8524532 0.010412218 0.8320457 0.8728608  0.95    LR    3
# 4  0.8577002 0.010425682 0.8372662 0.8781342  0.95    LR    4
# 5  0.8572182 0.010394050 0.8368463 0.8775902  0.95    LR    5
# 6  0.8584201 0.010422712 0.8379920 0.8788483  0.95    LR    6
# 7  0.8580585 0.010349747 0.8377734 0.8783436  0.95    LR    7
# 8  0.8606733 0.010458541 0.8401749 0.8811717  0.95    LR    8
# 9  0.8603688 0.010603705 0.8395859 0.8811517  0.95    LR    9
# 10 0.8648231 0.010455832 0.8443300 0.8853161  0.95    LR   10
# 11 0.8782162 0.010165830 0.8582916 0.8981409  0.95    LR   15
# 12 0.8881306 0.009216247 0.8700670 0.9061941  0.95    LR   20
# 13 0.8785544 0.010790768 0.8574049 0.8997040  0.95    LR   30
# 14 0.8732973 0.011070888 0.8515987 0.8949958  0.95    LR   40
# 15 0.8722955 0.011077408 0.8505842 0.8940068  0.95    LR   50
# 16 0.6561734 0.025096826 0.6069845 0.7053623  0.95    RF    1
# 17 0.7933909 0.012507896 0.7688758 0.8179059  0.95    RF    2
# 18 0.8173961 0.011106414 0.7956279 0.8391643  0.95    RF    3
# 19 0.8395647 0.011807364 0.8164227 0.8627067  0.95    RF    4
# 20 0.8390515 0.011126691 0.8172436 0.8608594  0.95    RF    5
# 21 0.8416746 0.010476016 0.8211420 0.8622072  0.95    RF    6
# 22 0.8456661 0.010105338 0.8258600 0.8654722  0.95    RF    7
# 23 0.8559542 0.009529574 0.8372765 0.8746318  0.95    RF    8
# 24 0.8550524 0.010210078 0.8350410 0.8750637  0.95    RF    9
# 25 0.8576928 0.009791566 0.8385017 0.8768839  0.95    RF   10
# 26 0.8747388 0.009411204 0.8562932 0.8931844  0.95    RF   15
# 27 0.8769555 0.010218494 0.8569276 0.8969834  0.95    RF   20
# 28 0.8913088 0.009352626 0.8729780 0.9096396  0.95    RF   30
# 29 0.8970984 0.009103690 0.8792555 0.9149413  0.95    RF   40
# 30 0.8989612 0.008754256 0.8818032 0.9161193  0.95    RF   50

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 -0.00843   -0.392    0.337 1.01      0.738     1.30 
# 2     2 -0.00416   -0.389    0.343 1.01      0.737     1.29 
# 3     3 -0.00511   -0.391    0.343 1.01      0.745     1.29 
# 4     4 -0.0141    -0.401    0.335 1.01      0.757     1.29 
# 5     5 -0.0141    -0.401    0.335 1.01      0.753     1.28 
# 6     6 -0.0209    -0.409    0.330 0.994     0.745     1.26 
# 7     7 -0.0146    -0.404    0.338 0.967     0.726     1.22 
# 8     8 -0.0130    -0.404    0.340 0.971     0.731     1.23 
# 9     9 -0.0139    -0.405    0.340 0.966     0.728     1.22 
# 10    10 -0.00928   -0.401    0.345 0.983     0.747     1.24 
# 11    15 -0.00339   -0.403    0.359 0.927     0.710     1.16 
# 12    20 -0.00883   -0.415    0.361 0.856     0.655     1.08 
# 13    30 -0.0219    -0.437    0.357 0.736     0.553     0.937
# 14    40 -0.0358    -0.454    0.347 0.689     0.515     0.881
# 15    50 -0.0576    -0.481    0.330 0.663     0.495     0.847
#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates too moderate
#slightly positive slope indicates slight underestimation

temp <- main[["decilesCC"]][c("1","2","3","4","5","6","7","8","9","10")]
names(temp) <- c("1-var","2-var","3-var","4-var","5-var","6-var","7-var","8-var","9-var","10-var")  #renaming
#jpeg("/DeathAll_CC_GEMS059_10iter.jpg",width=600,height=480,quality=400)
plot(x=seq(0,0.15,by=0.05),y=seq(0,0.15,by=0.05),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: death in cases 0-59mo in GEMS")),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(temp$`5-var`$`mean(pred_glm)`,temp$`5-var`$`mean(true)`,col="red",pch=1,cex=2,lwd=2)
points(temp$`10-var`$`mean(pred_glm)`,temp$`10-var`$`mean(true)`,col="blue",pch=2,cex=2,lwd=2)
legend("topleft",col=c("red","blue"),c("5-variable","10-variable"),pch=c(1,2),cex=1.5)
dev.off()


AUC_df <- main[["AUC_df"]]
#jpeg("/DeathAll_AUCs_GEMS059_100iter.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(main[["nvars_opts"]])[1]],AUC_df$AUC[1:length(main[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death in cases 0-59mo in GEMS",
     #main=expression(paste("">=0.5,Delta,"HAZ in cases 0-59mo in GEMS")),
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(main[["nvars_opts"]])[1]],AUC_df$AUC[(length(main[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()

glm_DeathAll <- glm(death_all~f4b_muac+f4b_resp+f4b_temp+base_age+f4a_ppl_house+
                      f4a_drh_days+f4a_offr_drink+f4a_yng_children+f4b_abn_hair+f4a_slp_rooms,
              data=complete_anydeath,family="binomial",control=glm.control(maxit=50))
summary(glm_DeathAll)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -12.251007   3.121788  -3.924 8.70e-05 ***
#   f4b_muac          -0.726858   0.057902 -12.553  < 2e-16 ***
#   f4b_resp           0.027355   0.007773   3.519 0.000432 ***
#   f4b_temp           0.412406   0.084312   4.891 1.00e-06 ***
#   base_age           0.016197   0.009438   1.716 0.086130 .  
# f4a_ppl_house     -0.001278   0.011781  -0.108 0.913635    
# f4a_drh_days       0.070985   0.060545   1.172 0.241026    
# f4a_offr_drink     0.298881   0.077334   3.865 0.000111 ***
#   f4a_yng_children  -0.017197   0.055232  -0.311 0.755524    
# f4b_abn_hair       1.393657   0.218244   6.386 1.71e-10 ***
#   f4a_slp_rooms      0.022111   0.031852   0.694 0.487566    
round(exp(coef(glm_DeathAll)),2)
# (Intercept)         f4b_muac         f4b_resp         f4b_temp         base_age 
# 0.00             0.48             1.03             1.51             1.02 
# f4a_ppl_house     f4a_drh_days   f4a_offr_drink f4a_yng_children     f4b_abn_hair 
# 1.00             1.07             1.35             0.98             4.03 
# f4a_slp_rooms 
# 1.02 
round(exp(confint(glm_DeathAll)),2)
# 2.5 % 97.5 %
#   (Intercept)       0.00   0.00
# f4b_muac          0.43   0.54
# f4b_resp          1.01   1.04
# f4b_temp          1.28   1.78
# base_age          1.00   1.03
# f4a_ppl_house     0.97   1.02
# f4a_drh_days      0.95   1.21
# f4a_offr_drink    1.16   1.57
# f4a_yng_children  0.88   1.09
# f4b_abn_hair      2.61   6.14
# f4a_slp_rooms     0.96   1.09

################### main (0-59mo) death ROC curve ####
#want a single ROC surve for CV in test datasets. use "result"
result<-main[["result"]]

roc.data=result %>% split(.,list(.$nvar),drop=TRUE) %>% .$"10" %>% #swith this .$"x" for number of var; #subset to all the iter's (v-fold cross validations) for a given nvar (number of predictor variables in CPR)
  split(.,list(.$iter),drop=TRUE) #now have a list for each iter

#iter are the folds
# table($iter)
#>>> now 1612 per iteration

#reformatting the data so in the same layout as example for 
predict_combo=list(NA)
true_combo=list(NA)

for (i in 1:main[["iter"]]){ #this iter from the main RF/LR loop
  temp.list <- list(roc.data[[i]]$pred_glm)
  predict_combo <- c(predict_combo,temp.list)
  
  temp.list2 <- list(roc.data[[i]]$true)
  true_combo <- c(true_combo,temp.list2)
  
}
predict_combo=predict_combo[-1]
true_combo=true_combo[-1]
str(predict_combo)
str(true_combo)
combo <- list(predict_combo,true_combo)
names(combo) <- c("pred_glm","true")  #renaming
str(combo)

CV.roc.data <- cvAUC(predictions=combo$pred_glm, labels=combo$true) #perf is an object of class "performance" from the ROCR pckg

CV.roc.data.DeathAll.2 <- CV.roc.data
CV.roc.data.DeathAll.5 <- CV.roc.data
CV.roc.data.DeathAll.10 <- CV.roc.data
#to save
CV.roc.data.DeathAll <- list(CV.roc.data.DeathAll.2,CV.roc.data.DeathAll.5,CV.roc.data.DeathAll.10)
names(CV.roc.data.DeathAll) <- c("DA.2","DA.5","DA.10")  #renaming
str(CV.roc.data.DeathAll)
#save(CV.roc.data.DeathAll, file = "/CV.roc.data.DeathAll.Rdata")

# #find the desired points for obtained Sp for given Se
# #options(max.print=2000)
# Se.Sp<-data.frame(Se=round(roc.obj.5var$sensitivities,2),Sp=round(roc.obj.5var$specificities,2))
# Se.Sp[which(Se.Sp$Se==0.85),]
# #between the point (x0[i], y0[i]) and the point (x1[i], y1[i])
# #x0, y0, x1 = x0, y1 = y0,

#load(file = "/CV.roc.data.DeathAll.Rdata")
#(there's a better way to do this as a list...)
CV.roc.data.DA.2 <- CV.roc.data.DeathAll$DA.2
CV.roc.data.DA.5 <- CV.roc.data.DeathAll$DA.5
CV.roc.data.DA.10 <- CV.roc.data.DeathAll$DA.10

#Plot CV AUC >>> SMA: a single line of the averaged cross-validated ROC curve
#tiff("/roc_DeathAll.tif",units="px",width=1500,height=1500,res=300)
plot(CV.roc.data.DA.2$perf, avg="vertical", main="Cross-validated ROC Curves for All Deaths",
     col="#1c61b6", lwd=2, lty=1,
     xlab="False positive rate (1-specificity)",ylab="Average true positive rate (sensitivity)") 
plot(CV.roc.data.DA.5$perf, avg="vertical", col="#1c61b6", lwd=2, lty=2, add=TRUE) 
plot(CV.roc.data.DA.10$perf, avg="vertical", col="#1c61b6", lwd=2, lty=3, add=TRUE) 
legend("bottomright", 
       legend = c("2-var","5-var","10-var"), 
       col = c("#1c61b6"),
       lty = c(1,2,3),
       lwd = 2)
segments(x0=0,y0=0.8,x1=0.25,y1=0.8,lty=2,col="gray")
segments(x0=0.25,y0=0.8,x1=0.25,y1=0.0,lty=2,col="gray")
segments(x0=0,y0=0.9,x1=0.38,y1=0.9,lty=2,col="gray")
segments(x0=0.38,y0=0.9,x1=0.38,y1=0.0,lty=2,col="gray")
dev.off()

################### main (0-59mo) death - what did low prob kids die of? ####
GEMS_glm_death <- glm(death_all~f4b_muac+f4b_resp,
                      data=complete_anydeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.560173   0.716579   6.364 1.97e-10 ***
#   f4b_muac    -0.749498   0.049352 -15.187  < 2e-16 ***
#   f4b_resp     0.030194   0.007255   4.162 3.15e-05 ***
round(exp(coef(GEMS_glm_death)),4)
# (Intercept)    f4b_muac    f4b_resp 
# 95.6000      0.4726      1.0307 

GEMS_2var<-complete_anydeath
GEMS_2var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_2var,type="response"))
hist(GEMS_2var$pred_glm)
summary(GEMS_2var$pred_glm)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.004756 0.009700 0.020472 0.019738 0.876977


#overlapping histograms comparing distributions of predicted prob of death
p1 <- hist(GEMS_2var[which(GEMS_2var$death_all==0),]$pred_glm,freq=F)
p2 <- hist(GEMS_2var[which(GEMS_2var$death_all==1),]$pred_glm,freq=F)
#jpeg("/overlap_hist_predictedprob_2var.jpg")
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,1),ylim=c(0,20),freq=F,xlab="predicted prob death",main="Predicted Probability of death any time in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,1),ylim=c(0,20),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()
#jpeg("/overlap_hist_predictedprob_2var_closeup.jpg")
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,0.2),ylim=c(0,20),freq=F,xlab="predicted prob death",main="Predicted Probability of death any time in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,0.2),ylim=c(0,20),freq=F, add=T)
legend('topright',c('survived','died'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()


#create a subset of observations of low predicted prob of death, see why they died
GEMS_2var_lower10 <- GEMS_2var %>% mutate(percentile_rank = ntile(temp$pred_glm,100)) %>% #lowest percentile has lowest predicted prob of death. other option: ntile(desc(temp$pred_glm),100)) 
  filter(percentile_rank<=25 & death_all==1)
summary(GEMS_2var_lower10$pred_glm)

#season, infection type, nutritional status
GEMS_2var_lower10_subset <- GEMS_2var_lower10 %>% select(caseid,death_all,death_hosp,death_home,base_age,
                                                         site,f4b_muac,f4b_haz,month,f4b_resp,f4b_temp)


################### main (0-59mo) death - what did high prob kids die of? ####
GEMS_glm_death <- glm(death_all~f4b_muac+f4b_resp,
                      data=complete_anydeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.560173   0.716579   6.364 1.97e-10 ***
#   f4b_muac    -0.749498   0.049352 -15.187  < 2e-16 ***
#   f4b_resp     0.030194   0.007255   4.162 3.15e-05 ***
round(exp(coef(GEMS_glm_death)),4)
# (Intercept)    f4b_muac    f4b_resp 
# 95.6000      0.4726      1.0307 

GEMS_2var<-complete_anydeath
GEMS_2var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_2var,type="response"))
hist(GEMS_2var$pred_glm)
summary(GEMS_2var$pred_glm)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.004756 0.009700 0.020472 0.019738 0.876977

#create a subset of observations of high predicted prob of death, see why they died
GEMS_2var_upper10 <- GEMS_2var %>% mutate(percentile_rank = ntile(GEMS_2var$pred_glm,100)) %>% #highest percentile has highest predicted prob of death. other option: ntile(desc(temp$pred_glm),100)) 
  #filter(percentile_rank>=90)
  filter(percentile_rank>=90 & death_all==1)
summary(GEMS_2var_upper10$pred_glm)

#season, infection type, nutritional status
GEMS_2var_upper10_subset <- GEMS_2var_upper10 %>% select(caseid,death_all,death_hosp,death_home,base_age,
                                                         site,f4b_muac,f4b_haz,month,f4b_resp,f4b_temp)
summary(GEMS_2var_upper10_subset)

complete_anydeath_subset <- complete_anydeath %>%  select(caseid,death_all,death_hosp,death_home,base_age,
                                                          site,f4b_muac,f4b_haz,month,f4b_resp,f4b_temp)
summary(complete_anydeath_subset)

################### main + 0-11(age_1) any death all cases -HAZ +MUAC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

main.011 <- CPR.funct(data=complete_anydeath_age1,outcome="death_all",iter=100,nvars_opts=c(1:10,15,20,30,40,50))
main.011[["df_imps"]]
main.011[["AUC_df"]]
main.011[["calib"]]

# names     var_red
# 1               f4b_muac 7.004748198
# 2               f4b_temp 3.809604854
# 3               f4b_resp 3.370398631
# 4               base_age 2.418585462
# 5          f4a_ppl_house 1.876994254
# 6           f4a_drh_days 1.609955791
# 7        f4b_chest_indrw 1.440756405
# 8           f4a_dad_live 1.416890693
# 9         f4a_offr_drink 1.413222993
# 10      f4a_yng_children 1.310763221
# 11         f4a_slp_rooms 1.230513818
# 12          f4a_ms_water 1.200520724
# 13        f4a_disp_feces 1.128384103
# 14         f4a_prim_schl 1.043137176
# 15         f4a_share_fac 1.028085604
# 16        f4b_under_nutr 1.005788928
# 17                  site 0.970187283
# 18    f4a_cur_fastbreath 0.966182483
# 19         f4b_recommend 0.957391444
# 20          f4b_abn_hair 0.944790627
# 21      f4a_drh_prolapse 0.937048701
# 22      f4a_relationship 0.933351523
# 23         f4a_breastfed 0.933104913
# 24       f4a_water_avail 0.925903007
# 25              f4b_skin 0.898059316
# 26           f3_drh_hosp 0.839250831
# 27            f4b_mental 0.834001117
# 28             f3_drh_iv 0.748908429
# 29             f4b_admit 0.741524089
# 30         f3_drh_turgor 0.725170209
# 31        f4a_max_stools 0.711175019
# 32        f4b_skin_flaky 0.682047688
# 33        f4a_trt_method 0.632986343
# 34     f4a_drh_bellypain 0.632400659
# 35          f4a_drh_conv 0.600211501
# 36          f4a_cur_skin 0.598204992
# 37             f4b_mouth 0.536538072
# 38         f4a_drh_vomit 0.536332219
# 39          f4a_wash_def 0.527819664
# 40          f4a_ani_goat 0.522132311
# 41        f4a_hometrt_ab 0.502491657
# 42        f4a_drh_thirst 0.501224416
# 43       f4a_cur_thirsty 0.500026597
# 44         f4a_fac_waste 0.497926550
# 45      f4a_seek_privdoc 0.497222476
# 46         f4a_drh_cough 0.496894911
# 47       f4a_ani_rodents 0.496378931
# 48  f4a_drh_lethrgy_miss 0.486586386
# 49           f4b_bipedal 0.484159619
# 50        f4a_house_elec 0.475550920
# 51     f4a_water_pubwell 0.475026084
# 52        f4a_wash_nurse 0.472724005
# 53         f4a_wash_cook 0.464874344
# 54           f4a_ani_cat 0.459367985
# 55             f3_gender 0.452965609
# 56       f4a_water_river 0.452799482
# 57      f4a_house_agland 0.445686752
# 58     f4a_water_covwell 0.440901755
# 59        f4a_water_bore 0.433059895
# 60      f4a_hometrt_none 0.426760640
# 61       f4a_hometrt_ors 0.420693719
# 62         f4a_ani_sheep 0.419515103
# 63      f4a_seek_outside 0.418253817
# 64       f4a_house_radio 0.417244083
# 65        f4a_house_bike 0.412964907
# 66           f4a_ani_cow 0.411565199
# 67          f4a_wash_eat 0.407472986
# 68        f4a_drh_strain 0.401064789
# 69          f4a_wash_use 0.392136604
# 70        f4a_wash_child 0.383841726
# 71      f4a_hometrt_herb 0.381762423
# 72        f4a_house_tele 0.376008068
# 73     f4a_fuel_charcoal 0.372625119
# 74       f4a_house_phone 0.372358981
# 75       f4a_seek_healer 0.371824616
# 76         f4a_trt_water 0.371317057
# 77          f4a_ani_fowl 0.367232801
# 78     f4a_drh_lessdrink 0.361355519
# 79       f4a_store_water 0.348258383
# 80             f4a_floor 0.340462260
# 81      f4a_cur_restless 0.340239250
# 82      f4a_drh_restless 0.339657743
# 83           f4a_ani_dog 0.338193354
# 84         f4a_ani_other 0.334838457
# 85      f4a_cur_drymouth 0.332349052
# 86        f4a_house_cart 0.325312034
# 87              f4b_eyes 0.323518282
# 88      f4a_water_pubtap 0.321693298
# 89     f4a_hometrt_othr1 0.319501941
# 90   f4a_water_prospring 0.311826870
# 91         f4a_drh_consc 0.309868618
# 92         f4a_fuel_wood 0.309713362
# 93        f4a_water_pond 0.306572552
# 94       f4a_house_scoot 0.304086044
# 95        f4a_water_yard 0.301622109
# 96         f4a_drh_blood 0.299430410
# 97      f4a_water_bought 0.297752821
# 98    f4a_water_deepwell 0.292558743
# 99         f4a_wash_othr 0.283913681
# 100         f4a_seek_doc 0.279067867
# 101        f4a_fuel_kero 0.276091985
# 102  f4a_hometrt_othrliq 0.268967635
# 103      f4a_fuel_biogas 0.256611488
# 104    f4a_hometrt_maize 0.256188230
# 105        f4a_fuel_elec 0.250089385
# 106     f4a_house_fridge 0.245623947
# 107       f4a_fuel_grass 0.241916450
# 108       f4a_seek_other 0.236025865
# 109           f4a_ani_no 0.225583730
# 110           f4b_rectal 0.220871451
# 111       f4a_fuel_other 0.215037647
# 112   f4a_water_covpwell 0.209790603
# 113    f4a_hometrt_othr2 0.206976844
# 114        f4a_house_car 0.182981423
# 115      f4a_fuel_natgas 0.177062474
# 116        f4a_fuel_crop 0.165937422
# 117      f4a_wash_animal 0.152496755
# 118       f4a_water_rain 0.132896787
# 119     f4a_hometrt_zinc 0.127617796
# 120  f4a_water_shallwell 0.113004040
# 121       f4a_water_othr 0.109965685
# 122       f4a_seek_pharm 0.105581647
# 123      f4a_water_house 0.094764403
# 124       f4a_water_well 0.093813875
# 125     f4b_nature_stool 0.088594996
# 126        f4a_fuel_dung 0.078113337
# 127     f4a_hometrt_milk 0.057817574
# 128      f4a_seek_friend 0.057486744
# 129       f4a_house_none 0.005166849
# 130       f4a_seek_remdy 0.005163565
# 131       f4a_house_boat 0.004151401
# 132        f4a_fuel_coal 0.002953594
# 133   f4a_water_unspring 0.001394135
# 134     f4a_fuel_propane 0.001248259
# 135    f4b_observe_stool 0.000000000

# AUC         SE     lower     upper level Model nvar
# 1  0.7650485 0.01672844 0.7322614 0.7978357  0.95    LR    1
# 2  0.8012750 0.01566039 0.7705811 0.8319688  0.95    LR    2
# 3  0.8188817 0.01460026 0.7902657 0.8474977  0.95    LR    3
# 4  0.8186197 0.01481497 0.7895829 0.8476565  0.95    LR    4
# 5  0.8161583 0.01493517 0.7868859 0.8454307  0.95    LR    5
# 6  0.8089597 0.01533762 0.7788985 0.8390208  0.95    LR    6
# 7  0.8062341 0.01571015 0.7754428 0.8370255  0.95    LR    7
# 8  0.8123044 0.01544725 0.7820284 0.8425805  0.95    LR    8
# 9  0.8129850 0.01587721 0.7818663 0.8441038  0.95    LR    9
# 10 0.8243153 0.01500431 0.7949074 0.8537232  0.95    LR   10
# 11 0.6426416 0.03163097 0.5806460 0.7046371  0.95    RF    1
# 12 0.7485024 0.01694027 0.7153001 0.7817047  0.95    RF    2
# 13 0.7761941 0.01663614 0.7435879 0.8088003  0.95    RF    3
# 14 0.7773573 0.01746014 0.7431360 0.8115785  0.95    RF    4
# 15 0.7689372 0.01665862 0.7362869 0.8015875  0.95    RF    5
# 16 0.7674749 0.01769223 0.7327987 0.8021510  0.95    RF    6
# 17 0.7662475 0.01720834 0.7325197 0.7999752  0.95    RF    7
# 18 0.7735256 0.01672133 0.7407524 0.8062988  0.95    RF    8
# 19 0.7833026 0.01751341 0.7489770 0.8176283  0.95    RF    9
# 20 0.7961453 0.01705809 0.7627121 0.8295786  0.95    RF   10

# nvar      intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>     <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 -0.000534   -0.491    0.429 0.945     0.531      1.37
# 2     2  0.000533   -0.497    0.439 0.956     0.594      1.34
# 3     3  0.0165     -0.482    0.456 0.984     0.631      1.36
# 4     4  0.0169     -0.481    0.456 0.985     0.632      1.37
# 5     5  0.0146     -0.485    0.455 0.969     0.620      1.34
# 6     6  0.0118     -0.494    0.458 0.904     0.578      1.25
# 7     7  0.0143     -0.494    0.463 0.862     0.546      1.20
# 8     8  0.0288     -0.482    0.480 0.865     0.554      1.20
# 9     9  0.0394     -0.474    0.493 0.850     0.547      1.18
# 10    10  0.0378     -0.479    0.495 0.829     0.534      1.15

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- main.011[["AUC_df"]]
#jpeg("/DeathAll_AUCs_GEMS011.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(main.011[["nvars_opts"]])[1]],AUC_df$AUC[1:length(main.011[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="All Deaths in cases 0-11mo in GEMS",
     #main=expression(paste("">=0.5,Delta,"HAZ in cases 0-59mo in GEMS")),
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(main.011[["nvars_opts"]])[1]],AUC_df$AUC[(length(main.011[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()



################### main + 12-23(age_2) any death all cases -HAZ +MUAC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

main.1223 <- CPR.funct(data=complete_anydeath_age2,outcome="death_all",iter=100,nvars_opts=c(1:10))
main.1223[["df_imps"]]
main.1223[["AUC_df"]]
main.1223[["calib"]]

# names      var_red
# 1               f4b_muac 3.484889e+00
# 2               f4b_temp 1.350539e+00
# 3          f4a_ppl_house 1.308293e+00
# 4               f4b_resp 1.156675e+00
# 5       f4a_yng_children 1.027906e+00
# 6           f4a_drh_days 9.668063e-01
# 7               base_age 8.981684e-01
# 8           f4b_abn_hair 8.303451e-01
# 9          f4a_slp_rooms 7.428259e-01
# 10      f4a_hometrt_herb 7.394404e-01
# 11      f4a_hometrt_milk 6.946596e-01
# 12             f3_drh_iv 6.653311e-01
# 13        f4a_offr_drink 6.293029e-01
# 14        f4b_under_nutr 6.169620e-01
# 15        f4a_disp_feces 5.633124e-01
# 16          f4a_dad_live 5.611145e-01
# 17              f4b_skin 5.391561e-01
# 18         f4a_drh_consc 4.972478e-01
# 19       f4a_water_avail 4.689189e-01
# 20                  site 4.598465e-01
# 21         f3_drh_turgor 4.471309e-01
# 22           f3_drh_hosp 4.261278e-01
# 23         f4a_breastfed 4.258462e-01
# 24         f4b_recommend 3.921551e-01
# 25         f4a_prim_schl 3.878498e-01
# 26        f4a_water_bore 3.874292e-01
# 27      f4a_relationship 3.573349e-01
# 28     f4a_drh_bellypain 3.558800e-01
# 29        f4a_trt_method 3.363019e-01
# 30            f4b_mental 3.342809e-01
# 31          f4a_ms_water 3.241314e-01
# 32        f4a_house_elec 3.205432e-01
# 33             f4b_admit 3.071714e-01
# 34  f4a_drh_lethrgy_miss 2.981460e-01
# 35     f4a_water_pubwell 2.933921e-01
# 36          f4a_cur_skin 2.924930e-01
# 37          f4a_wash_def 2.914903e-01
# 38             f4b_mouth 2.875545e-01
# 39     f4a_water_covwell 2.833083e-01
# 40     f4a_drh_lessdrink 2.711555e-01
# 41          f4a_ani_goat 2.564612e-01
# 42         f4a_drh_vomit 2.526112e-01
# 43       f4a_house_radio 2.481159e-01
# 44        f4a_water_well 2.426211e-01
# 45     f4a_hometrt_othr1 2.366691e-01
# 46      f4b_nature_stool 2.364490e-01
# 47       f4a_hometrt_ors 2.361095e-01
# 48        f4b_skin_flaky 2.347880e-01
# 49           f4b_bipedal 2.328379e-01
# 50       f4a_ani_rodents 2.314766e-01
# 51      f4a_hometrt_none 2.281126e-01
# 52        f4a_seek_pharm 2.267068e-01
# 53      f4a_cur_restless 2.250559e-01
# 54         f4a_wash_cook 2.235351e-01
# 55         f4a_drh_cough 2.230514e-01
# 56         f4a_ani_other 2.148119e-01
# 57          f4a_wash_use 2.035864e-01
# 58       f4a_seek_healer 2.029641e-01
# 59      f4a_house_agland 2.009529e-01
# 60             f4a_floor 1.976905e-01
# 61      f4a_seek_outside 1.946568e-01
# 62         f4a_share_fac 1.925302e-01
# 63         f4a_drh_blood 1.924575e-01
# 64             f3_gender 1.914262e-01
# 65    f4a_cur_fastbreath 1.865130e-01
# 66        f4a_wash_nurse 1.825435e-01
# 67        f4a_wash_child 1.819245e-01
# 68           f4a_ani_cow 1.802721e-01
# 69      f4a_water_pubtap 1.799619e-01
# 70    f4a_water_covpwell 1.784325e-01
# 71        f4a_max_stools 1.780303e-01
# 72         f4a_trt_water 1.772148e-01
# 73        f4a_house_tele 1.738361e-01
# 74        f4a_house_cart 1.703556e-01
# 75      f4a_drh_prolapse 1.665377e-01
# 76         f4a_ani_sheep 1.626639e-01
# 77      f4a_drh_restless 1.618437e-01
# 78       f4a_house_phone 1.583771e-01
# 79     f4a_fuel_charcoal 1.565299e-01
# 80       f4a_house_scoot 1.517810e-01
# 81       f4b_chest_indrw 1.511391e-01
# 82        f4a_house_bike 1.479113e-01
# 83           f4a_ani_dog 1.472380e-01
# 84        f4a_drh_thirst 1.441173e-01
# 85      f4a_cur_drymouth 1.422792e-01
# 86              f4b_eyes 1.416347e-01
# 87          f4a_ani_fowl 1.413788e-01
# 88       f4a_cur_thirsty 1.410548e-01
# 89           f4a_ani_cat 1.378857e-01
# 90            f4b_rectal 1.324595e-01
# 91       f4a_fuel_biogas 1.324168e-01
# 92        f4a_water_rain 1.268766e-01
# 93    f4a_water_deepwell 1.227461e-01
# 94     f4a_hometrt_maize 1.219033e-01
# 95          f4a_wash_eat 1.163537e-01
# 96    f4a_water_unspring 1.081212e-01
# 97        f4a_water_yard 9.494818e-02
# 98            f4a_ani_no 8.971647e-02
# 99       f4a_store_water 8.964127e-02
# 100       f4a_house_none 8.865105e-02
# 101  f4a_water_prospring 8.561354e-02
# 102        f4a_wash_othr 7.995833e-02
# 103        f4a_fuel_wood 7.903940e-02
# 104     f4a_house_fridge 7.663167e-02
# 105      f4a_water_river 7.306922e-02
# 106        f4a_fac_waste 6.536313e-02
# 107  f4a_hometrt_othrliq 6.287067e-02
# 108         f4a_drh_conv 6.076698e-02
# 109         f4a_seek_doc 5.967438e-02
# 110      f4a_fuel_natgas 5.877370e-02
# 111     f4a_water_bought 4.651101e-02
# 112        f4a_house_car 4.650688e-02
# 113     f4a_hometrt_zinc 4.412310e-02
# 114     f4a_seek_privdoc 4.352621e-02
# 115       f4a_water_pond 4.087204e-02
# 116       f4a_drh_strain 3.454505e-02
# 117      f4a_wash_animal 2.794433e-02
# 118       f4a_hometrt_ab 2.785004e-02
# 119      f4a_water_house 8.081658e-03
# 120       f4a_fuel_grass 7.905708e-03
# 121        f4a_fuel_crop 5.683387e-03
# 122        f4a_fuel_kero 4.954045e-03
# 123    f4a_hometrt_othr2 3.977254e-03
# 124       f4a_seek_other 3.127273e-03
# 125        f4a_fuel_dung 2.067997e-03
# 126  f4a_water_shallwell 1.438105e-03
# 127       f4a_house_boat 1.428571e-03
# 128       f4a_fuel_other 1.108802e-03
# 129       f4a_seek_remdy 1.035409e-03
# 130      f4a_seek_friend 6.523810e-04
# 131       f4a_water_othr 1.000000e-04
# 132        f4a_fuel_coal 3.333333e-05
# 133        f4a_fuel_elec 0.000000e+00
# 134     f4a_fuel_propane 0.000000e+00
# 135    f4b_observe_stool 0.000000e+00

# AUC         SE     lower     upper level Model nvar
# 1  0.7984606 0.02265193 0.7540637 0.8428576  0.95    LR    1
# 2  0.8019627 0.02222991 0.7583929 0.8455325  0.95    LR    2
# 3  0.8024847 0.02134362 0.7606519 0.8443174  0.95    LR    3
# 4  0.8068298 0.02101184 0.7656474 0.8480123  0.95    LR    4
# 5  0.8126056 0.02073666 0.7719625 0.8532487  0.95    LR    5
# 6  0.8158530 0.02059867 0.7754804 0.8562257  0.95    LR    6
# 7  0.8215597 0.02007274 0.7822178 0.8609015  0.95    LR    7
# 8  0.8259953 0.02019397 0.7864158 0.8655747  0.95    LR    8
# 9  0.8347924 0.01971703 0.7961478 0.8734371  0.95    LR    9
# 10 0.8453740 0.01916663 0.8078081 0.8829399  0.95    LR   10
# 11 0.6749678 0.05074265 0.5755140 0.7744215  0.95    RF    1
# 12 0.7146366 0.03325164 0.6494645 0.7798086  0.95    RF    2
# 13 0.7552448 0.02590932 0.7044635 0.8060262  0.95    RF    3
# 14 0.7700174 0.02803282 0.7150741 0.8249607  0.95    RF    4
# 15 0.7762663 0.02476468 0.7277284 0.8248042  0.95    RF    5
# 16 0.7975605 0.02180239 0.7548286 0.8402924  0.95    RF    6
# 17 0.8077531 0.02003207 0.7684910 0.8470152  0.95    RF    7
# 18 0.8270457 0.01749157 0.7927629 0.8613286  0.95    RF    8
# 19 0.8219994 0.02039546 0.7820250 0.8619738  0.95    RF    9
# 20 0.8352913 0.02069690 0.7947262 0.8758565  0.95    RF   10

# nvar   intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>  <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 0.0530   -0.690    0.665 1.01      0.506      1.56
# 2     2 0.0524   -0.691    0.664 1.02      0.517      1.58
# 3     3 0.0659   -0.677    0.677 1.02      0.512      1.57
# 4     4 0.0733   -0.672    0.687 0.978     0.484      1.51
# 5     5 0.0553   -0.696    0.676 0.895     0.438      1.38
# 6     6 0.0452   -0.709    0.670 0.880     0.435      1.35
# 7     7 0.0928   -0.661    0.718 0.846     0.416      1.30
# 8     8 0.0849   -0.678    0.720 0.806     0.404      1.23
# 9     9 0.0825   -0.687    0.725 0.772     0.394      1.17
# 10    10 0.0907   -0.680    0.738 0.751     0.390      1.14

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- main.1223[["AUC_df"]]
#jpeg("/DeathAll_AUCs_GEMS1223.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(main.1223[["nvars_opts"]])[1]],AUC_df$AUC[1:length(main.1223[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="All Deaths in cases 12-23mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(main.1223[["nvars_opts"]])[1]],AUC_df$AUC[(length(main.1223[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()

################### main + 24-59(age_3) any death all cases -HAZ +MUAC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

#trying to find why sometimes has convergence issues
#dataset of only predictor variables considering
temp<-complete_anydeath_age3[names]
#table all of these
temp %>% map( table )

#these all have small cell sizes
# $f4a_relationship
# $f4a_prim_schl
# $f4a_house_boat
# $f4a_water_covpwell
# $f4a_water_unspring
# $f4a_water_avail
# $f4a_trt_method
# $f4a_breastfed
# $f4a_drh_blood
# $f4a_drh_vomit
# $f4a_drh_lessdrink
# $f4a_drh_strain
# $f4a_drh_prolapse
# $f4a_cur_drymouth
# $f4a_cur_fastbreath
# $f4a_hometrt_milk
# $f4a_hometrt_othrliq
# $f4a_offr_drink
# $f4a_seek_friend
# $f4b_bipedal
# $f4b_recommend
# $f4a_fuel_elec
# $f4a_house_none
# $f4a_water_prospring
# $f4a_water_well
# $f4a_water_pond
# $f4a_water_othr
# $f4a_seek_other
# $f4b_chest_indrw
# $f4b_rectal
# $f4b_skin_flaky

complete_anydeath_age3_complete <- complete_anydeath_age3
summary(complete_anydeath_age3_complete$f4a_disp_feces)
#combine 1, 4, 6
complete_anydeath_age3_complete$f4a_disp_feces <- 
  as.factor(ifelse((complete_anydeath_age3_complete$f4a_disp_feces!=2 & complete_anydeath_age3_complete$f4a_disp_feces!=3),
  1,complete_anydeath_age3_complete$f4a_disp_feces))
summary(complete_anydeath_age3_complete$f4a_disp_feces)

summary(complete_anydeath_age3_complete$f4a_fac_waste)
#drop 5,6
complete_anydeath_age3_complete <- complete_anydeath_age3_complete %>% filter(f4a_fac_waste!=6)
summary(complete_anydeath_age3_complete$f4a_fac_waste)

#combine 3,4
summary(complete_anydeath_age3_complete$f4a_wash_use)
complete_anydeath_age3_complete$f4a_wash_use <- 
  as.factor(ifelse((complete_anydeath_age3_complete$f4a_wash_use==4),
                   3,complete_anydeath_age3_complete$f4a_wash_use))
summary(complete_anydeath_age3_complete$f4a_wash_use)

#make 9 into 0
summary(complete_anydeath_age3_complete$f4a_drh_thirst)
complete_anydeath_age3_complete$f4a_drh_thirst <- 
  (ifelse((complete_anydeath_age3_complete$f4a_drh_thirst==9 | complete_anydeath_age3_complete$f4a_drh_thirst==0),
                   0,complete_anydeath_age3_complete$f4a_drh_thirst))
complete_anydeath_age3_complete$f4a_drh_thirst <- 
  (ifelse((complete_anydeath_age3_complete$f4a_drh_thirst==2),
          1,complete_anydeath_age3_complete$f4a_drh_thirst))
complete_anydeath_age3_complete$f4a_drh_thirst <- as.factor(complete_anydeath_age3_complete$f4a_drh_thirst)
table(complete_anydeath_age3_complete$f4a_drh_thirst)
summary(complete_anydeath_age3_complete$f4a_drh_thirst)

#make 9 into 0
summary(complete_anydeath_age3_complete$f4a_cur_thirsty)
complete_anydeath_age3_complete$f4a_cur_thirsty <- 
  (ifelse((complete_anydeath_age3_complete$f4a_cur_thirsty==9 | complete_anydeath_age3_complete$f4a_cur_thirsty==0),
          0,complete_anydeath_age3_complete$f4a_cur_thirsty))
complete_anydeath_age3_complete$f4a_cur_thirsty <- 
  (ifelse((complete_anydeath_age3_complete$f4a_cur_thirsty==2),
          1,complete_anydeath_age3_complete$f4a_cur_thirsty))
complete_anydeath_age3_complete$f4a_cur_thirsty <- as.factor(complete_anydeath_age3_complete$f4a_cur_thirsty)
table(complete_anydeath_age3_complete$f4a_cur_thirsty)
summary(complete_anydeath_age3_complete$f4a_cur_thirsty)

#make 9 into 0
summary(complete_anydeath_age3_complete$f4a_cur_skin)
complete_anydeath_age3_complete$f4a_cur_skin <- 
  (ifelse((complete_anydeath_age3_complete$f4a_cur_skin==9 | complete_anydeath_age3_complete$f4a_cur_skin==0),
          0,complete_anydeath_age3_complete$f4a_cur_skin))
complete_anydeath_age3_complete$f4a_cur_skin <- 
  (ifelse((complete_anydeath_age3_complete$f4a_cur_skin==2),
          1,complete_anydeath_age3_complete$f4a_cur_skin))
complete_anydeath_age3_complete$f4a_cur_skin <- as.factor(complete_anydeath_age3_complete$f4a_cur_skin)
table(complete_anydeath_age3_complete$f4a_cur_skin)
summary(complete_anydeath_age3_complete$f4a_cur_skin)

#combine 4 and 5
summary(complete_anydeath_age3_complete$f4a_dad_live)
complete_anydeath_age3_complete$f4a_dad_live <- 
  as.factor(ifelse((complete_anydeath_age3_complete$f4a_dad_live==4 | complete_anydeath_age3_complete$f4a_dad_live==5),
                   4,complete_anydeath_age3_complete$f4a_dad_live))
summary(complete_anydeath_age3_complete$f4a_dad_live)

summary(complete_anydeath_age3_complete$f4a_ms_water)
#drop 4
complete_anydeath_age3_complete <- complete_anydeath_age3_complete %>% filter(f4a_ms_water!=4)
summary(complete_anydeath_age3_complete$f4a_ms_water)


names <- names[!names %in% c("f4b_haz","f4a_relationship","f4a_prim_schl",
                             "f4a_house_boat","f4a_water_covpwell","f4a_water_unspring",
                             "f4a_water_avail","f4a_trt_method","f4a_breastfed",
                             "f4a_drh_blood","f4a_drh_vomit","f4a_drh_lessdrink",
                             "f4a_drh_strain","f4a_drh_prolapse","f4a_cur_drymouth",
                             "f4a_cur_fastbreath","f4a_hometrt_milk","f4a_hometrt_othrliq",
                             "f4a_offr_drink","f4a_seek_friend","f4b_bipedal",
                             "f4b_recommend","f4a_fuel_elec","f4a_house_none",
                             "f4a_water_prospring","f4a_water_well","f4a_water_pond",
                             "f4a_water_othr","f4a_seek_other","f4b_chest_indrw",
                             "f4b_rectal","f4b_skin_flaky")]

main.2459 <- CPR.funct(data=complete_anydeath_age3_complete,outcome="death_all",iter=100,nvars_opts=c(1:10))
main.2459[["df_imps"]]
main.2459[["AUC_df"]]
main.2459[["calib"]]

# names      var_red
# 1               f4b_muac 2.4777221706
# 2           f4b_abn_hair 0.8780525373
# 3         f4b_skin_flaky 0.6803326255
# 4               f4b_resp 0.6055015762
# 5               f4b_temp 0.5370934615
# 6               base_age 0.4955140551
# 7           f4a_drh_days 0.4138605580
# 8         f4b_under_nutr 0.3682672931
# 9         f4a_offr_drink 0.3599162930
# 10      f4a_relationship 0.2732409195
# 11              f4b_skin 0.2650669219
# 12         f4a_ppl_house 0.2627140450
# 13          f4a_dad_live 0.2434376172
# 14                  site 0.2344008452
# 15         f4a_slp_rooms 0.2118071011
# 16          f4a_cur_skin 0.2092618900
# 17             f3_drh_iv 0.2055406425
# 18           f4b_bipedal 0.2016680343
# 19           f3_drh_hosp 0.1989540700
# 20            f4b_mental 0.1885224896
# 21             f4b_admit 0.1864698024
# 22         f4b_recommend 0.1859587374
# 23      f4a_yng_children 0.1799212985
# 24           f4a_ani_cow 0.1763802409
# 25         f4a_share_fac 0.1705844018
# 26        f4a_disp_feces 0.1672792630
# 27          f4a_ms_water 0.1536872683
# 28         f4a_prim_schl 0.1460652244
# 29         f3_drh_turgor 0.1393366912
# 30     f4a_drh_bellypain 0.1287018423
# 31        f4a_max_stools 0.1191471440
# 32       f4a_seek_healer 0.1091327080
# 33             f4a_floor 0.1042286592
# 34      f4a_hometrt_herb 0.1040882085
# 35         f4a_drh_cough 0.1027488026
# 36        f4a_house_bike 0.1012880908
# 37     f4a_fuel_charcoal 0.0980392762
# 38             f4b_mouth 0.0979956652
# 39       f4a_water_avail 0.0957441438
# 40          f4a_drh_conv 0.0945939869
# 41        f4a_house_elec 0.0915987842
# 42  f4a_drh_lethrgy_miss 0.0910443700
# 43         f4a_breastfed 0.0897802392
# 44      f4a_house_agland 0.0896648543
# 45      f4a_drh_restless 0.0896208668
# 46          f4a_wash_use 0.0884476675
# 47         f4a_drh_consc 0.0882725445
# 48         f4a_drh_vomit 0.0866442666
# 49       f4a_water_river 0.0857878974
# 50        f4a_water_bore 0.0844787275
# 51         f4a_ani_other 0.0840272171
# 52       f4a_house_phone 0.0837937719
# 53      f4b_nature_stool 0.0823551183
# 54        f4a_trt_method 0.0798910429
# 55       f4a_ani_rodents 0.0781300922
# 56      f4a_water_pubtap 0.0776773581
# 57      f4a_hometrt_none 0.0756083481
# 58     f4a_hometrt_maize 0.0741852227
# 59          f4a_ani_goat 0.0726360661
# 60         f4a_wash_cook 0.0703901823
# 61        f4a_wash_child 0.0691522080
# 62           f4a_ani_dog 0.0687295160
# 63          f4a_ani_fowl 0.0667777882
# 64          f4a_wash_def 0.0644943497
# 65      f4a_cur_restless 0.0635321687
# 66    f4a_cur_fastbreath 0.0630549973
# 67        f4a_wash_nurse 0.0630487882
# 68        f4a_water_rain 0.0620847755
# 69        f4a_house_cart 0.0617467452
# 70      f4a_seek_outside 0.0614896248
# 71       f4a_house_radio 0.0614866860
# 72     f4a_drh_lessdrink 0.0595886101
# 73         f4a_ani_sheep 0.0563202391
# 74        f4a_house_tele 0.0544851352
# 75       f4a_hometrt_ors 0.0517929521
# 76        f4a_drh_thirst 0.0470860914
# 77           f4a_ani_cat 0.0466522688
# 78       f4a_cur_thirsty 0.0460660298
# 79         f4a_drh_blood 0.0445790605
# 80         f4a_fuel_wood 0.0439105634
# 81         f4a_fac_waste 0.0431246386
# 82         f4a_trt_water 0.0424670628
# 83        f4a_seek_remdy 0.0407578368
# 84             f3_gender 0.0406126545
# 85    f4a_water_deepwell 0.0391951691
# 86      f4a_cur_drymouth 0.0362570106
# 87          f4a_wash_eat 0.0359330467
# 88       f4a_store_water 0.0349518060
# 89       f4a_fuel_biogas 0.0345365269
# 90       f4a_wash_animal 0.0341591362
# 91      f4a_house_fridge 0.0339051589
# 92         f4a_house_car 0.0332062509
# 93       f4a_house_scoot 0.0305823581
# 94              f4b_eyes 0.0301001068
# 95        f4a_seek_other 0.0280304141
# 96        f4a_seek_pharm 0.0267322518
# 97     f4a_hometrt_othr1 0.0266411098
# 98            f4a_ani_no 0.0223482216
# 99        f4a_drh_strain 0.0218127885
# 100       f4a_water_pond 0.0214383190
# 101        f4a_wash_othr 0.0207039233
# 102    f4a_hometrt_othr2 0.0204739817
# 103        f4a_fuel_dung 0.0182982641
# 104     f4a_water_bought 0.0169839361
# 105       f4a_hometrt_ab 0.0143403869
# 106      f4a_water_house 0.0126439834
# 107       f4a_fuel_grass 0.0119650676
# 108        f4a_fuel_crop 0.0116820931
# 109      f4a_fuel_natgas 0.0071356218
# 110       f4a_water_yard 0.0066327935
# 111  f4a_water_shallwell 0.0062957234
# 112    f4a_water_pubwell 0.0056433913
# 113        f4a_fuel_coal 0.0040000000
# 114        f4a_fuel_kero 0.0033405918
# 115     f4a_seek_privdoc 0.0032170635
# 116     f4a_drh_prolapse 0.0027014733
# 117      f4b_chest_indrw 0.0026269841
# 118         f4a_seek_doc 0.0023554573
# 119   f4a_water_covpwell 0.0023225957
# 120     f4a_hometrt_zinc 0.0018547619
# 121     f4a_fuel_propane 0.0017476190
# 122       f4a_fuel_other 0.0017083333
# 123       f4a_water_well 0.0015952381
# 124    f4a_water_covwell 0.0012168498
# 125       f4a_water_othr 0.0008181818
# 126           f4b_rectal 0.0003992923
# 127       f4a_house_boat 0.0000000000
# 128       f4a_house_none 0.0000000000
# 129        f4a_fuel_elec 0.0000000000
# 130  f4a_water_prospring 0.0000000000
# 131   f4a_water_unspring 0.0000000000
# 132     f4a_hometrt_milk 0.0000000000
# 133  f4a_hometrt_othrliq 0.0000000000
# 134      f4a_seek_friend 0.0000000000
# 135    f4b_observe_stool 0.0000000000

# AUC         SE     lower     upper level Model nvar
# 1  0.9132498 0.02017935 0.8736990 0.9528006  0.95    LR    1
# 2  0.9125029 0.02031271 0.8726908 0.9523151  0.95    LR    2
# 3  0.9128997 0.02014681 0.8734127 0.9523867  0.95    LR    3
# 4  0.9258366 0.01815123 0.8902608 0.9614123  0.95    LR    4
# 5  0.9262544 0.01804874 0.8908795 0.9616293  0.95    LR    5
# 6  0.9221260 0.01902097 0.8848456 0.9594064  0.95    LR    6
# 7  0.9203438 0.01902220 0.8830610 0.9576266  0.95    LR    7
# 8  0.9248501 0.01913535 0.8873455 0.9623547  0.95    LR    8
# 9  0.9275648 0.01908076 0.8901672 0.9649624  0.95    LR    9
# 10 0.9046080 0.02420601 0.8571651 0.9520509  0.95    LR   10
# 11 0.7141175 0.08115685 0.5550530 0.8731820  0.95    RF    1
# 12 0.8468977 0.05058806 0.7477470 0.9460485  0.95    RF    2
# 13 0.8278816 0.04983165 0.7302134 0.9255498  0.95    RF    3
# 14 0.8528594 0.04073112 0.7730278 0.9326909  0.95    RF    4
# 15 0.8683921 0.03051571 0.8085824 0.9282018  0.95    RF    5
# 16 0.8578396 0.03177267 0.7955663 0.9201128  0.95    RF    6
# 17 0.8572680 0.03321491 0.7921680 0.9223680  0.95    RF    7
# 18 0.8473881 0.03463057 0.7795134 0.9152628  0.95    RF    8
# 19 0.8575808 0.04676132 0.7659303 0.9492313  0.95    RF    9
# 20 0.8463527 0.04434418 0.7594397 0.9332657  0.95    RF   10

# nvar   intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>  <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 0.0564    -1.46     1.12 0.981     0.412      1.80
# 2     2 0.137     -1.47     1.27 0.854     0.357      1.54
# 3     3 0.122     -1.50     1.28 0.841     0.356      1.49
# 4     4 0.118     -1.51     1.29 0.837     0.365      1.48
# 5     5 0.0953    -1.54     1.27 0.827     0.363      1.48
# 6     6 0.0369    -1.63     1.24 0.757     0.333      1.36
# 7     7 0.120     -1.56     1.32 0.754     0.325      1.37
# 8     8 0.141     -1.56     1.36 0.776     0.334      1.50
# 9     9 0.140     -1.57     1.37 0.772     0.340      1.48
# 10    10 0.187     -1.50     1.42 0.688     0.280      1.31

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- main.2459[["AUC_df"]]
#jpeg("/DeathAll_AUCs_GEMS2459.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(main.2459[["nvars_opts"]])[1]],AUC_df$AUC[1:length(main.2459[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="All Deaths in cases 24-59mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(main.2459[["nvars_opts"]])[1]],AUC_df$AUC[(length(main.2459[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()

################### main (0-59mo) any death -HAZ + MUAC by country - Gambia ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz","site")]

complete_anydeath_Gambia_complete <- complete_anydeath_Gambia
#f4a_relationship recode everything that's not 1 to 2
summary(complete_anydeath_Gambia$f4a_relationship)
complete_anydeath_Gambia_complete$f4a_relationship <- as.factor(ifelse(complete_anydeath_Gambia_complete$f4a_relationship!=1,2,complete_anydeath_Gambia_complete$f4a_relationship))
summary(complete_anydeath_Gambia_complete$f4a_relationship)


# #f4a_wash_use recode 3 into 2
# summary(complete_anydeath_Gambia$f4a_wash_use)
# complete_anydeath_Gambia_complete$f4a_wash_use <- as.factor(ifelse(complete_anydeath_Gambia_complete$f4a_wash_use==3,2,complete_anydeath_Gambia_complete$f4a_wash_use))
# summary(complete_anydeath_Gambia_complete$f4a_wash_use)
# 
# #f4a_drh_blood drop 9 for missing
# summary(complete_anydeath_Gambia$f4a_drh_blood)
# complete_anydeath_Gambia_complete <- complete_anydeath_Gambia_complete %>% filter(f4a_drh_blood!=9) 
# summary(complete_anydeath_Gambia_complete$f4a_drh_blood)
# 
# #f4b_recommend recode 2 and 3 to 0
# summary(complete_anydeath_Gambia_complete$f4b_recommend)
# complete_anydeath_Gambia_complete$f4b_recommend <- as.factor(ifelse(complete_anydeath_Gambia_complete$f4b_recommend==0,0,
#                                 ifelse(complete_anydeath_Gambia_complete$f4b_recommend==2,0,
#                                 ifelse(complete_anydeath_Gambia_complete$f4b_recommend==3,0,
#                                 ifelse(complete_anydeath_Gambia_complete$f4b_recommend==1,1,4)))))
# table(complete_anydeath_Gambia_complete$f4b_recommend)
# 
# #f4a_disp_feces
# summary(complete_anydeath_Gambia_complete$f4a_disp_feces )


main.Gambia <- CPR.funct(data=complete_anydeath_Gambia_complete,outcome="death_all",iter=100,nvars_opts=c(1:10))
main.Gambia[["df_imps"]]
main.Gambia[["AUC_df"]]
main.Gambia[["calib"]]

#save(main.Gambia, file = "/main.Gambia.100iter.callib.Rdata")

# names     var_red
# 1               f4b_muac 3.452427229
# 2               f4b_resp 1.491006537
# 3               base_age 1.168131414
# 4           f4a_drh_days 0.964646725
# 5          f4a_ppl_house 0.884802460
# 6               f4b_temp 0.812596600
# 7          f4a_slp_rooms 0.782441630
# 8       f4a_yng_children 0.759018596
# 9             f4b_mental 0.685251261
# 10       f4b_chest_indrw 0.579453748
# 11        f4b_under_nutr 0.478801708
# 12      f4a_hometrt_herb 0.455162713
# 13              f4b_skin 0.414883187
# 14       f4a_water_avail 0.392993196
# 15        f4b_skin_flaky 0.379180580
# 16       f4a_seek_healer 0.371049667
# 17        f4a_offr_drink 0.370178937
# 18             f4b_admit 0.353601216
# 19          f4a_cur_skin 0.351573473
# 20             f3_drh_iv 0.347991142
# 21      f4a_drh_prolapse 0.341035336
# 22         f3_drh_turgor 0.340058740
# 23         f4b_recommend 0.329835471
# 24         f4a_drh_blood 0.328166057
# 25          f4b_abn_hair 0.327239169
# 26     f4a_drh_bellypain 0.319821151
# 27      f4a_hometrt_none 0.312723275
# 28           f3_drh_hosp 0.302252943
# 29          f4a_dad_live 0.298259881
# 30         f4a_drh_cough 0.294845724
# 31        f4a_disp_feces 0.283176089
# 32       f4a_hometrt_ors 0.264153133
# 33         f4a_wash_cook 0.259481082
# 34        f4a_seek_pharm 0.253141712
# 35         f4a_prim_schl 0.247731212
# 36             f4b_mouth 0.244317309
# 37     f4a_hometrt_othr1 0.240722934
# 38      f4a_seek_outside 0.238319063
# 39           f4a_ani_cow 0.226314715
# 40          f4a_wash_use 0.225028063
# 41     f4a_hometrt_maize 0.215424939
# 42     f4a_drh_lessdrink 0.210332083
# 43        f4a_wash_nurse 0.208346557
# 44         f4a_breastfed 0.206987574
# 45    f4a_water_deepwell 0.200957130
# 46      f4a_relationship 0.197883070
# 47          f4a_ms_water 0.197703420
# 48       f4a_cur_thirsty 0.195644507
# 49       f4a_house_scoot 0.193436859
# 50     f4a_water_pubwell 0.190201737
# 51      f4a_cur_restless 0.189003788
# 52          f4a_seek_doc 0.187005725
# 53    f4a_cur_fastbreath 0.183247933
# 54         f4a_drh_vomit 0.182172947
# 55    f4a_water_covpwell 0.176771710
# 56        f4a_house_elec 0.174785866
# 57      f4a_drh_restless 0.172566215
# 58       f4a_ani_rodents 0.170704207
# 59        f4a_drh_thirst 0.170554562
# 60           f4a_ani_cat 0.166373964
# 61       f4a_house_radio 0.159247180
# 62             f4a_floor 0.158997103
# 63         f4a_ani_other 0.158400614
# 64  f4a_drh_lethrgy_miss 0.155961679
# 65      f4a_cur_drymouth 0.149041535
# 66         f4a_trt_water 0.147088394
# 67         f4a_ani_sheep 0.146252054
# 68             f3_gender 0.145911496
# 69       f4a_store_water 0.145538810
# 70        f4a_house_bike 0.137549815
# 71        f4a_trt_method 0.135143415
# 72          f4a_wash_def 0.130129320
# 73      f4a_water_pubtap 0.129945324
# 74        f4a_water_othr 0.129250903
# 75          f4a_ani_goat 0.127193805
# 76        f4a_wash_child 0.124012912
# 77     f4a_hometrt_othr2 0.122835135
# 78     f4a_water_covwell 0.118137456
# 79        f4a_max_stools 0.118013802
# 80              f4b_eyes 0.102407389
# 81           f4a_ani_dog 0.101347366
# 82      f4a_house_fridge 0.100981434
# 83        f4a_house_tele 0.094942246
# 84        f4a_house_cart 0.094923487
# 85       f4a_house_phone 0.078396732
# 86        f4a_water_yard 0.077054909
# 87         f4a_house_car 0.070674303
# 88          f4a_drh_conv 0.070140336
# 89           f4b_bipedal 0.068178138
# 90        f4a_seek_other 0.067788630
# 91   f4a_hometrt_othrliq 0.060300220
# 92          f4a_ani_fowl 0.059757341
# 93       f4a_seek_friend 0.053621771
# 94            f4b_rectal 0.048413081
# 95      f4a_house_agland 0.044842569
# 96         f4a_share_fac 0.042901875
# 97        f4a_water_well 0.036135307
# 98          f4a_wash_eat 0.035820572
# 99       f4a_wash_animal 0.034781242
# 100       f4a_drh_strain 0.024361294
# 101       f4a_hometrt_ab 0.020498741
# 102        f4a_fac_waste 0.013965460
# 103        f4a_wash_othr 0.009401854
# 104     f4a_hometrt_zinc 0.004074603
# 105       f4a_seek_remdy 0.002791667
# 106     f4a_seek_privdoc 0.002497619
# 107           f4a_ani_no 0.002245238
# 108       f4a_water_bore 0.002211111
# 109        f4a_drh_consc 0.001777778
# 110    f4a_fuel_charcoal 0.001166667
# 111        f4a_fuel_coal 0.000400000
# 112                 site 0.000000000
# 113       f4a_house_boat 0.000000000
# 114       f4a_house_none 0.000000000
# 115        f4a_fuel_elec 0.000000000
# 116      f4a_fuel_biogas 0.000000000
# 117       f4a_fuel_grass 0.000000000
# 118     f4a_fuel_propane 0.000000000
# 119        f4a_fuel_dung 0.000000000
# 120      f4a_fuel_natgas 0.000000000
# 121        f4a_fuel_crop 0.000000000
# 122        f4a_fuel_kero 0.000000000
# 123        f4a_fuel_wood 0.000000000
# 124       f4a_fuel_other 0.000000000
# 125      f4a_water_house 0.000000000
# 126  f4a_water_prospring 0.000000000
# 127   f4a_water_unspring 0.000000000
# 128      f4a_water_river 0.000000000
# 129       f4a_water_pond 0.000000000
# 130       f4a_water_rain 0.000000000
# 131  f4a_water_shallwell 0.000000000
# 132     f4a_water_bought 0.000000000
# 133     f4a_hometrt_milk 0.000000000
# > main.Gambia[["AUC_df"]]
# AUC          SE     lower     upper level Model nvar
# 1 0.7849190 0.008941676 0.7673937 0.8024444  0.95    LR    5
# 2 0.7776139 0.009215502 0.7595518 0.7956759  0.95    LR   10
# 3 0.7344339 0.010049168 0.7147379 0.7541299  0.95    RF    5
# 4 0.7458098 0.009860614 0.7264834 0.7651363  0.95    RF   10
# > main.Gambia[["calib"]]
# # A tibble: 2 ? 7
# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     5  0.0211    -0.888    0.767 0.898     0.343      1.55
# 2    10 -0.00208   -0.943    0.774 0.740     0.267      1.28

# AUC          SE     lower     upper level Model nvar
# 1  0.7976624 0.009452858 0.7791352 0.8161897  0.95    LR    1
# 2  0.8090085 0.008843450 0.7916757 0.8263413  0.95    LR    2
# 3  0.8083496 0.008882485 0.7909402 0.8257589  0.95    LR    3
# 4  0.8076695 0.008894526 0.7902365 0.8251024  0.95    LR    4
# 5  0.8057025 0.008915963 0.7882275 0.8231775  0.95    LR    5
# 6  0.8057919 0.008936545 0.7882766 0.8233072  0.95    LR    6
# 7  0.8055826 0.008937386 0.7880656 0.8230995  0.95    LR    7
# 8  0.8048670 0.008951456 0.7873225 0.8224116  0.95    LR    8
# 9  0.8062706 0.009115996 0.7884036 0.8241377  0.95    LR    9
# 10 0.7996615 0.009269369 0.7814938 0.8178291  0.95    LR   10
# 11 0.6939819 0.016526400 0.6615907 0.7263730  0.95    RF    1
# 12 0.7743176 0.010823484 0.7531040 0.7955312  0.95    RF    2
# 13 0.7766464 0.009393620 0.7582352 0.7950575  0.95    RF    3
# 14 0.7562166 0.010710072 0.7352253 0.7772080  0.95    RF    4
# 15 0.7531385 0.010146721 0.7332513 0.7730257  0.95    RF    5
# 16 0.7521292 0.010084955 0.7323630 0.7718953  0.95    RF    6
# 17 0.7558029 0.010079991 0.7360465 0.7755593  0.95    RF    7
# 18 0.7590646 0.009988616 0.7394872 0.7786419  0.95    RF    8
# 19 0.7750424 0.009674056 0.7560816 0.7940032  0.95    RF    9
# 20 0.7763614 0.009652583 0.7574427 0.7952801  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 -0.0726   -0.989    0.662 1.24      0.498      2.28
# 2     2 -0.0613   -0.982    0.679 1.26      0.467     NA   
# 3     3 -0.0582   -0.984    0.688 1.27      0.478      3.24
# 4     4 -0.0455   -0.975    0.705 1.49    -23.1      100.  
# 5     5 -0.0459   -0.979    0.709 1.46    -23.1      100.  
# 6     6 -0.0489   -0.984    0.708 1.43    -23.0       99.5 
# 7     7 -0.0446   -0.986    0.718 1.39    -22.7       98.1 
# 8     8 -0.0434   -0.990    0.724 1.38    -22.3       96.7 
# 9     9 -0.0540   -1.01     0.725 0.881     0.354      1.53
# 10    10 -0.0625   -1.03     0.727 0.841     0.324      1.50

################### main (0-59mo) any death -HAZ + MUAC by country - Mali ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz","site")]

#trying to find why sometimes has convergence issues
#dataset of only predictor variables considering
temp<-complete_anydeath_Mali[names]
#table all of these
temp %>% map( table )

main.Mali <- CPR.funct(data=complete_anydeath_Mali,outcome="death_all",iter=100,nvars_opts=c(1:10))
main.Mali[["df_imps"]]
main.Mali[["AUC_df"]]
main.Mali[["calib"]]

#save(main.Mali, file = "/main.Mali.100iter.callib.Rdata")


# names      var_red
# 1               f4b_muac 1.309507e+00
# 2               f4b_temp 1.165726e+00
# 3          f4a_ppl_house 8.137959e-01
# 4               f4b_resp 8.043884e-01
# 5               base_age 7.556542e-01
# 6        f4b_chest_indrw 4.497847e-01
# 7          f4a_slp_rooms 4.048785e-01
# 8         f4a_disp_feces 3.634668e-01
# 9       f4a_yng_children 3.478204e-01
# 10      f4a_hometrt_milk 3.362221e-01
# 11         f4a_share_fac 3.351933e-01
# 12     f4a_drh_bellypain 3.194060e-01
# 13          f4a_drh_days 3.129439e-01
# 14    f4a_cur_fastbreath 2.654661e-01
# 15        f4a_offr_drink 2.526784e-01
# 16        f4b_under_nutr 2.444512e-01
# 17         f4b_recommend 2.310980e-01
# 18         f4a_prim_schl 2.203433e-01
# 19        f4a_wash_nurse 2.168079e-01
# 20          f4a_ms_water 2.085009e-01
# 21              f4b_skin 2.030052e-01
# 22             f3_drh_iv 2.019421e-01
# 23            f4b_mental 2.004085e-01
# 24          f4b_abn_hair 1.982176e-01
# 25           f4a_ani_cat 1.812205e-01
# 26         f4a_fuel_elec 1.808878e-01
# 27          f4a_drh_conv 1.802555e-01
# 28      f4a_seek_privdoc 1.744954e-01
# 29        f4a_house_tele 1.744877e-01
# 30          f4a_wash_eat 1.740617e-01
# 31  f4a_drh_lethrgy_miss 1.725112e-01
# 32         f3_drh_turgor 1.675558e-01
# 33        f4a_wash_child 1.674243e-01
# 34      f4a_house_fridge 1.656109e-01
# 35           f4b_bipedal 1.601646e-01
# 36           f3_drh_hosp 1.593495e-01
# 37        f4a_house_bike 1.528639e-01
# 38             f4b_admit 1.495389e-01
# 39          f4a_cur_skin 1.478966e-01
# 40         f4a_drh_blood 1.432905e-01
# 41       f4a_seek_healer 1.423711e-01
# 42        f4a_seek_other 1.402546e-01
# 43        f4a_house_elec 1.395136e-01
# 44        f4a_water_well 1.388837e-01
# 45         f4a_drh_vomit 1.380875e-01
# 46         f4a_house_car 1.367608e-01
# 47         f4a_drh_cough 1.356716e-01
# 48     f4a_fuel_charcoal 1.294105e-01
# 49      f4a_water_pubtap 1.281103e-01
# 50         f4a_breastfed 1.273416e-01
# 51       f4a_house_scoot 1.251916e-01
# 52        f4a_max_stools 1.244628e-01
# 53         f4a_fuel_wood 1.187240e-01
# 54       f4a_ani_rodents 1.180949e-01
# 55          f4a_dad_live 1.129790e-01
# 56        f4a_drh_strain 1.117518e-01
# 57          f4a_wash_def 1.114232e-01
# 58             f3_gender 1.079901e-01
# 59       f4a_hometrt_ors 1.030138e-01
# 60      f4a_water_bought 1.028207e-01
# 61          f4a_ani_fowl 1.022469e-01
# 62      f4a_hometrt_herb 9.740232e-02
# 63         f4a_wash_cook 9.686748e-02
# 64      f4a_hometrt_none 9.331540e-02
# 65          f4a_wash_use 9.293439e-02
# 66      f4a_cur_restless 8.788425e-02
# 67        f4a_drh_thirst 8.691953e-02
# 68      f4a_seek_outside 8.340415e-02
# 69       f4a_cur_thirsty 7.741769e-02
# 70      f4a_drh_restless 7.405962e-02
# 71       f4a_water_avail 7.198771e-02
# 72      f4a_drh_prolapse 7.175057e-02
# 73       f4a_house_radio 6.968366e-02
# 74      f4a_cur_drymouth 6.966165e-02
# 75        f4a_fuel_other 6.574266e-02
# 76     f4a_hometrt_othr1 6.319770e-02
# 77         f4a_drh_consc 5.986525e-02
# 78        f4a_hometrt_ab 5.912580e-02
# 79     f4a_hometrt_othr2 5.870748e-02
# 80         f4a_ani_sheep 5.259282e-02
# 81             f4b_mouth 5.252014e-02
# 82            f4a_ani_no 5.161265e-02
# 83        f4a_water_yard 4.992945e-02
# 84     f4a_water_covwell 4.710739e-02
# 85     f4a_drh_lessdrink 4.704268e-02
# 86        f4a_trt_method 4.687936e-02
# 87        f4a_seek_pharm 4.287958e-02
# 88        f4b_skin_flaky 4.215358e-02
# 89         f4a_trt_water 3.660354e-02
# 90        f4a_house_cart 3.024091e-02
# 91           f4a_ani_dog 2.777491e-02
# 92       f4a_store_water 2.568369e-02
# 93   f4a_hometrt_othrliq 2.421966e-02
# 94      f4a_relationship 5.244593e-03
# 95        f4a_seek_remdy 4.656746e-03
# 96         f4a_fac_waste 3.933671e-03
# 97       f4a_house_phone 2.282227e-03
# 98     f4a_hometrt_maize 1.762966e-03
# 99     f4a_water_pubwell 1.485218e-03
# 100       f4a_house_boat 1.237605e-03
# 101            f4a_floor 9.572990e-04
# 102   f4a_water_covpwell 9.362836e-04
# 103        f4a_ani_other 8.988095e-04
# 104        f4a_wash_othr 5.483766e-04
# 105         f4a_ani_goat 5.316258e-04
# 106      f4a_wash_animal 1.777778e-04
# 107             f4b_eyes 1.333333e-04
# 108      f4a_seek_friend 1.114996e-04
# 109     f4a_house_agland 9.742063e-05
# 110      f4a_water_house 6.666667e-05
# 111     f4a_fuel_propane 5.194805e-05
# 112   f4a_water_deepwell 7.575758e-06
# 113                 site 0.000000e+00
# 114       f4a_house_none 0.000000e+00
# 115      f4a_fuel_biogas 0.000000e+00
# 116       f4a_fuel_grass 0.000000e+00
# 117        f4a_fuel_coal 0.000000e+00
# 118        f4a_fuel_dung 0.000000e+00
# 119      f4a_fuel_natgas 0.000000e+00
# 120        f4a_fuel_crop 0.000000e+00
# 121        f4a_fuel_kero 0.000000e+00
# 122          f4a_ani_cow 0.000000e+00
# 123  f4a_water_prospring 0.000000e+00
# 124   f4a_water_unspring 0.000000e+00
# 125      f4a_water_river 0.000000e+00
# 126       f4a_water_pond 0.000000e+00
# 127       f4a_water_rain 0.000000e+00
# 128  f4a_water_shallwell 0.000000e+00
# 129       f4a_water_othr 0.000000e+00
# 130       f4a_water_bore 0.000000e+00
# 131     f4a_hometrt_zinc 0.000000e+00
# 132         f4a_seek_doc 0.000000e+00
# 133           f4b_rectal 0.000000e+00
# > main.Mali[["AUC_df"]]
# AUC          SE     lower     upper level Model nvar
# 1 0.8805149 0.005639497 0.8694617 0.8915681  0.95    LR    5
# 2 0.8658233 0.005916459 0.8542272 0.8774193  0.95    LR   10
# 3 0.7926260 0.011143771 0.7707846 0.8144674  0.95    RF    5
# 4 0.7833254 0.011030389 0.7617063 0.8049446  0.95    RF   10
# > main.Mali[["calib"]]
# # A tibble: 2 ? 7
# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     5 -0.0616    -1.30    0.845 0.930     0.315      1.69
# 2    10 -0.0225    -1.32    0.944 0.731     0.240      1.35

# AUC          SE     lower     upper level Model nvar
# 1  0.8069986 0.007891507 0.7915315 0.8224656  0.95    LR    1
# 2  0.8955410 0.005487912 0.8847849 0.9062971  0.95    LR    2
# 3  0.8978454 0.005062182 0.8879237 0.9077670  0.95    LR    3
# 4  0.8934389 0.005135882 0.8833728 0.9035051  0.95    LR    4
# 5  0.8890300 0.005422446 0.8784022 0.8996578  0.95    LR    5
# 6  0.8836620 0.005456767 0.8729670 0.8943571  0.95    LR    6
# 7  0.8768986 0.005606090 0.8659109 0.8878863  0.95    LR    7
# 8  0.8731268 0.005795757 0.8617673 0.8844863  0.95    LR    8
# 9  0.8720127 0.005815350 0.8606148 0.8834105  0.95    LR    9
# 10 0.8722369 0.005827473 0.8608153 0.8836586  0.95    LR   10
# 11 0.6031108 0.022701802 0.5586161 0.6476055  0.95    RF    1
# 12 0.7323794 0.017684406 0.6977186 0.7670402  0.95    RF    2
# 13 0.7809822 0.012299924 0.7568747 0.8050896  0.95    RF    3
# 14 0.7924053 0.013237561 0.7664602 0.8183505  0.95    RF    4
# 15 0.8099229 0.010753764 0.7888459 0.8309999  0.95    RF    5
# 16 0.8068283 0.010587859 0.7860765 0.8275801  0.95    RF    6
# 17 0.8014055 0.009792318 0.7822129 0.8205981  0.95    RF    7
# 18 0.7977577 0.009381799 0.7793698 0.8161457  0.95    RF    8
# 19 0.7900897 0.011114791 0.7683051 0.8118742  0.95    RF    9
# 20 0.7896290 0.010928077 0.7682104 0.8110477  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 -0.0276     -1.24    0.838 0.913     0.106      1.82
# 2     2 -0.0211     -1.25    0.870 4.46    -27.4      103.  
# 3     3 -0.0143     -1.24    0.880 1.05      0.384      1.92
# 4     4 -0.00566    -1.24    0.891 1.01      0.356      1.85
# 5     5  0.00379    -1.23    0.904 0.969     0.330      1.77
# 6     6  0.0100     -1.24    0.925 0.915     0.293      1.76
# 7     7  0.0254     -1.23    0.947 0.839     0.257      1.56
# 8     8  0.0313     -1.23    0.966 0.771     0.231      1.41
# 9     9  0.0333     -1.24    0.979 0.734     0.219      1.34
# 10    10  0.0302     -1.26    0.988 0.719     0.218      1.32

################### main (0-59mo) any death -HAZ + MUAC by country - Mozam ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz","site")]

complete_anydeath_Mozam_complete <- complete_anydeath_Mozam
#f4a_dad_live recode 4 to 5 
summary(complete_anydeath_Mozam$f4a_dad_live)
complete_anydeath_Mozam_complete$f4a_dad_live <- as.factor(ifelse(complete_anydeath_Mozam_complete$f4a_dad_live==4,5,complete_anydeath_Mozam_complete$f4a_dad_live))
summary(complete_anydeath_Mozam_complete$f4a_dad_live)

#f4a_relationship recode everything that's not 1 to 2
summary(complete_anydeath_Mozam_complete$f4a_relationship)
complete_anydeath_Mozam_complete$f4a_relationship <- as.factor(ifelse(complete_anydeath_Mozam_complete$f4a_relationship!=1,2,complete_anydeath_Mozam_complete$f4a_relationship))
summary(complete_anydeath_Mozam_complete$f4a_relationship)

main.Mozam <- CPR.funct(data=complete_anydeath_Mozam_complete,outcome="death_all",iter=100,nvars_opts=c(1:10))
main.Mozam[["df_imps"]]
main.Mozam[["AUC_df"]]
main.Mozam[["calib"]]

#save(main.Mozam, file = "/main.Mozam.100iter.callib.Rdata")

# names      var_red
# 1               f4b_muac 3.9423174460
# 2               f4b_temp 1.3465010274
# 3               base_age 1.3244289351
# 4               f4b_resp 0.9673045820
# 5         f4b_under_nutr 0.8907643477
# 6           f4b_abn_hair 0.8776171334
# 7         f4a_offr_drink 0.7993387027
# 8          f4a_ppl_house 0.7765707441
# 9           f4a_dad_live 0.7013743099
# 10      f4a_yng_children 0.6164699805
# 11        f4a_disp_feces 0.5794455780
# 12         f4a_breastfed 0.5271639565
# 13       f4a_water_avail 0.4915452216
# 14          f4a_drh_days 0.4795655772
# 15         f4a_slp_rooms 0.4231601357
# 16     f4a_water_covwell 0.3994874381
# 17      f4a_drh_prolapse 0.3800546978
# 18         f4a_prim_schl 0.3503315865
# 19      f4a_house_agland 0.3304544978
# 20          f4a_ms_water 0.3124354379
# 21  f4a_drh_lethrgy_miss 0.2950181934
# 22         f4b_recommend 0.2923198379
# 23             f4b_mouth 0.2894289025
# 24        f4b_skin_flaky 0.2882473476
# 25          f4a_ani_goat 0.2808771146
# 26        f4a_max_stools 0.2789924615
# 27      f4a_relationship 0.2772932432
# 28         f4a_wash_cook 0.2724654917
# 29     f4a_water_pubwell 0.2691852828
# 30         f4a_drh_vomit 0.2671014795
# 31             f4a_floor 0.2615359831
# 32        f4a_water_bore 0.2577382332
# 33       f4a_ani_rodents 0.2556189174
# 34             f3_drh_iv 0.2462055089
# 35          f4a_wash_use 0.2420770482
# 36             f3_gender 0.2409847813
# 37              f4b_skin 0.2355233445
# 38         f4a_drh_consc 0.2300012811
# 39        f4a_trt_method 0.2293723709
# 40          f4a_wash_def 0.2289319233
# 41       f4b_chest_indrw 0.2288548208
# 42         f4a_drh_cough 0.2281227763
# 43      f4a_seek_outside 0.2264555569
# 44          f4a_cur_skin 0.2253112488
# 45      f4a_hometrt_herb 0.2223120675
# 46       f4a_house_phone 0.2205059716
# 47         f3_drh_turgor 0.2201666497
# 48            f4b_mental 0.2189859721
# 49       f4a_house_radio 0.2164079029
# 50    f4a_cur_fastbreath 0.2108544120
# 51           f3_drh_hosp 0.1978723212
# 52     f4a_drh_lessdrink 0.1959818490
# 53      f4a_cur_drymouth 0.1930289450
# 54              f4b_eyes 0.1877048251
# 55         f4a_trt_water 0.1858598294
# 56      f4a_hometrt_none 0.1822496895
# 57       f4a_cur_thirsty 0.1776964308
# 58          f4a_ani_fowl 0.1753973942
# 59             f4b_admit 0.1737434203
# 60        f4a_house_none 0.1712662066
# 61        f4a_house_bike 0.1681275034
# 62         f4a_share_fac 0.1651327188
# 63        f4a_drh_thirst 0.1644714003
# 64     f4a_drh_bellypain 0.1621371740
# 65       f4a_store_water 0.1587015309
# 66        f4a_wash_child 0.1561779898
# 67      f4a_water_pubtap 0.1517663720
# 68           f4a_ani_cat 0.1489855139
# 69           f4a_ani_dog 0.1471903302
# 70            f4b_rectal 0.1412033652
# 71       f4a_hometrt_ors 0.1394437881
# 72     f4a_fuel_charcoal 0.1378705889
# 73      f4a_drh_restless 0.1324951066
# 74          f4a_wash_eat 0.1296438524
# 75            f4a_ani_no 0.1268345431
# 76        f4a_house_elec 0.1267178205
# 77        f4a_wash_nurse 0.1214857673
# 78      f4a_cur_restless 0.1045921841
# 79           f4a_ani_cow 0.1020213938
# 80         f4a_ani_sheep 0.0975879722
# 81         f4a_fuel_wood 0.0962101067
# 82       f4a_fuel_natgas 0.0956158349
# 83    f4a_water_covpwell 0.0922356948
# 84         f4a_drh_blood 0.0775398048
# 85        f4a_water_well 0.0769001657
# 86        f4a_house_tele 0.0723936056
# 87         f4a_ani_other 0.0641535023
# 88         f4a_fac_waste 0.0522966989
# 89        f4a_water_yard 0.0464238010
# 90        f4a_drh_strain 0.0385717698
# 91         f4a_house_car 0.0271801198
# 92      f4a_house_fridge 0.0160377435
# 93           f4b_bipedal 0.0084902858
# 94     f4a_hometrt_othr2 0.0067131583
# 95          f4a_drh_conv 0.0038047711
# 96     f4a_hometrt_othr1 0.0028285044
# 97       f4a_water_house 0.0026877594
# 98       f4a_water_river 0.0016875000
# 99       f4a_wash_animal 0.0015060107
# 100       f4a_hometrt_ab 0.0012019536
# 101       f4a_water_rain 0.0010666667
# 102        f4a_fuel_elec 0.0010552381
# 103       f4a_house_boat 0.0009000000
# 104    f4a_hometrt_maize 0.0008666667
# 105       f4a_house_cart 0.0006214176
# 106  f4a_water_shallwell 0.0001578947
# 107                 site 0.0000000000
# 108      f4a_house_scoot 0.0000000000
# 109      f4a_fuel_biogas 0.0000000000
# 110       f4a_fuel_grass 0.0000000000
# 111     f4a_fuel_propane 0.0000000000
# 112        f4a_fuel_coal 0.0000000000
# 113        f4a_fuel_dung 0.0000000000
# 114        f4a_fuel_crop 0.0000000000
# 115        f4a_fuel_kero 0.0000000000
# 116       f4a_fuel_other 0.0000000000
# 117  f4a_water_prospring 0.0000000000
# 118   f4a_water_unspring 0.0000000000
# 119       f4a_water_pond 0.0000000000
# 120   f4a_water_deepwell 0.0000000000
# 121     f4a_water_bought 0.0000000000
# 122       f4a_water_othr 0.0000000000
# 123        f4a_wash_othr 0.0000000000
# 124     f4a_hometrt_milk 0.0000000000
# 125     f4a_hometrt_zinc 0.0000000000
# 126  f4a_hometrt_othrliq 0.0000000000
# 127       f4a_seek_pharm 0.0000000000
# 128      f4a_seek_friend 0.0000000000
# 129      f4a_seek_healer 0.0000000000
# 130         f4a_seek_doc 0.0000000000
# 131     f4a_seek_privdoc 0.0000000000
# 132       f4a_seek_remdy 0.0000000000
# 133       f4a_seek_other 0.0000000000
# > main.Mozam[["AUC_df"]]
# AUC          SE     lower     upper level Model nvar
# 1 0.7860013 0.008690669 0.7689679 0.8030347  0.95    LR    5
# 2 0.7666795 0.009578473 0.7479061 0.7854530  0.95    LR   10
# 3 0.7424936 0.010166174 0.7225683 0.7624189  0.95    RF    5
# 4 0.7307943 0.010505919 0.7102031 0.7513855  0.95    RF   10
# > main.Mozam[["calib"]]
# # A tibble: 2 ? 7
# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     5 -0.0557   -0.993    0.716 0.894     0.282      1.67
# 2    10 -0.0641   -1.02     0.729 0.720     0.199      1.37

# AUC          SE     lower     upper level Model nvar
# 1  0.8110966 0.007634148 0.7961340 0.8260593  0.95    LR    1
# 2  0.8029129 0.007940776 0.7873493 0.8184765  0.95    LR    2
# 3  0.8122990 0.007970979 0.7966762 0.8279218  0.95    LR    3
# 4  0.8045354 0.008066750 0.7887249 0.8203460  0.95    LR    4
# 5  0.7979895 0.008196978 0.7819237 0.8140553  0.95    LR    5
# 6  0.7949972 0.008326855 0.7786769 0.8113175  0.95    LR    6
# 7  0.7927221 0.008420936 0.7762174 0.8092269  0.95    LR    7
# 8  0.7893091 0.008634358 0.7723861 0.8062321  0.95    LR    8
# 9  0.7865634 0.008800841 0.7693140 0.8038127  0.95    LR    9
# 10 0.7821700 0.009089561 0.7643548 0.7999852  0.95    LR   10
# 11 0.7233348 0.013964125 0.6959656 0.7507040  0.95    RF    1
# 12 0.7300486 0.010105983 0.7102413 0.7498560  0.95    RF    2
# 13 0.7480020 0.009997949 0.7284064 0.7675977  0.95    RF    3
# 14 0.7676085 0.009553986 0.7488830 0.7863339  0.95    RF    4
# 15 0.7646466 0.009440086 0.7461443 0.7831488  0.95    RF    5
# 16 0.7584579 0.009432713 0.7399701 0.7769457  0.95    RF    6
# 17 0.7566420 0.009520588 0.7379820 0.7753021  0.95    RF    7
# 18 0.7487476 0.009719122 0.7296984 0.7677967  0.95    RF    8
# 19 0.7541711 0.009658383 0.7352410 0.7731012  0.95    RF    9
# 20 0.7517044 0.009767378 0.7325607 0.7708481  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 0.00626   -0.873    0.736 1.14      0.421      2.06
# 2     2 0.0239    -0.862    0.761 1.10      0.411      2.00
# 3     3 0.0131    -0.882    0.760 1.10      0.435      2.00
# 4     4 0.0138    -0.883    0.763 1.03      0.388      1.88
# 5     5 0.0228    -0.878    0.775 0.966     0.350      1.75
# 6     6 0.0327    -0.870    0.788 0.931     0.332      1.69
# 7     7 0.0344    -0.874    0.795 0.911     0.325      1.66
# 8     8 0.0350    -0.878    0.800 0.882     0.314      1.61
# 9     9 0.0304    -0.887    0.800 0.849     0.299      1.55
# 10    10 0.0292    -0.894    0.804 0.822     0.277      1.56

################### main (0-59mo) any death -HAZ + MUAC by country - Kenya ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz","site")]

complete_anydeath_Kenya_complete <- complete_anydeath_Kenya
summary(complete_anydeath_Kenya$f4a_relationship)
#everything that not 1 recode to 2
complete_anydeath_Kenya$f4a_relationship <- as.factor(ifelse(complete_anydeath_Kenya$f4a_relationship!=1,2,complete_anydeath_Kenya$f4a_relationship))
summary(complete_anydeath_Kenya$f4a_relationship)

summary(complete_anydeath_Kenya_complete$f4a_trt_method )
#5 and 7 recode to 2
complete_anydeath_Kenya_complete$f4a_trt_method <- as.factor(ifelse((complete_anydeath_Kenya_complete$f4a_trt_method==5|complete_anydeath_Kenya_complete$f4a_trt_method==7)
                                                           ,2,complete_anydeath_Kenya_complete$f4a_trt_method))
complete_anydeath_Kenya_complete$f4a_trt_method <- as.factor(ifelse((complete_anydeath_Kenya_complete$f4a_trt_method==1)
                                                           ,0,complete_anydeath_Kenya_complete$f4a_trt_method))
summary(complete_anydeath_Kenya_complete$f4a_trt_method)

summary(complete_anydeath_Kenya$f4a_prim_schl)
#recode 5 into 4
#drop 1 and 7
complete_anydeath_Kenya_complete$f4a_prim_schl <- as.factor(ifelse((complete_anydeath_Kenya_complete$f4a_prim_schl==5)
                                                                    ,4,complete_anydeath_Kenya_complete$f4a_prim_schl))
complete_anydeath_Kenya_complete <- complete_anydeath_Kenya_complete %>% filter(f4a_prim_schl!=1 & f4a_prim_schl!=7)
summary(complete_anydeath_Kenya_complete$f4a_prim_schl)

summary(complete_anydeath_Kenya_complete$f4a_ms_water)
#3 recode into 2, drop 4; based on recode did in var creation
complete_anydeath_Kenya_complete$f4a_ms_water <- 
  (ifelse((complete_anydeath_Kenya_complete$f4a_ms_water==0),
          0,complete_anydeath_Kenya_complete$f4a_ms_water))
#table(complete_anydeath_Kenya_complete$f4a_ms_water)
complete_anydeath_Kenya_complete$f4a_ms_water <- 
  (ifelse((complete_anydeath_Kenya_complete$f4a_ms_water==2),
          1,complete_anydeath_Kenya_complete$f4a_ms_water))
#table(complete_anydeath_Kenya_complete$f4a_ms_water)
complete_anydeath_Kenya_complete$f4a_ms_water <- 
  (ifelse((complete_anydeath_Kenya_complete$f4a_ms_water==3 | complete_anydeath_Kenya_complete$f4a_ms_water==4),
          2,complete_anydeath_Kenya_complete$f4a_ms_water))
#table(complete_anydeath_Kenya_complete$f4a_ms_water)
complete_anydeath_Kenya_complete <- complete_anydeath_Kenya_complete %>% filter(f4a_ms_water!=5)
#table(complete_anydeath_Kenya_complete$f4a_ms_water)
complete_anydeath_Kenya_complete$f4a_ms_water <- as.factor(complete_anydeath_Kenya_complete$f4a_ms_water)
summary(complete_anydeath_Kenya_complete$f4a_ms_water)

main.Kenya <- CPR.funct(data=complete_anydeath_Kenya_complete,outcome="death_all",iter=100,nvars_opts=c(1:10))
main.Kenya[["df_imps"]]
main.Kenya[["AUC_df"]]
main.Kenya[["calib"]]

#save(main.Kenya, file = "/main.Kenya.100iter.callib.Rdata")


# names      var_red
# 1               f4b_muac 3.6300150238
# 2               base_age 1.5420481838
# 3               f4b_resp 1.4165820330
# 4               f4b_temp 1.2505191414
# 5         f4b_under_nutr 0.8142535797
# 6         f4a_offr_drink 0.7457057436
# 7            f4b_bipedal 0.7214852392
# 8           f4a_drh_days 0.7044880368
# 9          f4a_ppl_house 0.6696617802
# 10          f4a_ms_water 0.6447271820
# 11          f4b_abn_hair 0.6216284682
# 12        f4b_skin_flaky 0.6087990126
# 13         f4a_share_fac 0.5964709153
# 14         f4a_prim_schl 0.5316649146
# 15          f4a_dad_live 0.4903319380
# 16              f4b_skin 0.4789223457
# 17      f4a_yng_children 0.4422889717
# 18      f4a_relationship 0.4309762244
# 19             f4b_admit 0.4233000967
# 20            f4b_mental 0.4100040947
# 21         f4a_breastfed 0.4096661008
# 22             f3_drh_iv 0.3864831106
# 23        f4a_disp_feces 0.3706591817
# 24          f4a_cur_skin 0.3534054425
# 25           f3_drh_hosp 0.3468565768
# 26         f3_drh_turgor 0.3299958785
# 27        f4a_hometrt_ab 0.3274131519
# 28   f4a_water_prospring 0.3097770848
# 29        f4a_trt_method 0.3042766205
# 30        f4a_water_yard 0.2944001471
# 31       f4a_fuel_biogas 0.2685983061
# 32      f4a_hometrt_herb 0.2510876574
# 33          f4a_wash_def 0.2442959111
# 34    f4a_cur_fastbreath 0.2399406910
# 35         f4a_wash_othr 0.2370044020
# 36           f4a_ani_dog 0.2353677157
# 37       f4a_seek_healer 0.2351448262
# 38      f4a_seek_privdoc 0.2275279709
# 39         f4a_drh_vomit 0.2251699685
# 40        f4a_drh_strain 0.2242735172
# 41        f4a_wash_nurse 0.2240162360
# 42        f4a_max_stools 0.2224927261
# 43    f4a_water_deepwell 0.2165787545
# 44       f4a_water_river 0.2160689104
# 45     f4a_drh_bellypain 0.2096600885
# 46             f4b_mouth 0.2094900246
# 47        f4a_water_rain 0.2077960698
# 48       f4a_cur_thirsty 0.2059337837
# 49           f4a_ani_cow 0.2055028951
# 50        f4a_water_well 0.2045570076
# 51  f4a_drh_lethrgy_miss 0.2040009573
# 52       f4a_water_avail 0.2033323993
# 53      f4a_drh_restless 0.2022997191
# 54     f4a_fuel_charcoal 0.2010926936
# 55      f4a_drh_prolapse 0.2007027052
# 56         f4a_wash_cook 0.1988257554
# 57         f4a_slp_rooms 0.1985952555
# 58           f4a_ani_cat 0.1933994249
# 59       f4a_ani_rodents 0.1913947035
# 60       f4b_chest_indrw 0.1893588681
# 61        f4a_water_pond 0.1883345494
# 62         f4a_drh_consc 0.1865372138
# 63         f4a_ani_sheep 0.1862714805
# 64      f4a_cur_restless 0.1850134127
# 65             f3_gender 0.1808521144
# 66        f4a_house_bike 0.1795955499
# 67       f4a_wash_animal 0.1776855636
# 68      f4a_hometrt_none 0.1749553808
# 69     f4a_hometrt_maize 0.1741211830
# 70             f4a_floor 0.1712186083
# 71         f4a_trt_water 0.1695668077
# 72          f4a_wash_eat 0.1641302747
# 73       f4a_house_scoot 0.1627854077
# 74          f4a_ani_goat 0.1616670800
# 75         f4a_ani_other 0.1606390998
# 76       f4a_hometrt_ors 0.1592508880
# 77        f4a_water_bore 0.1581863043
# 78       f4a_house_phone 0.1544660070
# 79         f4a_drh_cough 0.1541343218
# 80       f4a_store_water 0.1504175417
# 81        f4a_drh_thirst 0.1454031005
# 82      f4a_cur_drymouth 0.1366389081
# 83     f4a_drh_lessdrink 0.1359335051
# 84      f4a_seek_outside 0.1352374822
# 85         f4a_fuel_kero 0.1338256855
# 86     f4a_hometrt_othr1 0.1283247335
# 87     f4a_water_pubwell 0.1282023519
# 88        f4a_wash_child 0.1254894056
# 89      f4a_hometrt_milk 0.1238826575
# 90       f4a_house_radio 0.1185762272
# 91        f4a_house_cart 0.1084427009
# 92        f4a_house_tele 0.1083008082
# 93          f4a_seek_doc 0.1031247790
# 94      f4a_house_agland 0.0976464484
# 95     f4a_water_covwell 0.0728429304
# 96            f4b_rectal 0.0704388855
# 97         f4a_fac_waste 0.0692980996
# 98          f4a_drh_conv 0.0641456126
# 99        f4a_seek_other 0.0618095055
# 100        f4a_drh_blood 0.0610952740
# 101   f4a_water_covpwell 0.0605194172
# 102         f4a_ani_fowl 0.0586056463
# 103       f4a_seek_pharm 0.0535353495
# 104         f4a_wash_use 0.0526855453
# 105             f4b_eyes 0.0504221945
# 106  f4a_hometrt_othrliq 0.0497063976
# 107        f4a_fuel_crop 0.0438216033
# 108       f4a_fuel_grass 0.0315729112
# 109        f4a_fuel_wood 0.0275305598
# 110     f4a_hometrt_zinc 0.0266203961
# 111   f4a_water_unspring 0.0250355844
# 112     f4a_water_pubtap 0.0091780957
# 113    f4a_hometrt_othr2 0.0073944256
# 114       f4a_house_elec 0.0030151261
# 115     f4b_nature_stool 0.0021288371
# 116        f4b_recommend 0.0020215949
# 117      f4a_seek_friend 0.0014285714
# 118       f4a_seek_remdy 0.0011796245
# 119        f4a_fuel_dung 0.0008571429
# 120  f4a_water_shallwell 0.0007783443
# 121     f4a_house_fridge 0.0006666667
# 122           f4a_ani_no 0.0003000000
# 123     f4a_water_bought 0.0003000000
# 124      f4a_water_house 0.0001269841
# 125                 site 0.0000000000
# 126        f4a_house_car 0.0000000000
# 127       f4a_house_boat 0.0000000000
# 128       f4a_house_none 0.0000000000
# 129        f4a_fuel_elec 0.0000000000
# 130     f4a_fuel_propane 0.0000000000
# 131        f4a_fuel_coal 0.0000000000
# 132      f4a_fuel_natgas 0.0000000000
# 133       f4a_fuel_other 0.0000000000
# 134       f4a_water_othr 0.0000000000
# 135    f4b_observe_stool 0.0000000000

# AUC         SE     lower     upper level Model nvar
# 1  0.8880011 0.01440306 0.8597716 0.9162306  0.95    LR    1
# 2  0.8860816 0.01533585 0.8560239 0.9161393  0.95    LR    2
# 3  0.8825447 0.01547045 0.8522231 0.9128662  0.95    LR    3
# 4  0.8778981 0.01544935 0.8476179 0.9081783  0.95    LR    4
# 5  0.8775933 0.01500334 0.8481873 0.9069993  0.95    LR    5
# 6  0.8641745 0.01701747 0.8308209 0.8975282  0.95    LR    6
# 7  0.8546275 0.01689173 0.8215203 0.8877347  0.95    LR    7
# 8  0.8512858 0.01718611 0.8176017 0.8849700  0.95    LR    8
# 9  0.8482877 0.01635094 0.8162404 0.8803349  0.95    LR    9
# 10 0.8393910 0.01636074 0.8073245 0.8714574  0.95    LR   10
# 11 0.7952438 0.03887227 0.7190555 0.8714320  0.95    RF    1
# 12 0.7990024 0.02080479 0.7582257 0.8397790  0.95    RF    2
# 13 0.8253050 0.01630276 0.7933522 0.8572578  0.95    RF    3
# 14 0.8280208 0.01846568 0.7918287 0.8642129  0.95    RF    4
# 15 0.8325981 0.01598887 0.8012605 0.8639357  0.95    RF    5
# 16 0.8424515 0.01473989 0.8135618 0.8713411  0.95    RF    6
# 17 0.8467445 0.01387406 0.8195518 0.8739372  0.95    RF    7
# 18 0.8469671 0.01375321 0.8200113 0.8739229  0.95    RF    8
# 19 0.8364854 0.01454151 0.8079846 0.8649863  0.95    RF    9
# 20 0.8271754 0.01588630 0.7960388 0.8583119  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1  0.00450   -0.821    0.686 1.23      0.640      1.96
# 2     2 -0.00269   -0.830    0.682 1.19      0.624      1.90
# 3     3  0.00192   -0.829    0.690 1.14      0.583      1.83
# 4     4  0.00399   -0.833    0.697 1.06      0.533      1.69
# 5     5  0.0248    -0.816    0.722 0.977     0.479      1.57
# 6     6  0.0261    -0.820    0.728 0.906     0.425      1.46
# 7     7 -0.00239   -0.867    0.715 0.806     0.364      1.31
# 8     8  0.0124    -0.851    0.729 0.807     0.366      1.31
# 9     9  0.0200    -0.849    0.743 0.751     0.330      1.23
# 10    10 -0.00378   -0.884    0.727 0.662     0.273      1.11

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation


################### main (0-59mo) any death -HAZ + MUAC by country - India ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz","site")]

#too few outcomes to converge
main.India <- CPR.funct(data=complete_anydeath_India,outcome="death_all",iter=10,nvars_opts=c(5,10))
main.India[["df_imps"]]
main.India[["AUC_df"]]
main.India[["calib"]]

################### main (0-59mo) any death -HAZ + MUAC by country - Bang ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz","site")]

#too few outcomes to converge
main.Bang <- CPR.funct(data=complete_anydeath_Bang,outcome="death_all",iter=10,nvars_opts=c(5,10))
main.Bang[["df_imps"]]
main.Bang[["AUC_df"]]
main.Bang[["calib"]]

################### main (0-59mo) any death -HAZ + MUAC by country - Pak ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz","site")]

#too few outcomes to converge reliably
main.Pak <- CPR.funct(data=complete_anydeath_Pak_complete,outcome="death_all",iter=100,nvars_opts=c(5,10))
main.Pak[["df_imps"]]
main.Pak[["AUC_df"]]
main.Pak[["calib"]]

################### by country any death - fit Africa, val in Asia ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

death.Afr <- CPR.funct(data=cases_anydeath_Afr,outcome="death_all",iter=100,nvars_opts=c(1:10))
death.Afr[["df_imps"]]
death.Afr[["AUC_df"]]
death.Afr[["calib"]]

#save(death.Afr, file = "/death.Afr.100iter.callib.Rdata")

# names      var_red
# 1               f4b_muac 1.208133e+01
# 2               base_age 4.766129e+00
# 3               f4b_resp 4.468878e+00
# 4               f4b_temp 4.368273e+00
# 5          f4a_ppl_house 3.240320e+00
# 6           f4a_drh_days 2.422367e+00
# 7       f4a_yng_children 2.239802e+00
# 8         f4a_offr_drink 2.216444e+00
# 9         f4b_under_nutr 2.135212e+00
# 10         f4a_slp_rooms 2.033796e+00
# 11          f4b_abn_hair 1.999719e+00
# 12        f4a_disp_feces 1.910797e+00
# 13          f4a_dad_live 1.785467e+00
# 14             f4b_admit 1.586280e+00
# 15           f3_drh_hosp 1.507089e+00
# 16          f4a_ms_water 1.450406e+00
# 17      f4a_relationship 1.366418e+00
# 18       f4b_chest_indrw 1.330267e+00
# 19         f4a_prim_schl 1.318111e+00
# 20             f3_drh_iv 1.285555e+00
# 21        f4b_skin_flaky 1.282919e+00
# 22         f4b_recommend 1.251537e+00
# 23         f4a_share_fac 1.240815e+00
# 24            f4b_mental 1.193433e+00
# 25              f4b_skin 1.177000e+00
# 26         f4a_breastfed 1.174320e+00
# 27     f4a_drh_bellypain 1.130650e+00
# 28       f4a_water_avail 1.127484e+00
# 29      f4a_drh_prolapse 1.121456e+00
# 30         f3_drh_turgor 1.028193e+00
# 31          f4a_cur_skin 8.958434e-01
# 32         f4a_drh_vomit 8.926731e-01
# 33    f4a_cur_fastbreath 8.874783e-01
# 34                  site 8.856284e-01
# 35        f4a_trt_method 8.316952e-01
# 36      f4a_hometrt_herb 8.221795e-01
# 37  f4a_drh_lethrgy_miss 8.008671e-01
# 38         f4a_drh_cough 7.516670e-01
# 39          f4a_wash_def 7.458866e-01
# 40           f4b_bipedal 7.340997e-01
# 41         f4a_wash_cook 7.337364e-01
# 42             f4b_mouth 7.299777e-01
# 43       f4a_ani_rodents 7.289018e-01
# 44        f4a_wash_nurse 7.147566e-01
# 45       f4a_cur_thirsty 6.998188e-01
# 46     f4a_water_pubwell 6.906526e-01
# 47        f4a_max_stools 6.902576e-01
# 48             f3_gender 6.805344e-01
# 49      f4a_hometrt_none 6.790375e-01
# 50          f4a_ani_goat 6.703381e-01
# 51          f4a_wash_use 6.653447e-01
# 52       f4a_hometrt_ors 6.399245e-01
# 53           f4a_ani_cow 6.325786e-01
# 54     f4a_drh_lessdrink 6.308545e-01
# 55     f4a_water_covwell 6.300238e-01
# 56        f4a_drh_thirst 6.294038e-01
# 57         f4a_drh_blood 6.292492e-01
# 58        f4a_house_bike 6.196760e-01
# 59       f4a_house_radio 6.071345e-01
# 60      f4a_house_agland 6.033673e-01
# 61        f4a_water_bore 6.016059e-01
# 62           f4a_ani_cat 6.007088e-01
# 63        f4a_house_elec 5.847746e-01
# 64          f4a_wash_eat 5.830405e-01
# 65      f4a_seek_outside 5.772039e-01
# 66     f4a_fuel_charcoal 5.693190e-01
# 67      f4a_cur_restless 5.618570e-01
# 68         f4a_ani_sheep 5.601446e-01
# 69        f4a_wash_child 5.556206e-01
# 70             f4a_floor 5.518843e-01
# 71      f4a_water_pubtap 5.443113e-01
# 72      f4a_drh_restless 5.434560e-01
# 73       f4a_seek_healer 5.409351e-01
# 74       f4a_house_scoot 5.146325e-01
# 75           f4a_ani_dog 5.067833e-01
# 76         f4a_drh_consc 5.004762e-01
# 77      f4a_cur_drymouth 4.859299e-01
# 78        f4a_house_tele 4.849843e-01
# 79         f4a_ani_other 4.758770e-01
# 80    f4a_water_deepwell 4.756361e-01
# 81              f4b_eyes 4.712630e-01
# 82       f4a_house_phone 4.602386e-01
# 83          f4a_ani_fowl 4.580999e-01
# 84         f4a_trt_water 4.549548e-01
# 85      f4a_hometrt_milk 4.350684e-01
# 86     f4a_hometrt_othr1 4.262835e-01
# 87     f4a_hometrt_maize 4.105692e-01
# 88       f4a_store_water 4.090584e-01
# 89      f4a_seek_privdoc 4.058048e-01
# 90        f4a_hometrt_ab 3.939340e-01
# 91       f4a_fuel_biogas 3.922620e-01
# 92        f4a_house_cart 3.900320e-01
# 93        f4a_drh_strain 3.876290e-01
# 94       f4a_water_river 3.853945e-01
# 95        f4a_water_pond 3.519267e-01
# 96        f4a_water_well 3.487685e-01
# 97   f4a_water_prospring 3.483217e-01
# 98        f4a_water_yard 3.403592e-01
# 99         f4a_wash_othr 3.334590e-01
# 100       f4a_water_rain 3.288109e-01
# 101       f4a_seek_pharm 3.232331e-01
# 102        f4a_fuel_wood 3.169618e-01
# 103   f4a_water_covpwell 3.001586e-01
# 104     f4a_house_fridge 2.997953e-01
# 105         f4a_seek_doc 2.907745e-01
# 106           f4a_ani_no 2.751980e-01
# 107           f4b_rectal 2.640636e-01
# 108        f4a_house_car 2.617291e-01
# 109         f4a_drh_conv 2.585679e-01
# 110       f4a_house_none 2.552049e-01
# 111        f4a_fuel_kero 2.381138e-01
# 112    f4a_hometrt_othr2 2.293001e-01
# 113     f4a_water_bought 2.274598e-01
# 114        f4a_fuel_elec 2.269001e-01
# 115      f4a_wash_animal 2.241672e-01
# 116       f4a_seek_other 1.903858e-01
# 117        f4a_fac_waste 1.820588e-01
# 118       f4a_fuel_other 1.796922e-01
# 119       f4a_water_othr 1.590669e-01
# 120  f4a_hometrt_othrliq 1.590493e-01
# 121      f4a_fuel_natgas 1.440999e-01
# 122        f4a_fuel_crop 8.171315e-02
# 123      f4a_seek_friend 5.956406e-02
# 124   f4a_water_unspring 5.866763e-02
# 125       f4a_fuel_grass 5.046869e-02
# 126     f4a_hometrt_zinc 4.134958e-02
# 127       f4a_seek_remdy 7.249047e-03
# 128  f4a_water_shallwell 5.534127e-03
# 129      f4a_water_house 2.764775e-03
# 130     f4a_fuel_propane 7.484848e-04
# 131       f4a_house_boat 0.000000e+00
# 132        f4a_fuel_coal 0.000000e+00
# 133        f4a_fuel_dung 0.000000e+00

# AUC          SE     lower     upper level Model nvar
# 1  0.8249597 0.004022214 0.8170763 0.8328431  0.95    LR    1
# 2  0.8273294 0.004000347 0.8194889 0.8351699  0.95    LR    2
# 3  0.8341331 0.003916417 0.8264570 0.8418091  0.95    LR    3
# 4  0.8410170 0.003824326 0.8335214 0.8485125  0.95    LR    4
# 5  0.8407945 0.003842039 0.8332642 0.8483247  0.95    LR    5
# 6  0.8374457 0.003896816 0.8298080 0.8450833  0.95    LR    6
# 7  0.8357007 0.003945839 0.8279670 0.8434344  0.95    LR    7
# 8  0.8345683 0.003984271 0.8267593 0.8423773  0.95    LR    8
# 9  0.8334682 0.004027660 0.8255741 0.8413623  0.95    LR    9
# 10 0.8322169 0.004061115 0.8242573 0.8401766  0.95    LR   10
# 11 0.6978250 0.008216251 0.6817215 0.7139286  0.95    RF    1
# 12 0.7815949 0.004652133 0.7724769 0.7907129  0.95    RF    2
# 13 0.8070756 0.004056064 0.7991258 0.8150253  0.95    RF    3
# 14 0.8277236 0.003973383 0.8199359 0.8355113  0.95    RF    4
# 15 0.8181875 0.004094178 0.8101631 0.8262119  0.95    RF    5
# 16 0.8178248 0.004084637 0.8098190 0.8258305  0.95    RF    6
# 17 0.8169484 0.004087663 0.8089368 0.8249601  0.95    RF    7
# 18 0.8195290 0.004047153 0.8115967 0.8274612  0.95    RF    8
# 19 0.8246133 0.004008987 0.8167558 0.8324707  0.95    RF    9
# 20 0.8277448 0.003965218 0.8199731 0.8355165  0.95    RF   10

# # A tibble: 10 ? 7
# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 -0.0371   -0.463    0.346 1.00      0.697      1.33
# 2     2 -0.0342   -0.461    0.349 1.00      0.699      1.33
# 3     3 -0.0341   -0.463    0.352 0.987     0.695      1.31
# 4     4 -0.0333   -0.464    0.354 0.993     0.705      1.31
# 5     5 -0.0333   -0.465    0.355 0.981     0.697      1.29
# 6     6 -0.0352   -0.467    0.354 0.963     0.683      1.27
# 7     7 -0.0370   -0.470    0.353 0.949     0.672      1.25
# 8     8 -0.0372   -0.471    0.354 0.938     0.665      1.24
# 9     9 -0.0376   -0.473    0.355 0.929     0.660      1.22
# 10    10 -0.0342   -0.470    0.359 0.922     0.656      1.21


GEMS_glm_death_Afr_2var <- glm(death_all~f4b_muac+f4b_resp,
                            data=cases_anydeath_Afr,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death_Afr_2var)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  6.573367   0.849457   7.738 1.01e-14 ***
#   f4b_muac    -0.847385   0.061360 -13.810  < 2e-16 ***
#   f4b_resp     0.025337   0.007887   3.213  0.00132 ** 

round(exp(coef(GEMS_glm_death_Afr_2var)),2)
# (Intercept)    f4b_muac    f4b_resp 
# 715.78        0.43        1.03  
round(exp(confint(GEMS_glm_death_Afr_2var)),2)
# 2.5 %  97.5 %
#   (Intercept) 137.46 3852.76
# f4b_muac      0.38    0.48
# f4b_resp      1.01    1.04

Afr_2var_AfrFit<-cases_anydeath_Afr %>% select(death_all,f4b_muac,f4b_resp)
Afr_2var_AfrFit$Afr_pred_glm <- as.numeric(predict(GEMS_glm_death_Afr_2var,newdata=Afr_2var_AfrFit,type="response"))
Afr_2var_AfrFit_AUC <- roc(response=Afr_2var_AfrFit$death_all,predictor=Afr_2var_AfrFit$Afr_pred_glm)
paste(round(Afr_2var_AfrFit_AUC$auc,2)," (",
      round(ci.auc(Afr_2var_AfrFit_AUC)[1],2),", ",
      round(ci.auc(Afr_2var_AfrFit_AUC)[3],2),")",sep="")
# "0.84 (0.81, 0.87)"

Asia_2var_AfrFit<-cases_anydeath_SEAsia %>% select(death_all,f4b_muac,f4b_resp)
Asia_2var_AfrFit$Afr_pred_glm <- as.numeric(predict(GEMS_glm_death_Afr_2var,newdata=Asia_2var_AfrFit,type="response"))
Asia_2var_AfrFit_AUC <- roc(response=Asia_2var_AfrFit$death_all,predictor=Asia_2var_AfrFit$Afr_pred_glm)
paste(round(Asia_2var_AfrFit_AUC$auc,2)," (",
      round(ci.auc(Asia_2var_AfrFit_AUC)[1],2),", ",
      round(ci.auc(Asia_2var_AfrFit_AUC)[3],2),")",sep="")
# "0.93 (0.9, 0.96)"

################### by country any death - Asia, val in Africa ####
death.Asia <- CPR.funct(data=cases_anydeath_SEAsia,outcome="death_all",iter=20,nvars_opts=c(1:10))
death.Asia[["df_imps"]]
death.Asia[["AUC_df"]]
death.Asia[["calib"]]

#too few outcomes to converge reliably
table(cases_anydeath_SEAsia$death_all)
# 0    1 
# 3755   23 

GEMS_glm_gf_Asia_2var <- glm(death_all~f4b_muac+f4b_resp,
                             data=cases_anydeath_SEAsia,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_gf_Asia_2var)


round(exp(coef(GEMS_glm_gf_Asia_2var)),2)
round(exp(confint(GEMS_glm_gf_Asia_2var)),2)

Asia_2var_AsiaFit<-cases_anydeath_SEAsia %>% select(death_all,f4b_muac,f4b_resp)
Asia_2var_AsiaFit$Asia_pred_glm <- as.numeric(predict(GEMS_glm_gf_Asia_2var,newdata=Asia_2var_AsiaFit,type="response"))
Asia_2var_AsiaFit_AUC <- roc(response=Asia_2var_AsiaFit$death_all,predictor=Asia_2var_AsiaFit$Asia_pred_glm)
paste(round(Asia_2var_AsiaFit_AUC$auc,2)," (",
      round(ci.auc(Asia_2var_AsiaFit_AUC)[1],2),", ",
      round(ci.auc(Asia_2var_AsiaFit_AUC)[3],2),")",sep="")

Afr_2var_AsiaFit<-cases_gf_Afr %>% select(death_all,f4b_muac,f4b_resp)
Afr_2var_AsiaFit$Asia_pred_glm <- as.numeric(predict(GEMS_glm_gf_Asia_2var,newdata=Afr_2var_AsiaFit,type="response"))
Afr_2var_AsiaFit_AUC <- roc(response=Afr_2var_AsiaFit$death_all,predictor=Afr_2var_AsiaFit$Asia_pred_glm)
paste(round(Afr_2var_AsiaFit_AUC$auc,2)," (",
      round(ci.auc(Afr_2var_AsiaFit_AUC)[1],2),", ",
      round(ci.auc(Afr_2var_AsiaFit_AUC)[3],2),")",sep="")


################### DEATH in hospital -HAZ +MUAC all cases ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

hosp <- CPR.funct(data=complete_hospdeath,outcome="death_hosp",iter=100,nvars_opts=c(1:10))
hosp[["df_imps"]]
hosp[["AUC_df"]]
hosp[["calib"]]

# names      var_red
# 1               f4b_muac 2.359422e+00
# 2               base_age 1.449808e+00
# 3               f4b_resp 1.431673e+00
# 4               f4b_temp 1.275080e+00
# 5          f4a_ppl_house 1.026777e+00
# 6           f4a_drh_days 8.808421e-01
# 7       f4a_yng_children 8.210972e-01
# 8           f4a_dad_live 7.423585e-01
# 9          f4a_slp_rooms 7.064064e-01
# 10        f4a_offr_drink 6.318173e-01
# 11        f4a_disp_feces 6.287027e-01
# 12       f4b_chest_indrw 6.109841e-01
# 13        f4b_skin_flaky 5.837760e-01
# 14             f4b_admit 5.734705e-01
# 15            f4b_mental 4.743845e-01
# 16          f4a_drh_conv 4.606078e-01
# 17                  site 4.282524e-01
# 18         f4b_recommend 4.247296e-01
# 19             f3_drh_iv 4.157831e-01
# 20           f3_drh_hosp 4.147997e-01
# 21       f4a_water_avail 4.140720e-01
# 22     f4a_water_pubwell 3.832993e-01
# 23         f4a_prim_schl 3.823774e-01
# 24          f4b_abn_hair 3.695923e-01
# 25        f4a_trt_method 3.530651e-01
# 26              f4b_skin 3.487412e-01
# 27         f4a_breastfed 3.462531e-01
# 28      f4a_hometrt_herb 3.385176e-01
# 29          f4a_ani_goat 3.311919e-01
# 30          f4a_cur_skin 3.271508e-01
# 31             f4b_mouth 3.224194e-01
# 32          f4a_ms_water 2.852305e-01
# 33        f4b_under_nutr 2.814764e-01
# 34          f4a_wash_use 2.772963e-01
# 35    f4a_water_deepwell 2.762303e-01
# 36          f4a_wash_def 2.749478e-01
# 37     f4a_drh_bellypain 2.708571e-01
# 38              f4b_eyes 2.703970e-01
# 39        f4a_wash_nurse 2.661499e-01
# 40         f4a_drh_vomit 2.649437e-01
# 41      f4a_relationship 2.579353e-01
# 42      f4b_nature_stool 2.504159e-01
# 43     f4a_hometrt_othr2 2.413186e-01
# 44         f3_drh_turgor 2.347041e-01
# 45         f4a_drh_consc 2.342211e-01
# 46        f4a_max_stools 2.288839e-01
# 47      f4a_drh_restless 2.163418e-01
# 48  f4a_drh_lethrgy_miss 2.130413e-01
# 49      f4a_house_agland 2.092414e-01
# 50    f4a_cur_fastbreath 2.055265e-01
# 51           f4a_ani_cow 2.017568e-01
# 52        f4a_water_bore 2.013551e-01
# 53    f4a_water_covpwell 1.976130e-01
# 54        f4a_seek_other 1.967482e-01
# 55        f4a_seek_pharm 1.912280e-01
# 56         f4a_drh_blood 1.908814e-01
# 57      f4a_seek_outside 1.908306e-01
# 58        f4a_house_elec 1.899830e-01
# 59         f4a_wash_cook 1.884276e-01
# 60          f4a_seek_doc 1.845546e-01
# 61       f4a_ani_rodents 1.798194e-01
# 62         f4a_drh_cough 1.763008e-01
# 63         f4a_share_fac 1.758363e-01
# 64             f3_gender 1.758117e-01
# 65        f4a_house_bike 1.751114e-01
# 66      f4a_cur_restless 1.719541e-01
# 67        f4a_house_cart 1.702138e-01
# 68        f4a_wash_child 1.665559e-01
# 69           f4a_ani_cat 1.662554e-01
# 70       f4a_house_radio 1.654136e-01
# 71        f4a_drh_thirst 1.632911e-01
# 72         f4a_ani_sheep 1.631489e-01
# 73       f4a_cur_thirsty 1.617125e-01
# 74      f4a_hometrt_none 1.613194e-01
# 75       f4a_seek_healer 1.572890e-01
# 76      f4a_water_pubtap 1.540111e-01
# 77       f4a_store_water 1.539255e-01
# 78     f4a_drh_lessdrink 1.526432e-01
# 79        f4a_fuel_other 1.515315e-01
# 80     f4a_hometrt_othr1 1.443591e-01
# 81         f4a_trt_water 1.441148e-01
# 82       f4a_hometrt_ors 1.437605e-01
# 83       f4a_house_phone 1.421512e-01
# 84   f4a_hometrt_othrliq 1.411213e-01
# 85         f4a_ani_other 1.366141e-01
# 86        f4a_water_pond 1.343972e-01
# 87             f4a_floor 1.320060e-01
# 88      f4a_cur_drymouth 1.315304e-01
# 89           f4a_ani_dog 1.294208e-01
# 90          f4a_ani_fowl 1.264290e-01
# 91        f4a_hometrt_ab 1.260365e-01
# 92     f4a_hometrt_maize 1.190129e-01
# 93            f4a_ani_no 1.183673e-01
# 94         f4a_fac_waste 1.180682e-01
# 95          f4a_wash_eat 1.159276e-01
# 96         f4a_fuel_wood 1.158754e-01
# 97        f4a_house_tele 1.148117e-01
# 98     f4a_fuel_charcoal 1.117984e-01
# 99         f4a_house_car 1.082377e-01
# 100      f4a_house_scoot 1.062280e-01
# 101           f4b_rectal 9.790963e-02
# 102       f4a_fuel_grass 9.397157e-02
# 103        f4a_fuel_crop 9.321779e-02
# 104     f4a_hometrt_zinc 9.189436e-02
# 105      f4a_seek_friend 9.088906e-02
# 106     f4a_seek_privdoc 8.569828e-02
# 107     f4a_drh_prolapse 8.393596e-02
# 108    f4a_water_covwell 8.292877e-02
# 109       f4a_water_yard 7.911209e-02
# 110          f4b_bipedal 7.515665e-02
# 111      f4a_wash_animal 7.339355e-02
# 112        f4a_fuel_dung 6.880912e-02
# 113     f4a_house_fridge 5.874588e-02
# 114       f4a_seek_remdy 5.495564e-02
# 115       f4a_water_well 5.189860e-02
# 116       f4a_drh_strain 4.763830e-02
# 117     f4a_water_bought 2.881449e-02
# 118  f4a_water_shallwell 2.795182e-02
# 119       f4a_water_rain 2.137451e-02
# 120      f4a_water_river 6.871763e-03
# 121      f4a_fuel_natgas 4.407559e-03
# 122       f4a_house_none 3.583333e-03
# 123      f4a_water_house 2.552372e-03
# 124      f4a_fuel_biogas 1.500000e-03
# 125        f4a_fuel_elec 1.461905e-03
# 126        f4a_fuel_kero 1.293241e-03
# 127  f4a_water_prospring 1.163961e-03
# 128       f4a_water_othr 8.333333e-04
# 129        f4a_wash_othr 8.180620e-04
# 130     f4a_hometrt_milk 3.000000e-04
# 131     f4a_fuel_propane 3.415992e-06
# 132       f4a_house_boat 0.000000e+00
# 133        f4a_fuel_coal 0.000000e+00
# 134   f4a_water_unspring 0.000000e+00
# 135    f4b_observe_stool 0.000000e+00

# AUC         SE     lower     upper level Model nvar
# 1  0.8416327 0.01611213 0.8100535 0.8732119  0.95    LR    1
# 2  0.8191332 0.01751326 0.7848078 0.8534585  0.95    LR    2
# 3  0.8428039 0.01565498 0.8121207 0.8734871  0.95    LR    3
# 4  0.8670808 0.01398672 0.8396673 0.8944942  0.95    LR    4
# 5  0.8686742 0.01329123 0.8426238 0.8947245  0.95    LR    5
# 6  0.8569549 0.01533526 0.8268984 0.8870115  0.95    LR    6
# 7  0.8523482 0.01570061 0.8215755 0.8831208  0.95    LR    7
# 8  0.8553485 0.01519677 0.8255634 0.8851336  0.95    LR    8
# 9  0.8593645 0.01576266 0.8284703 0.8902588  0.95    LR    9
# 10 0.8507035 0.01685880 0.8176609 0.8837461  0.95    LR   10
# 11 0.6297456 0.05597339 0.5200398 0.7394515  0.95    RF    1
# 12 0.7009844 0.04114936 0.6203332 0.7816357  0.95    RF    2
# 13 0.7793255 0.02465544 0.7310017 0.8276493  0.95    RF    3
# 14 0.7897292 0.04217323 0.7070712 0.8723872  0.95    RF    4
# 15 0.7947952 0.03707206 0.7221353 0.8674551  0.95    RF    5
# 16 0.8076227 0.02862618 0.7515164 0.8637290  0.95    RF    6
# 17 0.8063233 0.02636610 0.7546467 0.8579999  0.95    RF    7
# 18 0.7947707 0.02300438 0.7496829 0.8398584  0.95    RF    8
# 19 0.8286172 0.02731528 0.7750802 0.8821542  0.95    RF    9
# 20 0.8431419 0.02018580 0.8035785 0.8827053  0.95    RF   10

# nvar   intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>  <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 -0.139   -0.996    0.526 0.998     0.436      1.55
# 2     2 -0.145   -1.01     0.523 0.914     0.375      1.45
# 3     3 -0.217   -1.10     0.474 0.827     0.354      1.28
# 4     4 -0.187   -1.06     0.497 0.881     0.427      1.35
# 5     5 -0.209   -1.10     0.486 0.846     0.410      1.29
# 6     6 -0.201   -1.10     0.499 0.772     0.352      1.19
# 7     7 -0.195   -1.08     0.501 0.755     0.335      1.17
# 8     8 -0.178   -1.07     0.519 0.751     0.345      1.16
# 9     9 -0.164   -1.05     0.534 0.756     0.353      1.16
# 10    10 -0.219   -1.12     0.492 0.720     0.346      1.10

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

temp <- hosp[["decilesCC"]][c("1","2","3","4","5","6","7","8","9","10")]
names(temp) <- c("1-var","2-var","3-var","4-var","5-var","6-var","7-var","8-var","9-var","10-var")  #renaming
#jpeg("/DeathHosp_CC_GEMS059.jpg",width=600,height=480,quality=400)
plot(x=seq(0,0.15,by=0.05),y=seq(0,0.15,by=0.05),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve:","death in hospital in cases 0-59mo in GEMS")),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(temp$`5-var`$`mean(pred_glm)`,temp$`5-var`$`mean(true)`,col="red",pch=1,cex=2,lwd=2)
points(temp$`10-var`$`mean(pred_glm)`,temp$`10-var`$`mean(true)`,col="blue",pch=2,cex=2,lwd=2)
legend("topleft",col=c("red","blue"),c("5-variable","10-variable"),pch=c(1,2),cex=1.5)
dev.off()

AUC_df <- hosp[["AUC_df"]]
#jpeg("/DeathHosp_AUCs_GEMS059.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(hosp[["nvars_opts"]])[1]],AUC_df$AUC[1:length(hosp[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death in hospital in cases 0-59mo in GEMS",
     #main=expression(paste("">=0.5,Delta,"HAZ in cases 0-59mo in GEMS")),
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(hosp[["nvars_opts"]])[1]],AUC_df$AUC[(length(hosp[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()

################### hospital death ROC curve ####
#want a single ROC surve for CV in test datasets. use "result"
result<-hosp[["result"]]

roc.data=result %>% split(.,list(.$nvar),drop=TRUE) %>% .$"10" %>% #swith this .$"x" for number of var; #subset to all the iter's (v-fold cross validations) for a given nvar (number of predictor variables in CPR)
  split(.,list(.$iter),drop=TRUE) #now have a list for each iter

#iter are the folds
# table($iter)
#>>> now 1612 per iteration

#reformatting the data so in the same layout as example for 
predict_combo=list(NA)
true_combo=list(NA)

for (i in 1:hosp[["iter"]]){ #this iter from the main RF/LR loop
  temp.list <- list(roc.data[[i]]$pred_glm)
  predict_combo <- c(predict_combo,temp.list)
  
  temp.list2 <- list(roc.data[[i]]$true)
  true_combo <- c(true_combo,temp.list2)
  
}
predict_combo=predict_combo[-1]
true_combo=true_combo[-1]
str(predict_combo)
str(true_combo)
combo <- list(predict_combo,true_combo)
names(combo) <- c("pred_glm","true")  #renaming
str(combo)

CV.roc.data <- cvAUC(predictions=combo$pred_glm, labels=combo$true) #perf is an object of class "performance" from the ROCR pckg

CV.roc.data.DeathAll.2 <- CV.roc.data
CV.roc.data.DeathAll.5 <- CV.roc.data
CV.roc.data.DeathAll.10 <- CV.roc.data
#to save
CV.roc.data.DeathAll <- list(CV.roc.data.DeathAll.2,CV.roc.data.DeathAll.5,CV.roc.data.DeathAll.10)
names(CV.roc.data.DeathAll) <- c("DA.2","DA.5","DA.10")  #renaming
str(CV.roc.data.DeathAll)
#save(CV.roc.data.DeathAll, file = "/CV.roc.data.hospDeath.Rdata")

# #find the desired points for obtained Sp for given Se
# #options(max.print=2000)
# Se.Sp<-data.frame(Se=round(roc.obj.5var$sensitivities,2),Sp=round(roc.obj.5var$specificities,2))
# Se.Sp[which(Se.Sp$Se==0.85),]
# #between the point (x0[i], y0[i]) and the point (x1[i], y1[i])
# #x0, y0, x1 = x0, y1 = y0,

#load(file = "/CV.roc.data.hospDeath.Rdata")
#(there's a better way to do this as a list...)
CV.roc.data.DA.2 <- CV.roc.data.DeathAll$DA.2
CV.roc.data.DA.5 <- CV.roc.data.DeathAll$DA.5
CV.roc.data.DA.10 <- CV.roc.data.DeathAll$DA.10

#Plot CV AUC >>> SMA: a single line of the averaged cross-validated ROC curve
#jpeg("/roc_hospDeath.jpg",width=480,height=480,quality=400)
plot(CV.roc.data.DA.2$perf, avg="vertical", main="Cross-validated ROC Curves for Deaths in Hospital",
     col="#1c61b6", lwd=2, lty=1) 
plot(CV.roc.data.DA.5$perf, avg="vertical", col="#1c61b6", lwd=2, lty=2, add=TRUE) 
plot(CV.roc.data.DA.10$perf, avg="vertical", col="#1c61b6", lwd=2, lty=3, add=TRUE) 
legend("bottomright", 
       legend = c("2-var","5-var","10-var"), 
       col = c("#1c61b6"),
       lty = c(1,2,3),
       lwd = 2)
segments(x0=0,y0=0.8,x1=0.37,y1=0.8,lty=2,col="gray")
segments(x0=0.38,y0=0.8,x1=0.38,y1=0.0,lty=2,col="gray")
segments(x0=0.27,y0=0.8,x1=0.27,y1=0.0,lty=2,col="gray")
dev.off()


################### hosp death = 0-11(age_1) in hospital -HAZ +MUAC all cases ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

hosp.011 <- CPR.funct(data=complete_hospdeath_age1,outcome="death_hosp",iter=10,nvars_opts=c(1:10))
hosp.011[["df_imps"]]
hosp.011[["AUC_df"]]
hosp.011[["calib"]]

# names      var_red
# 1               f4b_muac 1.192304e+00
# 2               f4b_temp 1.119629e+00
# 3               f4b_resp 8.482389e-01
# 4          f4a_ppl_house 5.763893e-01
# 5         f4b_skin_flaky 4.244325e-01
# 6           f4a_dad_live 4.138023e-01
# 7               base_age 4.098997e-01
# 8        f4b_chest_indrw 4.077839e-01
# 9       f4a_yng_children 4.002948e-01
# 10         f4a_slp_rooms 3.854238e-01
# 11        f4a_disp_feces 3.817142e-01
# 12          f4a_drh_days 3.545632e-01
# 13        f4a_offr_drink 3.057430e-01
# 14          f4a_drh_conv 2.896198e-01
# 15             f4b_admit 2.832462e-01
# 16         f4a_prim_schl 2.814193e-01
# 17         f4b_recommend 2.678906e-01
# 18           f3_drh_hosp 2.566975e-01
# 19            f4b_mental 2.411711e-01
# 20    f4a_cur_fastbreath 2.250333e-01
# 21                  site 2.108935e-01
# 22             f3_drh_iv 2.069962e-01
# 23        f4a_seek_other 1.944788e-01
# 24     f4a_drh_bellypain 1.858302e-01
# 25     f4a_water_pubwell 1.842851e-01
# 26              f4b_skin 1.775662e-01
# 27          f4a_ms_water 1.762162e-01
# 28        f4a_trt_method 1.716333e-01
# 29          f4a_wash_use 1.702427e-01
# 30          f4a_ani_goat 1.637029e-01
# 31     f4a_hometrt_othr2 1.593938e-01
# 32             f4b_mouth 1.489969e-01
# 33    f4a_water_deepwell 1.465117e-01
# 34        f4a_fuel_other 1.389968e-01
# 35    f4a_water_covpwell 1.364984e-01
# 36        f4a_hometrt_ab 1.318237e-01
# 37         f3_drh_turgor 1.300694e-01
# 38        f4b_under_nutr 1.297937e-01
# 39         f4a_ani_sheep 1.287116e-01
# 40        f4a_fuel_grass 1.279745e-01
# 41              f4b_eyes 1.244308e-01
# 42           f4a_ani_cat 1.228866e-01
# 43          f4a_cur_skin 1.216778e-01
# 44        f4a_wash_nurse 1.215776e-01
# 45       f4a_cur_thirsty 1.213382e-01
# 46      f4a_hometrt_herb 1.193382e-01
# 47          f4b_abn_hair 1.186646e-01
# 48     f4a_hometrt_othr1 1.167599e-01
# 49      f4a_house_agland 1.161662e-01
# 50          f4a_wash_def 1.152700e-01
# 51          f4a_seek_doc 1.132089e-01
# 52        f4a_drh_thirst 1.114536e-01
# 53      f4a_hometrt_zinc 1.100224e-01
# 54        f4a_water_pond 1.041681e-01
# 55  f4a_drh_lethrgy_miss 1.033970e-01
# 56         f4a_breastfed 1.015990e-01
# 57            f4b_rectal 1.009829e-01
# 58        f4a_house_cart 1.009660e-01
# 59         f4a_drh_consc 9.870855e-02
# 60         f4a_drh_vomit 9.621112e-02
# 61        f4a_max_stools 9.308513e-02
# 62      f4a_seek_outside 9.196845e-02
# 63       f4a_house_radio 9.127746e-02
# 64             f4a_floor 9.052821e-02
# 65         f4a_drh_blood 9.010822e-02
# 66        f4a_house_elec 8.930037e-02
# 67      f4a_relationship 8.885581e-02
# 68          f4a_wash_eat 8.879315e-02
# 69             f3_gender 8.856053e-02
# 70         f4a_drh_cough 8.830177e-02
# 71           f4a_ani_dog 8.748794e-02
# 72       f4a_ani_rodents 8.564045e-02
# 73       f4a_seek_friend 8.396621e-02
# 74      f4a_drh_restless 8.158190e-02
# 75        f4a_water_yard 7.864465e-02
# 76        f4a_seek_pharm 7.805252e-02
# 77      f4a_cur_restless 7.567377e-02
# 78            f4a_ani_no 7.537246e-02
# 79         f4a_wash_cook 7.531233e-02
# 80         f4a_fuel_crop 7.514571e-02
# 81           f4b_bipedal 7.508501e-02
# 82         f4a_share_fac 7.507718e-02
# 83        f4a_water_bore 7.492383e-02
# 84       f4a_store_water 7.381578e-02
# 85      f4a_cur_drymouth 7.365547e-02
# 86       f4a_house_scoot 7.279823e-02
# 87        f4a_wash_child 7.243928e-02
# 88           f4a_ani_cow 7.186819e-02
# 89         f4a_fac_waste 7.108978e-02
# 90       f4a_water_avail 7.046788e-02
# 91      f4a_hometrt_none 6.928337e-02
# 92         f4a_fuel_dung 6.901568e-02
# 93       f4a_hometrt_ors 6.898348e-02
# 94        f4a_house_bike 6.869670e-02
# 95         f4a_trt_water 6.261464e-02
# 96          f4a_ani_fowl 6.255819e-02
# 97      f4b_nature_stool 5.820578e-02
# 98         f4a_ani_other 5.805340e-02
# 99      f4a_drh_prolapse 5.723608e-02
# 100     f4a_seek_privdoc 5.629211e-02
# 101      f4a_house_phone 5.451345e-02
# 102    f4a_drh_lessdrink 5.414386e-02
# 103    f4a_hometrt_maize 4.934920e-02
# 104       f4a_house_tele 4.819192e-02
# 105     f4a_water_pubtap 4.777033e-02
# 106        f4a_fuel_wood 4.774081e-02
# 107    f4a_fuel_charcoal 4.726974e-02
# 108     f4a_house_fridge 4.696707e-02
# 109  f4a_hometrt_othrliq 4.328158e-02
# 110  f4a_water_shallwell 3.841524e-02
# 111        f4a_house_car 3.623178e-02
# 112     f4a_water_bought 3.225066e-02
# 113       f4a_water_rain 3.116985e-02
# 114      f4a_seek_healer 2.796215e-02
# 115      f4a_wash_animal 2.507574e-02
# 116       f4a_drh_strain 7.725491e-03
# 117      f4a_fuel_natgas 4.982927e-03
# 118      f4a_water_river 3.843918e-03
# 119        f4a_fuel_kero 3.483992e-03
# 120    f4a_water_covwell 2.798104e-03
# 121        f4a_wash_othr 2.402399e-03
# 122       f4a_house_none 1.714286e-03
# 123  f4a_water_prospring 1.333333e-03
# 124       f4a_water_well 1.114286e-03
# 125      f4a_water_house 1.064935e-03
# 126       f4a_water_othr 9.523810e-04
# 127       f4a_house_boat 5.000000e-04
# 128        f4a_fuel_coal 8.888889e-05
# 129        f4a_fuel_elec 1.696833e-05
# 130      f4a_fuel_biogas 0.000000e+00
# 131     f4a_fuel_propane 0.000000e+00
# 132   f4a_water_unspring 0.000000e+00
# 133     f4a_hometrt_milk 0.000000e+00
# 134       f4a_seek_remdy 0.000000e+00
# 135    f4b_observe_stool 0.000000e+00

# AUC          SE     lower     upper level Model nvar
# 1  0.8215838 0.025986888 0.7706505 0.8725172  0.95    LR    1
# 2  0.9024594 0.017579971 0.8680033 0.9369155  0.95    LR    2
# 3  0.9321467 0.011808049 0.9090034 0.9552901  0.95    LR    3
# 4  0.9336894 0.010420381 0.9132659 0.9541130  0.95    LR    4
# 5  0.9186926 0.013451439 0.8923282 0.9450569  0.95    LR    5
# 6  0.9218959 0.014094939 0.8942704 0.9495215  0.95    LR    6
# 7  0.9274749 0.014414492 0.8992230 0.9557268  0.95    LR    7
# 8  0.9312875 0.014836216 0.9022091 0.9603660  0.95    LR    8
# 9  0.9422815 0.009441684 0.9237762 0.9607869  0.95    LR    9
# 10 0.9250224 0.023499370 0.8789644 0.9710803  0.95    LR   10
# 11 0.5947989 0.075808321 0.4462173 0.7433805  0.95    RF    1
# 12 0.8281877 0.040062171 0.7496673 0.9067081  0.95    RF    2
# 13 0.8515065 0.033252875 0.7863321 0.9166809  0.95    RF    3
# 14 0.8667761 0.041374430 0.7856837 0.9478685  0.95    RF    4
# 15 0.8835642 0.026961117 0.8307214 0.9364070  0.95    RF    5
# 16 0.8756421 0.029051440 0.8187023 0.9325819  0.95    RF    6
# 17 0.8716640 0.022007562 0.8285300 0.9147981  0.95    RF    7
# 18 0.8710306 0.020704656 0.8304502 0.9116110  0.95    RF    8
# 19 0.8459887 0.044224964 0.7593094 0.9326680  0.95    RF    9
# 20 0.8402104 0.045138507 0.7517406 0.9286803  0.95    RF   10

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- hosp.011[["AUC_df"]]
#jpeg("/DeathHosp_AUCs_GEMS011.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[1:length(hosp.011[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death in hospital in cases 0-11mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[(length(hosp.011[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()


################### hosp death = 12-23(age_2) in hospital -HAZ +MUAC all cases ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

hosp.1223 <- CPR.funct(data=complete_hospdeath_age2,outcome="death_hosp",iter=10,nvars_opts=c(1:10))
hosp.1223[["df_imps"]]
hosp.1223[["AUC_df"]]
hosp.1223[["calib"]]

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- hosp.011[["AUC_df"]]
#jpeg("/DeathHosp_AUCs_GEMS011.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[1:length(hosp.011[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death in hospital in cases 0-11mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[(length(hosp.011[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()


################### hosp death = 24-59(age_3) in hospital -HAZ +MUAC all cases ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

hosp.2459 <- CPR.funct(data=complete_hospdeath_age3,outcome="death_hosp",iter=10,nvars_opts=c(1:10))
hosp.2459[["df_imps"]]
hosp.2459[["AUC_df"]]
hosp.2459[["calib"]]

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- hosp.011[["AUC_df"]]
#jpeg("/DeathHosp_AUCs_GEMS011.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[1:length(hosp.011[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death in hospital in cases 0-11mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[(length(hosp.011[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()



################### DEATH at home -HAZ +MUAC all cases ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

home <- CPR.funct(data=complete_homedeath,outcome="death_home",iter=100,nvars_opts=c(1:10))
home[["df_imps"]]
home[["AUC_df"]]
home[["calib"]]

# names     var_red
# 1               f4b_muac 9.339578119
# 2               f4b_resp 4.218702840
# 3               f4b_temp 4.106628137
# 4               base_age 3.838894583
# 5          f4a_ppl_house 2.859858852
# 6           f4a_drh_days 2.065860626
# 7         f4a_offr_drink 1.808649112
# 8       f4a_yng_children 1.784661774
# 9           f4b_abn_hair 1.673441863
# 10        f4b_under_nutr 1.624329022
# 11         f4a_slp_rooms 1.593789469
# 12          f4a_dad_live 1.514076414
# 13          f4a_ms_water 1.394930057
# 14        f4a_disp_feces 1.359394144
# 15         f4a_share_fac 1.344068424
# 16      f4a_relationship 1.329173155
# 17         f4a_prim_schl 1.230695796
# 18       f4a_water_avail 1.134305304
# 19              f4b_skin 1.103986775
# 20         f4a_breastfed 1.070210649
# 21      f4a_drh_prolapse 1.013580438
# 22                  site 1.001842149
# 23            f4b_mental 0.984243420
# 24         f4b_recommend 0.927394685
# 25     f4a_drh_bellypain 0.902477361
# 26             f3_drh_iv 0.886232950
# 27        f4a_trt_method 0.791444813
# 28        f4a_max_stools 0.789125406
# 29         f3_drh_turgor 0.755455374
# 30          f4a_cur_skin 0.722429108
# 31           f4b_bipedal 0.717596389
# 32    f4a_cur_fastbreath 0.713222926
# 33         f4a_drh_vomit 0.684540900
# 34        f4a_water_bore 0.672779361
# 35          f4a_wash_def 0.672039883
# 36       f4b_chest_indrw 0.653101750
# 37      f4a_hometrt_herb 0.648823023
# 38  f4a_drh_lethrgy_miss 0.635532398
# 39     f4a_water_covwell 0.632553486
# 40           f3_drh_hosp 0.630077916
# 41       f4a_ani_rodents 0.629130171
# 42       f4a_seek_healer 0.615092036
# 43       f4a_hometrt_ors 0.614216339
# 44      f4a_hometrt_none 0.609956308
# 45       f4a_cur_thirsty 0.599798167
# 46             f4b_mouth 0.594486480
# 47         f4a_wash_cook 0.589823131
# 48        f4a_house_elec 0.585865811
# 49     f4a_drh_lessdrink 0.585101713
# 50      f4a_hometrt_milk 0.579910071
# 51        f4b_skin_flaky 0.575205622
# 52           f4a_ani_cow 0.574297503
# 53             f4b_admit 0.571900870
# 54             f3_gender 0.559616025
# 55        f4a_wash_nurse 0.542300668
# 56        f4a_drh_thirst 0.531298702
# 57        f4a_wash_child 0.530268222
# 58       f4a_house_radio 0.529469332
# 59         f4a_drh_cough 0.526308859
# 60           f4a_ani_cat 0.526180355
# 61         f4a_ani_sheep 0.523201399
# 62        f4a_house_bike 0.512081262
# 63          f4a_wash_eat 0.491906173
# 64       f4a_house_phone 0.490641880
# 65      f4a_seek_outside 0.489751191
# 66         f4a_ani_other 0.485737031
# 67          f4a_ani_goat 0.481127765
# 68          f4a_wash_use 0.478974434
# 69             f4a_floor 0.466086309
# 70        f4a_house_tele 0.464674293
# 71      f4a_house_agland 0.463839415
# 72     f4a_water_pubwell 0.461124852
# 73      f4a_cur_restless 0.460075817
# 74     f4a_fuel_charcoal 0.457152264
# 75         f4a_fac_waste 0.453088075
# 76      f4a_drh_restless 0.451946696
# 77       f4a_water_river 0.450543940
# 78         f4a_drh_consc 0.442148934
# 79      f4a_seek_privdoc 0.436488687
# 80       f4a_fuel_biogas 0.429884237
# 81           f4a_ani_dog 0.423753629
# 82     f4a_hometrt_othr1 0.409807304
# 83       f4a_house_scoot 0.409218359
# 84   f4a_water_prospring 0.407717634
# 85          f4a_ani_fowl 0.405770583
# 86        f4a_drh_strain 0.401044322
# 87         f4a_trt_water 0.396600131
# 88         f4a_wash_othr 0.393137998
# 89         f4a_drh_blood 0.387794489
# 90        f4a_water_yard 0.381355014
# 91      f4a_cur_drymouth 0.376614438
# 92     f4a_hometrt_maize 0.371205713
# 93      f4a_water_bought 0.360296905
# 94      f4a_water_pubtap 0.359381433
# 95       f4a_store_water 0.359270073
# 96        f4a_house_cart 0.352548205
# 97        f4a_hometrt_ab 0.349568382
# 98        f4a_water_rain 0.328136803
# 99        f4a_water_pond 0.326455968
# 100        f4a_fuel_wood 0.316140439
# 101   f4a_water_deepwell 0.294383647
# 102        f4a_fuel_kero 0.292617919
# 103       f4a_water_well 0.290839808
# 104     f4a_house_fridge 0.274588039
# 105             f4b_eyes 0.261610816
# 106           f4a_ani_no 0.261597124
# 107      f4a_fuel_natgas 0.250543167
# 108      f4a_wash_animal 0.234354156
# 109     f4b_nature_stool 0.225287972
# 110         f4a_drh_conv 0.212539044
# 111   f4a_water_covpwell 0.211889681
# 112  f4a_hometrt_othrliq 0.206861127
# 113        f4a_fuel_elec 0.204610705
# 114         f4a_seek_doc 0.197333223
# 115           f4b_rectal 0.194753806
# 116       f4a_seek_pharm 0.193671427
# 117        f4a_house_car 0.181322161
# 118       f4a_seek_other 0.166042022
# 119       f4a_water_othr 0.163955528
# 120       f4a_house_none 0.139367080
# 121      f4a_water_house 0.137276823
# 122       f4a_fuel_grass 0.131270721
# 123  f4a_water_shallwell 0.099018187
# 124       f4a_fuel_other 0.088846330
# 125    f4a_hometrt_othr2 0.084292015
# 126        f4a_fuel_crop 0.084288883
# 127   f4a_water_unspring 0.081950420
# 128     f4a_hometrt_zinc 0.047669653
# 129        f4a_fuel_dung 0.015187219
# 130       f4a_house_boat 0.007800052
# 131       f4a_seek_remdy 0.004529112
# 132        f4a_fuel_coal 0.004480000
# 133      f4a_seek_friend 0.002166667
# 134     f4a_fuel_propane 0.001557540
# 135    f4b_observe_stool 0.000000000

# AUC          SE     lower     upper level Model nvar
# 1  0.8256319 0.004306643 0.8171910 0.8340727  0.95    LR    1
# 2  0.8395731 0.003961456 0.8318088 0.8473374  0.95    LR    2
# 3  0.8408466 0.003955121 0.8330947 0.8485985  0.95    LR    3
# 4  0.8417601 0.003985705 0.8339483 0.8495719  0.95    LR    4
# 5  0.8407724 0.004013378 0.8329064 0.8486385  0.95    LR    5
# 6  0.8384213 0.004094221 0.8303968 0.8464458  0.95    LR    6
# 7  0.8387904 0.004117773 0.8307197 0.8468610  0.95    LR    7
# 8  0.8406014 0.004105251 0.8325553 0.8486476  0.95    LR    8
# 9  0.8421800 0.004097315 0.8341494 0.8502106  0.95    LR    9
# 10 0.8439231 0.004105757 0.8358760 0.8519703  0.95    LR   10
# 11 0.8633373 0.003803389 0.8558828 0.8707918  0.95    LR   15
# 12 0.8575804 0.004010400 0.8497202 0.8654406  0.95    LR   20
# 13 0.8513378 0.004229318 0.8430484 0.8596271  0.95    LR   30
# 14 0.8447645 0.004396648 0.8361472 0.8533817  0.95    LR   40
# 15 0.8415850 0.004461649 0.8328403 0.8503297  0.95    LR   50
# 16 0.6286751 0.009579060 0.6099005 0.6474498  0.95    RF    1
# 17 0.7469393 0.005606331 0.7359511 0.7579275  0.95    RF    2
# 18 0.7914656 0.004821522 0.7820156 0.8009156  0.95    RF    3
# 19 0.8105647 0.005159962 0.8004514 0.8206780  0.95    RF    4
# 20 0.8089315 0.004830747 0.7994634 0.8183996  0.95    RF    5
# 21 0.8066731 0.004626494 0.7976053 0.8157408  0.95    RF    6
# 22 0.8137748 0.004439853 0.8050728 0.8224767  0.95    RF    7
# 23 0.8188211 0.004225256 0.8105398 0.8271025  0.95    RF    8
# 24 0.8234764 0.004376115 0.8148993 0.8320534  0.95    RF    9
# 25 0.8273637 0.004224487 0.8190839 0.8356436  0.95    RF   10
# 26 0.8531264 0.003786257 0.8457055 0.8605473  0.95    RF   15
# 27 0.8658069 0.003464014 0.8590175 0.8725962  0.95    RF   20
# 28 0.8731369 0.003289706 0.8666892 0.8795846  0.95    RF   30
# 29 0.8743391 0.003375291 0.8677236 0.8809545  0.95    RF   40
# 30 0.8740122 0.003429614 0.8672903 0.8807341  0.95    RF   50

# nvar   intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <int>  <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     1 0.0395   -0.394    0.424 1.06      0.754      1.38
# 2     2 0.0502   -0.384    0.436 1.04      0.733      1.35
# 3     3 0.0457   -0.390    0.433 1.04      0.743      1.35
# 4     4 0.0420   -0.394    0.430 1.02      0.734      1.33
# 5     5 0.0414   -0.395    0.429 1.01      0.725      1.32
# 6     6 0.0512   -0.386    0.440 1.00      0.719      1.31
# 7     7 0.0498   -0.390    0.441 0.940     0.669      1.23
# 8     8 0.0523   -0.388    0.445 0.930     0.663      1.21
# 9     9 0.0541   -0.387    0.447 0.931     0.668      1.21
# 10    10 0.0560   -0.387    0.451 0.932     0.675      1.21

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

temp <- home[["decilesCC"]][c("1","2","3","4","5","6","7","8","9","10")]
names(temp) <- c("1-var","2-var","3-var","4-var","5-var","6-var","7-var","8-var","9-var","10-var")  #renaming
#jpeg("/DeathHome_CC_GEMS059.jpg",width=600,height=480,quality=400)
plot(x=seq(0,0.15,by=0.05),y=seq(0,0.15,by=0.05),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve:","death after discharge in cases 0-59mo in GEMS")),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(temp$`5-var`$`mean(pred_glm)`,temp$`5-var`$`mean(true)`,col="red",pch=1,cex=2,lwd=2)
points(temp$`10-var`$`mean(pred_glm)`,temp$`10-var`$`mean(true)`,col="blue",pch=2,cex=2,lwd=2)
legend("topleft",col=c("red","blue"),c("5-variable","10-variable"),pch=c(1,2),cex=1.5)
dev.off()

AUC_df <- home[["AUC_df"]]
#jpeg("/DeathHome_AUCs_GEMS059.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(home[["nvars_opts"]])[1]],AUC_df$AUC[1:length(home[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death after discharge in cases 0-59mo in GEMS",
     #main=expression(paste("">=0.5,Delta,"HAZ in cases 0-59mo in GEMS")),
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(home[["nvars_opts"]])[1]],AUC_df$AUC[(length(home[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()



################### home death ROC curve ####
#want a single ROC surve for CV in test datasets. use "result"
result<-home[["result"]]

roc.data=result %>% split(.,list(.$nvar),drop=TRUE) %>% .$"2" %>% #swith this .$"x" for number of var; #subset to all the iter's (v-fold cross validations) for a given nvar (number of predictor variables in CPR)
  split(.,list(.$iter),drop=TRUE) #now have a list for each iter

#iter are the folds
# table($iter)
#>>> now 1612 per iteration

#reformatting the data so in the same layout as example for 
predict_combo=list(NA)
true_combo=list(NA)

for (i in 1:home[["iter"]]){ #this iter from the main RF/LR loop
  temp.list <- list(roc.data[[i]]$pred_glm)
  predict_combo <- c(predict_combo,temp.list)
  
  temp.list2 <- list(roc.data[[i]]$true)
  true_combo <- c(true_combo,temp.list2)
  
}
predict_combo=predict_combo[-1]
true_combo=true_combo[-1]
str(predict_combo)
str(true_combo)
combo <- list(predict_combo,true_combo)
names(combo) <- c("pred_glm","true")  #renaming
str(combo)

CV.roc.data <- cvAUC(predictions=combo$pred_glm, labels=combo$true) #perf is an object of class "performance" from the ROCR pckg

CV.roc.data.DeathAll.2 <- CV.roc.data
CV.roc.data.DeathAll.5 <- CV.roc.data
CV.roc.data.DeathAll.10 <- CV.roc.data
#to save
CV.roc.data.DeathAll <- list(CV.roc.data.DeathAll.2,CV.roc.data.DeathAll.5,CV.roc.data.DeathAll.10)
names(CV.roc.data.DeathAll) <- c("DA.2","DA.5","DA.10")  #renaming
str(CV.roc.data.DeathAll)
#save(CV.roc.data.DeathAll, file = "/CV.roc.data.homeDeath.Rdata")

# #find the desired points for obtained Sp for given Se
# #options(max.print=2000)
# Se.Sp<-data.frame(Se=round(roc.obj.5var$sensitivities,2),Sp=round(roc.obj.5var$specificities,2))
# Se.Sp[which(Se.Sp$Se==0.85),]
# #between the point (x0[i], y0[i]) and the point (x1[i], y1[i])
# #x0, y0, x1 = x0, y1 = y0,

#load(file = "/CV.roc.data.homeDeath.Rdata")
#(there's a better way to do this as a list...)
CV.roc.data.DA.2 <- CV.roc.data.DeathAll$DA.2
CV.roc.data.DA.5 <- CV.roc.data.DeathAll$DA.5
CV.roc.data.DA.10 <- CV.roc.data.DeathAll$DA.10

#Plot CV AUC >>> SMA: a single line of the averaged cross-validated ROC curve
#jpeg("/roc_homeDeath.jpg",width=480,height=480,quality=400)
plot(CV.roc.data.DA.2$perf, avg="vertical", main="Cross-validated ROC Curves for Deaths after Discharge",
     col="#1c61b6", lwd=2, lty=1) 
plot(CV.roc.data.DA.5$perf, avg="vertical", col="#1c61b6", lwd=2, lty=2, add=TRUE) 
plot(CV.roc.data.DA.10$perf, avg="vertical", col="#1c61b6", lwd=2, lty=3, add=TRUE) 
legend("bottomright", 
       legend = c("2-var","5-var","10-var"), 
       col = c("#1c61b6"),
       lty = c(1,2,3),
       lwd = 2)
segments(x0=0,y0=0.8,x1=0.23,y1=0.8,lty=2,col="gray")
segments(x0=0.23,y0=0.8,x1=0.23,y1=0.0,lty=2,col="gray")
dev.off()


################### home death = 0-11(age_1) in hospital -HAZ +MUAC all cases ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

home.011 <- CPR.funct(data=complete_homedeath_age1,outcome="death_home",iter=10,nvars_opts=c(1:10))
home.011[["df_imps"]]
home.011[["AUC_df"]]
home.011[["calib"]]

# names      var_red
# 1               f4b_muac 4.9617261596
# 2               f4b_resp 2.7970324838
# 3               f4b_temp 2.6252133362
# 4               base_age 1.9572748110
# 5          f4a_ppl_house 1.4989540989
# 6           f4a_drh_days 1.2685373238
# 7         f4a_offr_drink 1.1763190882
# 8           f4a_ms_water 1.0564576701
# 9          f4a_share_fac 0.9707633057
# 10          f4a_dad_live 0.9275068893
# 11      f4a_drh_prolapse 0.8979684617
# 12      f4a_yng_children 0.8975568071
# 13      f4a_relationship 0.8499393216
# 14         f4a_slp_rooms 0.8426360612
# 15         f4a_breastfed 0.8259140692
# 16       f4a_water_avail 0.8193183176
# 17        f4b_under_nutr 0.8077552354
# 18        f4a_disp_feces 0.8065822687
# 19       f4b_chest_indrw 0.7761822264
# 20          f4b_abn_hair 0.7623581448
# 21         f4a_prim_schl 0.7497950490
# 22         f4b_recommend 0.6925278519
# 23              f4b_skin 0.6901192451
# 24            f4b_mental 0.6629180575
# 25                  site 0.6594180700
# 26        f4a_trt_method 0.5694687197
# 27    f4a_cur_fastbreath 0.5602814483
# 28          f4a_cur_skin 0.5181000297
# 29        f4a_max_stools 0.5092666199
# 30     f4a_drh_bellypain 0.4881738072
# 31      f4a_seek_privdoc 0.4783929834
# 32         f3_drh_turgor 0.4680922896
# 33     f4a_water_covwell 0.4586249320
# 34          f4a_wash_def 0.4249380047
# 35         f4a_drh_vomit 0.4194339329
# 36             f3_drh_iv 0.4182115558
# 37       f4a_water_river 0.4142265214
# 38           f3_drh_hosp 0.4087813349
# 39             f4b_mouth 0.4053962877
# 40             f4b_admit 0.4022480744
# 41        f4a_drh_strain 0.3987001344
# 42      f4a_hometrt_none 0.3951990071
# 43        f4a_wash_nurse 0.3937404922
# 44         f4a_wash_cook 0.3839771272
# 45          f4a_ani_goat 0.3808342355
# 46         f4a_fac_waste 0.3616600939
# 47        f4a_water_bore 0.3600033092
# 48  f4a_drh_lethrgy_miss 0.3585251271
# 49       f4a_cur_thirsty 0.3566736399
# 50       f4a_ani_rodents 0.3561425200
# 51        f4a_wash_child 0.3552348147
# 52         f4a_drh_cough 0.3551152139
# 53           f4b_bipedal 0.3549976002
# 54        f4a_hometrt_ab 0.3542353899
# 55        f4a_drh_thirst 0.3500666297
# 56          f4a_wash_eat 0.3492763347
# 57           f4a_ani_cat 0.3469033156
# 58     f4a_drh_lessdrink 0.3410071817
# 59       f4a_hometrt_ors 0.3390974496
# 60   f4a_water_prospring 0.3370537930
# 61         f4a_ani_sheep 0.3271240440
# 62     f4a_fuel_charcoal 0.3239484787
# 63        f4a_house_elec 0.3229890899
# 64      f4a_seek_outside 0.3224281256
# 65             f3_gender 0.3211353025
# 66       f4a_fuel_biogas 0.3104863715
# 67        f4a_house_bike 0.3094878402
# 68           f4a_ani_cow 0.3076073773
# 69        f4a_house_tele 0.3057475641
# 70       f4a_house_radio 0.3045795293
# 71         f4a_trt_water 0.3014596806
# 72         f4a_ani_other 0.3014050715
# 73       f4a_seek_healer 0.2976559928
# 74      f4a_house_agland 0.2940441463
# 75      f4a_drh_restless 0.2918547832
# 76          f4a_drh_conv 0.2904422875
# 77      f4a_water_bought 0.2880284939
# 78             f4a_floor 0.2848877148
# 79           f4a_ani_dog 0.2826418399
# 80         f4a_wash_othr 0.2826119497
# 81     f4a_water_pubwell 0.2794335873
# 82        f4a_house_cart 0.2787365174
# 83          f4a_wash_use 0.2734787269
# 84      f4a_cur_drymouth 0.2708454122
# 85         f4a_fuel_kero 0.2703194379
# 86      f4a_hometrt_herb 0.2696133700
# 87       f4a_store_water 0.2613928320
# 88       f4a_house_phone 0.2573307034
# 89              f4b_eyes 0.2560480956
# 90      f4a_cur_restless 0.2553074557
# 91        f4a_water_yard 0.2544962860
# 92          f4a_ani_fowl 0.2489249810
# 93      f4a_water_pubtap 0.2480013859
# 94         f4a_fuel_wood 0.2479442133
# 95         f4a_drh_blood 0.2477542054
# 96        f4a_water_pond 0.2469461203
# 97   f4a_hometrt_othrliq 0.2418759754
# 98       f4a_house_scoot 0.2346706449
# 99         f4a_fuel_elec 0.2309813376
# 100    f4a_hometrt_maize 0.2250578594
# 101     f4a_house_fridge 0.2013881225
# 102    f4a_hometrt_othr1 0.2005447623
# 103       f4b_skin_flaky 0.1994301105
# 104   f4a_water_deepwell 0.1994207945
# 105        f4a_drh_consc 0.1936250133
# 106      f4a_fuel_natgas 0.1799369965
# 107         f4a_seek_doc 0.1781214181
# 108       f4a_water_rain 0.1726939234
# 109           f4a_ani_no 0.1562195982
# 110        f4a_house_car 0.1511944278
# 111       f4a_water_othr 0.1384996690
# 112       f4a_fuel_grass 0.1304015260
# 113           f4b_rectal 0.1251860808
# 114      f4a_water_house 0.1064378854
# 115       f4a_seek_other 0.1024499697
# 116      f4a_wash_animal 0.1017038151
# 117  f4a_water_shallwell 0.1004825519
# 118       f4a_water_well 0.1003172488
# 119   f4a_water_covpwell 0.0995560904
# 120        f4a_fuel_crop 0.0984716727
# 121       f4a_fuel_other 0.0888910157
# 122     f4a_hometrt_milk 0.0614902132
# 123    f4a_hometrt_othr2 0.0561623887
# 124       f4a_seek_pharm 0.0327041555
# 125        f4a_fuel_dung 0.0320065945
# 126     f4b_nature_stool 0.0119388621
# 127       f4a_house_boat 0.0086495169
# 128     f4a_hometrt_zinc 0.0057221807
# 129       f4a_seek_remdy 0.0042085863
# 130        f4a_fuel_coal 0.0023300232
# 131     f4a_fuel_propane 0.0017846443
# 132   f4a_water_unspring 0.0012642857
# 133       f4a_house_none 0.0009142857
# 134      f4a_seek_friend 0.0007533613
# 135    f4b_observe_stool 0.0000000000

# AUC         SE     lower     upper level Model nvar
# 1  0.7863790 0.01780670 0.7514785 0.8212795  0.95    LR    1
# 2  0.8029893 0.01636613 0.7709123 0.8350664  0.95    LR    2
# 3  0.8134824 0.01603557 0.7820532 0.8449115  0.95    LR    3
# 4  0.8126710 0.01643594 0.7804571 0.8448848  0.95    LR    4
# 5  0.8129537 0.01657895 0.7804595 0.8454478  0.95    LR    5
# 6  0.8082282 0.01662093 0.7756517 0.8408046  0.95    LR    6
# 7  0.8085930 0.01677517 0.7757143 0.8414717  0.95    LR    7
# 8  0.8043300 0.01685500 0.7712948 0.8373652  0.95    LR    8
# 9  0.8053984 0.01703653 0.7720074 0.8387894  0.95    LR    9
# 10 0.8016220 0.01724304 0.7678262 0.8354177  0.95    LR   10
# 11 0.5896979 0.03736423 0.5164654 0.6629305  0.95    RF    1
# 12 0.6833089 0.02289927 0.6384272 0.7281907  0.95    RF    2
# 13 0.7658223 0.01887424 0.7288295 0.8028152  0.95    RF    3
# 14 0.7545564 0.02162988 0.7121626 0.7969502  0.95    RF    4
# 15 0.7339397 0.02170245 0.6914036 0.7764757  0.95    RF    5
# 16 0.7239649 0.02091937 0.6829637 0.7649662  0.95    RF    6
# 17 0.7395141 0.01949319 0.7013081 0.7777200  0.95    RF    7
# 18 0.7404153 0.02054218 0.7001533 0.7806772  0.95    RF    8
# 19 0.7441537 0.02218767 0.7006667 0.7876407  0.95    RF    9
# 20 0.7501279 0.02077481 0.7094100 0.7908458  0.95    RF   10

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- home.011[["AUC_df"]]
#jpeg("/DeathHome_AUCs_GEMS011.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(home.011[["nvars_opts"]])[1]],AUC_df$AUC[1:length(home.011[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death in hospital in cases 0-11mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(home.011[["nvars_opts"]])[1]],AUC_df$AUC[(length(home.011[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()


################### home death = 12-23(age_2) in hospital -HAZ +MUAC all cases ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

hosp.1223 <- CPR.funct(data=complete_homedeath_age2,outcome="death_home",iter=10,nvars_opts=c(1:10))
hosp.1223[["df_imps"]]
hosp.1223[["AUC_df"]]
hosp.1223[["calib"]]

# names      var_red
# 1               f4b_muac 3.0880691878
# 2          f4a_ppl_house 1.0751010201
# 3               f4b_temp 0.9866620478
# 4               f4b_resp 0.9097460027
# 5       f4a_hometrt_milk 0.6976084622
# 6       f4a_yng_children 0.6425276923
# 7               base_age 0.6222009032
# 8          f4a_slp_rooms 0.6194469925
# 9         f4a_offr_drink 0.5549150882
# 10          f4a_drh_days 0.5434642905
# 11          f4b_abn_hair 0.4618218460
# 12        f4b_under_nutr 0.4526143914
# 13              f4b_skin 0.4478621012
# 14      f4a_hometrt_herb 0.3723895831
# 15             f3_drh_iv 0.3666149853
# 16          f4a_dad_live 0.3578184474
# 17                  site 0.3410808899
# 18        f4a_water_bore 0.3266594932
# 19         f4a_drh_consc 0.3116178985
# 20         f4a_breastfed 0.2984341140
# 21         f3_drh_turgor 0.2982737755
# 22            f4b_mental 0.2889058066
# 23          f4a_ms_water 0.2762351711
# 24       f4a_water_avail 0.2730937597
# 25        f4a_disp_feces 0.2726407950
# 26         f4a_prim_schl 0.2623039505
# 27     f4a_drh_bellypain 0.2600509169
# 28           f4b_bipedal 0.2423398104
# 29     f4a_drh_lessdrink 0.2362374529
# 30        f4a_water_well 0.2280899912
# 31          f4a_wash_def 0.2203218106
# 32     f4a_hometrt_othr1 0.2150086185
# 33  f4a_drh_lethrgy_miss 0.2120423754
# 34       f4a_house_radio 0.2052537043
# 35         f4b_recommend 0.2046017536
# 36        f4a_trt_method 0.2013881556
# 37           f3_drh_hosp 0.1961225530
# 38         f4a_drh_vomit 0.1916175916
# 39      f4a_relationship 0.1903206765
# 40       f4a_hometrt_ors 0.1897323515
# 41          f4a_ani_goat 0.1879214428
# 42      f4a_drh_prolapse 0.1832589104
# 43        f4a_seek_pharm 0.1823220026
# 44         f4a_wash_cook 0.1809639523
# 45         f4a_share_fac 0.1793943902
# 46          f4a_cur_skin 0.1747585695
# 47     f4a_water_covwell 0.1714907944
# 48        f4a_house_elec 0.1662442307
# 49      f4b_nature_stool 0.1644220130
# 50         f4a_ani_other 0.1614864272
# 51      f4a_hometrt_none 0.1614469735
# 52      f4a_seek_outside 0.1593924067
# 53             f4b_mouth 0.1575346582
# 54        f4a_wash_child 0.1565603190
# 55    f4a_cur_fastbreath 0.1556305349
# 56        f4a_drh_thirst 0.1552467352
# 57         f4a_drh_cough 0.1543993601
# 58     f4a_water_pubwell 0.1523796859
# 59             f3_gender 0.1507700846
# 60       f4a_ani_rodents 0.1485977450
# 61             f4b_admit 0.1483400798
# 62           f4a_ani_cat 0.1466429360
# 63      f4a_cur_restless 0.1458333443
# 64      f4a_house_agland 0.1426328712
# 65             f4a_floor 0.1401507643
# 66        f4a_max_stools 0.1397701029
# 67       f4a_fuel_biogas 0.1393444475
# 68       f4a_seek_healer 0.1368373716
# 69       f4a_cur_thirsty 0.1340800322
# 70         f4a_drh_blood 0.1314293871
# 71        f4a_water_rain 0.1310326009
# 72        f4a_wash_nurse 0.1305763441
# 73         f4a_ani_sheep 0.1293390380
# 74      f4a_water_pubtap 0.1276046022
# 75           f4a_ani_cow 0.1265670314
# 76     f4a_fuel_charcoal 0.1215481986
# 77        f4a_house_tele 0.1212672696
# 78       f4a_house_phone 0.1209036640
# 79   f4a_water_prospring 0.1181037775
# 80            f4b_rectal 0.1160635558
# 81        f4a_house_bike 0.1156698919
# 82        f4a_house_none 0.1081358213
# 83         f4a_trt_water 0.1067192432
# 84      f4a_cur_drymouth 0.1017674391
# 85      f4a_drh_restless 0.1000136096
# 86    f4a_water_covpwell 0.0994026723
# 87         f4a_wash_othr 0.0986018535
# 88       f4a_house_scoot 0.0982316191
# 89          f4a_wash_use 0.0966868460
# 90    f4a_water_unspring 0.0962856265
# 91        f4b_skin_flaky 0.0960224657
# 92        f4a_house_cart 0.0918416000
# 93        f4a_water_yard 0.0878168280
# 94           f4a_ani_dog 0.0877379780
# 95    f4a_water_deepwell 0.0822290820
# 96       f4a_store_water 0.0777804952
# 97     f4a_hometrt_maize 0.0777299592
# 98        f4a_water_pond 0.0775189915
# 99            f4a_ani_no 0.0746760037
# 100     f4a_house_fridge 0.0733104767
# 101         f4a_wash_eat 0.0723385935
# 102        f4a_fuel_wood 0.0716040191
# 103         f4a_ani_fowl 0.0687196771
# 104      f4a_water_river 0.0632657430
# 105     f4a_hometrt_zinc 0.0529706612
# 106     f4a_water_bought 0.0506740534
# 107     f4a_seek_privdoc 0.0488055862
# 108        f4a_fac_waste 0.0450546197
# 109      f4a_fuel_natgas 0.0349581237
# 110             f4b_eyes 0.0301505867
# 111       f4a_hometrt_ab 0.0261388421
# 112        f4a_house_car 0.0252769886
# 113       f4a_drh_strain 0.0147849651
# 114      f4b_chest_indrw 0.0065048937
# 115       f4a_fuel_grass 0.0060126713
# 116      f4a_water_house 0.0056557008
# 117        f4a_fuel_crop 0.0035392768
# 118         f4a_drh_conv 0.0034779249
# 119    f4a_hometrt_othr2 0.0032600148
# 120       f4a_seek_remdy 0.0031363636
# 121      f4a_wash_animal 0.0025681451
# 122       f4a_seek_other 0.0025500000
# 123        f4a_fuel_kero 0.0017861324
# 124  f4a_hometrt_othrliq 0.0013184926
# 125        f4a_fuel_dung 0.0011703617
# 126      f4a_seek_friend 0.0008366013
# 127       f4a_fuel_other 0.0007285714
# 128        f4a_fuel_coal 0.0005952381
# 129     f4a_fuel_propane 0.0003333333
# 130  f4a_water_shallwell 0.0003050649
# 131         f4a_seek_doc 0.0001134454
# 132       f4a_house_boat 0.0000000000
# 133        f4a_fuel_elec 0.0000000000
# 134       f4a_water_othr 0.0000000000
# 135    f4b_observe_stool 0.0000000000

# AUC         SE     lower     upper level Model nvar
# 1  0.8486663 0.02424976 0.8011376 0.8961949  0.95    LR    1
# 2  0.8507111 0.02347867 0.8046937 0.8967284  0.95    LR    2
# 3  0.8530747 0.02364577 0.8067298 0.8994196  0.95    LR    3
# 4  0.8533942 0.02123654 0.8117713 0.8950171  0.95    LR    4
# 5  0.8542274 0.02074035 0.8135770 0.8948777  0.95    LR    5
# 6  0.8533611 0.01933277 0.8154696 0.8912526  0.95    LR    6
# 7  0.8566357 0.01894184 0.8195104 0.8937610  0.95    LR    7
# 8  0.8552155 0.01985324 0.8163039 0.8941272  0.95    LR    8
# 9  0.8486091 0.02018134 0.8090544 0.8881638  0.95    LR    9
# 10 0.8461550 0.02149103 0.8040334 0.8882766  0.95    LR   10
# 11 0.7334413 0.05889402 0.6180111 0.8488714  0.95    RF    1
# 12 0.7635018 0.03874078 0.6875713 0.8394324  0.95    RF    2
# 13 0.7712061 0.03516738 0.7022793 0.8401329  0.95    RF    3
# 14 0.7907399 0.03697952 0.7182614 0.8632184  0.95    RF    4
# 15 0.7794093 0.03800737 0.7049163 0.8539024  0.95    RF    5
# 16 0.8061115 0.03247258 0.7424664 0.8697566  0.95    RF    6
# 17 0.8138204 0.02729177 0.7603295 0.8673113  0.95    RF    7
# 18 0.8253956 0.02551381 0.7753894 0.8754017  0.95    RF    8
# 19 0.8500862 0.02626645 0.7986049 0.9015675  0.95    RF    9
# 20 0.8711154 0.02007705 0.8317651 0.9104657  0.95    RF   10

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- hosp.011[["AUC_df"]]
#jpeg("/DeathHosp_AUCs_GEMS011.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[1:length(hosp.011[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death in hospital in cases 0-11mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[(length(hosp.011[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()


################### home death = 24-59(age_3) in hospital -HAZ +MUAC all cases ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

hosp.2459 <- CPR.funct(data=complete_homedeath_age3,outcome="death_home",iter=10,nvars_opts=c(1:10))
hosp.2459[["df_imps"]]
hosp.2459[["AUC_df"]]
hosp.2459[["calib"]]

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

AUC_df <- hosp.011[["AUC_df"]]
#jpeg("/DeathHosp_AUCs_GEMS011.jpg",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[1:length(hosp.011[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Death in hospital in cases 0-11mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(hosp.011[["nvars_opts"]])[1]],AUC_df$AUC[(length(hosp.011[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()




####################
#KIDMS for rederive growth faltering prediction
################### import, org Kilifi data ####
kilifi_orig <- read.csv("/Kilifi.csv", header=T)

kilifi <- kilifi_orig %>% rename(f4b_muac=muac,
                                 f4b_resp=resp_rate,
                                 f4b_temp=temp_axilla,
                                 base_age=agemths) %>%
                          mutate(death_all=ifelse(is.na(days_death),0,1),
                                 death_hosp=ifelse((!is.na(days_death)&days_death<=hosp_instay),1,0),
                                 death_home=ifelse((!is.na(days_death)&days_death>hosp_instay),1,0),
                                 agegroup=ifelse(base_age<12,1,
                                                 ifelse((base_age>=12 & base_age<24),2,3))) %>%
                          #select(death_all,death_hosp,death_home,f4b_muac,f4b_resp,f4b_temp,base_age,days_death,hosp_instay) %>%
                          #filter(days_death<0) %>% #dropping kids w/ negative death days, some of whom still have positive host_instay, don't know what that means    
                            #checked w/ Moses, email Aug 19 2022. days_death calc from admission and hosp_instay. have manually checked which observations dropped, ok as is
                            #I don't have specific dates to calc non-neg days_death, but create binary outcome for all and ok how calculates out
                          select(death_all,death_hosp,death_home,f4b_muac,f4b_resp,f4b_temp,base_age,agegroup) %>%
                          filter(base_age>=0 & f4b_resp>=18) %>% #ask Daniel if like this resp rate criteria
                          na.omit() #drop any missing
#3009 to 2901
#most dropped are missing resp or MUAC

table(is.na(kilifi_orig$days_death))
table(kilifi$death_all)
summary(complete_anydeath$f4b_muac)
summary(kilifi$f4b_muac)
summary(complete_anydeath$f4b_resp)
summary(kilifi$f4b_resp) #how can there by 0 here?
summary(complete_anydeath$f4b_temp)
summary(kilifi$f4b_temp)
summary(complete_anydeath$base_age)
summary(kilifi$base_age) #one of these is negative >>> stillbirth?

################### descriptive for pub ####
table(kilifi$death_all)
# 0    1 
# 2642  259 

table(kilifi$agegroup)
table(kilifi$death_hosp,kilifi$agegroup)
table(kilifi$death_home,kilifi$agegroup)
table(kilifi$death_all,kilifi$agegroup)

table(kilifi$death_hosp)
table(kilifi$death_home)

################### external validation 4var: MAIN any death ####
GEMS_glm_death <- glm(death_all~f4b_muac+f4b_resp+f4b_temp+base_age,
                   data=complete_anydeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -11.709354   2.946540  -3.974 7.07e-05 ***
#   f4b_muac     -0.837129   0.055226 -15.158  < 2e-16 ***
#   f4b_resp      0.023904   0.007842   3.048   0.0023 ** 
#   f4b_temp      0.461924   0.081235   5.686 1.30e-08 ***
#   base_age      0.022321   0.008879   2.514   0.0119 * 
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathAll_4var.Rdata")
round(exp(coef(glm(death_all~f4b_muac+f4b_resp+f4b_temp+base_age,
                   data=complete_anydeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac    f4b_resp    f4b_temp    base_age 
# 0.0000      0.4330      1.0242      1.5871      1.0226 

GEMS_4var<-complete_anydeath %>% select(death_all,f4b_muac,f4b_resp,f4b_temp,base_age)
GEMS_4var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_4var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_4var$pred_glm,labels=GEMS_4var$death_all)
GEMS_AUC
# [1] 0.8558919
#using pROC package since need CI
GEMS_AUC <- roc(response=GEMS_4var$death_all,predictor=GEMS_4var$pred_glm)
paste(round(GEMS_AUC$auc,2)," (",
      round(ci.auc(GEMS_AUC)[1],2),", ",
      round(ci.auc(GEMS_AUC)[3],2),")",sep="")
# "0.86 (0.83, 0.89)"


GEMS_decilesCC <- GEMS_4var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_all),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathAll_4var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_all)
kilifi_AUC
# [1] 0.7374956
kilifi_AUC <- roc(response=kilifi$death_all,predictor=kilifi$GEMS_pred_glm)
paste(round(kilifi_AUC$auc,2)," (",
      round(ci.auc(kilifi_AUC)[1],2),", ",
      round(ci.auc(kilifi_AUC)[3],2),")",sep="")
# "0.73 (0.7, 0.77)"

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_all),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathAll_4var.Rdata")

#jpeg("/DeathAll_CC_ExternalVal_GEMS059.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: All Deaths")),
     xlim=c(0,0.4),ylim=c(0,0.4))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_all)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_all)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.86","Kilifi data (0-59mo), AUC=0.74"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_all~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_all~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.62388    0.07462    8.36   <2e-16 ***
#   2.5 %    97.5 % 
#   0.4755069 0.7681226 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -0.46235    0.13838  -3.341 0.000834 ***
#   log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.57404    0.04347  13.206  < 2e-16 ***
#   2.5 %     97.5 %
#   (Intercept)                            -0.7350533 -0.1920384
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.4897328  0.6602851


################### external validation 2var: MAIN any death ####
GEMS_glm_death <- glm(death_all~f4b_muac+f4b_resp,
                      data=complete_anydeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.560173   0.716579   6.364 1.97e-10 ***
#   f4b_muac    -0.749498   0.049352 -15.187  < 2e-16 ***
#   f4b_resp     0.030194   0.007255   4.162 3.15e-05 ***
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathAll_2var.Rdata")
round(exp(coef(glm(death_all~f4b_muac+f4b_resp,
                   data=complete_anydeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac    f4b_resp 
# 95.6000      0.4726      1.0307 

GEMS_2var<-complete_anydeath %>% select(death_all,f4b_muac,f4b_resp)
GEMS_2var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_2var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_2var$pred_glm,labels=GEMS_2var$death_all)
GEMS_AUC
# [1] 0.8505391
#using pROC package since need CI
GEMS_AUC <- roc(response=GEMS_2var$death_all,predictor=GEMS_2var$pred_glm)
paste(round(GEMS_AUC$auc,2)," (",
      round(ci.auc(GEMS_AUC)[1],2),", ",
      round(ci.auc(GEMS_AUC)[3],2),")",sep="")
# "0.85 (0.82, 0.88)"


GEMS_decilesCC <- GEMS_2var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_all),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathAll_2var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_all)
kilifi_AUC
# [1] 0.7401151
kilifi_AUC <- roc(response=kilifi$death_all,predictor=kilifi$GEMS_pred_glm)
paste(round(kilifi_AUC$auc,2)," (",
      round(ci.auc(kilifi_AUC)[1],2),", ",
      round(ci.auc(kilifi_AUC)[3],2),")",sep="")
# "0.74 (0.71, 0.77)"

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_all),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathAll_2var.Rdata")

# load(file = "/GEMS_decilesCC_deathAll_2var.Rdata")
# load(file = "/kilifi_decilesCC_deathAll_2var.Rdata")

#tiff("/DeathAll_CC_ExternalVal_GEMS059_2var.tif",units="px",width=1980,height=1650,res=300)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: All Deaths")),
     xlim=c(0,0.4),ylim=c(0,0.4))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_all)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_all)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.85 (95% CI: 0.82, 0.88)","Kilifi data (0-59mo), AUC=0.74 (95% CI: 0.71, 0.77"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_all~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_all~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.82305    0.07389   11.14   <2e-16 ***

# Waiting for profiling to be done...
# 2.5 %    97.5 % 
#   0.6760895 0.9658204 

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -0.26483    0.15123  -1.751   0.0799 .  
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.60762    0.04574  13.286   <2e-16 ***

# Waiting for profiling to be done...
# 2.5 %     97.5 %
#   (Intercept)                            -0.5624790 0.03098968
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.5188765 0.69832576


################### external validation 1var: MAIN any death ####
GEMS_glm_death <- glm(death_all~f4b_muac,
                      data=complete_anydeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   6.2958     0.5839   10.78   <2e-16 ***
#   f4b_muac     -0.7910     0.0483  -16.38   <2e-16 ***
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathAll_1var.Rdata")
round(exp(coef(glm(death_all~f4b_muac,
                   data=complete_anydeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac 
# 542.2858      0.4534 

GEMS_1var<-complete_anydeath %>% select(death_all,f4b_muac)
GEMS_1var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_1var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_1var$pred_glm,labels=GEMS_1var$death_all)
GEMS_AUC
# [1] 0.8307989
#using pROC package since need CI
GEMS_AUC <- roc(response=GEMS_1var$death_all,predictor=GEMS_1var$pred_glm)
paste(round(GEMS_AUC$auc,2)," (",
      round(ci.auc(GEMS_AUC)[1],2),", ",
      round(ci.auc(GEMS_AUC)[3],2),")",sep="")
# "0.83 (0.8, 0.86)"


GEMS_decilesCC <- GEMS_1var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_all),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathAll_1var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_all)
kilifi_AUC
# [1] 0.7179465
kilifi_AUC <- roc(response=kilifi$death_all,predictor=kilifi$GEMS_pred_glm)
paste(round(kilifi_AUC$auc,2)," (",
      round(ci.auc(kilifi_AUC)[1],2),", ",
      round(ci.auc(kilifi_AUC)[3],2),")",sep="")
# "0.72 (0.68, 0.75)"

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_all),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathAll_1var.Rdata")

#jpeg("/DeathAll_CC_ExternalVal_GEMS059_1var.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: All Deaths")),
     xlim=c(0,0.4),ylim=c(0,0.4))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_all)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_all)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.83","Kilifi data (0-59mo), AUC=0.72"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_all~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_all~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.89631    0.07338   12.21   <2e-16 ***
#   2.5 %    97.5 % 
#   0.7503625 1.0380892 
# 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -0.33633    0.15754  -2.135   0.0328 *  
#   log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.57094    0.04627  12.338   <2e-16 ***
#   2.5 %      97.5 %
#   (Intercept)                            -0.6470885 -0.02885205
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.4809503  0.66250778


################### external validation 2var: MAIN any death - Kenya only ####
GEMS_glm_death <- glm(death_all~f4b_muac+f4b_resp,
                      data=complete_anydeath_Kenya,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  8.515605   1.609084   5.292 1.21e-07 ***
#   f4b_muac    -0.917609   0.114183  -8.036 9.26e-16 ***
#   f4b_resp     0.008297   0.014459   0.574    0.566  
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathAll_2var_Kenya.Rdata")
round(exp(coef(glm(death_all~f4b_muac+f4b_resp,
                   data=complete_anydeath_Kenya,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac    f4b_resp 
# 4992.0634      0.3995      1.0083  

GEMS_2var_Kenya<-complete_anydeath_Kenya %>% select(death_all,f4b_muac,f4b_resp)
GEMS_2var_Kenya$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_2var_Kenya,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_2var_Kenya$pred_glm,labels=GEMS_2var_Kenya$death_all)
GEMS_AUC
# [1] 0.8585684
#using pROC package since need CI
GEMS_AUC <- roc(response=GEMS_2var_Kenya$death_all,predictor=GEMS_2var_Kenya$pred_glm)
paste(round(GEMS_AUC$auc,2)," (",
      round(ci.auc(GEMS_AUC)[1],2),", ",
      round(ci.auc(GEMS_AUC)[3],2),")",sep="")
# "0.86 (0.8, 0.92)"


GEMS_decilesCC <- GEMS_2var_Kenya %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_all),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathAll_2var_Kenya.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_all)
kilifi_AUC
# [1] 0.7249678
kilifi_AUC <- roc(response=kilifi$death_all,predictor=kilifi$GEMS_pred_glm)
paste(round(kilifi_AUC$auc,2)," (",
      round(ci.auc(kilifi_AUC)[1],2),", ",
      round(ci.auc(kilifi_AUC)[3],2),")",sep="")
# "0.72 (0.69, 0.76)"

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_all),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathAll_2var_Kenya.Rdata")

#jpeg("/DeathAll_CC_ExternalVal_GEMS059_2var_Kenya.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: All Deaths in Kenya")),
     xlim=c(0,0.4),ylim=c(0,0.4))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_all)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_all)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.86 (0.80, 0.92)","Kilifi data (0-59mo), AUC=0.72 (0.69, 0.76)"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_all~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_all~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)   
# (Intercept) -0.24638    0.07642  -3.224  0.00126 **
# Waiting for profiling to be done...
# 2.5 %      97.5 % 
#   -0.39829671 -0.09865555 

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -1.08926    0.10362  -10.51   <2e-16 ***
#   log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.50069    0.03957   12.65   <2e-16 ***
  
# Waiting for profiling to be done...
# 2.5 %     97.5 %
#   (Intercept)                            -1.294439 -0.8878263
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.423787  0.5790449


################### external validation 4var: DEATH in hospital ####
GEMS_glm_death <- glm(death_hosp~f4b_muac+f4b_resp+f4b_temp+base_age,
                      data=complete_hospdeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -18.61016    5.31338  -3.503 0.000461 ***
#   f4b_muac     -0.86213    0.09791  -8.805  < 2e-16 ***
#   f4b_resp      0.02879    0.01245   2.313 0.020718 *  
#   f4b_temp      0.59891    0.14650   4.088 4.35e-05 ***
#   base_age      0.04694    0.01475   3.182 0.001465 **
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathHosp_4var.Rdata")
round(exp(coef(glm(death_hosp~f4b_muac+f4b_resp+f4b_temp+base_age,
                   data=complete_hospdeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac    f4b_resp    f4b_temp    base_age 
# 0.0000      0.4223      1.0292      1.8201      1.0481

GEMS_4var<-complete_anydeath %>% select(death_hosp,f4b_muac,f4b_resp,f4b_temp,base_age)
GEMS_4var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_4var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_4var$pred_glm,labels=GEMS_4var$death_hosp)
GEMS_AUC
# [1] 0.8869814

GEMS_decilesCC <- GEMS_4var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_hosp),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathHosp_4var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_hosp)
kilifi_AUC
# [1] 0.7236402

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_hosp),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathHosp_4var.Rdata")

#jpeg("/DeathHosp_CC_ExternalVal_GEMS059_4var.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Deaths in Hospital")),
     xlim=c(0,0.2),ylim=c(0,0.2))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_hosp)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_hosp)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.87","Kilifi data (0-59mo), AUC=0.72"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_hosp~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_hosp~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.0380     0.1018    10.2   <2e-16 ***
#   2.5 %    97.5 % 
#   0.8338884 1.2332007 
# 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -0.80730    0.23518  -3.433 0.000598 ***
#   log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.48898    0.05311   9.207  < 2e-16 ***
#   2.5 %     97.5 %
#   (Intercept)                            -1.2772765 -0.3535589
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.3852787  0.5937959


################### external validation 2var: DEATH in hospital ####
GEMS_glm_death <- glm(death_hosp~f4b_muac+f4b_resp,
                      data=complete_hospdeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  2.33125    1.26216   1.847   0.0647 .  
# f4b_muac    -0.69050    0.08639  -7.993 1.32e-15 ***
#   f4b_resp     0.03212    0.01167   2.753   0.0059 ** 
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathHosp_2var.Rdata")
round(exp(coef(glm(death_hosp~f4b_muac+f4b_resp,
                   data=complete_hospdeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac    f4b_resp 
# 10.2908      0.5013      1.0326

GEMS_2var<-complete_anydeath %>% select(death_hosp,f4b_muac,f4b_resp)
GEMS_2var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_2var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_2var$pred_glm,labels=GEMS_2var$death_hosp)
GEMS_AUC
# [1] 0.8660985

GEMS_decilesCC <- GEMS_2var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_hosp),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathHosp_2var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_hosp)
kilifi_AUC
# [1] 0.7492435

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_hosp),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathHosp_2var.Rdata")

#jpeg("/DeathHosp_CC_ExternalVal_GEMS059_2var.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Deaths in Hospital")),
     xlim=c(0,0.2),ylim=c(0,0.2))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_hosp)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_hosp)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.87","Kilifi data (0-59mo), AUC=0.75"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_hosp~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_hosp~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.35374    0.09928   13.63   <2e-16 ***
#   2.5 %   97.5 % 
#   1.154535 1.543982 
# 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -0.06859    0.27467   -0.25    0.803    
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.63891    0.06181   10.34   <2e-16 ***
#   2.5 %    97.5 %
#   (Intercept)                            -0.6147587 0.4641589
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.5183292 0.7610292



################### external validation 1var: DEATH in hospital ####
GEMS_glm_death <- glm(death_hosp~f4b_muac,
                      data=complete_hospdeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.35087    0.98844   4.402 1.07e-05 ***
#   f4b_muac    -0.74754    0.08303  -9.003  < 2e-16 ***
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathHosp_1var.Rdata")
round(exp(coef(glm(death_hosp~f4b_muac,
                   data=complete_hospdeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac 
# 77.5456      0.4735 

GEMS_1var<-complete_anydeath %>% select(death_hosp,f4b_muac)
GEMS_1var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_1var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_1var$pred_glm,labels=GEMS_1var$death_hosp)
GEMS_AUC
# [1] 0.8403117

GEMS_decilesCC <- GEMS_1var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_hosp),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathHosp_1var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_hosp)
kilifi_AUC
# [1] 0.7106909

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_hosp),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathHosp_1var.Rdata")

#jpeg("/DeathHosp_CC_ExternalVal_GEMS059_1var.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Deaths in Hospital")),
     xlim=c(0,0.2),ylim=c(0,0.2))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_hosp)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_hosp)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.84","Kilifi data (0-59mo), AUC=0.71"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_hosp~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_hosp~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.42333    0.09868   14.42   <2e-16 ***
#   2.5 %   97.5 % 
#   1.225313 1.612396 
# 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -0.31007    0.28955  -1.071    0.284    
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.57374    0.06294   9.116   <2e-16 ***
#   2.5 %    97.5 %
#   (Intercept)                            -0.8881754 0.2490529
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.4504611 0.6975548



################### external validation 4var: DEATH at home ####
GEMS_glm_death <- glm(death_home~f4b_muac+f4b_resp+f4b_temp+base_age,
                      data=complete_homedeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -10.440602   3.418743  -3.054  0.00226 ** 
#   f4b_muac     -0.790555   0.061745 -12.804  < 2e-16 ***
#   f4b_resp      0.017926   0.008823   2.032  0.04217 *  
#   f4b_temp      0.415972   0.093757   4.437 9.14e-06 ***
#   base_age      0.010675   0.010809   0.988  0.32335  
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathHome_4var.Rdata")
round(exp(coef(glm(death_home~f4b_muac+f4b_resp+f4b_temp+base_age,
                   data=complete_homedeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac    f4b_resp    f4b_temp    base_age 
# 0.0000      0.4536      1.0181      1.5158      1.0107 

GEMS_4var<-complete_anydeath %>% select(death_home,f4b_muac,f4b_resp,f4b_temp,base_age)
GEMS_4var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_4var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_4var$pred_glm,labels=GEMS_4var$death_home)
GEMS_AUC
# [1] 0.8440134

GEMS_decilesCC <- GEMS_4var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_home),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathHome_4var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_all)
kilifi_AUC
# [1] 0.7332254

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_home),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathHome_4var.Rdata")

#jpeg("/DeathHome_CC_ExternalVal_GEMS059_4var.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Deaths at Home")),
     xlim=c(0,0.4),ylim=c(0,0.4))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_home)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_home)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.84","Kilifi data (0-59mo), AUC=0.73"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_home~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_home~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.08997    0.09844   0.914    0.361
# 2.5 %     97.5 % 
#   -0.1073781  0.2787509 
# 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -1.20129    0.19700  -6.098 1.08e-09 ***
#   log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.51830    0.05732   9.043  < 2e-16 ***
#   2.5 %     97.5 %
#   (Intercept)                            -1.596569 -0.8228016
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.406129  0.6311350


################### external validation 2var: DEATH at home ####
GEMS_glm_death <- glm(death_home~f4b_muac+f4b_resp,
                      data=complete_homedeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.325252   0.816932   5.295 1.19e-07 ***
#   f4b_muac    -0.736566   0.055535 -13.263  < 2e-16 ***
#   f4b_resp     0.024780   0.008285   2.991  0.00278 ** 
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathHome_2var.Rdata")
round(exp(coef(glm(death_home~f4b_muac+f4b_resp,
                   data=complete_homedeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac    f4b_resp 
# 75.5846      0.4788      1.0251

GEMS_2var<-complete_anydeath %>% select(death_home,f4b_muac,f4b_resp)
GEMS_2var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_2var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_2var$pred_glm,labels=GEMS_2var$death_home)
GEMS_AUC
# [1] 0.8402765

GEMS_decilesCC <- GEMS_2var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_home),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathHome_2var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_all)
kilifi_AUC
# [1] 0.7381043

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_home),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathHome_2var.Rdata")

#jpeg("/DeathHome_CC_ExternalVal_GEMS059_2var.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Deaths at Home")),
     xlim=c(0,0.4),ylim=c(0,0.4))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_home)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_home)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.84","Kilifi data (0-59mo), AUC=0.74"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_home~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_home~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  0.26169    0.09817   2.666  0.00768 **
#   2.5 %     97.5 % 
#   0.06486085 0.44993471 
# 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                             -1.1509     0.2111  -5.452 4.97e-08 ***
#   log(GEMS_pred_glm/(1 - GEMS_pred_glm))   0.5099     0.0587   8.686  < 2e-16 ***
#   2.5 %     97.5 %
#   (Intercept)                            -1.5748309 -0.7457446
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.3947467  0.6252127



################### external validation 1var: DEATH at home ####
GEMS_glm_death <- glm(death_home~f4b_muac,
                      data=complete_homedeath,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_death)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   5.7778     0.6533   8.845   <2e-16 ***
#   f4b_muac     -0.7733     0.0541 -14.293   <2e-16 ***
round(coef(GEMS_glm_death),4)
#save(GEMS_glm_death, file = "/GEMS_glm_deathHome_1var.Rdata")
round(exp(coef(glm(death_home~f4b_muac,
                   data=complete_homedeath,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)    f4b_muac 
# 323.0630      0.4615 

GEMS_1var<-complete_anydeath %>% select(death_home,f4b_muac)
GEMS_1var$pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=GEMS_1var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_1var$pred_glm,labels=GEMS_1var$death_home)
GEMS_AUC
# [1] 0.8238288

GEMS_decilesCC <- GEMS_1var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_home),mean(pred_glm))
#save(GEMS_decilesCC, file = "/GEMS_decilesCC_deathHome_1var.Rdata")

kilifi$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_death,newdata=kilifi,type="response"))
kilifi_AUC <- AUC(predictions=kilifi$GEMS_pred_glm,labels=kilifi$death_all)
kilifi_AUC
# [1] 0.7179465

kilifi_decilesCC <- kilifi %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(death_home),mean(GEMS_pred_glm))
#save(kilifi_decilesCC, file = "/kilifi_decilesCC_deathHome_1var.Rdata")

#jpeg("/DeathHome_CC_ExternalVal_GEMS059_1var.jpg",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Deaths at Home")),
     xlim=c(0,0.4),ylim=c(0,0.4))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(death_home)`,col="red",pch=1)
points(kilifi_decilesCC$`mean(GEMS_pred_glm)`,kilifi_decilesCC$`mean(death_home)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.82","Kilifi data (0-59mo), AUC=0.72"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(death_home~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(death_home~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.32582    0.09768   3.336 0.000851 ***
#   2.5 %    97.5 % 
#   0.1299591 0.5131159 
# 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            -1.15113    0.21883   -5.26 1.44e-07 ***
#   log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.50137    0.05983    8.38  < 2e-16 ***
#   2.5 %     97.5 %
#   (Intercept)                            -1.5907744 -0.7313277
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.3839185  0.6187939


