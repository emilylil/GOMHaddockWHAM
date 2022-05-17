# Gulf of Maine Haddock in WHAM
# Emily Liljestrand
# Created: May 2, 2022
# Updates: May 4, 2022

# Remove Old Objects and set working directory to file location:
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load Packages:
library(wham)
library(TMB)
library(dplyr)

##################################### DATA #####################################

#Read in Gulf of Maine Haddock Data from ASAP Dat File
GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB.DAT")
#Read in Gulf of Maine Haddock Data from ASAP Dat File, version for selectivity random effects
# GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB_EVERYYEARSEL.DAT")

##################################### MODEL(s) #####################################

#------------------------------------ Test 3: Compare models with different Rec and NAA structure ------------------------------------
# Make a folder to put the results of Test 3:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
TestName <- "Test3_Comparison_RecNAA"
if(!dir.exists(TestName)) dir.create(TestName)
setwd(file.path(getwd(),TestName))

# Read in ASAP model fitted values
ASAP <- read_asap3_fit(wd="C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM/GOMHaddockASAP",asap.name="GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB")

agecomp <- c('dirichlet','dir-mult','logistic-normal-pool0','logistic-normal-miss0')
recruitment<-c(1,2)
NAA_re_sigma<-c('rec','rec+1')
NAA_re_cor<-c('iid','ar1_y','ar1_a','2dar1')
# selectivity_re<-c('none','iid','ar1','ar1_y','2dar1')


# Set up a data frame for information of model options:
df.mods <- as.data.frame(tidyr::crossing(agecomp,recruitment,NAA_re_sigma,NAA_re_cor))

n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

for(m in 1:n.mods){
  # input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment")
  input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=df.mods[m,"recruitment"], model_name="GOMHaddock_RW_Recruitment",
                              selectivity=list(model=rep("age-specific",5),
                                               re=rep("none",5),
                                               initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                                 c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                                 c(0.01,0.1,0.3,0.5,0.8,0.9,1,1,0.999),
                                                                 c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                 c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                               fix_pars=list(8:9,8:9,7:8,4:9,6:9)),
                              NAA_re=list(sigma=df.mods[m,"NAA_re_sigma"],cor=df.mods[m,"NAA_re_cor"]),
                              # age_comp=df.mods[m,"agecomp"])
                              age_comp=df.mods[m,"agecomp"])
  
  if(df.mods[m,"agecomp"]!='dirichlet')
  {
    input$data$catch_Neff <- input$data$catch_Neff*100
    input$data$index_Neff <- input$data$index_Neff*100
  }
  
  # Mapping the estimation of numbers at age in the first year:
  # input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
  input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
  # input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8
  
  # Fit model:
  mod <- fit_wham(input, do.osa = F, do.check=F) # turn off OSA residuals to save time
  # m1 <- fit_wham(input, do.osa = F) # turn off OSA residuals to save time
  
  saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))
}

# List Models:
mod.list <- paste0(df.mods$Model,".rds")
mods <- lapply(mod.list, readRDS)

opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)

df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
not_conv <- !df.mods$conv | !df.mods$pdHess
mods2 <- mods
mods2[not_conv] <- NULL
df.aic.tmp <- as.data.frame(compare_wham_models(mods2, table.opts=list(sort=FALSE, calc.rho=TRUE))$tab)
df.aic <- df.aic.tmp[FALSE,]
ct = 1
for(i in 1:n.mods){
  if(not_conv[i]){
    df.aic[i,] <- rep(NA,5)
  } else {
    df.aic[i,] <- df.aic.tmp[ct,]
    ct <- ct + 1
  }
}
df.aic[,1:2] <- format(round(df.aic[,1:2], 1), nsmall=1)
df.aic[,3:5] <- format(round(df.aic[,3:5], 3), nsmall=3)
df.aic[grep("NA",df.aic$dAIC),] <- "---"
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL

# Compare across models
write.csv(df.mods, file="df.mods.csv",quote=F, row.names=F)

# plot_wham_output(mod=mods[[27]],out.type="html")

models <- list(ASAP=ASAP,iid=mods[[24]],AR1A=mods[[22]],AR1Y=mods[[23]])

# Compare across models
# compare_wham_models(models,fdir=getwd())

#------------------------------------ Test 4: Compare models with different Rec, NAA, sel structure ------------------------------------

# White noise: recruit_model-1, cor-iid
# Random walk: recruit_model-2, cor-iid
# AR1: recruit_model-1, cor-ar1_y
# 
# if(df.mods[m,'recruitment']=='AR(1)'){corid="ar1_y"}
# else{corid='iid'}

GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_1BLOCK_NEWCALIB.DAT")

# Make a folder to put the results of Test 4:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
TestName <- "Test5_Comparison_RecNAASel"
if(!dir.exists(TestName)) dir.create(TestName)
setwd(file.path(getwd(),TestName))

# Read in ASAP model fitted values
ASAP <- read_asap3_fit(wd="C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM/GOMHaddockASAP",asap.name="GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB")

agecomp <- c('dir-mult','logistic-normal-pool0','logistic-normal-miss0')
recruitment<-1
NAA_re_sigma<-c('rec','rec+1')
NAA_re_cor<-c('ar1_y','2dar1')
selectivity_re<-c('iid','ar1_y','2dar1')


# Set up a data frame for information of model options:
df.mods <- as.data.frame(tidyr::crossing(agecomp,recruitment,NAA_re_sigma,NAA_re_cor,selectivity_re))

n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

for(m in 1:n.mods){
  # initial_pars<- list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999))
  # input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment")
  input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=df.mods[m,"recruitment"], model_name="GOMHaddock_RW_Recruitment",
                              selectivity=list(model=rep("age-specific",3),
                                               re=c(rep(df.mods[m,"selectivity_re"],1),rep("none",2)),
                                               initial_pars <- list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                                    c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                    c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                               fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                              NAA_re=list(sigma=df.mods[m,"NAA_re_sigma"],cor=df.mods[m,"NAA_re_cor"]),
                              age_comp=df.mods[m,"agecomp"])
                              # age_comp="dirichlet")
  
  input$par$log_F1 <- log(5)
  input$map$log_F1 <- as.factor(NA)
  input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))
  
  # Mapping the estimation of numbers at age in the first year:
  # input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
  input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
  # input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8
  
  if(df.mods[m,"agecomp"]!='dirichlet')
  {
    # input$data$catch_Neff <- input$data$catch_Neff*100
    # input$data$index_Neff <- input$data$index_Neff*100
    
    input$data$catch_Neff <- input$data$catch_Neff*0+150.0
    input$data$index_Neff <- input$data$index_Neff*0+150.0
  }
  
  # Fit model:
  mod <- try(fit_wham(input, do.osa = F, do.check=F)) # turn off OSA residuals to save time
  # m1 <- fit_wham(input, do.osa = F) # turn off OSA residuals to save time
  
  saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))
}


##################################### Epilogue/Footer, junk code that might be useful #####################################
# How to parameterize the 3 block model into two blocks:
input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                            selectivity=list(model=rep("age-specific",5),
                                             re=rep("none",5),
                                             initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                               c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                               c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                               c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                               c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                             fix_pars=list(8:9,8:9,1:9,4:9,6:9)))


# Change mapping so the third catch selectivity block necessarily matches the second
input$map$logit_selpars[1:45] <- as.factor(c(1,2,2,3,4,5,6,6,7,8,9,10,10,11,12,13,14,14,NA,15,16,17,17,NA,18,19,20,20,NA,NA,21,22,22,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))

#Basic fit model:
# Fit model:
mod <- fit_wham(input, do.osa = F, do.check=T) # turn off OSA residuals to save time
# m1 <- fit_wham(input, do.osa = F) # turn off OSA residuals to save time

check_convergence(mod)

plot_wham_output(mod=mod,out.type="html")

saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))

models <- list(ASAP=ASAP,iid=mod,iid=mod)

# Compare across models
# compare_wham_models(models,fdir=getwd())