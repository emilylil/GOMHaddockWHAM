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
# GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB.DAT")
#Read in Gulf of Maine Haddock Data from ASAP Dat File, version for selectivity random effects (one block for all catch)
GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_1BLOCK_NEWCALIB.DAT")

##################################### MODEL(s) #####################################

# Make a folder to put the results:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
TestName <- "OperatingModel"
if(!dir.exists(TestName)) dir.create(TestName)
setwd(file.path(getwd(),TestName))

input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                            selectivity=list(model=rep("age-specific",3),
                                             re=c(rep('ar1_y',1),rep("none",2)),
                                             initial_pars <- list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999)*0.2,
                                                                  c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                  c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                             fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                            NAA_re=list(sigma='rec+1',cor='ar1_y'),
                            age_comp='dir-mult')

input$par$log_F1 <- log(2)
input$map$log_F1 <- as.factor(NA)
input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))

# Mapping the estimation of numbers at age in the first year:
# input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
# input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8

input$data$catch_Neff <- input$data$catch_Neff*100
input$data$index_Neff <- input$data$index_Neff*100

# Fit model:
mod <- try(fit_wham(input, do.osa = F, do.check=F,do.retro=T)) # turn off OSA residuals to save time
saveRDS(mod, file=paste0(TestName,".rds"))

plot_wham_output(mod=mod,out.type="html")


##################################### Epilogue/Footer, junk code that might be useful #####################################
#
#------------------------------------ Simulation tests without selectivity random effects ------------------------------------
# # Make a folder to put the results:
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# TestName <- "Test_3Block_Sel"
# if(!dir.exists(TestName)) dir.create(TestName)
# setwd(file.path(getwd(),TestName))
# 
# # Read in ASAP model fitted values
# ASAP <- read_asap3_fit(wd="C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM/GOMHaddockASAP",asap.name="GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB")
# 
# # Some modeling options
# agecomp <- c('dirichlet','dir-mult','logistic-normal-pool0','logistic-normal-miss0') #Likelihood of age compositions
# recruitment<-c(1,2) # Recruitment is 1- random walk or 2- random about the mean
# NAA_re_sigma<-c('rec','rec+1') # Process error on Recruitment or Recruitment AND survival
# NAA_re_cor<-c('iid','ar1_y','ar1_a','2dar1') # Correlation structure on process errors of Recruitment (+survival)
# 
# # Set up a data frame that contains all combinations of the modeling options:
# df.mods <- as.data.frame(tidyr::crossing(agecomp,recruitment,NAA_re_sigma,NAA_re_cor))
# 
# # Count number of models and add column with model labels
# n.mods <- dim(df.mods)[1]
# df.mods$Model <- paste0("m",1:n.mods)
# df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
# 
# for(m in 1:n.mods){
#   # input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment")
#   input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=df.mods[m,"recruitment"], model_name="GOMHaddock_RW_Recruitment",
#                               selectivity=list(model=rep("age-specific",5),
#                                                re=rep("none",5),
#                                                initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
#                                                                  c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
#                                                                  c(0.01,0.1,0.3,0.5,0.8,0.9,1,1,0.999),
#                                                                  c(0.2,0.4,0.8,1,1,1,1,1,1),
#                                                                  c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
#                                                fix_pars=list(8:9,8:9,7:8,4:9,6:9)),
#                               NAA_re=list(sigma=df.mods[m,"NAA_re_sigma"],cor=df.mods[m,"NAA_re_cor"]),
#                               age_comp=df.mods[m,"agecomp"])
#   
#   # When the composition model involves weighting the effective sample size, make the input sample size large
#   if(df.mods[m,"agecomp"]!='dirichlet')
#   {
#     input$data$catch_Neff <- input$data$catch_Neff*10
#     input$data$index_Neff <- input$data$index_Neff*10
#   }
#   
#   # Mapping the estimation of numbers at age in the first year:
#   # input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
#   input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
#   # input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8
#   
#   # Fit model:
#   mod <- fit_wham(input, do.osa = F, do.check=F,do.retro=F) # turn off OSA residuals to save time
#   # m1 <- fit_wham(input, do.osa = F) # Default
#   
#   # Save modeling object to folder:
#   saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))
# }
# 
# # List Models:
# mod.list <- paste0(df.mods$Model,".rds")
# mods <- lapply(mod.list, readRDS)
# 
# # Add columns for convergence information and if hessian was produced:
# opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
# ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
# df.mods$conv <- as.logical(opt_conv)
# df.mods$pdHess <- as.logical(ok_sdrep)
# 
# # Additional columns added to data frame
# df.mods$runtime <- sapply(mods, function(x) x$runtime)
# df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
# not_conv <- !df.mods$conv | !df.mods$pdHess
# mods2 <- mods
# mods2[not_conv] <- NULL
# df.aic.tmp <- as.data.frame(compare_wham_models(mods2, table.opts=list(sort=FALSE, calc.rho=TRUE))$tab)
# df.aic <- df.aic.tmp[FALSE,]
# ct = 1
# for(i in 1:n.mods){
#   if(not_conv[i]){
#     df.aic[i,] <- rep(NA,5)
#   } else {
#     df.aic[i,] <- df.aic.tmp[ct,]
#     ct <- ct + 1
#   }
# }
# df.aic[,1:2] <- format(round(df.aic[,1:2], 1), nsmall=1)
# df.aic[,3:5] <- format(round(df.aic[,3:5], 3), nsmall=3)
# df.aic[grep("NA",df.aic$dAIC),] <- "---"
# df.mods <- cbind(df.mods, df.aic)
# rownames(df.mods) <- NULL
# 
# # Compare across models
# write.csv(df.mods, file="df.mods.csv",quote=F, row.names=F)
# 
# # Plot a particular model
# # plot_wham_output(mod=mods[[1]],out.type="html")
# 
# # Specify models to compare across. If you want to include asap model, must have more than 1 WHAM model in list
# # models <- list(ASAP=ASAP,model1=mods[[1]],model2=mods[[2]],model3=mods[[3]])
# 
# # Compare across models
# # compare_wham_models(models,fdir=getwd())
# 
# #------------------------------------ Simulation tests WITH selectivity random effects ------------------------------------
# 
# GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_1BLOCK_NEWCALIB.DAT")
# 
# # Make a folder to put the results of Test 4:
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# # TestName <- "Test6_Comparison_RecNAASel_NoScale"
# TestName <- "Test7_finalmodel"
# if(!dir.exists(TestName)) dir.create(TestName)
# setwd(file.path(getwd(),TestName))
# 
# # Read in ASAP model fitted values
# ASAP <- read_asap3_fit(wd="C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM/GOMHaddockASAP",asap.name="GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB")
# 
# agecomp <- c('dir-mult')
# recruitment<-c(2)
# NAA_re_sigma<-c('rec+1')
# NAA_re_cor<-c('ar1_y')
# selectivity_re<-c('ar1_y')
# ESSoptions<-c(500)
# 
# # Set up a data frame for information of model options:
# df.mods <- as.data.frame(tidyr::crossing(agecomp,recruitment,NAA_re_sigma,NAA_re_cor,selectivity_re,ESSoptions))
# 
# n.mods <- dim(df.mods)[1]
# df.mods$Model <- paste0("m",1:n.mods)
# df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
# 
# for(m in 1:n.mods){
#   # initial_pars<- list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999))
#   # input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment")
#   input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=df.mods[m,"recruitment"], model_name="GOMHaddock_RW_Recruitment",
#                               selectivity=list(model=rep("age-specific",3),
#                                                re=c(rep(df.mods[m,"selectivity_re"],1),rep("none",2)),
#                                                initial_pars <- list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
#                                                                     c(0.2,0.4,0.8,1,1,1,1,1,1),
#                                                                     c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
#                                                fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
#                               NAA_re=list(sigma=df.mods[m,"NAA_re_sigma"],cor=df.mods[m,"NAA_re_cor"]),
#                               age_comp=df.mods[m,"agecomp"])
#                               # age_comp="dirichlet")
#   
#   input$par$log_F1 <- log(2)
#   input$map$log_F1 <- as.factor(NA)
#   input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))
#   
#   # Mapping the estimation of numbers at age in the first year:
#   # input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
#   input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
#   # input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8
#   
#   # if(df.mods[m,"ESSoptions"]<200)
#   # {
#     input$data$catch_Neff <- input$data$catch_Neff*df.mods[m,"ESSoptions"]
#     input$data$index_Neff <- input$data$index_Neff*df.mods[m,"ESSoptions"]
#   # }
#   # else
#   # {
#   #   input$data$catch_Neff <- input$data$catch_Neff*0.0+500.0
#   #   input$data$index_Neff <- input$data$index_Neff*0.0+500.0
#   # }
#   
#   # Fit model:
#   mod <- try(fit_wham(input, do.osa = F, do.check=F,do.retro=T,n.peels=1)) # turn off OSA residuals to save time
#   # m1 <- fit_wham(input, do.osa = F) # turn off OSA residuals to save time
#   
#   saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))
# }
# 
# # List Models:
# mod.list <- paste0(df.mods$Model,".rds")
# mods <- lapply(mod.list, readRDS)
# 
# opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
# ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
# df.mods$conv <- as.logical(opt_conv)
# df.mods$pdHess <- as.logical(ok_sdrep)
# 
# df.mods$runtime <- sapply(mods, function(x) x$runtime)
# df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
# not_conv <- !df.mods$conv | !df.mods$pdHess
# mods2 <- mods
# mods2[not_conv] <- NULL
# df.aic.tmp <- as.data.frame(compare_wham_models(mods2, table.opts=list(sort=FALSE, calc.rho=TRUE))$tab)
# df.aic <- df.aic.tmp[FALSE,]
# ct = 1
# for(i in 1:n.mods){
#   if(not_conv[i]){
#     df.aic[i,] <- rep(NA,5)
#   } else {
#     df.aic[i,] <- df.aic.tmp[ct,]
#     ct <- ct + 1
#   }
# }
# df.aic[,1:2] <- format(round(df.aic[,1:2], 1), nsmall=1)
# df.aic[,3:5] <- format(round(df.aic[,3:5], 3), nsmall=3)
# df.aic[grep("NA",df.aic$dAIC),] <- "---"
# df.mods <- cbind(df.mods, df.aic)
# rownames(df.mods) <- NULL
# 
# # Compare across models
# write.csv(df.mods, file="df.mods.csv",quote=F, row.names=F)
# 
# fittedmodel <- mods[[1]]
# betaC <- exp(fittedmodel$parList$catch_paa_pars)
# betaI1 <- exp(fittedmodel$parList$index_paa_pars[1])
# betaI2 <- exp(fittedmodel$parList$index_paa_pars[2])
# 
# ISSC <- fittedmodel$input$data$catch_Neff
# ISSI1 <- fittedmodel$input$data$index_Neff[,1]
# ISSI2 <- fittedmodel$input$data$index_Neff[,2]
# 
# ESSC <- (ISSC+ISSC*betaC)/(ISSC+betaC)
# ESSI1 <- (ISSI1+ISSI1*betaI1)/(ISSI1+betaI1)
# ESSI2 <- (ISSI2+ISSI2*betaI2)/(ISSI2+betaI2)
# 
# ESSC
# 
# # plot_wham_output(mod=mods[[1]],out.type="html")
# 
# # models <- list(ASAP=ASAP,model1=mods[[1]],model2=mods[[2]],model3=mods[[3]])
# models <- list(first=mods[[1]],fourth=mods[[4]])
# 
# # Compare across models
# # compare_wham_models(models,fdir=getwd())
# 
# #------------------------------------ Parameterizing three sel blocks as two blocks ------------------------------------
# 
# # How to parameterize the 3 block model into two blocks:
# input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
#                             selectivity=list(model=rep("age-specific",5),
#                                              re=rep("none",5),
#                                              initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
#                                                                c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
#                                                                c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
#                                                                c(0.2,0.4,0.8,1,1,1,1,1,1),
#                                                                c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
#                                              fix_pars=list(8:9,8:9,1:9,4:9,6:9)))
# 
# 
# # Change mapping so the third catch selectivity block necessarily matches the second
# input$map$logit_selpars[1:45] <- as.factor(c(1,2,2,3,4,5,6,6,7,8,9,10,10,11,12,13,14,14,NA,15,16,17,17,NA,18,19,20,20,NA,NA,21,22,22,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
# 
# #Basic fit model:
# # Fit model:
# mod <- fit_wham(input, do.osa = F, do.check=T,do.retro=F) # turn off OSA residuals to save time
# # m1 <- fit_wham(input, do.osa = F) # turn off OSA residuals to save time
# 
# check_convergence(mod)
# 
# plot_wham_output(mod=mods[[3]],out.type="html")
# 
# saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))
# 
# models <- list(Best=mods[[11]],Second=mods[[13]],Third=mods[[15]])
# 
# # Compare across models
# # compare_wham_models(models,fdir=getwd())