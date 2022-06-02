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
GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_1BLOCK_NEWCALIB.DAT")

##################################### MODEL(s) #####################################

#------------------------------------ Operating Model ------------------------------------

# Make a folder to put the results:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
TestName <- "OperatingModel"
if(!dir.exists(TestName)) dir.create(TestName)
setwd(file.path(getwd(),TestName))

input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                            selectivity=list(model=rep("age-specific",3),
                                             re=c(rep('ar1_y',1),rep("none",2)),
                                             initial_pars <- list(c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                                                                  c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                  c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                             fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                            NAA_re=list(sigma='rec+1',cor='ar1_y'),
                            age_comp='multinomial')
                            # age_comp="dirichlet")
  
input$par$log_F1 <- log(2)
input$map$log_F1 <- as.factor(NA)
input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))
  
# Mapping the estimation of numbers at age in the first year:
# input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
# input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8
  
# Fit model:
mod <- try(fit_wham(input, do.osa = F, do.check=F,do.retro=T)) # turn off OSA residuals to save time

# plot_wham_output(mod=mod,out.type="html")
