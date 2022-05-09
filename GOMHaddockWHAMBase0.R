# Gulf of Maine Haddock in WHAM
# Just a copy of the Base0 model, for posterity
# Emily Liljestrand
# Created: May 2, 2022
# Updated: May 9, 2022

# Remove Old Objects and set working directory to file location:
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load Packages:
library(wham)
library(TMB)

##################################### DATA #####################################

#Read in Gulf of Maine Haddock Data from ASAP Dat File
GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB.DAT")

##################################### MODEL(s) #####################################
# Read in ASAP model fitted values
#This is the version that sets lamba=0 for recruitment deviations
ASAP <- read_asap3_fit(wd="C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM/GOMHaddockASAP",asap.name="GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB_RECDEV0")

# Base WHAM Model
# input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment")
input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                            selectivity=list(model=rep("age-specific",5),
                                             re=rep("none",5),
                                             initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                               c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                               c(0.01,0.1,0.3,0.5,0.8,0.9,1,1,0.999),
                                                               c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                               c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                             fix_pars=list(8:9,8:9,7:8,4:9,6:9)))

# Fix age 7-8 in the first year to 500 individuals
input$par$log_N1_pars[7:8] <- 0.5

# Mapping the estimation of numbers at age in the first year:
# input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
# input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8

# Fit model:
m0 <- fit_wham(input, do.osa = F, do.check=T) # turn off OSA residuals to save time
# m0 <- fit_wham(input, do.osa = F) # turn off OSA residuals to save time

# Check that m0 converged (m1$opt$convergence should be 0, and the maximum gradient should be < 1e-06)
# check_convergence(m0)
m0$rep
m0$rep$SSB

# Plot results
# plot_wham_output(mod=m0,out.type="html")

# Make list of models
models <- list(ASAP=ASAP,WHAMBase0=m0,WHAMBase0=m0)

# Compare across models
compare_wham_models(models,fdir=getwd())
