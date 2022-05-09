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

##################################### DATA #####################################

#Read in Gulf of Maine Haddock Data from ASAP Dat File
GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB.DAT")
GOM_HADDOCK_DAT_2BLOCK <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_2BLOCK_NEWCALIB.DAT")

##################################### MODEL(s) #####################################

#------------------------------------ Test 0: Compare WHAM and ASAP ------------------------------------
# Make a folder to put the results of Test 0:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
TestName <- "Test0_Comparison"
if(!dir.exists(TestName)) dir.create(TestName)
setwd(file.path(getwd(),TestName))

# Read in ASAP model fitted values
ASAP <- read_asap3_fit(wd="C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM/GOMHaddockASAP",asap.name="GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB")

# Model 1: Base WHAM Model, 3 selectivity blocks, all asymptote
input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment")
input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                            selectivity=list(model=rep("age-specific",5),
                                             re=rep("none",5),
                                             initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                               c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                               c(0.01,0.1,0.3,0.5,0.8,0.9,1,1,0.8),
                                                               c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                               c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                             fix_pars=list(8:9,8:9,7:8,4:9,6:9)))

# Mapping the estimation of numbers at age in the first year:
# input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
# input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8

# Fit model:
m1 <- fit_wham(input, do.osa = F, do.check=T) # turn off OSA residuals to save time
# m1 <- fit_wham(input, do.osa = F) # turn off OSA residuals to save time

# Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradient should be < 1e-06)
# check_convergence(m1)
m1$rep
m1$rep$SSB

# Plot results
# plot_wham_output(mod=m1,out.type="html")

# Model 2: Modified WHAM Model, 2 selectivity blocks, all asymptote

input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                            selectivity=list(model=rep("age-specific",5),
                                             re=rep("none",5),
                                             initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                               c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                               c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                               c(0.2,0.4,0.8,0.999,0.999,0.999,0.999,0.999,0.999),
                                                               c(0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999,0.999)),
                                             fix_pars=list(8:9,8:9,1:9,4:9,6:9)))

# Mapping the estimation of numbers at age in the first year:
# input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
# input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8

# Change mapping so the third catch selectivity block necessarily matches the second
input$map$logit_selpars[1:45] <- as.factor(c(1,2,2,3,4,5,6,6,7,8,9,10,10,11,12,13,14,14,NA,15,16,17,17,NA,18,19,20,20,NA,NA,21,22,22,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))

# Fit model:
# m2 <- fit_wham(input, do.osa = F, do.check=T) # turn off OSA residuals to save time
m2 <- fit_wham(input, do.osa = F) # turn off OSA residuals to save time

# Check that model converged
# check_convergence(m2)

# Plot results
# plot_wham_output(mod=m2,out.type="html")

# Make list of models
models <- list(ASAP=ASAP,Sel3Block=m1,Sel2Block=m2)
models <- list(ASAP=ASAP,WHAM=m1,WHAM=m1)

# Compare across models
compare_wham_models(models,fdir=getwd())

#------------------------------------ Test 1: Sel at age 9 ------------------------------------

# create directory for analysis, based on TestName,
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
TestName <- "Test1_Sel9"
if(!dir.exists(TestName)) dir.create(TestName)
setwd(file.path(getwd(),TestName))

# Specify possible fixed values of catch selectivity in final age in last time block
sel9values <- c(0.85,0.90,0.925,0.95,0.99)
models <- vector("list",length(sel9values))

for(i in 1:length(sel9values))
{
  input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                               selectivity=list(model=rep("age-specific",5),
                                                re=rep("none",5),
                                                initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                                  c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                                  c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,sel9vals[i]),
                                                                  c(0.2,0.4,0.8,0.999,0.999,0.999,0.999,0.999,0.999),
                                                                  c(0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999,0.999)),
                                                fix_pars=list(8:9,8:9,7:9,4:9,6:9)))
  
  # Mapping the estimation of numbers at age in the first year:
  # input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
  input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
  # input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8
  
  # Fit model:
  m1 <- fit_wham(input, do.osa = F, do.check=T) # turn off OSA residuals to save time
  models[[i]] <- m1
  #Look at parameter estimates:
  # models[i]$sdrep
  
  # Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradient should be < 1e-06)
  # check_convergence(models[i])
  
  # Plot results
  # plot_wham_output(mod=m1, out.type="html")
}

# Compare models
compare_wham_models(models,table.opts=list(fname="",sort=T))



##################################### Epilogue/Footer #####################################

#Example WHAM code:
library(wham)
wham.dir <- find.package("wham")
file.path(wham.dir,"example_scripts")

setwd("C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM")
source(file.path(wham.dir,"example_scripts","ex1_basics.R"))



