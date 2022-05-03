# Gulf of Maine Haddock in WHAM
# Emily Liljestrand
# Created: May 2, 2022

# Remove Old Objects and set working directory to file location:
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load Packages:
library(wham)
library(TMB)

##################################### DATA #####################################

#Read in Gulf of Maine Haddock Data from ASAP Dat File
GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB.DAT")

#Look at dat file:
GOM_HADDOCK_DAT$dat

##################################### MODEL #####################################
# create directory for analysis, E.g.,
#write.dir <- "/path/to/save/output"
if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)
file.path(getwd(),"Model1")


input1 <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                             selectivity=list(model=rep("age-specific",5),
                                              re=rep("none",5),
                                              initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                                c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999),
                                                                c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.9),
                                                                c(0.2,0.4,0.8,0.999,0.999,0.999,0.999,0.999,0.999),
                                                                c(0.1,0.3,0.5,0.8,0.9,0.999,0.999,0.999,0.999)),
                                              fix_pars=list(8:9,8:9,7:8,4:9,6:9)))

#Mapping the estimation of numbers at age in the first year:
input1$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
input1$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
input1$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8

m1 <- fit_wham(input1, do.osa = F, do.check=T) # turn off OSA residuals to save time
#Look at parameter estimates:
m1$sdrep

# Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradient should be < 1e-06)
check_convergence(m1)

plot_wham_output(mod=m1, out.type="html")

# Helpful code for later:
# mods <- list(m1=m1,...)
# compare_wham_models(mods,table.opts=list(fname="",sort=T))
# m1_proj <- project_wham
# plot_wham_output(mod=m1_proj, out.type="html")


#Example WHAM code:
library(wham)
wham.dir <- find.package("wham")
file.path(wham.dir,"example_scripts")

setwd("C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM")
source(file.path(wham.dir,"example_scripts","ex1_basics.R"))



