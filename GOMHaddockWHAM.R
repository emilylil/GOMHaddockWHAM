# Gulf of Maine Haddock in WHAM
# Emily Liljestrand
# Created: May 2, 2022

# Remove Old Objects and set working directory to file location:
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load Packages:
library(wham)

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

input1 <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=1, model_name="GOMHaddock RW Recruitment")


input1 <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=1, model_name="GOMHaddock_RW_Recruitment",
                             selectivity=list(model=rep("age-specific",5),
                                              re=rep("none",5),
                                              initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,1,1,1),
                                                                c(0.01,0.1,0.3,0.5,0.8,0.9,1,1,1),
                                                                c(0.01,0.1,0.3,0.5,0.8,0.9,1,1,1),
                                                                c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                              fix_pars=list(8:9,8:9,7:9,4:9,6:9)))
m1 <- fit_wham(input1, do.osa = F) # turn off OSA residuals to save time

# Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradiet should be < 1e-06)
check_convergence(m1)



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



