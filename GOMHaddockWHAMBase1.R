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

##################################### MODEL(s) #####################################

#------------------------------------ Test 0: Compare WHAM and ASAP ------------------------------------
# Make a folder to put the results of Test 0:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
TestName <- "Test0_Comparison"
if(!dir.exists(TestName)) dir.create(TestName)
setwd(file.path(getwd(),TestName))

# Read in ASAP model fitted values
ASAP <- read_asap3_fit(wd="C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM/GOMHaddockASAP",asap.name="GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB")

compopts <- c('multinomial','logistic-normal-miss0','logistic-normal-pool0')
# Set up a data frame for information of model options:
df.mods <- data.frame(fleets=rep(compopts,length(compopts)^2),
                      index1=rep(compopts,each=length(compopts)),
                      index2=rep(compopts,each=length(compopts)^2))
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

for(m in 1:n.mods){
  # input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment")
  input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                              selectivity=list(model=rep("age-specific",5),
                                               re=rep("none",5),
                                               initial_pars=list(c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                                 c(0.01,0.1,0.3,0.5,0.8,0.9,0.999,1,1),
                                                                 c(0.01,0.1,0.3,0.5,0.8,0.9,1,1,0.999),
                                                                 c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                 c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                               fix_pars=list(8:9,8:9,7:8,4:9,6:9)),
                              age_comp=list(fleets=df.mods[m,"fleets"],indices=c(df.mods[m,"index1"],df.mods[m,"index2"])))
  
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
