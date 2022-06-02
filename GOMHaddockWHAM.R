# Gulf of Maine Haddock in WHAM
# Emily Liljestrand
# Created: May 2, 2022
# Updates: June 2, 2022

# Remove Old Objects and set working directory to file location:
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
originalwd = getwd()

# uninstall("wham")
#Load Packages:
# devtools::install_github("emilylil/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
# devtools::install_github("timjmiller/wham", dependencies=TRUE,ref='dvel')
# install.packages("TMB")

library(wham)
library(TMB)
library(dplyr)
library(ggplot2)

##################################### DATA #####################################

#Read in Gulf of Maine Haddock Data from ASAP Dat File
# GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB.DAT")
#Read in Gulf of Maine Haddock Data from ASAP Dat File, version for selectivity random effects (one block for all catch)
GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_1BLOCK_NEWCALIB.DAT")

#Read in existing seeds file, or make one if it doesn't exist
if(!file.exists("SIM.seeds.csv"))
{
  set.seed(799291)
  nseeds <- 1000
  r.seed.set <- trunc(1e7*runif(n=nseeds),7) + trunc(1e3*runif(n=nseeds),3)
  write.table(r.seed.set, file="SIM.seeds.csv", quote=F, row.names=F, col.names=F, sep=",")
} else(r.seed.set <- read.table("SIM.seeds.csv"))

# ##################################### MODEL(s) #####################################
#
# ------------------------------------ OM with selectivity random effects ------------------------------------

# Make a folder to put the particular OM results:
# OM.Var.Name <- "OperatingModelRecM.SurM.SelM"
# if(!dir.exists(OM.Var.Name)) dir.create(OM.Var.Name)
# setwd(file.path(getwd(),OM.Var.Name))

# Set and fit the OM, if the file doesn't exist
input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                            selectivity=list(model=rep("age-specific",3),
                                             re=c(rep('ar1_y',1),rep("none",2)),
                                             initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15,0.15),
                                                                  c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                  c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                             fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                            NAA_re=list(sigma='rec+1',cor='ar1_y'),
                            age_comp='multinomial')

input$par$log_F1 <- log(3)
input$map$log_F1 <- as.factor(NA)
input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))

# Mapping the estimation of numbers at age in the first year:
# input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
# input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,NA),nrow=1,ncol=9)) # Fix ages 7-8
input$par$log_N1_pars[7:8] <- c(log(5),log(1))
input$par$log_N1_pars[9] <- log(20)

# Fit model:
if(!file.exists("mod.OM.rds"))
{
  mod.OM <- try(fit_wham(input, do.osa = F, do.check=F,do.retro=T)) # turn off OSA residuals to save time
  saveRDS(mod.OM, file='mod.OM.rds')
  mod.OM.proj <- project_wham(mod.OM,proj.opts=list(use.FXSPR=T))
  saveRDS(mod.OM.proj, file='mod.OM.proj.rds')

  realwd <- getwd()
  Foldername <- "mod.OM.Output"
  if(!dir.exists(Foldername)) dir.create(Foldername)
  setwd(file.path(getwd(),Foldername))
  plot_wham_output(mod=mod.OM,out.type="html")
  setwd(realwd)

} else {mod.OM<-readRDS('mod.OM.rds')
        mod.OM.proj <-readRDS('mod.OM.proj.rds')}

# ------------------------------------ Make Simulation Folder and Establish Sims ------------------------------------

# For the operating model, if the process error variance is low (0.75), baseline (1), or high (1.25)
recruitment_RE_mag <- 1
survival_RE_mag <- 1
selectivity_RE_mag <- 1

# Set a name-specific estimation model folder for each combination of process errors
if(recruitment_RE_mag<1) recname="L"
if(recruitment_RE_mag>1) recname="H"
if(recruitment_RE_mag=1) recname="M"
if(survival_RE_mag<1) survivalname="L"
if(survival_RE_mag>1) survivalname="H"
if(survival_RE_mag=1) survivalname="M"
if(selectivity_RE_mag<1) selectivityname="L"
if(selectivity_RE_mag>1) selectivityname="H"
if(selectivity_RE_mag=1) selectivityname="M"

OM.Var.Name <- paste0("OperatingModelRec",recname,".Sur",survivalname,".Sel",selectivityname)
setwd(originalwd)
if(!dir.exists(OM.Var.Name)) dir.create(OM.Var.Name)
setwd(file.path(originalwd,OM.Var.Name))

# Set number of simulations and the age comps of the catch and indices
nsims <- 100
simsvec <- c(1:nsims)
agecomp <- 'multinomial'
catch_Neff_adjust <- 1
index_Neff_adjust <- 1

# Make a folder to put the particular OM results:
OM.CorLik.Name <- "SimResults_AR1_MULT"
if(!dir.exists(OM.CorLik.Name)) dir.create(OM.CorLik.Name)
setwd(file.path(getwd(),OM.CorLik.Name))

# Set the combinations of recruitment, survival, selectivity random effects in estimation models
recruitment_RE <- c(0,1)
survival_RE <- c(0,1)
selectivity_RE <- c(0,1)

df.scenarios <- as.data.frame(tidyr::crossing(recruitment_RE,survival_RE,selectivity_RE))
df.scenarios <- dplyr::mutate(df.scenarios,'scenario'=paste0(recruitment_RE,survival_RE,selectivity_RE))
# Remove some combinations of estimation models (always remove the survival without recruitment models)
df.scenarios <- df.scenarios[-which(df.scenarios$scenario=='010'),]
df.scenarios <- df.scenarios[-which(df.scenarios$scenario=='011'),]
# df.scenarios <- df.scenarios[-which(df.scenarios$scenario=='001'),]
# df.scenarios <- df.scenarios[-which(df.scenarios$scenario=='111'),]
# df.scenarios <- df.scenarios[-which(df.scenarios$scenario=='110'),]
# df.scenarios <- df.scenarios[-which(df.scenarios$scenario=='000'),]
# df.scenarios <- df.scenarios[-which(df.scenarios$scenario=='100'),]
# df.scenarios <- df.scenarios[-which(df.scenarios$scenario=='101'),]

# Initiate the simsRE csv file which will collect all results from all scenarios and simulations
n.scenarios <- dim(df.scenarios)[1]
simsRE <-  as.data.frame(tidyr::crossing(scenario=df.scenarios$scenario,simsvec))
simsRE$FAA_RE <-rep('',nsims*n.scenarios)
simsRE$FAA_ARE <- rep('',nsims*n.scenarios)
simsRE$NAA_RE <- rep('',nsims*n.scenarios)
simsRE$NAA_ARE <- rep('',nsims*n.scenarios)
simsRE$SSB_RE <- rep('',nsims*n.scenarios)
simsRE$SSB_ARE <- rep('',nsims*n.scenarios)
simsRE$Rec_Var <- rep('',nsims*n.scenarios)
simsRE$Sur_Var <- rep('',nsims*n.scenarios)
simsRE$Sur_Cor <- rep('',nsims*n.scenarios)
simsRE$Sel_Var <- rep('',nsims*n.scenarios)
simsRE$Sel_Cor <- rep('',nsims*n.scenarios)
simsRE$mohnsSSB <- rep('',nsims*n.scenarios)
simsRE$mohnsFbar <- rep('',nsims*n.scenarios)
simsRE$mohnsR<- rep('',nsims*n.scenarios)
simsRE$mohnsN2 <- rep('',nsims*n.scenarios)
simsRE$mohnsN3 <- rep('',nsims*n.scenarios)
simsRE$mohnsN4 <- rep('',nsims*n.scenarios)
simsRE$mohnsN5 <- rep('',nsims*n.scenarios)
simsRE$mohnsN6 <- rep('',nsims*n.scenarios)
simsRE$mohnsN7 <- rep('',nsims*n.scenarios)
simsRE$mohnsN8 <- rep('',nsims*n.scenarios)
simsRE$mohnsN9 <- rep('',nsims*n.scenarios)
simsRE$F40Catch1 <- rep('',nsims*n.scenarios)
simsRE$F40Catch2 <- rep('',nsims*n.scenarios)
simsRE$F40Catch3 <- rep('',nsims*n.scenarios)
simsRE$F40Catch1_RE <- rep('',nsims*n.scenarios)
simsRE$F40Catch2_RE <- rep('',nsims*n.scenarios)
simsRE$F40Catch3_RE <- rep('',nsims*n.scenarios)
simsRE$scenarioname <- rep('',nsims*n.scenarios)
simsRE$convergence <- rep(1,nsims*n.scenarios)
simsRE$hessian1 <- rep(1,nsims*n.scenarios)
simsRE$hessian2 <- rep(1,nsims*n.scenarios)

for(s in 1:n.scenarios)
{
  # Set a name-specific estimation model folder for each combination of process errors
  if(df.scenarios[s,'recruitment_RE']==1) recname="Rec" else recname=""
  if(df.scenarios[s,'survival_RE']==1) survivalname="Sur" else survivalname=""
  if(df.scenarios[s,'selectivity_RE']==1) selectivityname="Sel" else selectivityname=""
  if(df.scenarios[s,'recruitment_RE']==0&&df.scenarios[s,'survival_RE']==0&&df.scenarios[s,'selectivity_RE']==0) recname="None"
  if(df.scenarios[s,'survival_RE']==1) NAA_re_sigma="rec+1" else NAA_re_sigma="rec"

  # Make a folder to put the particular EM results:
  EM.Name <- paste0(recname,survivalname,selectivityname)
  if(!dir.exists(EM.Name)) dir.create(EM.Name)
  setwd(file.path(originalwd,OM.Var.Name,OM.CorLik.Name,EM.Name))
  
  df.mods.sim <- data.frame("Sim"=rep(0,nsims))
  df.mods.sim$Sim <- paste0("Sim",1:nsims)
  df.mods.sim <- df.mods.sim %>% select(Sim, everything()) # moves Model to first col
  
  # Change some process error variance components before simulating data:
  simpar <- mod.OM$env$last.par.best
  simpar[names(simpar)=="log_NAA_sigma"] = c(log(exp(simpar[names(simpar)=="log_NAA_sigma"][1])*recruitment_RE_mag),log(exp(simpar[names(simpar)=="log_NAA_sigma"][2])*survival_RE_mag))
  simpar[names(simpar)=='sel_repars'] = c(log(exp(simpar[names(simpar)=='sel_repars'][1])*selectivity_RE_mag),simpar[names(simpar)=='sel_repars'][2])

  # Change some observation error variance components before simulating data:
  # mod.OM$env$data$agg_catch_sigma <- mod.OM$env$data$agg_catch_sigma
  # mod.OM$env$data$agg_index_sigma <- mod.OM$env$data$agg_index_sigma
  mod.OM$env$data$catch_Neff <- mod.OM$env$data$catch_Neff * catch_Neff_adjust
  mod.OM$env$data$index_Neff <- mod.OM$env$data$index_Neff * index_Neff_adjust
  
  #------------------------------------ First Round: Regular Test ------------------------------------
  nohessian <- c()
  noconverge <- c()
  for(m in c(1:nsims))
  {
    simsRE$scenarioname[(s-1)*nsims+m] <- EM.Name
    set.seed(r.seed.set[m,1])
    simdata <- mod.OM$simulate(par=simpar,complete=T)
    
    if(df.scenarios[s,'recruitment_RE']==1)
    {
      if(df.scenarios[s,'selectivity_RE']==1)
      {
        input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                                    selectivity=list(model=rep("age-specific",3),
                                                     re=c(rep('ar1_y',1),rep("none",2)),
                                                     initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.33,0.33,0.33),
                                                                          c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                          c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                                     fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                                     # fix_pars=list(c(7:9),c(4:9),c(6:9))),
                                    NAA_re=list(sigma=NAA_re_sigma,cor='ar1_y'),
                                    age_comp=agecomp)
        
        input$par$log_F1 <- log(3)
        input$map$log_F1 <- as.factor(NA)
        input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))
      }
      if(df.scenarios[s,'selectivity_RE']==0)
      {
        input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                                    selectivity=list(model=rep("age-specific",3),
                                                     re=c(rep('none',1),rep("none",2)),
                                                     initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.33,0.33,0.33),
                                                                          c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                          c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                                     # fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                                     fix_pars=list(c(7:9),c(4:9),c(6:9))),
                                    NAA_re=list(sigma=NAA_re_sigma,cor='ar1_y'),
                                    age_comp=agecomp)
      }
    }
    if(df.scenarios[s,'recruitment_RE']==0)
    {
      if(df.scenarios[s,'selectivity_RE']==1)
      {
        input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                                    selectivity=list(model=rep("age-specific",3),
                                                     re=c(rep('ar1_y',1),rep("none",2)),
                                                     initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.33,0.33,0.33),
                                                                          c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                          c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                                     fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                    # fix_pars=list(c(7:9),c(4:9),c(6:9))),
                                    # NAA_re=list(sigma=NAA_re_sigma,cor='ar1_y'),
                                    age_comp=agecomp)
        
        input$par$log_F1 <- log(3)
        input$map$log_F1 <- as.factor(NA)
        input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))
      }
      if(df.scenarios[s,'selectivity_RE']==0)
      {
        input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                                    selectivity=list(model=rep("age-specific",3),
                                                     re=c(rep('none',1),rep("none",2)),
                                                     initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.33,0.33,0.33),
                                                                          c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                          c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                                     # fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                                     fix_pars=list(c(7:9),c(4:9),c(6:9))),
                                    # NAA_re=list(sigma=NAA_re_sigma,cor='ar1_y'),
                                    age_comp=agecomp)
      }
    }
    
    input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
    input$par$log_N1_pars[7:8] <- c(log(5),log(1))
    input$par$log_N1_pars[9] <- log(20)
    
    siminput <- input
    
    siminput$data$catch_paa <- simdata$catch_paa
    siminput$data$index_paa <- simdata$index_paa
    siminput$data$agg_catch <- simdata$agg_catch
    siminput$data$agg_indices <- simdata$agg_indices
    siminput$data$obsvec <- simdata$obsvec
    
    simmod.proj.real.input <-  mod.OM
    simmod.proj.real.input$input$data <- siminput$data
    simmod.proj.real.input$input$par <- siminput$par
    simmod.proj.real.input$input$map <- siminput$map
    simmod.proj.real <- project_wham(simmod.proj.real.input,proj.opts=list(use.FXSPR=T))
    
    skip_to_next<- FALSE
    simmod <- tryCatch(fit_wham(siminput, do.osa = F, do.check=F,do.retro=T),error=function(e){skip_to_next <<- T}) # turn off OSA residuals to save time
    if(skip_to_next) {next}
    
    if(!(simmod$na_sdrep==FALSE & !is.na(simmod$na_sdrep)))
    {
      nohessian <- c(nohessian,m)
      simsRE$hessian1[(s-1)*nsims+m] <- 0
    }
    if(simmod$opt$convergence==1)
    {
      noconverge <- c(noconverge,m)
      simsRE$convergence[(s-1)*nsims+m] <- 0
    }
    
    simmod.proj <- project_wham(simmod,proj.opts=list(use.FXSPR=T))
    
    # Collect all the relative and absolute relative errors between estimated and true values: 
    simmod$simfit$FAA_RE <- (simmod$rep$FAA_tot - simdata$FAA[,,])/ simdata$FAA[,,]
    simmod$simfit$FAA_ARE <- abs(simmod$rep$FAA_tot - simdata$FAA[,,])/ simdata$FAA[,,]
    simmod$simfit$NAA_RE <- (simmod$rep$NAA - simdata$NAA) / simdata$NAA
    simmod$simfit$NAA_ARE <- abs(simmod$rep$NAA - simdata$NAA) / simdata$NAA
    simmod$simfit$SSB_RE <- (simmod$rep$SSB - simdata$SSB) / simdata$SSB 
    simmod$simfit$SSB_ARE <- abs(simmod$rep$SSB - simdata$SSB) / simdata$SSB 
    simmod$simfit$Rec_Var <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='log_NAA_sigma')[1]])
    simmod$simfit$Sur_Var <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='log_NAA_sigma')[2]])
    simmod$simfit$Sur_Cor <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='trans_NAA_rho')[1]])/(1+exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='trans_NAA_rho')[1]]))
    simmod$simfit$Sel_Var <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='sel_repars')[1]])
    simmod$simfit$Sel_Cor <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='sel_repars')[2]])/(1+exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='sel_repars')[2]]))
    simmod$simfit$mohns <- try(mohns_rho(simmod))
    simmod$simfit$F40Catch <- simmod.proj$rep$catch_proj
    simmod$simfit$F40Catch1_RE <- (simmod.proj$rep$catch_proj[1]-simmod.proj.real$rep$catch_proj[1])/simmod.proj.real$rep$catch_proj[1]
    simmod$simfit$F40Catch2_RE <- (simmod.proj$rep$catch_proj[2]-simmod.proj.real$rep$catch_proj[2])/simmod.proj.real$rep$catch_proj[2]
    simmod$simfit$F40Catch3_RE <- (simmod.proj$rep$catch_proj[3]-simmod.proj.real$rep$catch_proj[3])/simmod.proj.real$rep$catch_proj[3]
    
    saveRDS(simmod, file=paste0(df.mods.sim$Sim[m],".rds"))
    saveRDS(simmod.proj, file=paste0(df.mods.sim$Sim[m],".rds",".proj"))
    
    simsRE$FAA_RE[(s-1)*nsims+m] <- mean(simmod$simfit$FAA_RE[,8])
    simsRE$FAA_ARE[(s-1)*nsims+m] <- mean(simmod$simfit$FAA_ARE[,8])
    simsRE$NAA_RE[(s-1)*nsims+m] <- mean(simmod$simfit$NAA_RE[,1])
    simsRE$NAA_ARE[(s-1)*nsims+m] <- mean(simmod$simfit$NAA_ARE[,1])
    simsRE$SSB_RE[(s-1)*nsims+m] <- mean(simmod$simfit$SSB_RE)
    simsRE$SSB_ARE[(s-1)*nsims+m] <- mean(simmod$simfit$SSB_ARE)
    simsRE$mohnsSSB[(s-1)*nsims+m] <- simmod$simfit$mohns[1]
    simsRE$mohnsFbar[(s-1)*nsims+m] <- simmod$simfit$mohns[2]
    simsRE$mohnsR[(s-1)*nsims+m] <- simmod$simfit$mohns[3]
    simsRE$mohnsN2[(s-1)*nsims+m] <- simmod$simfit$mohns[4]
    simsRE$mohnsN3[(s-1)*nsims+m] <- simmod$simfit$mohns[5]
    simsRE$mohnsN4[(s-1)*nsims+m] <- simmod$simfit$mohns[6]
    simsRE$mohnsN5[(s-1)*nsims+m] <- simmod$simfit$mohns[7]
    simsRE$mohnsN6[(s-1)*nsims+m] <- simmod$simfit$mohns[8]
    simsRE$mohnsN7[(s-1)*nsims+m] <- simmod$simfit$mohns[9]
    simsRE$mohnsN8[(s-1)*nsims+m] <- simmod$simfit$mohns[10]
    simsRE$mohnsN9[(s-1)*nsims+m] <- simmod$simfit$mohns[11]
    simsRE$F40Catch1[(s-1)*nsims+m] <- simmod$simfit$F40Catch[1]
    simsRE$F40Catch2[(s-1)*nsims+m] <- simmod$simfit$F40Catch[2]
    simsRE$F40Catch3[(s-1)*nsims+m] <- simmod$simfit$F40Catch[3]
    simsRE$F40Catch1_RE[(s-1)*nsims+m] <- simmod$simfit$F40Catch1_RE
    simsRE$F40Catch2_RE[(s-1)*nsims+m] <- simmod$simfit$F40Catch2_RE
    simsRE$F40Catch3_RE[(s-1)*nsims+m] <- simmod$simfit$F40Catch3_RE
    
    if(df.scenarios[s,'recruitment_RE']==1) simsRE$Rec_Var[(s-1)*nsims+m] <- simmod$simfit$Rec_Var
    if(df.scenarios[s,'survival_RE']==1) simsRE$Sur_Cor[(s-1)*nsims+m] <- simmod$simfit$Sur_Cor
    if(df.scenarios[s,'survival_RE']==1) simsRE$Sur_Var[(s-1)*nsims+m] <- simmod$simfit$Sur_Var
    if(df.scenarios[s,'selectivity_RE']==1) simsRE$Sel_Var[(s-1)*nsims+m] <- simmod$simfit$Sel_Var
    if(df.scenarios[s,'selectivity_RE']==1) simsRE$Sel_Cor[(s-1)*nsims+m] <- simmod$simfit$Sel_Cor
    
    write.csv(simsRE[((s-1)*nsims+1):((s-1)*nsims+100),],file='simsRE.csv',append=F)
    
    # Plot output in new subfolder
    # plot_wham_output(mod=simmod, dir.main=file.path(getwd(),df.mods.sim$Sim[m]), out.type='html')
  }
  
  #------------------------------------ Second Round: Test with fixed abundance age 9 first year------------------------------------
  # If the hessian wasn't made (NA's in the SD estimates), try again with fixing the abundance at age 9 in first year
  stillnohessian <- c()
  for(m in nohessian)
  {
    set.seed(r.seed.set[m,1])
    simdata <- mod.OM$simulate(par=simpar,complete=T)
    
    if(df.scenarios[s,'recruitment_RE']==1)
    {
      if(df.scenarios[s,'selectivity_RE']==1)
      {
        input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                                    selectivity=list(model=rep("age-specific",3),
                                                     re=c(rep('ar1_y',1),rep("none",2)),
                                                     initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.33,0.33,0.33),
                                                                          c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                          c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                                     fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                    # fix_pars=list(c(7:9),c(4:9),c(6:9))),
                                    NAA_re=list(sigma=NAA_re_sigma,cor='ar1_y'),
                                    age_comp=agecomp)
        
        input$par$log_F1 <- log(3)
        input$map$log_F1 <- as.factor(NA)
        input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))
      }
      if(df.scenarios[s,'selectivity_RE']==0)
      {
        input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                                    selectivity=list(model=rep("age-specific",3),
                                                     re=c(rep('none',1),rep("none",2)),
                                                     initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.33,0.33,0.33),
                                                                          c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                          c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                                     # fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                                     fix_pars=list(c(7:9),c(4:9),c(6:9))),
                                    NAA_re=list(sigma=NAA_re_sigma,cor='ar1_y'),
                                    age_comp=agecomp)
      }
    }
    if(df.scenarios[s,'recruitment_RE']==0)
    {
      if(df.scenarios[s,'selectivity_RE']==1)
      {
        input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                                    selectivity=list(model=rep("age-specific",3),
                                                     re=c(rep('ar1_y',1),rep("none",2)),
                                                     initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.33,0.33,0.33),
                                                                          c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                          c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                                     fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                    # fix_pars=list(c(7:9),c(4:9),c(6:9))),
                                    # NAA_re=list(sigma=NAA_re_sigma,cor='ar1_y'),
                                    age_comp=agecomp)
        
        input$par$log_F1 <- log(3)
        input$map$log_F1 <- as.factor(NA)
        input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))
      }
      if(df.scenarios[s,'selectivity_RE']==0)
      {
        input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                                    selectivity=list(model=rep("age-specific",3),
                                                     re=c(rep('none',1),rep("none",2)),
                                                     initial_pars <- list(c(0.1,0.1,0.1,0.1,0.15,0.15,0.33,0.33,0.33),
                                                                          c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                          c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                                     # fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                                     fix_pars=list(c(7:9),c(4:9),c(6:9))),
                                    # NAA_re=list(sigma=NAA_re_sigma,cor='ar1_y'),
                                    age_comp=agecomp)
      }
    }
    
    input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,NA),nrow=1,ncol=9)) # Fix ages 7-8
    input$par$log_N1_pars[7:8] <- c(log(5),log(1))
    input$par$log_N1_pars[9] <- log(20)
    
    siminput <- input
    
    # siminput$data <- simdata
    siminput$data$catch_paa <- simdata$catch_paa
    siminput$data$index_paa <- simdata$index_paa
    siminput$data$agg_catch <- simdata$agg_catch
    siminput$data$agg_indices <- simdata$agg_indices
    siminput$data$obsvec <- simdata$obsvec
    
    simmod.proj.real.input <-  mod.OM
    simmod.proj.real.input$input$data <- siminput$data
    simmod.proj.real.input$input$par <- siminput$par
    simmod.proj.real.input$input$map <- siminput$map
    simmod.proj.real <- project_wham(simmod.proj.real.input,proj.opts=list(use.FXSPR=T))
    
    skip_to_next<- FALSE
    simmod <- tryCatch(fit_wham(siminput, do.osa = F, do.check=F,do.retro=T),error=function(e){skip_to_next <<- T}) # turn off OSA residuals to save time
    if(skip_to_next) {next}
    
    simmod.proj <- project_wham(mod.OM,proj.opts=list(use.FXSPR=T))
    
    # Collect all the relative and absolute relative errors between estimated and true values: 
    if(simmod$na_sdrep==FALSE & !is.na(simmod$na_sdrep))
    {
      simmod$simfit$FAA_RE <- (simmod$rep$FAA_tot - simdata$FAA[,,])/ simdata$FAA[,,]
      simmod$simfit$FAA_ARE <- abs(simmod$rep$FAA_tot - simdata$FAA[,,])/ simdata$FAA[,,]
      simmod$simfit$NAA_RE <- (simmod$rep$NAA - simdata$NAA) / simdata$NAA
      simmod$simfit$NAA_ARE <- abs(simmod$rep$NAA - simdata$NAA) / simdata$NAA
      simmod$simfit$SSB_RE <- (simmod$rep$SSB - simdata$SSB) / simdata$SSB 
      simmod$simfit$SSB_ARE <- abs(simmod$rep$SSB - simdata$SSB) / simdata$SSB 
      simmod$simfit$Rec_Var <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='log_NAA_sigma')[1]])
      simmod$simfit$Sur_Var <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='log_NAA_sigma')[2]])
      simmod$simfit$Sur_Cor <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='trans_NAA_rho')[1]])/(1+exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='trans_NAA_rho')[1]]))
      simmod$simfit$Sel_Var <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='sel_repars')[1]])
      simmod$simfit$Sel_Cor <- exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='sel_repars')[2]])/(1+exp(simmod$env$last.par.best[which(names(simmod$env$last.par.best)=='sel_repars')[2]]))
      simmod$simfit$mohns <- mohns_rho(simmod)
      simmod$simfit$F40Catch <- simmod.proj$rep$catch_proj
      simmod$simfit$F40Catch1_RE <- (simmod.proj$rep$catch_proj[1]-simmod.proj.real$rep$catch_proj[1])/simmod.proj.real$rep$catch_proj[1]
      simmod$simfit$F40Catch2_RE <- (simmod.proj$rep$catch_proj[2]-simmod.proj.real$rep$catch_proj[2])/simmod.proj.real$rep$catch_proj[2]
      simmod$simfit$F40Catch3_RE <- (simmod.proj$rep$catch_proj[3]-simmod.proj.real$rep$catch_proj[3])/simmod.proj.real$rep$catch_proj[3]
      
    }
    
    if(!(simmod$na_sdrep==FALSE & !is.na(simmod$na_sdrep)))
    {
      stillnohessian <- c(stillnohessian,m)
      simsRE$hessian2[(s-1)*nsims+m] <- 0
    }
    else{
      saveRDS(simmod, file=paste0(df.mods.sim$Sim[m],".rds"))
      saveRDS(simmod.proj, file=paste0(df.mods.sim$Sim[m],".rds",".proj"))
      
      simsRE$FAA_RE[(s-1)*nsims+m] <- mean(simmod$simfit$FAA_RE[,8])
      simsRE$FAA_ARE[(s-1)*nsims+m] <- mean(simmod$simfit$FAA_ARE[,8])
      simsRE$NAA_RE[(s-1)*nsims+m] <- mean(simmod$simfit$NAA_RE[,1])
      simsRE$NAA_ARE[(s-1)*nsims+m] <- mean(simmod$simfit$NAA_ARE[,1])
      simsRE$SSB_RE[(s-1)*nsims+m] <- mean(simmod$simfit$SSB_RE)
      simsRE$SSB_ARE[(s-1)*nsims+m] <- mean(simmod$simfit$SSB_ARE)
      simsRE$mohnsSSB[(s-1)*nsims+m] <- simmod$simfit$mohns[1]
      simsRE$mohnsFbar[(s-1)*nsims+m] <- simmod$simfit$mohns[2]
      simsRE$mohnsR[(s-1)*nsims+m] <- simmod$simfit$mohns[3]
      simsRE$mohnsN2[(s-1)*nsims+m] <- simmod$simfit$mohns[4]
      simsRE$mohnsN3[(s-1)*nsims+m] <- simmod$simfit$mohns[5]
      simsRE$mohnsN4[(s-1)*nsims+m] <- simmod$simfit$mohns[6]
      simsRE$mohnsN5[(s-1)*nsims+m] <- simmod$simfit$mohns[7]
      simsRE$mohnsN6[(s-1)*nsims+m] <- simmod$simfit$mohns[8]
      simsRE$mohnsN7[(s-1)*nsims+m] <- simmod$simfit$mohns[9]
      simsRE$mohnsN8[(s-1)*nsims+m] <- simmod$simfit$mohns[10]
      simsRE$mohnsN9[(s-1)*nsims+m] <- simmod$simfit$mohns[11]
      simsRE$F40Catch1[(s-1)*nsims+m] <- simmod$simfit$F40Catch[1]
      simsRE$F40Catch2[(s-1)*nsims+m] <- simmod$simfit$F40Catch[2]
      simsRE$F40Catch3[(s-1)*nsims+m] <- simmod$simfit$F40Catch[3]
      simsRE$hessian2[(s-1)*nsims+m] <- 1
      simsRE$F40Catch1_RE[(s-1)*nsims+m] <- simmod$simfit$F40Catch1_RE
      simsRE$F40Catch2_RE[(s-1)*nsims+m] <- simmod$simfit$F40Catch2_RE
      simsRE$F40Catch3_RE[(s-1)*nsims+m] <- simmod$simfit$F40Catch3_RE
      
      if(df.scenarios[s,'recruitment_RE']==1) simsRE$Rec_Var[(s-1)*nsims+m] <- simmod$simfit$Rec_Var
      if(df.scenarios[s,'survival_RE']==1) simsRE$Sur_Cor[(s-1)*nsims+m] <- simmod$simfit$Sur_Cor
      if(df.scenarios[s,'survival_RE']==1) simsRE$Sur_Var[(s-1)*nsims+m] <- simmod$simfit$Sur_Var
      if(df.scenarios[s,'selectivity_RE']==1) simsRE$Sel_Var[(s-1)*nsims+m] <- simmod$simfit$Sel_Var
      if(df.scenarios[s,'selectivity_RE']==1) simsRE$Sel_Cor[(s-1)*nsims+m] <- simmod$simfit$Sel_Cor
      
      write.csv(simsRE[((s-1)*nsims+1):((s-1)*nsims+100),],file='simsRE.csv',append=F)
    }
    # Plot output in new subfolder
    # plot_wham_output(mod=simmod, dir.main=file.path(getwd(),df.mods.sim$Sim[m]), out.type='html')
  }
  write.table(noconverge,file='noconverge.txt')
  write.table(nohessian,file='nohessian.txt')
  write.table(stillnohessian,file='stillnohessian.txt')

  #------------------------------------ Summarize and Visualize Results within Sims------------------------------------
  n.fail <- length(stillnohessian)
  goodmodels <- c(1:nsims)
  if(n.fail>0) goodmodels <- goodmodels[-stillnohessian]
  if(n.fail>0)  df.mods.sim <- data.frame(Sim=df.mods.sim[goodmodels,])
  
  mod.list <- paste0(df.mods.sim$Sim,".rds")
  mods <- lapply(mod.list, readRDS)
  names(mods) <- goodmodels
  # mods<-append(mods,mod.OM)
  
  # Add columns for convergence information and if hessian was produced:
  opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
  ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
  df.mods.sim$conv <- as.logical(opt_conv)
  df.mods.sim$pdHess <- as.logical(ok_sdrep)
  
  # Compare across models
  write.csv(df.mods.sim, file="df.mods.sim.csv",quote=F, row.names=F)
  
  # Plot error across years for F at age 8, R, and SSB
  years<- c(mod.OM$years)
  sims <- goodmodels
  FAA_RE <- c()
  NAA_RE <- c()
  SSB_RE <- c()
  for(m in 1:length(goodmodels)) FAA_RE <- c(FAA_RE,array(mods[[m]]$simfit$FAA_RE[,8]))
  for(m in 1:length(goodmodels)) NAA_RE <- c(NAA_RE,array(mods[[m]]$simfit$NAA_RE[,1]))
  for(m in 1:length(goodmodels)) SSB_RE <- c(SSB_RE,array(mods[[m]]$simfit$SSB_RE))
  ModsRE.df <- as.data.frame(tidyr::crossing(sims,years))
  ModsRE.df$FAA_RE <- FAA_RE
  ModsRE.df$NAA_RE <- NAA_RE
  ModsRE.df$SSB_RE <- SSB_RE
  
  jpeg('FAA_RE.jpg',width=800,height=400)
    plot(ggplot(data=ModsRE.df,aes(x=as.factor(years),y=FAA_RE)) + geom_violin() + theme_classic()+ geom_hline(yintercept=0, color = "red", size=1)
    + stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black"))
  dev.off()
  
  jpeg('NAA_RE.jpg',width=800,height=400)
    plot(ggplot(data=ModsRE.df,aes(x=as.factor(years),y=NAA_RE)) + geom_violin() + theme_classic()+ geom_hline(yintercept=0, color = "red", size=1)
    + stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black"))
  dev.off()
  
  jpeg('SSB_RE.jpg',width=800,height=400)
    plot(ggplot(data=ModsRE.df,aes(x=as.factor(years),y=SSB_RE)) + geom_violin() + theme_classic()+ geom_hline(yintercept=0, color = "red", size=1)
    + stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black"))
  dev.off()
  
}

#------------------------------------ Summarize and Visualize Results across Sims------------------------------------

setwd(file.path(originalwd,OM.Var.Name,OM.CorLik.Name))
write.csv(simsRE,file='simsREFull.csv')

# Read in existing file if it needed manipulating outside of R
# simsRE <- read.csv(file="simsREFull.csv",head=T)

# Extract the number of successful runs, as defined by not having a failed second hessian
successes <- data.frame(simsRE %>% group_by(scenarioname) %>% summarize(sum(hessian2)))
successes <- paste0("n = ",as.character(successes[,2]))

jpeg('FAA_RE.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=FAA_RE)) + geom_violin() + theme_classic()+ 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    xlab("Estimation Model") + ylab ("F age 8 Relative Error") +
    stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
    annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(FAA_RE,na.rm=T)))[,2]*1.1,label=successes) +
    geom_hline(yintercept=0, color = "red", size=1))
dev.off()
jpeg('R_RE.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=NAA_RE)) + geom_violin() + theme_classic()+ 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    xlab("Estimation Model") + ylab ("Recruitment Relative Error") +
      stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(NAA_RE,na.rm=T)))[,2]*1.1,label=successes) +
    geom_hline(yintercept=0, color = "red", size=1))
dev.off()
jpeg('SSB_RE.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=SSB_RE)) + geom_violin() + theme_classic()+ 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    xlab("Estimation Model") + ylab ("SSB Relative Error") +
      stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(SSB_RE,na.rm=T)))[,2]*1.1,label=successes) +
    geom_hline(yintercept=0, color = "red", size=1))
dev.off()

jpeg('mohnsSSB.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=mohnsSSB)) + geom_violin() + theme_classic()+
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
      stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(mohnsSSB,na.rm=T)))[,2]*1.1,label=successes) +
    xlab("Estimation Model") + ylab ("Mohns Rho SSB"))
dev.off()
jpeg('mohnsFbar.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=mohnsFbar)) + geom_violin() + theme_classic()+
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(mohnsFbar,na.rm=T)))[,2]*1.1,label=successes) +
    xlab("Estimation Model") + ylab ("Mohns Rho Fbar"))
dev.off()
jpeg('mohnsR.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=mohnsR)) + geom_violin() + theme_classic()+
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(mohnsR,na.rm=T)))[,2]*1.1,label=successes) +
    xlab("Estimation Model") + ylab ("Mohns Rho Recruitment"))
dev.off()

jpeg('Rec_Var.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=Rec_Var)) + geom_violin() + theme_classic()+ 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    xlab("Estimation Model") + ylab ("Recruitment Variance") +
      stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(Rec_Var,na.rm=T)))[,2]*1.1,label=successes) +
    geom_hline(yintercept=exp(mod.OM$env$last.par.best[which(names(mod.OM$env$last.par.best)=='log_NAA_sigma')][1]), color = "red", size=1))
dev.off()
jpeg('Sur_Var.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=Sur_Var)) + geom_violin() + theme_classic()+ 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    xlab("Estimation Model") + ylab ("Survival Variance") +
      stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(Surv_Var,na.rm=T)))[,2]*1.1,label=successes) +
    geom_hline(yintercept=exp(mod.OM$env$last.par.best[which(names(mod.OM$env$last.par.best)=='log_NAA_sigma')][2]), color = "red", size=1))
dev.off()
jpeg('Sel_Var.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=Sel_Var)) + geom_violin() + theme_classic()+ 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    xlab("Estimation Model") + ylab ("Selectivity Variance") +
      stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(Sel_Var,na.rm=T)))[,2]*1.1,label=successes) +
    geom_hline(yintercept=exp(mod.OM$env$last.par.best[which(names(mod.OM$env$last.par.best)=='sel_repars')][1]), color = "red", size=1))
dev.off()
jpeg('Sur_Cor.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=Sur_Cor)) + geom_violin() + theme_classic()+ 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    xlab("Estimation Model") + ylab ("Survival Correlation") +
      stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(Sur_Cor,na.rm=T)))[,2]*1.1,label=successes) +
    geom_hline(yintercept=exp(mod.OM$env$last.par.best[which(names(mod.OM$env$last.par.best)=='trans_NAA_rho')])/(1+exp(mod.OM$env$last.par.best[which(names(mod.OM$env$last.par.best)=='trans_NAA_rho')])), color = "red", size=1))
dev.off()
jpeg('Sel_Cor.jpg',width=800,height=400)
  plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=Sel_Cor)) + geom_violin() + theme_classic()+ 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
    xlab("Estimation Model") + ylab ("Selectivity Correlation") +
      stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
      annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(Sel_Cor,na.rm=T)))[,2]*1.1,label=successes) +
    geom_hline(yintercept=exp(mod.OM$env$last.par.best[which(names(mod.OM$env$last.par.best)=='sel_repars')][2])/(1+exp(mod.OM$env$last.par.best[which(names(mod.OM$env$last.par.best)=='sel_repars')][2])), color = "red", size=1))
dev.off()
jpeg('SSB_RE.jpg',width=800,height=400)
plot(ggplot(data=simsRE,aes(x=as.factor(scenarioname),y=SSB_RE)) + geom_violin() + theme_classic()+ 
       theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
       xlab("Estimation Model") + ylab ("SSB Relative Error") +
       stat_summary(fun.y=mean, geom="point", size=2, color="red") + stat_summary(fun.y=median, geom="point", size=2, color="black") +
       annotate("text",x=c(1:length(successes)),y=data.frame(simsRE %>% group_by(scenarioname) %>% summarize(max(SSB_RE,na.rm=T)))[,2]*1.1,label=successes) +
       geom_hline(yintercept=0, color = "red", size=1))
dev.off()


##################################### EPILOGUE/JUNK CODE #####################################
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
# TestName <- "Test8_finalmodel"
# if(!dir.exists(TestName)) dir.create(TestName)
# setwd(file.path(getwd(),TestName))
# 
# # Read in ASAP model fitted values
# ASAP <- read_asap3_fit(wd="C:/Users/Emily/Documents/WoodsHoleSabbotical/GOMHaddockWHAM/GOMHaddockASAP",asap.name="GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB")
# 
# # agecomp <- c('dir-mult')
# # recruitment<-c(2)
# # NAA_re_sigma<-c('rec+1')
# NAA_re_cor<-c('ar1_y','ar1_a','2dar1')
# selectivity_re<-c('iid','ar1','ar1_y','2dar1')
# # ESSoptions<-c(500)
# 
# # Set up a data frame for information of model options:
# df.mods <- as.data.frame(tidyr::crossing(NAA_re_cor,selectivity_re))
# 
# n.mods <- dim(df.mods)[1]
# df.mods$Model <- paste0("m",1:n.mods)
# df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
# 
# for(m in 1:n.mods){
#   # input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment")
#   input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
#                               selectivity=list(model=rep("age-specific",3),
#                                                re=c(rep(df.mods$selectivity_re[m],1),rep("none",2)),
#                                                initial_pars <- list(c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
#                                                                     c(0.2,0.4,0.8,1,1,1,1,1,1),
#                                                                     c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
#                                                fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
#                               NAA_re=list(sigma='rec+1',cor=df.mods$NAA_re_cor[m]),
#                               age_comp='multinomial')
#   
#   input$par$log_F1 <- log(3)
#   input$map$log_F1 <- as.factor(NA)
#   input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))
#   
#   # Mapping the estimation of numbers at age in the first year:
#   # input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
#   input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
#   # input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,7,8),nrow=1,ncol=9)) # Fix ages 8
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