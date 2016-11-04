# ------------------------------------------------------------------------------------
# This file compiles the steady-state approximation for an accelerated model spinup
# at individual points (this will need to be modified to work efficiently with spatially 
# files)
#
# References: 
#   1. Xia, J.Y., Y.Q. Luo, Y.-P. Wang, E.S. Weng, and O. Hararuk. 2012. A semi-analytical 
#      solution to accelerate spin-up of a coupled carbon and nitrogen land model to 
#      steady state. Geoscientific Model Development 5:1259-1271.
#
#   2. Xia, J., Y. Luo, Y.-P. Wang, and O. Hararuk. 2013. Traceable components of terrestrial 
#      carbon storage capacity in biogeochemical models.  Global Change Biology 19:2104-2116
#
#
# Original ED SAS solution Script at PalEON modeling HIPS sites:
# Jaclyn Hatala Matthes, 2/18/14
# jaclyn.hatala.matthes@gmail.com
#
# Modifications for greater site flexibility & updated ED
# Christine Rollinson, Aug 2015
# crollinson@gmail.com
#
# Adaptation for regional-scale runs (single-cells run independently, but executed in batches)
# Christine Rollinson, Jan 2016
# crollinson@gmail.com
# 
# Adaptation to generate an ensemble of initial conditions
# Christine Rollinson, November 2016
# crollinson@gmail.com
# ------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------
# NOTES ON THE SAS SPINUP:
# ------------------------------------------------------------------------------------
# The SAS (semi-analytical solution) should be perfomed on ED runs 
#          *******WITH DISTURBANCE TURNED OFF*******
# Turning off the disturbance (both treefall & fire) means the model will run with a 
# single patch AND we have a robust patch saying what a theoretical old growth looks like
#
# Abreviations:
# FSC = Fast soil C 
# SSC = Structural soil C
# SSL = structural soil L
# MSN = Mineralized Soil N
# FSN = Fast soil N
#
# Note: there are a number of constants and parameters that should be 
#       pulled from the code and/or XML
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Setting things up to run equations, etc
# ------------------------------------------------------------------------------------
#---------------------------------------
#Load libraries
#---------------------------------------
library(chron)
library(ncdf4)
library(colorspace)
#---------------------------------------

#---------------------------------------
# 1. Define File Structures & steps
# Additional fixed constants and file paths that don't depend on the site
#---------------------------------------
# ----------
# 1.a. Directory file structure
# ----------
#Setup analysis file structure
in.base  <- "~/Desktop/init_test_data/lat42.75lon-77.25/"

path.ED2IN <- file.path(in.base, "ED2IN")
path.histo <- file.path(in.base, "histo")
path.analy <- file.path(in.base, "analy")
site.ID    <- "HARVARD"

path.out <- file.path(in.base, "SAS_out")
if(!dir.exists(path.out)) dir.create(path.out)
# ----------

# ----------
# 1.b. List of files
# ----------
ann.files  <- dir(path.analy, "-Y-") #yearly files only  
# ----------

#---------------------------------------

#---------------------------------------
# 2. Set up constants & parameters taken from ED structure
# These should be constant across all sites for the regional runs
#---------------------------------------
# ----------
# 2.a. SAS conditions & paramters
# ----------
# Basic Settings
blckyr  <- 50 #number of years to chunk data by
dpm     <- c(31,28,31,30,31,30,31,31,30,31,30,31) # days per month cheat list
sufx    <- "g01.h5"
month.begin <- 1
month.end   <- 12
restart.year <- 1850 # Year the runs using the .css & .pss files will start in 

# Settings that can be perturbed to help generate an ensemble of veg & carbon initial conditions
# *** PLACE TO GET VALUES TO GENERATE ENSEMBLE ****
disturb <- 0.005 # the treefall disturbance rate you will prescribe in the actual runs (or close to it)

# These are calculated from the output, but could be tweaked to get the distribution of soil carb
# # soil_tempk <- 0:35+273
# # rel_soil_moist <- seq(0, 1, by=0.1)
# ----------

# ----------
# 2.b. Parameters to extract from ED2IN
# ----------
site.lat      <- 42.75            # NL%POI_LAT
site.lon      <- -72.25           # NL%POI_LON
spin.start    <- 1850             # NL%IYEARA
spin.end      <- 2850             # NL%IYEARZ
yrs.met       <- 20               # NL%METCYC1 - NL$METCYF +1
pft           <- c(5,6,8,9,10,11) # NL%INCLUDE_THESE_PFT
decomp_scheme <- 2                # NL%DECOMP_SCHEME
xml.path      <- file.path(in.base, "PalEON_Phase2.v1.xml") # NL%IEDCNFGF
# ----------

# ----------
# 2.c. Paramters to extract from XML
# Note: The names are as they're specified in the xml
# ----------
rh_decay_low   <- 0.14 
rh_decay_high  <- 0.60
rh_low_temp    <- 291.15
rh_high_temp   <- 313.15
rh_decay_dry   <- 12.0
rh_decay_wet   <- 36.0
rh_dry_smoist  <- 0.20
rh_wet_smoist  <- 0.98
# ----------

# ----------
# 2.d. Parameters defined in ED2 source code
# ----------
kh_active_depth = -0.2

# # Note: these values should be per YEAR (i.e. do not divide these by yr_day like happens in ed_params.f90)
decay_rate_fsc <- 11
decay_rate_stsc <- 4.5
decay_rate_ssc <- 0.2 # New ED2 Defaults; old would do weird things
# decay_rate_ssc <- 100.2  # Old ED2

if(decomp_scheme!=2){
  # These are used if you're not working with decomp_scheme = 2 (i.e. LLoyd-Taylor)
  resp_opt_water            <- 0.8938
  resp_water_below_opt      <- 5.0786
  resp_water_above_opt	  <- 4.5139
  #resp_temperature_increase <- 0.0757 # Jackie had greatly modified this value
  resp_temperature_increase <- 0.23503 # Jackie
}

if(decomp_scheme==2){
  # Constants used in DECOMP_SCHEME=2
  Lc                        <- 0.049787 # in soil_respiration; =exp(-3.0)
  c2n_slow                  <- 10.0
  c2n_structural            <- 150.0
  r_stsc                    <- 0.3
}
# ----------
#---------------------------------------


# ------------------------------------------------------------------------------------
# Running the SAS Solution
# ------------------------------------------------------------------------------------
#---------------------------------------
# 3. Setting up some specifics that vary by site (like soil depth)
#---------------------------------------
# ----------
# 3.a. Derived values 
# ----------
yeara  <- spin.start+1 #first full year
yearz  <- spin.end-1 #last full year
yrs    <- seq(yeara, yearz, by=blckyr) # The years we're going to use as time steps for the demography
nsteps <- length(yrs) # The number of blocks = the number steps we'll have
# ----------

# ----------
# 3.b. Set up patch matrix
# ----------
#create an emtpy storage for the patch info
pss.big <- matrix(nrow=length(yrs),ncol=14) # save every X yrs according to chunks specified above
colnames(pss.big) <- c("site","year","patch","dst","age","area","water","fsc","stsc","stsl",
                       "ssc","psc","msn","fsn")

# Age-area matrix for size distribution
# Calculate area distribution based on geometric decay based loosely on your disturbance rates
# Note: Current configuration throws all leftover area in to the oldest age bin because it assumes 
#       you're going to do a transient run with disturbance on to equilibrate structure.
#       If you're not planning on doing this run, you may want to go through and evenly distribute
#       everything to the geometric distribution (what Jackie's original SAS looked like)
stand.age <- seq(yrs[1]-yeara,nrow(pss.big)*blckyr,by=blckyr)
area.dist <- vector(length=nrow(pss.big))
area.dist[1] <- sum(dgeom(0:(stand.age[2]-1), disturb))
for(i in 2:(length(area.dist)-1)){
  area.dist[i] <- sum(dgeom((stand.age[i]):(stand.age[i+1]-1),disturb))
}
area.dist[length(area.dist)] <- 1 - sum(area.dist[1:(length(area.dist)-1)])
pss.big[,"area"] <- area.dist
# ----------


# ----------
# 3.c. Soil layers
# ----------
# Need to get the layers being used for calculating temp & moist
# Note: In ED there's a pain in the butt way of doing this with the energy, but we're going to approximate
# slz  <- c(-5.50, -4.50, -2.17, -1.50, -1.10, -0.80, -0.60, -0.45, -0.30, -0.20, -0.12, -0.06)
# dslz <- c(1.00,   2.33,  0.67,  0.40,  0.30,  0.20,  0.15,  0.15,  0.10,  0.08,  0.06,  0.06)
nc.temp <- nc_open(file.path(dat.dir, ann.files[1]))
slz <- ncvar_get(nc.temp, "SLZ")
dslz <- vector(length=length(slz))
dslz[length(dslz)] <- 0-slz[length(dslz)]

for(i in 1:(length(dslz)-1)){
  dslz[i] <- slz[i+1] - slz[i]    
}

# nsoil=which(slz>= kh_active_depth)
nsoil=length(slz)
# ----------

#---------------------------------------

# First loop over analy files (faster than histo) to aggregate initial 
# 	 .css and .pss files for each site

#---------------------------------------
# 4 Finding the mean soil temp & moisture
# *** PLACE TO GET VALUES TO GENERATE ENSEMBLE ****
# NOTE:  I've been plyaing around with finding the best temp & soil moisture to initialize things
#        with; if using the means from the spin met cycle work best, insert them here
#---------------------------------------
tempk.air <- tempk.soil <- moist.soil <- moist.soil.mx <- vector()
for(y in yrs){
  air.temp.tmp <- soil.temp.tmp <- soil.moist.tmp <- soil.mmax.tmp <- vector()
  ind <- which(yrs == y)
  for(m in month.begin:month.end){
    #Make the file name. 
    year.now  <-sprintf("%4.4i",y)
    month.now <- sprintf("%2.2i",m)
    day.now   <- sprintf("%2.2i",0)
    hour.now  <- sprintf("%6.6i",0)
    
    file.now    <- paste(sites[s],"-E-",year.now,"-",month.now,"-",day.now,"-"
                         ,hour.now,"-",sufx,sep="")
    
    # cat(" - Reading file :",file.now,"...","\n")
    now <- nc_open(paste(dat.dir,file.now,sep=""))
    
    air.temp.tmp  [m] <- ncvar_get(now, "MMEAN_ATM_TEMP_PY")
    soil.temp.tmp [m] <- sum(ncvar_get(now, "MMEAN_SOIL_TEMP_PY")[nsoil]*dslz[nsoil]/sum(dslz[nsoil]))
    soil.moist.tmp[m] <- sum(ncvar_get(now, "MMEAN_SOIL_WATER_PY")[nsoil]*dslz[nsoil]/sum(dslz[nsoil]))
    soil.mmax.tmp [m] <- max(ncvar_get(now, "MMEAN_SOIL_WATER_PY"))
    
    nc_close(now)
  }
  # Finding yearly means
  tempk.air    [ind] <- mean(air.temp.tmp)
  tempk.soil   [ind] <- max(soil.temp.tmp) # means are typically making temps too low
  moist.soil   [ind] <- mean(soil.moist.tmp)
  moist.soil.mx[ind] <- max(soil.mmax.tmp)
}
# adjusting the soil moist to be relative based on the max observed
moist.soil <- moist.soil/max(moist.soil.mx)

soil_tempk     <- mean(tempk.soil)
# rel_soil_moist <- mean(moist.soil)+.2
rel_soil_moist <- 0.5
#---------------------------------------

#---------------------------------------  
# 5. Cohort Extraction Loop
# This loop does the following:
#  -- Extract cohort info from each age slice from *annual* *analy* files (these are annual means)
#  -- Write cohort info to the .css file as a new patch for each age slice
#  -- Dummy extractions of patch-level variables; all of the important variables here are place holders
#---------------------------------------  
for (y in yrs){
  cat(" - Reading file :",ann.files[y-yeara+1],"...","\n")
  now <- nc_open(paste(dat.dir,ann.files[y-yeara+1],sep=""))
  ind <- which(yrs == y)
  
  #Grab variable to see how many cohorts there are
  ipft      <- ncvar_get(now,'PFT')
  
  #---------------------------------------
  # organize into .css variables (Cohorts)
  # Note: all cohorts from a time slice are assigned to a single patch representing a stand of age X
  #---------------------------------------
  css.tmp <- matrix(nrow=length(ipft),ncol=10)
  css.tmp[,1] <- rep(restart.year,length(ipft))
  css.tmp[,2] <- rep(floor((y-yeara)/blckyr)+1,length(ipft))
  css.tmp[,3] <- 1:length(ipft)
  css.tmp[,4] <- ncvar_get(now,'DBH')
  css.tmp[,5] <- ncvar_get(now,'HITE')
  css.tmp[,6] <- ipft
  css.tmp[,7] <- ncvar_get(now,'NPLANT')
  css.tmp[,8] <- ncvar_get(now,'BDEAD')
  css.tmp[,9] <- ncvar_get(now,'BALIVE')
  css.tmp[,10] <- rep(-999,length(ipft))
  colnames(css.tmp) <- c("year","patch","cohort","dbh","ht","pft","n","bdead","balive","Avgrg")
  
  #save big .css matrix
  if(y==yrs[1]){
    css.big <- css.tmp
  } else{
    css.big <- rbind(css.big,css.tmp)
  }
  #---------------------------------------
  
  
  #---------------------------------------
  # save .pss variables (Patches)
  # NOTE: patch AREA needs to be adjusted to be equal to the probability of a stand of age x on the landscape
  #---------------------------------------
  pss.big[ind,1]  <- 1
  pss.big[ind,2]  <- restart.year
  pss.big[ind,3]  <- floor((y-yeara)/blckyr)+1
  pss.big[ind,4]  <- 1
  pss.big[ind,5]  <- y-yeara
  # Note: the following are just place holders that will be overwritten post-SAS
  # pss.big[ind,6]  <- ncvar_get(now,"AREA") # This is done earlier
  pss.big[ind,7]  <- 0.5 
  pss.big[ind,8]  <- ncvar_get(now,"FAST_SOIL_C")
  pss.big[ind,9]  <- ncvar_get(now,"STRUCTURAL_SOIL_C")
  pss.big[ind,10] <- ncvar_get(now,"STRUCTURAL_SOIL_L")
  pss.big[ind,11] <- ncvar_get(now,"SLOW_SOIL_C")
  pss.big[ind,12] <- 0
  pss.big[ind,13] <- ncvar_get(now,"MINERALIZED_SOIL_N")
  pss.big[ind,14] <- ncvar_get(now,"FAST_SOIL_N")
  
  nc_close(now)
}
#---------------------------------------  

#---------------------------------------  
# 6. Path Extraction Loop
# This loop does the following:
#  -- Extract age slice (new patch) soil carbon conditions from *monthly* *histo* files
#       -- Note: this is done because most of the necessary inputs for SAS are instantaneous values that 
#                are not currently tracked in analy files, let alone annual analy files; this could 
#                theoretically change in the future
#       -- Monthly data is then aggregated to a yearly value: sum for carbon inputs; mean for temp/moist 
#          (if not calculated above)
#---------------------------------------
pss.big <- pss.big[complete.cases(pss.big),]

# some empty vectors for storage etc
fsc_in_y <- ssc_in_y <- ssl_in_y <- fsn_in_y <- pln_up_y <- vector()
fsc_in_m <- ssc_in_m <- ssl_in_m <- fsn_in_m <- pln_up_m <-  vector()

# switch to the histo directory
dat.dir    <- paste(in.base,sites[s],"/histo/",sep="")
mon.files  <- dir(dat.dir, "-S-") # monthly files only  

for (y in yrs){      
  for(m in month.begin:month.end){
    #Make the file name. 
    year.now  <-sprintf("%4.4i",y)
    month.now <- sprintf("%2.2i",m)
    day.now   <- sprintf("%2.2i",1)
    hour.now  <- sprintf("%6.6i",0)
    
    dat.dir     <- paste(in.base,sites[s],"/histo/",sep="")
    file.now    <- paste(sites[s],"-S-",year.now,"-",month.now,"-",day.now,"-"
                         ,hour.now,"-",sufx,sep="")
    
    cat(" - Reading file :",file.now,"...","\n")
    now <- nc_open(paste(dat.dir,file.now,sep=""))
    
    # Note: we have to convert the daily value for 1 month by days per month to get a monthly estimate
    fsc_in_m[m-month.begin+1] <- ncvar_get(now,"FSC_IN")*dpm[m] #kg/(m2*day) --> kg/(m2*month)
    ssc_in_m[m-month.begin+1] <- ncvar_get(now,"SSC_IN")*dpm[m]
    ssl_in_m[m-month.begin+1] <- ncvar_get(now,"SSL_IN")*dpm[m]
    fsn_in_m[m-month.begin+1] <- ncvar_get(now,"FSN_IN")*dpm[m]
    pln_up_m[m-month.begin+1] <- ncvar_get(now,"TOTAL_PLANT_NITROGEN_UPTAKE")*dpm[m]
    # ssc_in_m[m-month.begin+1] <- ncvar_get(now,"SSC_IN")*dpm[m]
    
    # # NOTE: the following lines shoudl get removed if using 20-year means
    # soil_tempk_m[m-month.begin+1] <- ncvar_get(now,"SOIL_TEMPK_PA")[nsoil] # Surface soil temp
    # swc_max_m[m-month.begin+1] <- max(ncvar_get(now,"SOIL_WATER_PA")) # max soil moist to avoid digging through water capacity stuff
    # swc_m[m-month.begin+1] <- ncvar_get(now,"SOIL_WATER_PA")[nsoil] #Surface soil moist
    
    nc_close(now)
  }
  
  # Find which patch we're working in
  ind <- (y-yeara)/blckyr + 1
  
  # Sum monthly values to get a total estimated carbon input
  fsc_in_y[ind] <- sum(fsc_in_m,na.rm=TRUE)
  ssc_in_y[ind] <- sum(ssc_in_m,na.rm=TRUE)
  ssl_in_y[ind] <- sum(ssl_in_m,na.rm=TRUE)
  fsn_in_y[ind] <- sum(fsn_in_m,na.rm=TRUE)
  pln_up_y[ind] <- sum(pln_up_m,na.rm=TRUE)
  
  # # Soil temp & moisture here should get deleted if using the 20-year means
  # soil_tempk_y[ind] <- mean(soil_tempk_m,na.rm=TRUE) 
  # swc_y[ind] <- mean(swc_m,na.rm=TRUE)/max(swc_max_m,na.rm=TRUE) 
}
#---------------------------------------

#---------------------------------------  
# 7. Calculate steady-state soil pools!
#
# These are the equations from soil_respiration.f90 -- if this module has changed, these need 
# Note: We ignore the unit conversions here because we're now we're working with the yearly 
#       sum so that we end up with straight kgC/m2
# fast_C_loss <- kgCday_2_umols * A_decomp * decay_rate_fsc * fast_soil_C
# struc_C_loss <- kgCday_2_umols * A_decomp * Lc * decay_rate_stsc * struct_soil_C * f_decomp
# slow_C_loss <- kcCday_2_umols * A_decomp * decay_rate_ssc * slow_soil_C
#---------------------------------------

# -----------------------
# Calculate the annual carbon loss if things are stable
# -----------
fsc_loss <- decay_rate_fsc
ssc_loss <- decay_rate_ssc
ssl_loss <- decay_rate_stsc
# -----------


# *************************************
# 7.a. Calculate A_decomp according to your DECOMP_SCPEME
# A_decomp <- temperature_limitation * water_limitation # aka het_resp_weight
# *************************************
# ========================
# 7.a.1. Temperature Limitation 
# ========================
# soil_tempk <- sum(soil_tempo_y*area.dist)
if(decomp_scheme==0){
  temperature_limitation = exp(resp_temperature_increase * (soil_tempk-318.15))
}
if(decomp_scheme==1){
  temperature_limitation = resp_temperature_increase * exp(308.56 * (1./56.02 - 1./(soil_tempk-227.15))) # this is LloydTaylor
}
if(decomp_scheme==2){
  # Low Temp Limitation
  lnexplow <- rh_decay_low * (rh_low_temp - soil_tempk)
  tlow_fun <- 1 + exp(lnexplow)
  
  # High Temp Limitation
  lnexphigh <- rh_decay_high*(soil_tempk - rh_high_temp)
  thigh_fun <- 1 + exp(lnexphigh)
  
  temperature_limitation <- 1/(tlow_fun*thigh_fun)
}
# ========================

# ========================
# 7.a.2 Moisture Limitation 
# ========================
# rel_soil_moist <- sum(swc_y*area.dist)
# # -----------
# # Moyano (not implemented; what Jackie used)
# # -----------
# water_limitation <- rel_soil_moist*4.0893 + rel_soil_moist^2*-3.1681 - 0.3195897 # This is Jackie's Moyano et al equation
# # -----------
if(decomp_scheme %in% c(0,1)){
  water_limitation <- exp((rel_soil_moist - resp_opt_water) * resp_water_below_opt)
}

if(decomp_scheme == 2){
  # Dry soil Limitation
  lnexpdry <- rh_decay_dry * (rh_dry_smoist - rel_soil_moist)
  smdry_fun <- 1+exp(lnexpdry)
  
  # Wet Soil limitation
  lnexpwet <- rh_decay_wet * (rel_soil_moist - rh_wet_smoist)
  smwet_fun <- 1+exp(lnexpwet)
  
  water_limitation <- 1/(smdry_fun * smwet_fun)
}
# ========================

A_decomp <- temperature_limitation * water_limitation # aka het_resp_weight
# *************************************

# *************************************
# 7.b Calculate the steady-state pools
# NOTE: Current implementation weights carbon input by patch size rather than using the 
#       carbon balance from the oldest state (as was the first implementation)
# *************************************
# -------------------
# 7.b.1 Do the carbon and fast nitrogen pools
# -------------------
fsc_ss <- fsc_in_y[length(fsc_in_y)]/(fsc_loss * A_decomp)
ssl_ss <- ssl_in_y[length(ssl_in_y)]/(ssl_loss * A_decomp * Lc) # Structural soil C
ssc_ss <- ((ssl_loss * A_decomp * Lc * ssl_ss)*(1 - r_stsc))/(ssc_loss * A_decomp )
fsn_ss <- fsn_in_y[length(fsn_in_y)]/(fsc_loss * A_decomp)
# -------------------

# -------------------
# 7.b.2 Do the mineralized nitrogen calculation
# -------------------
#ED2: csite%mineralized_N_loss  = csite%total_plant_nitrogen_uptake(ipa)             
# + csite%today_Af_decomp(ipa) * Lc * K1 * csite%structural_soil_C(ipa)                     
# * ( (1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)
msn_loss <- pln_up_y[length(pln_up_y)] + 
  A_decomp*Lc*ssl_loss*ssl_in_y[length(ssl_in_y)]*
  ((1.0-r_stsc)/c2n_slow - 1.0/c2n_structural)

#fast_N_loss + slow_C_loss/c2n_slow
msn_med  <- fsc_loss*A_decomp*fsn_in_y[length(fsn_in_y)]+ (ssc_loss * A_decomp)/c2n_slow 

msn_ss   <- msn_med/msn_loss
# -------------------
# *************************************

# *************************************
# 7.c. Replace dummy values in patch matrix with the steady state calculations
# *************************************
# Figure out what steady state value index we shoudl use 
# Note: In the current implementaiton this should be 1 because we did the weighted averaging up front, 
#       but if something went wrong and dimensions are off, use this to pick the last (etc)
p.use <- 1

# write the values to file
pss.big[,3]  <- 1:nrow(pss.big)
pss.big[,6]  <- area.dist
pss.big[,8]  <- rep(fsc_ss[p.use],nrow(pss.big)) # fsc
pss.big[,9]  <- rep(ssl_ss[p.use],nrow(pss.big)) # stsc
pss.big[,10] <- rep(ssl_ss[p.use],nrow(pss.big)) # stsl (not used)
pss.big[,11] <- rep(ssc_ss[p.use],nrow(pss.big)) # ssc
pss.big[,13] <- rep(msn_ss[p.use],nrow(pss.big)) # msn
pss.big[,14] <- rep(fsn_ss[p.use],nrow(pss.big)) # fsn
# *************************************
#---------------------------------------

#---------------------------------------
# 8. Write everything to file!!
# NOTE: ED is silly and requires the site lat/lon to be part of the init file name
#---------------------------------------
sas.name <- paste0(site.ID, "_lat", site.lat, "lon", site.lon) 

write.table(css.big,file=paste(out,sas.name,".css",sep=""),row.names=FALSE,append=FALSE,
            col.names=TRUE,quote=FALSE)

write.table(pss.big,file=paste(out,sas.name,".pss",sep=""),row.names=FALSE,append=FALSE,
            col.names=TRUE,quote=FALSE)
#---------------------------------------
# -------------------------------------
