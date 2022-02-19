################################
################################
#### analyse_current.R

#### This script: 
# 1) Analyses current speeds in the Red Rocks and Longay MPA
# ... A maximum predicted current speed for the area over a one
# ... year period is derived. 

#### Steps preceding this script: 
# 1) Define global parameters     (define_global_param.R)
# 2) Obtain and process raw data  (process_data_raw.R)


################################
################################
#### Set up

#### Wipe workspace and source essential packages and variables
source("./R/define_global_param.R")

#### Load data 
mesh  <- readRDS("./data/spatial/mesh/mesh_around_elements_in_mpa.rds")


################################
################################
#### Derive current speed estimates 

#### Show that current speeds have been correctly estimated 
# Load files 
uv <- R.matlab::readMat(paste0(wc_root, "uvelocity/160301.mat"))$data
vv <- R.matlab::readMat(paste0(wc_root, "vvelocity/160301.mat"))$data
cs <- readRDS(paste0(wc_root, "current_speed/160301.rds"))
# Define random index of bottom velocities to check 
ind <- cbind(1:20, 10, 1:100)
identical(sqrt(uv[ind]^2 + vv[ind]^2), cs[ind]) # TRUE

#### Define data for extraction 
dates       <- seq(as.Date("2016-03-01"), as.Date("2017-02-28"), 1)
date_names  <- date_name(dates)
hours       <- 0:23
# Define layers (1 for surface or 10 for near seabed)
layers      <- 10 # 1
mesh_IDs    <- as.integer(as.character(mesh$ID))
wc <- 
  expand.grid(date_name = date_names, 
              hour = hours, 
              layer = layers, 
              mesh_ID = mesh_IDs) %>%
  dplyr::arrange(date_name, 
                 hour, 
                 mesh_ID)

#### Implement extraction (~49 mins)
run <- FALSE
if(run){
  t1 <- Sys.time()
  wc <- fvcom.tbx::extract(dat = wc, # [1:10, ], 
                           read_fvcom = readRDS,
                           dir2load = paste0(wc_root, "/current_speed/"), 
                           extension = ".rds", 
                           cl = parallel::makeCluster(10L)
  )
  t2 <- Sys.time()
  difftime(t2, t1)
  if(layers == 10){
    saveRDS(wc, "./data/wc/current/wc_bot.rds")
  } else if(layers == 1){
    saveRDS(wc, "./data/wc/current/wc_sur.rds")
  }
} else {
  wc_bot <- readRDS("./data/wc/current/wc_bot.rds")
  wc_sur <- readRDS("./data/wc/current/wc_sur.rds")
}

#### Estimate maximum current speeds 
hist(wc_bot$wc)
max(wc_bot$wc)
# 0.1974771
hist(wc_sur$wc)
max(wc_sur$wc)
# 0.7945567


#### End of code. 
################################
################################