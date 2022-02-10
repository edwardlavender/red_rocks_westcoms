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
uv <- R.matlab::readMat(paste0(wc_root, "uvelocity/160301.mat"))$data
vv <- R.matlab::readMat(paste0(wc_root, "vvelocity/160301.mat"))$data
cs <- readRDS(paste0(wc_root, "current_speed/160301.rds"))
sqrt(uv[1, 1, 1]^2 + vv[1, 1, 1]^2)
cs[1, 1, 1]

#### Define data for extraction 
dates       <- seq(as.Date("2016-03-01"), as.Date("2017-02-28"), 1)
date_names  <- date_name(dates)
hours       <- 0:23
layers      <- 10
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
  saveRDS(wc, "./data/wc/current/wc.rds")
} else {
  wc <- readRDS("./data/wc/current/wc.rds")
}

#### Estimate maximum current speeds 
hist(wc$wc)
max(wc$wc)
# 0.1974771


#### End of code. 
################################
################################