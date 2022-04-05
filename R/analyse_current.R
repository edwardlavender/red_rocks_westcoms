################################
################################
#### analyse_current.R

#### This script: 
# 1) Analyses current speeds in the Red Rocks and Longay MPA & egg laying sites

#### Steps preceding this script: 
# 1) Define global parameters     (define_global_param.R)
# 2) Obtain and process raw data  (process_data_raw.R)


################################
################################
#### Set up

#### Wipe workspace and source essential packages and variables
source("./R/define_global_param.R")

#### Load data 
mesh_mpa  <- readRDS("./data/spatial/mesh/mesh_around_elements_in_mpa.rds")
mesh_reg  <- readRDS("./data/spatial/mesh/mesh_around_elements_in_region.rds")
eggs      <- readRDS("./data/eggs/eggs.rds")


################################
################################
#### Derive current speed estimates 

#### Define analysis
# mpa - extracts conditions in cells in MPA
# eggs - extracts conditions with/without eggs (across a wider region)
analysis <- "eggs" # mpa 
stopifnot(analysis %in% c("mpa", "eggs"))
if(analysis == "mpa"){
  mesh <- mesh_mpa
} else if(analysis == "eggs"){
  mesh <- mesh_reg
}
raster::plot(mesh)

#### Show that current speeds have been correctly estimated 
check <- FALSE
if(check){
  # Load files 
  uv <- R.matlab::readMat(paste0(wc_root, "uvelocity/160301.mat"))$data
  vv <- R.matlab::readMat(paste0(wc_root, "vvelocity/160301.mat"))$data
  cs <- readRDS(paste0(wc_root, "current_speed/160301.rds"))
  # Define random index of bottom velocities to check 
  ind <- cbind(1:20, 10, 1:100)
  identical(sqrt(uv[ind]^2 + vv[ind]^2), cs[ind]) # TRUE
}

#### Define data for extraction 
dates       <- seq(as.Date("2016-03-01"), as.Date("2017-02-28"), 1)
date_names  <- date_name(dates)
hours       <- 0:23
# Define layers (1 for surface or 10 for near seabed)
layers      <- 10 # 1
if(analysis == "mpa"){ 
  mesh_IDs    <- as.integer(as.character(mesh$ID))
} else if(analysis == "eggs"){
  eggs$cell <-    
    find_cells(eggs$lat, 
               eggs$long, 
               proj = raster::crs(mesh), 
               mesh = mesh,
               return = 4)
  eggs$cell <- as.integer(as.character(eggs$cell))
  table(eggs$cell[eggs$present == 1] %in% mesh_mpa$ID)
  table(eggs$cell[eggs$present == 0] %in% mesh_mpa$ID)
  mesh_IDs <- unique(eggs$cell)
}
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
  if(analysis == "mpa"){
    if(layers == 10){
      wc_bot <- wc
      saveRDS(wc, "./data/wc/current/wc_bot.rds")
    } else if(layers == 1){
      wc_sur <- wc
      saveRDS(wc, "./data/wc/current/wc_sur.rds")
    }
  } else if(analysis == "eggs"){
    if(layers == 10){
      wc_bot <- wc
      saveRDS(wc, "./data/wc/current/wc_bot_eggs.rds")
    } else if(layers == 1){
      wc_sur <- wc
      saveRDS(wc, "./data/wc/current/wc_sur_eggs.rds")
    }
  }
} else {
  if(analysis == "mpa"){
    wc_bot <- readRDS("./data/wc/current/wc_bot.rds")
    # wc_sur <- readRDS("./data/wc/current/wc_sur.rds")
  } else if(analysis == "eggs"){
    wc_bot <- readRDS("./data/wc/current/wc_bot_eggs.rds")
    # wc_sur <- readRDS("./data/wc/current/wc_sur_eggs.rds")
  }
}

#### Estimate maximum current speeds 
utils.add::basic_stats(wc_bot$wc)
if(analysis == "eggs"){
  cells_with_eggs <- unique(eggs$cell[eggs$present == 1])
  utils.add::basic_stats(wc_bot$wc[wc_bot$mesh_ID %in% cells_with_eggs])
}

#### Examine annual trends
## Define timestamps
wc$date           <- date_name(wc$date_name, define = "date")
wc$hour_char      <- as.character(wc$hour)
pos               <- nchar(wc$hour_char) == 1
wc$hour_char[pos] <- paste0("0", wc$hour_char[pos])
wc$timestamp      <- fasttime::fastPOSIXct(paste0(wc$date, " ", wc$hour_char, ":00:00"), tz = "UTC")
## Visualise trends 
wc_by_hour <- 
  wc %>%
  dplyr::group_by(timestamp) %>%
  dplyr::summarise(wc = median(wc))
pretty_plot(wc_by_hour$timestamp, wc_by_hour$wc, type = "l")



#### End of code. 
################################
################################