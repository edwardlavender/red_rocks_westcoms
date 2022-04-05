################################
################################
#### analyse_temp.R

#### This script: 
# 1) Analyses temperatures in the Red Rocks and Longay MPA & egg-laying sites. 
# ... Across the MPA:
# ... ... The annual average and variability are explored
# ... ... Monthly averages/spatiotemporal variation are explored
# ... At egg-laying sites: 
# ... ... Temperature time series are related to trends in the MPA at large. 

#### Steps preceding this script: 
# 1) Define global parameters     (define_global_param.R)
# 2) Obtain and process raw data  (process_data_raw.R)


################################
################################
#### Set up

#### Wipe workspace and source essential packages and variables
source("./R/define_global_param.R")

#### Load data 
coast      <- readRDS("./data/spatial/coast/coast.rds")
mpa        <- readRDS("./data/spatial/mpa/mpa.rds")
mesh_mpa   <- readRDS("./data/spatial/mesh/mesh_around_nodes_in_mpa.rds")
mesh_reg   <- readRDS("./data/spatial/mesh/mesh_around_nodes_in_region.rds")
nodexy     <- readRDS("./data/spatial/mesh/mesh_nodexy.rds")
eggs       <- readRDS("./data/eggs/eggs.rds")

#### Essential packages
library(prettyGraphics)


################################
################################
#### Extract temperatures 

#### Examine example file 
# temp <- R.matlab::readMat(paste0(wc_root, "temp/160301.mat"))$data
# str(temp)

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

#### Define data for extraction 
dates       <- seq(as.Date("2016-03-01"), as.Date("2017-02-28"), 1)
date_names  <- date_name(dates)
hours       <- 0:23
layers      <- 10
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
str(wc)

#### Implement extraction (~5 mins)
run <- FALSE
if(run){
  t1 <- Sys.time()
  wc <- fvcom.tbx::extract(dat = wc, # [1:10, ], 
                           dir2load = paste0(wc_root, "/temp/"), 
                           extension = ".mat", 
                           cl = parallel::makeCluster(10L)
  )
  t2 <- Sys.time()
  difftime(t2, t1)
  if(analysis == "mpa"){
    saveRDS(wc, "./data/wc/temp/wc.rds")
  } else if(analysis == "eggs"){
    saveRDS(wc, "./data/wc/temp/wc_eggs.rds")
  }
} else {
  if(analysis == "mpa"){
    wc <- readRDS("./data/wc/temp/wc.rds")
  } else if(analysis == "eggs"){
    wc <- readRDS("./data/wc/temp/wc_eggs.rds")
  }
}


################################
################################
#### Data processing 

#### Define timestamps
# dates
wc$date           <- date_name(wc$date_name, define = "date")
# timestamps 
wc$hour_char      <- as.character(wc$hour)
pos               <- nchar(wc$hour_char) == 1
wc$hour_char[pos] <- paste0("0", wc$hour_char[pos])
wc$timestamp      <- fasttime::fastPOSIXct(paste0(wc$date, " ", wc$hour_char, ":00:00"), tz = "UTC")
# months
wc$month          <- lubridate::month(wc$date, label = TRUE)
# check
head(wc)

#### Define depths
wc$depth <- nodexy$z[match(wc$mesh_ID, nodexy$node_id)]
utils.add::basic_stats(wc$depth)
# 0.41 34.08  30.88 117.81 18.1 12.41 8.24

#### Temperature range
utils.add::basic_stats(wc$wc)


################################
################################
#### Analyse temperatures (MPA)

if(analysis == "mpa"){
  
  ################################
  #### Summary statistics 
  
  #### Annual average
  ## Bottom temps
  utils.add::basic_stats(wc$wc)
  # min  mean median   max   sd  IQR  MAD
  # 8.69 10.49  10.61 12.43 1.11 2.23 1.55
  
  #### Monthly averages
  wc_by_month <- 
    wc %>% 
    dplyr::group_by(month) %>%
    dplyr::summarise(wc_mean = mean(wc), 
                     wc_q1   = quantile(wc, 0.25), 
                     wc_q3   = quantile(wc, 0.75), 
                     wc_iqr  = IQR(wc), 
                     wc_sd   = sd(wc)) %>%
    dplyr::arrange(wc_mean)
  wc_by_month
  
  #### Monthly averages by depth bin
  wc_by_month_by_bin <- 
    wc %>% 
    dplyr::mutate(bin = cut(depth, 10)) %>%
    dplyr::group_by(bin, month) %>%
    dplyr::summarise(depth   = depth[1], 
                     wc_mean = mean(wc), 
                     # wc_q1   = quantile(wc, 0.25), 
                     # wc_q3   = quantile(wc, 0.75), 
                     wc_iqr  = IQR(wc)) %>%
    dplyr::arrange(month, wc_mean)
  wc_by_month_by_bin
  # View(wc_by_month_by_bin)
  
  
  ################################
  #### Monthly maps
  
  #### Summarise temperatures by month and mesh cell
  wc_by_month_and_cell <- 
    wc %>% 
    dplyr::mutate(mm_yy = Tools4ETS::mmyy(date)) %>%
    dplyr::group_by(mm_yy, mesh_ID) %>%
    dplyr::summarise(fvcom = mean(wc)) %>%
    dplyr::arrange(mm_yy) %>%
    dplyr::select(mm_yy = mm_yy, ID = mesh_ID, fvcom = fvcom) %>%
    data.frame()
  
  #### Define graphical parameters
  ext  <- raster::extent(mpa)
  xlim <- ext[1:2]
  ylim <- ext[3:4]
  coast_around_mpa <- raster::crop(coast, mpa)
  
  #### Make monthly maps 
  run <- FALSE
  if(run){
    png("./fig/temp_maps_by_month.png", 
        height = 12, width = 8, units = "in", res = 800)
    pp <- par(mfrow = c(5, 3), oma = c(2, 2, 2, 2), mar = c(0.5, 2, 0.5, 2))
    
    ## Plot WeStCOMS bathymetry 
    # Define mesh cell depths
    bathy <- nodexy %>% 
      dplyr::filter(node_id %in% mesh$ID) %>%
      dplyr::mutate(z = abs(z) * -1) %>%
      dplyr::select(ID = node_id, fvcom = z)
    # Plot depths 
    plot_field_2d(coastline = coast_around_mpa, 
                  mesh = mesh, 
                  data = bathy, 
                  col_fn = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues"))),
                  xlim = xlim, ylim = ylim, 
                  xlab = "", ylab = "", axes = FALSE, 
                  colour_bar_x = c(xlim[2], xlim[2] + 0.005), 
                  colour_bar_y = ylim)
    rect(xlim[1], ylim[1], xlim[2], ylim[2])
    mtext(side = 3, "Bathymetry (m)", line = -1, adj = 0.5, font = 2)
    
    ## Plot monthly temperature averages 
    # Loop over each month and make plot 
    wc_by_month_and_cell_ls <-
      split(wc_by_month_and_cell, wc_by_month_and_cell$mm_yy)
    flapper::cl_lapply(wc_by_month_and_cell_ls, function(wc_for_month){
      # wc_for_month <- wc_by_month_and_cell_ls[[1]]
      plot_field_2d(coastline = coast_around_mpa, 
                    mesh = mesh, 
                    data = wc_for_month[, c("ID", "fvcom")], 
                    col_fn = function(x) rev(grDevices::heat.colors(x)),
                    xlim = xlim, ylim = ylim, 
                    xlab = "", ylab = "", axes = FALSE, 
                    colour_bar_x = c(xlim[2], xlim[2] + 0.005),
                    colour_bar_y = ylim)
      rect(xlim[1], ylim[1], xlim[2], ylim[2])
      mtext(side = 3, wc_for_month$mm_yy[1], line = -1, adj = 0.5, font = 2)
    }) %>% invisible()
    mtext(side = 4, expression("Bottom temperature (" * degree * "C)"), line = 5)
    dev.off()
  }
  
  
  ################################
  #### Temperature time series
  
  #### Summarise mean and variability across the area through time 
  wc_by_hour <- 
    wc %>% 
    dplyr::group_by(timestamp) %>%
    dplyr::summarise(wc_mean = mean(wc), 
                     wc_q1   = quantile(wc, 0.25), 
                     wc_q3   = quantile(wc, 0.75), 
                     wc_min  = min(wc), 
                     wc_max  = max(wc), 
                     wc_depth_for_min = depth[which.min(wc)], 
                     wc_depth_for_max = depth[which.max(wc)]) %>%
    dplyr::arrange(timestamp)
  saveRDS(wc_by_hour, "./data/wc/temp/wc_by_hour.rds")
  
  # Explore variation across the area
  wc_by_hour %>% 
    dplyr::mutate(date = as.Date(timestamp), 
                  mm_yy = Tools4ETS::mmyy(date)) %>%
    dplyr::group_by(mm_yy) %>%
    dplyr::summarise(max(wc_mean) - min(wc_mean))
  
  # Explore variation with depth 
  wc_by_hour %>% 
    dplyr::mutate(date = as.Date(timestamp), 
                  mm_yy = Tools4ETS::mmyy(date)) %>%
    dplyr::group_by(mm_yy) %>%
    dplyr::summarise(min(wc_depth_for_max), 
                     max(wc_depth_for_max), 
                     min(wc_depth_for_min), 
                     max(wc_depth_for_min)
    )
  
  #### Visualise mean and variability across the area through time (~1 sec)
  
  ## Set up plot to save 
  t1 <- Sys.time()
  png("./fig/temp_time_series.png", 
      height = 5, width = 5, res = 800, units = "in")
  
  ## Define blank plot 
  pretty_plot(wc_by_hour$timestamp, wc_by_hour$wc_mean, 
              pretty_axis_args = list(x = list(x = range(wc_by_hour$timestamp), 
                                               y = range(c(wc_by_hour$wc_min, wc_by_hour$wc_max)))),
              xlab = "", ylab = "",
              type = "n")
  
  ## Add range in temperature 
  add_range_by_arrow <- TRUE
  s <- nrow(wc_by_hour)
  scale_lwd <- 25
  arrows(wc_by_hour$timestamp[1:(s-1)], 
         wc_by_hour$wc_min[1:(s-1)], 
         wc_by_hour$timestamp[2:s], 
         wc_by_hour$wc_min[2:s], 
         lwd = wc_by_hour$wc_depth_for_min/scale_lwd, 
         length = 0, 
         col = "royalblue")
  arrows(wc_by_hour$timestamp[1:(s-1)], 
         wc_by_hour$wc_max[1:(s-1)], 
         wc_by_hour$timestamp[2:s], 
         wc_by_hour$wc_max[2:s], 
         lwd = wc_by_hour$wc_depth_for_max/scale_lwd,
         length = 0,
         col = "red3")
  add_range_poly <- FALSE
  if(add_range_poly){
    polygon(c(wc_by_hour$timestamp, rev(wc_by_hour$timestamp)), 
            c(wc_by_hour$wc_min, wc_by_hour$wc_max), 
            col = scales::alpha("grey", 0.9), 
            border = FALSE)
  }
  
  ## Add IQR in temperature 
  add_iqr_poly <- FALSE
  if(add_iqr_poly){
    polygon(c(wc_by_hour$timestamp, rev(wc_by_hour$timestamp)), 
            c(wc_by_hour$wc_q1, wc_by_hour$wc_q3), 
            col = scales::alpha("dimgrey", 0.9), 
            border = FALSE)
  }
  
  ## Add mean temperature 
  # lines(wc_by_hour$timestamp, wc_by_hour$wc_mean)
  
  ## Add legend
  text(1485103235 + 24*60*60*22, 10.2, "max", col = "red3")
  text(1484027049 - 24*60*60*20, 9.6, "min", col = "royalblue")
  zl <- c(5, 25, 100)
  legend(min(wc_by_hour$timestamp) + 60*60*24, 13.2,
         lwd = zl/scale_lwd, 
         col = "red3",
         legend = zl, 
         title = "Depth (m)",
         bty = "n")
  
  ## Add titles
  mtext(side = 1, "Time (months)", line = 2)
  mtext(side = 2, expression("Bottom temperature (" * degree * "C)"), line = 2)
  dev.off()
  t2 <- Sys.time()
  difftime(t2, t1)
  
  #### Visualise difference in temperature through time
  wc_by_hour$range <- wc_by_hour$wc_max - wc_by_hour$wc_min
  wc_by_hour %>% 
    dplyr::mutate(date = as.Date(timestamp), 
                  mm_yy = Tools4ETS::mmyy(date)) %>%
    dplyr::group_by(mm_yy) %>%
    dplyr::summarise(max(range))
  pretty_plot(wc_by_hour$timestamp, wc_by_hour$range,
              pretty_axis_args = list(side = 1:2, control_digits = 1),
              xlab = "", ylab = "",
              type = "n")
  arrows(wc_by_hour$timestamp[1:(s-1)], 
         wc_by_hour$range[1:(s-1)], 
         wc_by_hour$timestamp[2:s], 
         wc_by_hour$range[2:s], 
         lwd = abs(wc_by_hour$wc_depth_for_max - wc_by_hour$wc_depth_for_min)/50, 
         length = 0)
  mtext(side = 1, "Time (months)", line = 2)
  mtext(side = 2, expression("Difference in bottom temperature (" * degree * "C)"), line = 2)
}


################################
################################
#### Analysis (eggs)

if(analysis == "eggs"){
  
  ################################
  #### Summary statistics 
  
  #### Egg station counts
  # Count the number of locations with eggs present/absent 
  length(eggs$cell[eggs$present == 1]); length(eggs$cell[eggs$present == 0])
  # Count the number of locations with eggs present/absent 
  length(unique(eggs$cell[eggs$present == 1])); length(unique(eggs$cell[eggs$present == 0]))
  
  #### Define cells with/without eggs
  cells_with_eggs <- unique(eggs$cell[eggs$present == 1])
  cells_wo_eggs   <- unique(eggs$cell[eggs$present == 0])
  wc_with_eggs <- wc[wc$mesh_ID %in% cells_with_eggs, ]
  
  #### Examine depth and temperature range in areas with/without eggs
  # Depth range
  utils.add::basic_stats(eggs$depth[eggs$present == 1], na.rm = TRUE)
  # Temperature ranges in areas with eggs
  wc_with_eggs[wc_with_eggs$wc == min(wc_with_eggs$wc), ]
  wc_with_eggs[wc_with_eggs$wc == max(wc_with_eggs$wc), ]
  # Average annual temperature in areas with/without eggs 
  wc %>%
    dplyr::mutate(present = eggs$present[match(mesh_ID, eggs$cell)]) %>%
    dplyr::group_by(present) %>%
    # dplyr::group_by(month) %>%
    dplyr::summarise(wc = mean(wc)) # %>%
  # tidyr::pivot_wider(names_from = eggs, values_from = wc)
  
  
  ################################
  #### Visualise time series
  
  ## Set up plot to save 
  t1 <- Sys.time()
  png("./fig/temp_time_series_eggs.png", 
      height = 5, width = 5, res = 800, units = "in")
  
  ## Define blank plot 
  pretty_plot(wc$timestamp, wc$wc, 
              xlab = "", ylab = "",
              type = "n")
  
  ## Add lines for min/max temperature (across the pMPA at large)
  wc_by_hour <- readRDS("./data/wc/temp/wc_by_hour.rds")
  lines(wc_by_hour$timestamp, wc_by_hour$wc_min, col = "royalblue", lwd = 0.5)
  lines(wc_by_hour$timestamp, wc_by_hour$wc_max, col = "red3", lwd = 0.5)
  
  ## Add time series for areas with/without eggs
  add_without <- FALSE
  if(add_without){
    pbapply::pblapply(cells_wo_eggs, function(cell){
      # cell <- cells_wo_eggs[1]
      wc_for_cell <- wc %>% dplyr::filter(mesh_ID == cell)
      lines(wc_for_cell$timestamp, wc_for_cell$wc, 
            col = scales::alpha("black", 0.25),
            lwd = 0.1, 
            lty = 1)
    }) %>% invisible()
  }
  add_with <- TRUE
  if(add_with){
    pbapply::pblapply(cells_with_eggs, function(cell){
      # cell <- cells_with_eggs[1]
      wc_for_cell <- wc %>% dplyr::filter(mesh_ID == cell)
      lines(wc_for_cell$timestamp, wc_for_cell$wc, 
            col = scales::alpha("grey", 0.3),
            lwd = 0.1, 
            lty = 1)
    }) %>% invisible()
  }
  
  ## Add legend
  text(1485103235 + 24*60*60*22, 10.2, "max", col = "red3")
  text(1484027049 - 24*60*60*20, 9.6, "min", col = "royalblue")
  text(1479498413 + 24*60*60*4, 11.4, "eggs", col = "dimgrey")
  add_legend <- FALSE
  if(add_legend){
    legend(min(wc$timestamp) + 60*60*24, 13.2,
           lty = c(1, 1, 1), 
           col = c(scales::alpha("grey", 0.3)), 
           legend = c("Eggs"),
           bty = "n")
  }

  ## Add titles
  mtext(side = 1, "Time (months)", line = 2)
  mtext(side = 2, expression("Bottom temperature (" * degree * "C)"), line = 2)
  dev.off()
  
}


#### End of code. 
################################
################################