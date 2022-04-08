################################
################################
#### process_data_raw.R

#### This script: 
# 1) Processes raw data for use in this project:
# ... Defines MPA boundaries
# ... Defines coastline in area
# ... Builds WeStCOMS meshes for the area 

#### Steps preceding this script: 
# 1) Obtain raw data          (see README)


################################
################################
#### Set up

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)
library(fvcom.tbx)

#### Load data
## MPA boundaries
mpa <- rgdal::readOGR("./data-raw/spatial/mpa/2022/original_MPA_boundary.shp")
## WeStCOMS mesh
# mesh coordinates (nodes)
nodexy <- read.csv("./data-raw/spatial/mesh/mesh_x.csv")
str(nodexy)
# trinodes
trinodes <- read.csv("./data-raw/spatial/mesh/mesh_trinodes.csv")
str(trinodes)

#### Define global parameters
wgs84 <- sp::CRS(as.character("+init=epsg:4326"))


################################
################################
#### Coastline and MPA boundaries

#### Define MPA boundaries
mpa <- sp::spTransform(mpa, wgs84)
raster::crs(mpa)
raster::plot(mpa)
saveRDS(mpa, "./data/spatial/mpa/mpa.rds")

#### Define coastline
# Download a SpatialPolygonsDataFrame defining the administrative areas of the UK:
download <- FALSE
if(download){
  download.file("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_GBR_0_sp.rds",
                destfile = "./data-raw/spatial/coast/GBR_adm0.rds")
}
coast <- readRDS("./data-raw/spatial/coast/GBR_adm0.rds")
coast <- sp::spTransform(coast, wgs84)
coast <- raster::crop(coast, flapper::update_extent(raster::extent(mpa), 0.5))
raster::plot(coast)
raster::lines(mpa, col = "blue")
saveRDS(coast, "./data/spatial/coast/coast.rds")


################################
################################
#### Build mesh 

#### Data processing
## Node coordinates 
colnames(nodexy) <- c("x", "y", "z")
nodexy$node_id   <- 1:nrow(nodexy)
nodexy           <- dplyr::select(nodexy, node_id, x, y, z)
## Node connections
colnames(trinodes)  <- c("node1", "node2", "node3")
trinodes$element_id <- 1:nrow(trinodes)
trinodes            <- dplyr::select(trinodes, element_id, node1, node2, node3)
## Save dataframes
saveRDS(nodexy, 
        "./data/spatial/mesh/mesh_nodexy.rds")
saveRDS(trinodes, 
        "./data/spatial/mesh/mesh_trinodes.rds")

#### Build meshs 
run <- FALSE
if(run){
  ## Mesh around nodes
  mesh_around_nodes <- build_mesh(nodexy = nodexy,
                                  trinodes = trinodes,
                                  mesh_type = "element")
  ## Mesh around elements 
  mesh_around_elements <- build_mesh(nodexy = nodexy,
                                     trinodes = trinodes,
                                     mesh_type = "node")
  
  ## Save meshes
  saveRDS(mesh_around_nodes, 
          "./data/spatial/mesh/mesh_around_nodes.rds")
  saveRDS(mesh_around_elements, 
          "./data/spatial/mesh/mesh_around_elements.rds")
} else {
  mesh_around_nodes <- 
    readRDS("./data/spatial/mesh/mesh_around_nodes.rds")
  mesh_around_elements <- 
    readRDS("./data/spatial/mesh/mesh_around_elements.rds")
}

#### Crop meshes within the MPA
run <- FALSE
if(run){
  ## Mesh around nodes
  mesh_around_nodes_in_mpa <- raster::crop(mesh_around_nodes, mpa)
  mesh_around_nodes_in_mpa$ID
  raster::plot(mesh_around_nodes_in_mpa, col = "royalblue")
  
  ## Mesh around elements
  mesh_around_elements_in_mpa <- raster::crop(mesh_around_elements, mpa)
  mesh_around_elements_in_mpa$ID
  raster::plot(mesh_around_elements_in_mpa, col = "royalblue")
  
  ## Save meshes
  saveRDS(mesh_around_nodes_in_mpa, 
          "./data/spatial/mesh/mesh_around_nodes_in_mpa.rds")
  saveRDS(mesh_around_elements_in_mpa, 
          "./data/spatial/mesh/mesh_around_elements_in_mpa.rds")
} else {
  mesh_around_nodes_in_mpa <- 
    readRDS("./data/spatial/mesh/mesh_around_nodes_in_mpa.rds")
  mesh_around_elements_in_mpa <- 
    readRDS("./data/spatial/mesh/mesh_around_elements_in_mpa.rds")
}

#### Crop meshes within region 
run <- FALSE
if(run){
  
  ## Define sensible region boundaries 
  ext <- raster::extent(c(-6.0352, -5.8239, 57.25837, 57.40612))
  
  ## Mesh around nodes
  mesh_around_nodes_in_region <- raster::crop(mesh_around_nodes, ext)
  mesh_around_nodes_in_region$ID
  raster::plot(mesh_around_nodes_in_region, col = "royalblue")
  
  ## Mesh around elements
  mesh_around_elements_in_region <- raster::crop(mesh_around_elements, ext)
  mesh_around_elements_in_region$ID
  raster::plot(mesh_around_elements_in_region, col = "royalblue")
  
  ## Save meshes
  saveRDS(mesh_around_nodes_in_region, 
          "./data/spatial/mesh/mesh_around_nodes_in_region.rds")
  saveRDS(mesh_around_elements_in_region, 
          "./data/spatial/mesh/mesh_around_elements_in_region.rds")
  
} else {
  mesh_around_nodes_in_region <- 
    readRDS("./data/spatial/mesh/mesh_around_nodes_in_region.rds")
  mesh_around_elements_in_region <- 
    readRDS("./data/spatial/mesh/mesh_around_elements_in_region.rds")
}

#### Checks
check <- FALSE
if(check){
  ## Check grid resolution in MPA (overestimate because some cells are 'cut' by MPA boundary)
  # Nodes
  length(unique(mesh_around_nodes_in_mpa$ID))
  utils.add::basic_stats(raster::area(mesh_around_nodes_in_mpa, byid = TRUE)/1e6)
  # Elements
  length(unique(mesh_around_elements_in_mpa$ID))
  utils.add::basic_stats(raster::area(mesh_around_elements_in_mpa, byid = TRUE)/1e6)
  
  ## Check that mesh cells in the MPA have been identified correctly, e.g., for elements:
  mesh_around_elements$col <- "black"
  mesh_around_elements$col[
    mesh_around_elements$ID %in% unique(mesh_around_elements_in_mpa$ID)] <- 
    "red"
  raster::plot(mesh_around_elements, col = mesh_around_elements$col, 
               xlim = raster::extent(mesh_around_elements_in_mpa)[1:2], 
               ylim = raster::extent(mesh_around_elements_in_mpa)[3:4])
}


################################
################################
#### Egg locations

#### Read data on egg locations (collected from recent DDV and ROV surveys)
eggs <- 
  readxl::read_excel("./data-raw/eggs/RR&L - locations where eggs were found and depth.xlsx") %>%
  data.frame()
head(eggs)

#### Define clean dataframe of egg locations 
eggs <- data.frame(survey  = "DDV/ROV",
                   station = eggs$Stn, 
                   date    = eggs$Date, 
                   time_1  = eggs$Time_Start, 
                   time_2  = eggs$Time_End_U, 
                   lat_1   = eggs$Lat_start_,
                   long_1  = eggs$Long_start,
                   lat_2   = eggs$Lat_end_DD,
                   long_2  = eggs$Long_end_D,
                   depth_1 = as.numeric(eggs$Depth_Star), # coerce n/a to NA
                   depth_2 = as.numeric(eggs$Depth_St_1), 
                   present = eggs$Skate_egg
                   )
eggs$present[which(is.na(eggs$present))]   <- 0
eggs$present[which(eggs$present == "YES")] <- 1
unique(eggs$present)
eggs$present <- factor(eggs$present)

#### Define average locations
eggs$lat   <- apply(eggs[, c("lat_1", "lat_2")], 1, mean)
eggs$long  <- apply(eggs[, c("long_1", "long_2")], 1, mean)
eggs$depth <- apply(eggs[, c("depth_1", "depth_2")], 1, mean)

#### Get mesh cells
eggs$nodes    <- find_cells(lat = eggs$lat, long = eggs$long, 
                            mesh = mesh_around_nodes_in_region,  
                            proj = raster::crs(mesh_around_nodes_in_region), 
                            f = function(x) as.integer(as.character(x)), 
                            return = 4)
eggs$elements <- find_cells(lat = eggs$lat, long = eggs$long, 
                            mesh = mesh_around_elements_in_region,  
                            proj = raster::crs(mesh_around_elements_in_region), 
                            f = function(x) as.integer(as.character(x)), 
                            return = 4)

#### Check the locations of eggs on the mesh
# Locations in and around the MPA were surveyed
# Eggs were almost exclusively found in the (updated) MPA's boundaries
# One site beyond the MPA was identified with eggs 
png("./fig/map_mesh.png", 
    height = 10, width = 10, units = "in", res = 600)
ext <- raster::extent(mesh_around_nodes_in_mpa)
ext <- flapper::update_extent(ext, 0.05)
xlim <- ext[1:2]
ylim <- ext[3:4]
raster::plot(coast, 
             col = "dimgrey",
             xlim = xlim, ylim = ylim)
raster::lines(mesh_around_nodes_in_mpa, col = "royalblue")
raster::lines(mpa, col = "royalblue")
box()
arrows(x0 = eggs$long_1, 
       y0 = eggs$lat_1, 
       x1 = eggs$long_2, 
       y1 = eggs$lat_2, 
       col = c("red", "darkgreen")[eggs$present],
       length = 0.01)
dev.off()

#### Check the locations of previous DDV and diver surveys (using data from JD)
## -> All of these locations within nodes represented above. 
# Define additional survey locations 
extras <- data.frame(survey = c("DDV", "DDV", "Diver"),
                     date = as.Date(c("2019-03-02", "2018-07-19", "2020-03-03")),
                     lat = c(57.31282, 57.32254, 57.33325), 
                     long = c(-5.87816, -5.89849, -5.92718))
extras$nodes    <- find_cells(lat = extras$lat, long = extras$long, 
                              mesh = mesh_around_nodes_in_region,  
                              proj = raster::crs(mesh_around_nodes_in_region), 
                              f = function(x) as.integer(as.character(x)), 
                              return = 4)
extras$elements <- find_cells(lat = extras$lat, long = extras$long, 
                              mesh = mesh_around_elements_in_region,  
                              proj = raster::crs(mesh_around_elements_in_region), 
                              f = function(x) as.integer(as.character(x)), 
                              return = 4)
# Check whether these occurred in the same WeStCOMS nodes/elements as above
extras$nodes %in% eggs$nodes
extras$elements %in% eggs$elements
# Examine visual check
raster::plot(mesh_around_nodes_in_mpa)
arrows(x0 = eggs$long_1, 
       y0 = eggs$lat_1, 
       x1 = eggs$long_2, 
       y1 = eggs$lat_2, 
       col = c("red", "darkgreen")[eggs$present],
       length = 0.01)
points(extras$long,extras$lat, cex = 2, pch = 4, lwd = 3, col = "orange")

#### Save dataframe
saveRDS(eggs, "./data/eggs/eggs.rds")




#### End of code. 
################################
################################