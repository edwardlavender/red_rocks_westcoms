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
library(fvcom.tbx)

#### Load data
## MPA boundaries
mpa <- rgdal::readOGR("./data-raw/spatial/mpa/NC_MPA.shp")
## WeStCOMS mesh
# mesh coordinates (nodes)
nodexy <- read.csv("./data-raw/spatial/mesh/mesh_x.csv")
str(nodexy)
# trinodes
trinodes <- read.csv("./data-raw/spatial/mesh/mesh_trinodes.csv")
str(trinodes)


################################
################################
#### Coastline and MPA boundaries

#### Define MPA boundaries
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
coast <- raster::crop(coast, flapper::update_extent(raster::extent(mpa), 0.5))
raster::plot(coast)
raster::crs(coast)
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

#### Build meshs 
## Mesh around nodes
mesh_around_nodes <- build_mesh(nodexy = nodexy,
                                trinodes = trinodes,
                                mesh_type = "element")
## Mesh around elements 
mesh_around_elements <- build_mesh(nodexy = nodexy,
                                   trinodes = trinodes,
                                   mesh_type = "node")

#### Crop meshes within the MPA
## Mesh around nodes
mesh_around_nodes_in_mpa <- raster::crop(mesh_around_nodes, mpa)
mesh_around_nodes_in_mpa$ID
raster::plot(mesh_around_nodes_in_mpa, col = "royalblue")
## Mesh around elements
mesh_around_elements_in_mpa <- raster::crop(mesh_around_elements, mpa)
mesh_around_elements_in_mpa$ID
raster::plot(mesh_around_elements_in_mpa, col = "royalblue")

#### Save meshes
save <- FALSE
if(save){
  saveRDS(nodexy, 
          "./data/spatial/mesh/mesh_nodexy.rds")
  saveRDS(trinodes, 
          "./data/spatial/mesh/mesh_trinodes.rds")
  saveRDS(mesh_around_nodes, 
          "./data/spatial/mesh/mesh_around_nodes.rds")
  saveRDS(mesh_around_elements, 
          "./data/spatial/mesh/mesh_around_elements.rds")
  saveRDS(mesh_around_nodes_in_mpa, 
          "./data/spatial/mesh/mesh_around_nodes_in_mpa.rds")
  saveRDS(mesh_around_elements_in_mpa, 
          "./data/spatial/mesh/mesh_around_elements_in_mpa.rds")
}


#### End of code. 
################################
################################