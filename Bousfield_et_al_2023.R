##################
##### Intro ######
##################

### R scripts for Bousfield et al., Nature Geoscience, 2023 ###
### "Substantial and increasing global losses of timber-producing forest due to wildfires" ###

## This script uses freely available datasets to estimate the spatial extent and temporal trends of stand replacing wildfires on timber-producing forests ##

## Datasets used ##
# Logging layers - Lesiv et al (2022) global forest management map (https://doi.org/10.5281/zenodo.5879022)
#                - Curtis et al (2018) map of forest loss by driver (https://data.globalforestwatch.org/documents/gfw::tree-cover-loss-by-dominant-driver-2022/about)
# Fire layers - Tyukavina et al (2022) global map of wildfire-driven forest loss, regular layer taken directly from (https://glad.umd.edu/dataset/Fire_GFL/)
#             - Burn estimate + 1SE layers (filename: annual_plus_SE...) (low certainty) taken using all burned pixels from Tyukavina with confidence codes 2,3 and 4, 
#             - Burn estimate - 1SE layers (filenmae: high_certainty...) (high certainty) taken using only burned pixels with confidence code 4
#             - No confidence levels supplied for Africa region by Tyukavina et al., hence no uncertainty estimates for this region

## Methodology
# The general methodology is as follows:
# create continental scale grid network for efficient analysis
# calculate the global extent of timber-producing forest for each logging layer
# Divide the globe into 5 different regions (following Tyukavina), each region divided into 0.25 degree grid cells for increased computational efficiency
# In each grid cell, calculate the total area of timber producing forest, the total area of stand-replacing wildfires, and the total overlap between the two
# record year of burn also for temporal trend analysis
# aggregate burned area data at national, regional and global scale as outputs
# repeat this for each logging layer

############################
#### Set up grid layers ####
############################


#### Create empty grid cells for analysis ####

library(sf)
library(raster)
library(sp)
library(dplyr)
library(RColorBrewer)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(s2)
library(terra)

### get overlap of skeleton raster with each fire raster - they are split by continent

#create new skeleton global raster grid
ref_grid<-extent(-180, 180, -90, 90)
ref_grid<-raster(ref_grid)
res(ref_grid)<-0.2498264
#dummy values
projection(ref_grid)<-CRS("+proj=longlat +datum=WGS84 +no_defs    +ellps=WGS84 +towgs84=0,0,0")
ref_grid<- rast(ref_grid)
ref_grid<- as.polygons(ref_grid)
global.skeleton.raster<- st_as_sf(ref_grid)
global.skeleton.raster<- mutate(global.skeleton.raster, global.id = row_number())
global.skeleton.raster<- st_transform(global.skeleton.raster, 4326 )

#country boundaries to retain only grids overlapping with land mass - here using GADM shapefiles (https://gadm.org/data.html)
new.sf.world<- st_read("gadm.countries.territories.shp")
new.sf.world<- st_make_valid(new.sf.world)


### Create empty grid objects for each continent - following fire rasters available from tyukavina ###

#Europe
sf_use_s2(F)
europe.fire.layer<-rast("EUR_fire_forest_loss_2001-21_annual.tif") # europe/eurasia fire raster
europe.bbox<- st_as_sfc(st_bbox(europe.fire.layer))
europe.bbox<- st_make_valid(europe.bbox)
#cut down global grid by bounding box for region
europe.grids<- global.skeleton.raster %>%
  st_filter(y=europe.bbox, .predicate = st_covered_by)
#filter for just grids over land mass
europe.grids<- europe.grids %>%
  st_filter(y=new.sf.world, .predicate = st_intersects) %>%
  mutate(cont = "Europe", grid_id = row_number())

st_write(europe.grids, "europe.grids.shp")

#North America
na.fire.layer<- rast("new.fire.rasters/NAM_fire_forest_loss_2001-21_annual.tif") # north america fire raster
na.bbox<- st_as_sfc(st_bbox(na.fire.layer))
na.bbox<- st_make_valid(na.bbox)
#cut down global grid by bounding box for region
na.grids<- global.skeleton.raster %>%
  st_filter(y=na.bbox, .predicate = st_intersects)
#filter for just grids over land mass
na.grids<-  na.grids %>%
  st_filter(y=new.sf.world, .predicate = st_intersects) %>%
  mutate(cont = "North America", grid_id = row_number())
st_write(na.grids, "north.america.grids.shp")

#South America
sa.fire.layer<- rast("new.fire.rasters/LAM_fire_forest_loss_2001-21_annual.tif") # latin america fire raster
sa.bbox<- st_as_sfc(st_bbox(sa.fire.layer))
sa.bbox<- st_make_valid(sa.bbox)
#cut down global grid by bounding box for region
sa.grids<- global.skeleton.raster %>%
  st_filter(y=sa.bbox, .predicate = st_intersects)
#filter for just grids over land mass
sa.grids<-  sa.grids %>%
  st_filter(y=new.sf.world, .predicate = st_intersects) %>%
  mutate(cont = "South America", grid_id = row_number())
st_write(sa.grids, "south.america.grids.shp")

#South-East Asia / Australasia
aus.fire.layer<- rast("new.fire.rasters/SEA_AUS_fire_forest_loss_2001-21_annual.tif") # SEA/Aus fire raster
aus.bbox<- st_as_sfc(st_bbox(aus.fire.layer))
aus.bbox<- st_make_valid(aus.bbox)
#cut down global grid by bounding box for region
aus.grids<- global.skeleton.raster %>%
  st_filter(y=aus.bbox, .predicate = st_intersects)
#filter for just grids over land mass
aus.grids<-  aus.grids %>%
  st_filter(y=new.sf.world, .predicate = st_intersects) %>%
  mutate(cont = "Australasia", grid_id = row_number())
st_write(aus.grids, "australasia.grids.shp")

#Africa
africa.fire.layer<- rast("new.fire.rasters/AFR_fire_forest_loss_2001-21_annual.tif") # africa fire raster
africa.bbox<- st_as_sfc(st_bbox(africa.fire.layer))
africa.bbox<- st_make_valid(africa.bbox)
#cut down global grid by bounding box for region
africa.grids<- global.skeleton.raster %>%
  st_filter(y=africa.bbox, .predicate = st_intersects)
#filter for just grids over land mass
africa.grids<-  africa.grids %>%
  st_filter(y=new.sf.world, .predicate = st_intersects) %>%
  mutate(cont = "Africa", grid_id = row_number())
st_write(africa.grids, "africa.grids.shp")

#merge all for global grids
all.global.grids<- rbind( aus.grids, europe.grids, na.grids, sa.grids, africa.grids)
#remove duplicate grids from contintent overlaps
all.global.grids<- all.global.grids %>%
  group_by(global.id) %>%
  slice(1) %>%
  dplyr::select(global.id)
st_write(all.global.grids, "all.global.grids.no.overlap.shp")


#######################################
#### Get logging and fire overlaps ####
#######################################

#### First logging layer - Lesiv et al ####

## Global extent of Lesiv et al. timber-producing forest ##
library(dplyr)
library(sf)
library(terra)
sf_use_s2(F)
#set wd
setwd("/shared/edwards_lab1/User/bop20cgb/Fire/Data")

# load in all grids #
all.global.grids.c<- st_read("smaller.grids/all.global.smaller.grids.no.overlap.shp") # this is just a global grid at 0.25 degree resolution, clipped to include only grids overlapping with land area


forest.management.raster<- rast("FML_v3-2_with-colorbar.tif") # this is the Lesiv forest management map


new.sf.world<- vect("new.sf.world.shp") # this is a map of the world based on rnaturalearth, to allocate countries to each grid cell


global.ids<- all.global.grids$global_id

for(i in 1:length(global.ids)){
  print(i)
  task.id<- global.ids[i]
  
  random.grid<- dplyr::filter(all.global.grids, global_id== task.id)
  
  
  #natural logged forest
  #clip by this grid
  grid.forest.management<- terra::crop(forest.management.raster, random.grid)
  grid.forest.management<- as.polygons(grid.forest.management)
  grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
  colnames(grid.forest.management.sf)[1] <- "forest.type"
  
  #filter for logged forests, plantations and planted forest
  grid.forest.management.sf<- grid.forest.management.sf %>%
    filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
    mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
  
  if(nrow(grid.forest.management.sf)>0){
    grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    
    #create dummy layer of no fire to retain info on total forest area of each type
    forest.area.dummy<- mutate(grid.forest.management.sf, global_id = random.grid$global_id)
    forest.area.dummy<- st_transform(forest.area.dummy, "ESRI:54009")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    forest.area.dummy<- st_intersection(forest.area.dummy, crop.sf.world)
    forest.area.dummy$forest.area<- as.numeric(st_area(forest.area.dummy))
    st_geometry(forest.area.dummy)<- NULL  
    save(forest.area.dummy, file =paste0("../Outputs/new.data/burnt.grids/just.raster/smaller.grids/forest.areas/", task.id, ".forest.area.dummy.moll.RData"))
    print("some logged forest")
  } else {print("no logged forest")}
}


#### Fire overlap with timber producing forest ####

## Lesiv North America ##

 
#na grids
  na.grids<- st_read("smaller.grids/north.america.grids.shp")
  
  #load in rasters
  forest.management.raster<- rast("FML_v3-2_with-colorbar.tif")
  new.fire.layer<- rast("new.fire.rasters/NAM_fire_forest_loss_2001-21_annual.tif")
  
  
  new.sf.world<- vect("new.sf.world.shp")
  #set out tasks - each job to do 25
  
  tasks<- na.grids$grid_id
  
  for(i in 1:nrow(na.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(na.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "North America", confidence = "avg")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".north.america.intersection.area.geom.avg.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".north.america.intersection.area.table.avg.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
  ## Estimate + 1SE ##
  
  #load in rasters
  new.fire.layer<- rast("new.fire.rasters/NAM_fire_forest_loss_2001-21_annual_plus_SE.tif")
  
  
  
  for(i in 1:nrow(na.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(na.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- subst(grid.fire, NA, 0)
    
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "North America", confidence = "plus")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".north.america.intersection.area.geom.plus.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".north.america.intersection.area.table.plus.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
  ## Estimate - 1SE ##
  
  #load in rasters
  new.fire.layer<- rast("new.fire.rasters/NAM_fire_forest_loss_2001-21_annual_high_certainty.tif")
  
  
  
  for(i in 1:nrow(na.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(na.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- subst(grid.fire, NA, 0)
    
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "North America", confidence = "minus")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".north.america.intersection.area.geom.minus.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".north.america.intersection.area.table.minus.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  

## Lesiv - Latin Amreica ##

  #Sa grids
  sa.grids<- st_read("smaller.grids/south.america.grids.shp")
  
  #load in rasters
  forest.management.raster<- rast("FML_v3-2_with-colorbar.tif")
  new.fire.layer<- rast("new.fire.rasters/LAM_fire_forest_loss_2001-21_annual.tif")
  
  
  new.sf.world<- vect("new.sf.world.shp")
  #set out tasks - each job to do 25
  
  tasks<- sa.grids$grid_id
  
  for(i in 1:nrow(sa.grids)){
    print(i)
    task.id<- i

        random.grid<- dplyr::filter(sa.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "South America", confidence = "avg")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".south.america.intersection.area.geom.avg.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".south.america.intersection.area.table.avg.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
  ## Estimate + 1SE ##
  
  #load in rasters
  new.fire.layer<- rast("new.fire.rasters/LAM_fire_forest_loss_2001-21_annual_plus_SE.tif")
  
  
  
  for(i in 1:nrow(sa.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(sa.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- subst(grid.fire, NA, 0)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "South America", confidence = "plus")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".south.america.intersection.area.geom.plus.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".south.america.intersection.area.table.plus.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
  ## Estimate - 1SE ##
  
  #load in rasters
  new.fire.layer<- rast("new.fire.rasters/LAM_fire_forest_loss_2001-21_annual_high_certainty.tif")
  
  
  
  for(i in 1:nrow(sa.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(sa.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- subst(grid.fire, NA, 0)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "South America", confidence = "minus")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".south.america.intersection.area.geom.minus.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".south.america.intersection.area.table.minus.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  

## Lesiv - Eurasia ##

  #eur grids
  eur.grids<- st_read("smaller.grids/europe.grids.shp")
  
  #load in rasters
  forest.management.raster<- rast("FML_v3-2_with-colorbar.tif")
  new.fire.layer<- rast("new.fire.rasters/EUR_fire_forest_loss_2001-21_annual.tif")
  
  
  new.sf.world<- vect("new.sf.world.shp")
  #set out tasks - each job to do 25
  
  tasks<- eur.grids$grid_id
  
  for(i in 1:nrow(eur.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(eur.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Eurasia", confidence = "avg")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".eurasia.intersection.area.geom.avg.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".eurasia.intersection.area.table.avg.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
  ## Estimate + 1SE ##
  
  #load in rasters
  new.fire.layer<- rast("new.fire.rasters/EUR_fire_forest_loss_2001-21_annual_plus_SE.tif")
  
  
  
  for(i in 1:nrow(eur.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(eur.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- subst(grid.fire, NA, 0)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Eurasia", confidence = "plus")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".eurasia.intersection.area.geom.plus.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".eurasia.intersection.area.table.plus.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
  ## Estimate - 1SE ##
  
  #load in rasters
  new.fire.layer<- rast("new.fire.rasters/EUR_fire_forest_loss_2001-21_annual_high_certainty.tif")
  
  
  for(i in 1:nrow(eur.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(eur.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- subst(grid.fire, NA, 0)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Eurasia", confidence = "minus")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".eurasia.intersection.area.geom.minus.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".eurasia.intersection.area.table.minus.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
## Lesiv - Africa ##

  #africa grids
  africa.grids<- st_read("smaller.grids/africa.grids.shp")
  
  #load in rasters
  forest.management.raster<- rast("FML_v3-2_with-colorbar.tif")
  new.fire.layer<- rast("new.fire.rasters/AFR_fire_forest_loss_2001-21_annual.tif")
  
  
  new.sf.world<- vect("new.sf.world.shp")
  #set out tasks - each job to do 25
  
  tasks<- africa.grids$grid_id
  
  for(i in 1:nrow(africa.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(africa.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        
        
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Africa", confidence = "avg")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".africa.intersection.area.geom.avg.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".africa.intersection.area.table.avg.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  

## Lesiv - Australasia/South-East Asia/Oceania

  #aus grids
  aus.grids<- st_read("smaller.grids/australasia.grids.shp")
  
  #load in rasters
  forest.management.raster<- rast("FML_v3-2_with-colorbar.tif")
  new.fire.layer<- rast("new.fire.rasters/SEA_AUS_fire_forest_loss_2001-21_annual.tif")
  
  
  new.sf.world<- vect("new.sf.world.shp")
  #set out tasks - each job to do 25
  
  tasks<- aus.grids$grid_id
  
  for(i in 1:nrow(aus.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(aus.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Australasia", confidence = "avg")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".australasia.intersection.area.geom.avg.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".australasia.intersection.area.table.avg.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
  ## Estimate + 1SE ##
  
  #load in rasters
  new.fire.layer<- rast("new.fire.rasters/SEA_AUS_fire_forest_loss_2001-21_annual_plus_SE.tif")
  
  
  
  for(i in 1:nrow(aus.grids)){
    print(i)
    task.id<- i    
    random.grid<- dplyr::filter(aus.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- subst(grid.fire, NA, 0)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Australasia", confidence = "plus")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".australasia.intersection.area.geom.plus.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".australasia.intersection.area.table.plus.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  
  ## Estimate - 1SE ##
  
  #load in rasters
  new.fire.layer<- rast("new.fire.rasters/SEA_AUS_fire_forest_loss_2001-21_annual_high_certainty.tif")
  
  
  
  for(i in 1:nrow(aus.grids)){
    print(i)
    task.id<- i
    
    random.grid<- dplyr::filter(aus.grids, grid_id == task.id)
    
    #check if any fires first
    ## fires ##
    grid.fire<- terra::crop(new.fire.layer, random.grid)
    grid.fire<- subst(grid.fire, NA, 0)
    grid.fire<- as.polygons(grid.fire)
    grid.fire.sf<- sf::st_as_sf(grid.fire)
    colnames(grid.fire.sf)[1] <- "fire.year"
    #filter for fires ( fire.year > 0)
    grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
    grid.fire.sf<- st_make_valid(grid.fire.sf)
    
    if (nrow(grid.fire.sf) > 0) {
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      
      #natural logged forest
      #clip by this grid
      grid.forest.management<- terra::crop(forest.management.raster, random.grid)
      grid.forest.management<- as.polygons(grid.forest.management)
      grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
      colnames(grid.forest.management.sf)[1] <- "forest.type"
      
      #filter for logged forests, plantations and planted forest
      grid.forest.management.sf<- grid.forest.management.sf %>%
        filter(forest.type == 20 | forest.type == 31 | forest.type == 32) %>%
        mutate(forest.type = case_when(forest.type == 20 ~ "Logged forest", forest.type == 31 ~ "Planted forest", forest.type == 32 ~ "Plantation"))
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      
      
      if (nrow(grid.forest.management.sf)>0){
        
        #get area of each intersection of fires over logging
        grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
        grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
        
        grid.fire.sf<- st_make_valid(grid.fire.sf)
        grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
        intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
        #mutate(fire.area = as.numeric(st_area(.))) %>% 
        #dplyr::select(forest.type,  fire.area, fire.year)
        if (nrow(intersection.area) > 0){
          intersection.area <- st_intersection(intersection.area, crop.sf.world)
          intersection.area$fire.area <- as.vector(st_area(intersection.area))
          intersection.area<- intersection.area %>%
            group_by(country, forest.type, fire.year) %>%
            summarise(fire.area = sum(fire.area))
          
          intersection.area<- st_make_valid(intersection.area)
          intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Australasia", confidence = "minus")
          
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".australasia.intersection.area.geom.minus.RData"))
          #put rest in a table
          st_geometry(intersection.area)<- NULL
          save(intersection.area, file =paste0("../Outputs/final.outputs/lesiv/2021/", task.id, ".australasia.intersection.area.table.minus.RData"))
          
          print(paste0(task.id, " - FOREST BURNED"))
        }else {print(paste0(task.id,"- No forest burned"))}
      } else {print(paste0(task.id,"- No forest"))}
    } else {print(paste0(task.id,"- No fires"))}
  }
  


#### Second logging layer - Curtis et al ####


## Global extent of Curtis et al. timber-producing forest ##

library(dplyr)
library(sf)
library(terra)
sf_use_s2(F)
#set wd
setwd("/shared/edwards_lab1/User/bop20cgb/Fire/Data")

# load in all grids #
all.global.grids<- st_read("smaller.grids/all.global.smaller.grids.no.overlap.shp")
forest.cover.raster<- rast("hansen.forest.cover/hansen.100m.global.2000.forest.cover.tif") # Hansen forest cover map for year 2000 at 100m resolution


logging.grids.sf<- st_read("all.curtis.logging.grids.shp") # Shapefile of Curtis et al. raster map of forest loss, where only category 3 (forest loss caused by forestry) has been retained

new.sf.world<- vect("new.sf.world.shp")

#to speed up - get grid id of grids that overlap with logging layer - since curtis files are shapefiles, we can do this easily to identify only grids that overlap with forestry land, could not do this with Lesiv et al. as converting the global 100m raster to shapefile takes too much memory/time
overlapping.grids<-  all.global.grids %>%
  st_filter(y=logging.grids.sf, .predicate = st_intersects)
overlapping.grids<- overlapping.grids$global_id


for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  
  random.grid<- dplyr::filter(all.global.grids, global_id== task.id)
  
  #natural logged forest
  #clip by this grid
  grid.forest.management<- terra::crop(logging.grids.terra, random.grid)
  grid.forest.management<- as.polygons(grid.forest.management)
  grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
  colnames(grid.forest.management.sf)[1] <- "forest.type"
  if (nrow(grid.forest.management.sf)>0){
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    #clip by hansen raster - 100m forest cover
    forest.cover.gfw<- crop(forest.cover.raster, random.grid)
    #keep only area where forest cover > 10%
    forest.cover.gfw[forest.cover.gfw < 10] <- NA
    forest.cover.gfw[forest.cover.gfw >= 10] <- 1
    forest.cover.gfw<- as.polygons(forest.cover.gfw)
    forest.cover.gfw<- st_as_sf(forest.cover.gfw)
    colnames(forest.cover.gfw)[1]<- "forest.cover"
    if(nrow(forest.cover.gfw)>0){
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #trim world map down to grid
      
      crop.sf.world<- crop(new.sf.world, random.grid)
      crop.sf.world<- st_as_sf(crop.sf.world)
      crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
      crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
      #creat dummy layer of no fire to retain info on total forest area of each type
      forest.area.dummy<- mutate(grid.forest.management.sf, global_id = random.grid$global_id,  forest.type = "Logged forest")
      forest.area.dummy<- st_intersection(forest.area.dummy, crop.sf.world)
      forest.area.dummy$forest.area<- as.numeric(st_area(forest.area.dummy))
      st_geometry(forest.area.dummy)<- NULL  
      save(forest.area.dummy, file =paste0("../Outputs/new.data/burnt.grids/curtis/smaller.grids/forest.areas/", task.id, ".forest.area.dummy.moll.RData"))
      print("some logged forest")
    } }else {print("no logged forest")}
}

## Fire overlap with timber producing forest ##

## Curtis - North America ##

#na grids
na.grids<- st_read("smaller.grids/north.america.grids.shp")
forest.cover.raster<- rast("hansen.forest.cover/hansen.100m.global.2000.forest.cover.tif")



#load in rasters
new.fire.layer<- rast("new.fire.rasters/NAM_fire_forest_loss_2001-21_annual.tif")
#curtis logging
logging.grids<- vect("all.curtis.logging.grids.shp")
logging.grids<- crop(logging.grids, new.fire.layer)
logging.grids.sf<- st_read("all.curtis.logging.grids.shp")
new.sf.world<- vect("new.sf.world.shp")
#to speed up - get grid id of grids that overlap with plantation layer
overlapping.grids<-  na.grids %>%
  st_filter(y=logging.grids.sf, .predicate = st_intersects)
overlapping.grids<- overlapping.grids$grid_id

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  
  print(task.id)
  random.grid<- dplyr::filter(na.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "North America", forest.type = "Logged forest", confidence = "avg")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".north.america.intersection.area.geom.avg.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".north.america.intersection.area.table.avg.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Estimate + 1SE ##

#load in rasters
new.fire.layer<- rast("new.fire.rasters/NAM_fire_forest_loss_2001-21_annual_plus_SE.tif")

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(na.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- subst(grid.fire, NA, 0)
  
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "North America", forest.type = "Logged forest", confidence = "plus")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".north.america.intersection.area.geom.plus.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".north.america.intersection.area.table.plus.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Estimate - 1SE ##

#load in rasters
new.fire.layer<- rast("new.fire.rasters/NAM_fire_forest_loss_2001-21_annual_high_certainty.tif")

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(na.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- subst(grid.fire, NA, 0)
  
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "North America", forest.type = "Logged forest", confidence = "minus")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".north.america.intersection.area.geom.minus.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".north.america.intersection.area.table.minus.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Curtis - Latin America ##

#latin america grids
sa.grids<- st_read("smaller.grids/south.america.grids.shp")
forest.cover.raster<- rast("hansen.forest.cover/hansen.100m.global.2000.forest.cover.tif")



#load in rasters
new.fire.layer<- rast("new.fire.rasters/LAM_fire_forest_loss_2001-21_annual.tif")
#curtis logging
logging.grids<- vect("all.curtis.logging.grids.shp")
logging.grids<- crop(logging.grids, new.fire.layer)
logging.grids.sf<- st_read("all.curtis.logging.grids.shp")
new.sf.world<- vect("new.sf.world.shp")
#to speed up - get grid id of grids that overlap with plantation layer
overlapping.grids<-  sa.grids %>%
  st_filter(y=logging.grids.sf, .predicate = st_intersects)
overlapping.grids<- overlapping.grids$grid_id

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(sa.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "South America", forest.type = "Logged forest", confidence = "avg")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".south.america.intersection.area.geom.avg.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".south.america.intersection.area.table.avg.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Estimate + 1SE ##

#load in rasters
new.fire.layer<- rast("new.fire.rasters/LAM_fire_forest_loss_2001-21_annual_plus_SE.tif")

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(sa.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- subst(grid.fire, NA, 0)
  
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "South America", forest.type = "Logged forest", confidence = "plus")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".south.america.intersection.area.geom.plus.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".south.america.intersection.area.table.plus.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Estimate - 1SE ##

#load in rasters
new.fire.layer<- rast("new.fire.rasters/LAM_fire_forest_loss_2001-21_annual_high_certainty.tif")

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(sa.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- subst(grid.fire, NA, 0)
  
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "South America", forest.type = "Logged forest", confidence = "minus")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".south.america.intersection.area.geom.minus.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".south.america.intersection.area.table.minus.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}


## Curtis - Eurasia ##

#Eurasia grids

eur.grids<- st_read("smaller.grids/europe.grids.shp")
forest.cover.raster<- rast("hansen.forest.cover/hansen.100m.global.2000.forest.cover.tif")



#load in rasters
new.fire.layer<- rast("new.fire.rasters/EUR_fire_forest_loss_2001-21_annual.tif")
#curtis logging
logging.grids<- vect("all.curtis.logging.grids.shp")
logging.grids<- crop(logging.grids, new.fire.layer)
logging.grids.sf<- st_read("all.curtis.logging.grids.shp")
new.sf.world<- vect("new.sf.world.shp")
#to speed up - get grid id of grids that overlap with plantation layer
overlapping.grids<-  eur.grids %>%
  st_filter(y=logging.grids.sf, .predicate = st_intersects)
overlapping.grids<- overlapping.grids$grid_id

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(eur.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Eurasia", forest.type = "Logged forest", confidence = "avg")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".eurasia.intersection.area.geom.avg.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".eurasia.intersection.area.table.avg.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Estimate + 1SE ##

#load in rasters
new.fire.layer<- rast("new.fire.rasters/EUR_fire_forest_loss_2001-21_annual_plus_SE.tif")

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(eur.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- subst(grid.fire, NA, 0)
  
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Eurasia", forest.type = "Logged forest", confidence = "plus")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".eurasia.intersection.area.geom.plus.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".eurasia.intersection.area.table.plus.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Estimate - 1SE ##

#load in rasters
new.fire.layer<- rast("new.fire.rasters/EUR_fire_forest_loss_2001-21_annual_high_certainty.tif")

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(eur.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- subst(grid.fire, NA, 0)
  
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Eurasia", forest.type = "Logged forest", confidence = "minus")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".eurasia.intersection.area.geom.minus.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".eurasia.intersection.area.table.minus.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Curtis - Africa ##

#africa grids
afr.grids<- st_read("smaller.grids/africa.grids.shp")
forest.cover.raster<- rast("hansen.forest.cover/hansen.100m.global.2000.forest.cover.tif")



#load in rasters
new.fire.layer<- rast("new.fire.rasters/AFR_fire_forest_loss_2001-21_annual.tif")
#curtis logging
logging.grids<- vect("all.curtis.logging.grids.shp")
logging.grids<- crop(logging.grids, new.fire.layer)
logging.grids.sf<- st_read("all.curtis.logging.grids.shp")
new.sf.world<- vect("new.sf.world.shp")
#to speed up - get grid id of grids that overlap with plantation layer
overlapping.grids<-  afr.grids %>%
  st_filter(y=logging.grids.sf, .predicate = st_intersects)
overlapping.grids<- overlapping.grids$grid_id

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  print(task.id)
  random.grid<- dplyr::filter(afr.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Africa", forest.type = "Logged forest", confidence = "avg")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".africa.intersection.area.geom.avg.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".africa.intersection.area.table.avg.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Curtis - Australasia/South-East Asia/Oceania ##


#Aus grids
aus.grids<- st_read("smaller.grids/australasia.grids.shp")
forest.cover.raster<- rast("hansen.forest.cover/hansen.100m.global.2000.forest.cover.tif")



#load in rasters
new.fire.layer<- rast("new.fire.rasters/SEA_AUS_fire_forest_loss_2001-21_annual.tif")
#curtis logging
logging.grids<- vect("all.curtis.logging.grids.shp")
logging.grids<- crop(logging.grids, new.fire.layer)
logging.grids.sf<- st_read("all.curtis.logging.grids.shp")
new.sf.world<- vect("new.sf.world.shp")
#to speed up - get grid id of grids that overlap with plantation layer
overlapping.grids<-  aus.grids %>%
  st_filter(y=logging.grids.sf, .predicate = st_intersects)
overlapping.grids<- overlapping.grids$grid_id

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  
  print(task.id)
  random.grid<- dplyr::filter(aus.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Australasia", forest.type = "Logged forest", confidence = "avg")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".australasia.intersection.area.geom.avg.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".australasia.intersection.area.table.avg.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Estimate + 1SE ##

#load in rasters
new.fire.layer<- rast("new.fire.rasters/SEA_AUS_fire_forest_loss_2001-21_annual_plus_SE.tif")

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  
  print(task.id)
  random.grid<- dplyr::filter(aus.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- subst(grid.fire, NA, 0)
  
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Australasia", forest.type = "Logged forest", confidence = "plus")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".australasia.intersection.area.geom.plus.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".australasia.intersection.area.table.plus.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}

## Estimate - 1SE ##

#load in rasters
new.fire.layer<- rast("new.fire.rasters/SEA_AUS_fire_forest_loss_2001-21_annual_high_certainty.tif")

for(i in 1:length(overlapping.grids)){
  print(i)
  task.id<- overlapping.grids[i]
  
  print(task.id)
  random.grid<- dplyr::filter(aus.grids, grid_id == task.id)
  
  #check if any fires first
  ## fires ##
  grid.fire<- terra::crop(new.fire.layer, random.grid)
  grid.fire<- subst(grid.fire, NA, 0)
  
  grid.fire<- as.polygons(grid.fire)
  grid.fire.sf<- sf::st_as_sf(grid.fire)
  colnames(grid.fire.sf)[1] <- "fire.year"
  #filter for fires ( fire.year > 0)
  grid.fire.sf<- filter(grid.fire.sf, fire.year > 0)
  grid.fire.sf<- st_make_valid(grid.fire.sf)
  
  if (nrow(grid.fire.sf) > 0) {
    #trim world map down to grid
    
    crop.sf.world<- crop(new.sf.world, random.grid)
    crop.sf.world<- st_as_sf(crop.sf.world)
    crop.sf.world<- rename(crop.sf.world, "country" = "SOVEREIGN")
    crop.sf.world<- st_transform(crop.sf.world, "ESRI:54009")
    
    #natural logged forest
    #clip by this grid
    grid.forest.management<- terra::crop(logging.grids, random.grid)
    grid.forest.management<- as.polygons(grid.forest.management)
    grid.forest.management.sf<- sf::st_as_sf(grid.forest.management)
    colnames(grid.forest.management.sf)[1] <- "forest.type"
    grid.forest.management.sf<- st_transform(grid.forest.management.sf, "ESRI:54009")
    
    if (nrow(grid.forest.management.sf)>0){
      #clip by lesiv raster - 100m forest cover
      forest.cover.gfw<- crop(forest.cover.raster, random.grid)
      #keep only area where forest cover > 10%
      forest.cover.gfw[forest.cover.gfw < 10] <- NA
      forest.cover.gfw[forest.cover.gfw >= 10] <- 1
      forest.cover.gfw<- as.polygons(forest.cover.gfw)
      forest.cover.gfw<- st_as_sf(forest.cover.gfw)
      colnames(forest.cover.gfw)[1]<- "forest.cover"
      forest.cover.gfw<- st_transform(forest.cover.gfw, "ESRI:54009")
      
      forest.cover.gfw<- st_make_valid(forest.cover.gfw)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      #only keep where is forest > 10%
      grid.forest.management.sf<- st_intersection(grid.forest.management.sf, forest.cover.gfw)
      grid.forest.management.sf<- grid.forest.management.sf %>%
        group_by(forest.type, global_id) %>%
        summarise()
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
    }
    if (nrow(grid.forest.management.sf)>0){
      
      
      #get area of each intersection of fires over logging
      grid.fire.sf<- st_transform(grid.fire.sf, "ESRI:54009")
      
      grid.fire.sf<- st_make_valid(grid.fire.sf)
      grid.forest.management.sf<- st_make_valid(grid.forest.management.sf)
      intersection.area <- st_intersection(grid.forest.management.sf, grid.fire.sf)  
      #mutate(fire.area = as.numeric(st_area(.))) %>% 
      #dplyr::select(forest.type,  fire.area, fire.year)
      if (nrow(intersection.area) > 0){
        intersection.area <- st_intersection(intersection.area, crop.sf.world)
        intersection.area$fire.area <- as.vector(st_area(intersection.area))
        intersection.area<- intersection.area %>%
          group_by(country, forest.type, fire.year) %>%
          summarise(fire.area = sum(fire.area))
        
        intersection.area<- st_make_valid(intersection.area)
        intersection.area<- mutate(intersection.area, global_id = random.grid$global_id, grid_id = task.id, region = "Australasia", forest.type = "Logged forest", confidence = "minus")
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".australasia.intersection.area.geom.minus.RData"))
        #put rest in a table
        st_geometry(intersection.area)<- NULL
        save(intersection.area, file =paste0("../Outputs/final.outputs/curtis/2021/", task.id, ".australasia.intersection.area.table.minus.RData"))
        print(paste0(task.id, " - FOREST BURNED"))
      }else {print(paste0(task.id,"- No forest burned"))}
    } else {print(paste0(task.id,"- No forest"))}
  } else {print(paste0(task.id,"- No fires"))}
}