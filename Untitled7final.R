#### set working directory ####
setwd("~/Desktop/MSc thesis/Data")

#### Load Packages ####
library(httr)
library(terra)
library(sf)
library(inlabru)
library(raster)
library(geodata)
library(osmdata)
library(dplyr)
library(readxl)
library(ggplot2)
library(stringr)
library(mapview)
library(INLA)
library(fmesher)
library(RColorBrewer)
library(concaveman)
library(sn)
library(spatstat.geom)
library(lubridate)
library(tidyterra)
library(RColorBrewer)
library(raster)
library(MASS)
library(nlme)
library(ncf)
library(stringr)
library(measurements)
library(nngeo)
library(kableExtra)
library(spatstat)
library(tmap)
#### Kruger Boundary ####
#kruger vector
kruger<-st_read("knp/knp/Kruger_National_Park.shp")
kruger<-st_transform(kruger,32735)
kruger <- sf::st_union(kruger)
krugerkm<- st_transform(kruger, "+proj=utm +zone=35 +south +units=km +datum=WGS84")
krugerkm<-st_as_sf(krugerkm)
kruger_vect <- vect(krugerkm)
krugersampler<- st_buffer(krugerkm, dist = 1)
names(krugersampler)[names(krugersampler) == attr(krugersampler, "sf_column")] <- "geometry"
attr(krugersampler, "sf_column") <- "geometry"

 
####Latest Sightings Data####
data1 <- "posts_all_years.xlsx"
sheet_names <- excel_sheets(data1)
print(sheet_names)
allsightings <- data.frame()

data2 <- lapply(sheet_names, function(sheet) {
  data <- read_excel(data1, sheet = sheet)
  data$year <- sheet
  # Force 'visibility' to numeric in all sheets
  data$visibility <- as.numeric(data$visibility)
  return(data)
})

allsightings <- bind_rows(data2)
#View(allsightings)
#head(allsightings)
#unique(allsightings$park)
#subset to just kruger and greater kruger
krugersightingsdf <- allsightings %>%
  filter(park %in% c("Kruger National Park", "Greater Kruger"))

#cleaned data - extract coordinates
krugersightingsdf <- krugersightingsdf %>%
  mutate(
    latitude  = as.numeric(str_extract(coordinates, "(?<=latitude': )-?\\d+\\.\\d+")),
    longitude = as.numeric(str_extract(coordinates, "(?<=longitude': )-?\\d+\\.\\d+"))
  ) %>%
  filter(!is.na(latitude) & !is.na(longitude))  # Remove any rows without coords

#convert to spatial object
krugersightingspoints <- st_as_sf(krugersightingsdf, coords = c("longitude", "latitude"), crs = 4326)
krugersightingspoints<- st_transform(krugersightingspoints,32735)
#Clean the data - restrict sightings to kruger
cleanedsf <- st_filter(krugersightingspoints, kruger)
#2024 sightings
kruger2024 <- cleanedsf %>%
  filter(year== 2024)
#mapview(cleanedsf)
#cleaneddf<-st_drop_geometry(cleanedsf)
#### Lion Sightings ####
#Extract only lion sightings (for now)
lion <- cleanedsf %>%
  filter(str_detect(animal, regex("lion", ignore_case = TRUE)))
#mapview(lion)
#Extract only 2024 sightings
lion2024 <- lion %>%
  filter(year== 2024)
# Add date column
lion2024$Date <- as.Date(lion2024$tingTime)
kruger2024$Date <- as.Date(kruger2024$tingTime)
kruger2024km <- st_transform(kruger2024, "+proj=utm +zone=35 +south +units=km +datum=WGS84")
### Preparing data for analysis ####
#Create Mesh

#max.edge is the max length
mesh <- fm_mesh_2d(
  boundary = krugerkm,
  max.edge = c(3,80),
  offset=c(10,20),
  cutoff=3
)
fm_crs(mesh)<-st_crs(krugerkm)
mesh$n
mesh
plot(mesh)

#Create the grid
triangles <- st_sf(geometry = fm_as_sfc(mesh, mask=krugerkm))
#Use triangles as grid
triangles$grid_id <- 1:nrow(triangles)
#### Only take 1 point/day ####
#First step is to find which  grid the points are in
# Join each lion sighting to a grid cell
lion2024km <- st_transform(lion2024, "+proj=utm +zone=35 +south +units=km +datum=WGS84")
liongrid2024 <- st_join(lion2024km, triangles, left = FALSE)
#Spatial and Temporal Bias
effortpergridday<-liongrid2024 %>%
  group_by(Date,grid_id) %>%
  summarise(usersingrid=n_distinct(username),.groups="drop")
totaleffortperday<-kruger2024km %>% ##not sure if I should put for kruger in general or just lion sightings
  group_by(Date) %>%
  summarise(totalusers=n_distinct(username),.groups="drop")
effort<-left_join(effortpergridday,
                  st_drop_geometry(totaleffortperday), 
                  by="Date") %>%
  mutate(releffort=usersingrid/totalusers)
effortdaily <- liongrid2024 %>%
  st_join(effort, by = c("grid_id"="grid_id", "Date"="Date"))

View(effortdaily)
# Keep one sighting per grid cell per day
cleanedlion2024km <- effortdaily %>%
  group_by(grid_id.x, Date.x) %>%
  slice(1) %>%
  ungroup()

 
#### Create boundary ####
# Extract node coordinates from mesh
nodes <- mesh$loc[, 1:2]
# Create convex hull
hull_indices <- chull(nodes)
hull_coords <- nodes[c(hull_indices, hull_indices[1]), ]  # close the polygon
# Create sf polygon
boundary_sf <- st_sfc(st_polygon(list(hull_coords)), crs = fm_crs(mesh))  # or use mesh CRS
boundary<-st_as_sf(boundary_sf)
boundaryvect<-vect(boundary_sf)
plot(boundary)

#### Roads Map ####
# Get the Kruger bounding box
boundary_wgs <- st_transform(boundary, crs = 4326)
boundbox<- st_bbox(boundary_wgs)
# Download roads
roadsosm <- opq(bbox = boundbox) %>%
  add_osm_feature(key = "highway") %>%
  osmdata_sf()
roadskruger<-roadsosm$osm_lines
roadskruger<-st_transform(roadskruger,st_crs(krugerkm))
roadskruger <- st_intersection(roadskruger, krugerkm)
roadlines <- st_as_sf(roadskruger$geometry)
#### Create weights for samplers ####
#Split road by tessellation
roadsegments <- st_intersection(roadlines, triangles)
roadsegments <- st_zm(roadsegments)  # Removes Z and M dimensions
roadsegments <- st_cast(st_collection_extract(roadsegments, "LINESTRING"), "LINESTRING")
#Calculate length of each segment
roadsegments$length <- st_length(roadsegments)
#Count sightinsg per tessellation
sightings_per_cell <- st_join(kruger2024km, triangles, join = st_within) %>%
  group_by(grid_id) %>%
  summarise(sightings = n(),.groups="drop")
#Calculate road length per tessellation
road_effort_per_cell <- roadsegments %>%
  group_by(grid_id) %>%
  summarise(roadlength = sum(length),.groups="drop")
#Get the weights
roadsampler <- left_join(sightings_per_cell, st_drop_geometry(road_effort_per_cell), by = "grid_id") %>%
  mutate(weight = sightings / as.numeric(roadlength))
# Geometry
roadsampler <- roadsegments %>%
  left_join(st_drop_geometry(roadsampler), by = "grid_id")
st_geometry(roadsampler) <- roadsampler$x
st_geometry(roadsampler) <- "geometry"
roadsampler$weight[is.na(roadsampler$weight)] <- 0
roadsampler$weight <- roadsampler$weight + 1e-6
View(roadsampler)

#### NDVI ####
ndvi<- rast("NDVI_Export.tif")
ndvi<-project(ndvi,crs(krugerkm))
crs(ndvi)
ndvi<-crop(ndvi,boundary)
ndvi<-mask(ndvi,boundary)
plot(ndvi)
lines(kruger_vect,col="red",lwd=2)
lines(boundary)
#Coarsen the raster
aggfactor<-round(1/res(ndvi)[1])
ndvi<-aggregate(ndvi,fact=aggfactor,fun=mean)
res(ndvi)
plot(ndvi)
plot(kruger_vect)
#### Canopy ####
# Load canopy height raster
canopy <- rast("canopy_kruger_2020.tif")
canopy<-project(canopy,crs(krugerkm))
# Clip to Kruger boundary
canopykruger <- crop(canopy, boundary)
canopykruger <- mask(canopy, boundary)
res(canopykruger)
plot(canopykruger)

#### Vegetation ####
#kmboundbox<-st_transform(boundbox,crs(krugerkm))
#kmboundbox<-st_as_sfc(kmboundbox)
vegmap <- st_read("vegtype/NVM2018_AEA_V22_7_16082019_final.shp")
vegmap<- st_transform(vegmap,crs(krugerkm))
#Clip to Kruger
vegmapkruger<-st_intersection(vegmap,boundary)
#Make sure it's a factor
vegmapkruger2<-vegmapkruger
vegmapkruger2$BIOREGION_<-as.factor(vegmapkruger2$BIOREGION_)
vegmapkruger2$vegetation<-as.integer(vegmapkruger2$BIOREGION_)
vegmapkruger2$vegetation<-as.factor(vegmapkruger2$vegetation)
#levels(vegmapkruger$BIOME_18)
#levels(vegmapkruger$BIOMEID_18)
#Convert to raster
vegnames <- rasterize(vegmapkruger2, ndvi, field = "BIOREGION_", touches = TRUE)
vegraster <- rasterize(vegmapkruger2, ndvi, field = "vegetation", touches = TRUE)
freq(vegraster)
freq(vegnames)
plot(vegraster)

#The problem might be that there are too few pixels of some vegetation categories. Reclassify:
rc<-as.matrix(data.frame(from=c(0,1,2,3,4,5,6),to=c(5,3,3,3,3,5,3)))
vegraster2<-classify(vegraster,rc)
freq(vegraster2)

#### rest camps ####
# Read CSV exported from Numbers
rest <- read.csv("restcamps.csv")

convert_dms <- function(x) {
  # Remove degree, minute, second symbols and clean whitespace
  cleaned <- gsub("[°º]", " ", x)          # degree symbols
  cleaned <- gsub("[’‘′']", " ", cleaned)  # minute symbols
  cleaned <- gsub("[”″\"]", "", cleaned)   # second symbols
  cleaned <- gsub("[[:space:]]+", " ", cleaned) # collapse whitespace
  
  direction <- toupper(gsub("[^NSEW]", "", cleaned))
  dms <- trimws(gsub("[NSEW]", "", cleaned))
  
  # Try converting
  dec <- suppressWarnings(conv_unit(dms, from = "deg_min_sec", to = "dec_deg"))
  
  # Convert to numeric safely
  dec <- suppressWarnings(as.numeric(dec))
  
  if (is.na(dec)) return(NA)
  if (direction %in% c("S", "W")) dec <- -dec
  
  return(dec)
}

rest <- rest %>%
  mutate(
    latitude = sapply(Latitude, convert_dms),
    longitude = sapply(Longitude, convert_dms)
  )

rest <- rest %>% dplyr::select(-Latitude, -Longitude)
View(rest)
#write.csv(rest, "rest.csv", row.names = FALSE)
restsf<-st_as_sf(rest,coords=c("longitude","latitude"),crs=4326)
restkm<- st_transform(restsf, "+proj=utm +zone=35 +south +units=km +datum=WGS84")
# 2. Create a template raster that covers your study area
raster_template <- rast(ext(boundary), resolution = 1, crs=boundary)
# 3. Rasterize: calculate distance from every cell to the nearest rest camp
# Convert sf to SpatVector
rest_vect <- vect(restkm)
# 4. Create distance raster
restdist <- distance(raster_template, y = rest_vect)
plot(restdist, main = "Distance to nearest rest camp")
#### waterbodies ####
# Download waterways
waterways <- opq(bbox = boundbox) %>%
  add_osm_feature(key = "waterway") %>%
  osmdata_sf()
waterways_sf<-waterways$osm_lines
waterways_sf<-st_transform(waterways_sf,st_crs(boundary))
waterways_sf <- st_intersection(waterways_sf, boundary)
waterways_sf$waterway<-as.factor(waterways_sf$waterway)

# distance raster
unique(waterways_sf$waterway)
riversf <- waterways_sf[waterways_sf$waterway %in% c("river","dam"), ]
riversfclean <- riversf[!st_is_empty(riversf), ]
rivers_vect <- terra::vect(riversfclean)

waterdist <- distance(raster_template, y = rivers_vect)
plot(waterdist)

#### Create a single raster stack ####
#Align them together
base_raster<-ndvi
layers <- list(canopykruger, vegraster2,waterdist,restdist)
aligned <- lapply(layers, function(r) {
  project(r, base_raster, method = "near") 
})

#Stack rasters
# Stack all aligned layers with aspect_raster
gcov <- c(
  ndvi,
  aligned[[1]],  # canopy
  aligned[[2]],  # vegetation
  aligned[[3]],  # waterways
  aligned[[4]]  # rest camps
)

names(gcov) <- c(
  "ndvi", "canopy", "vegetation","waterways","restcamp"
)

#Save file
writeRaster(gcov, "gcovimpala.tif", overwrite = TRUE)

#### Load raster ####
gcov<- rast("gcovimpala.tif")
 
#### Collinearity check ####
covars <- data.frame(
  canopy = values(gcov$canopy),
  waterways = values(gcov$waterways),
  ndvi = values(gcov$ndvi),
  restcamp = values(gcov$restcamp),
  impala = values(gcov$impala)
)

cor(covars, method = "pearson", use = "pairwise.complete.obs")

lm_fit <- lm(canopy ~ waterways + ndvi + restcamp + impala, data = covars)

# Calculate VIF
vif(lm_fit)
#### create matern ####
matern <- inla.spde2.pcmatern(
  mesh,
  prior.sigma = c(1, 0.01), #probability that true sigma is larger than specified
  prior.range = c(30, 0.5) #probability that true range is smaller than specified
)

#### Create prediction grid ####
pred.df <- fm_pixels(mesh, mask = krugerkm)
 
#### vegetation ####
#comp1 <- geometry ~ vegetation(gcov$vegetation, model = "factor_full") - 1
#fit1 <- lgcp(comp1, cleanedlion2024km, samplers = roadlines, domain = list(geometry = mesh),weights=weight_scaled,options = list(control.fixed = list(expand.factor.strategy = "inla")))
#int1 <- predict(fit1, pred.df, ~ exp(vegetation))

#### waterways ####
waterways <- gcov$waterways
waterways<-scale(waterways)
# Note: this method is usually not needed; the automatic invocation of
# `eval_spatial()` method by the component input evaluator is usually
# sufficient.
f.waterways <- function(where) {
  # Extract the values
  v <- eval_spatial(waterways, where, layer = "waterways")
  # Fill in missing values; this example would work for
  # SpatialPixelsDataFrame data
  if (any(is.na(v))) {
    v <- bru_fill_missing(waterways, where, v)
  }
  v
}
#### ndvi####
ndvi <- gcov$ndvi
ndvi<-scale(ndvi)
# Note: this method is usually not needed; the automatic invocation of
# `eval_spatial()` method by the component input evaluator is usually
# sufficient.
f.ndvi <- function(where) {
  # Extract the values
  v <- eval_spatial(ndvi, where, layer = "ndvi")
  # Fill in missing values; this example would work for
  # SpatialPixelsDataFrame data
  if (any(is.na(v))) {
    v <- bru_fill_missing(ndvi, where, v)
  }
  v
}
#### canopy ####
canopy <- gcov$canopy
hist(values(gcov$canopy[[1]]), main="Canopy (range_min)", breaks=50)
canopy<-scale(canopy)
# Note: this method is usually not needed; the automatic invocation of
# `eval_spatial()` method by the component input evaluator is usually
# sufficient.
f.canopy <- function(where) {
  # Extract the values
  v <- eval_spatial(canopy, where, layer = "canopy")
  # Fill in missing values; this example would work for
  # SpatialPixelsDataFrame data
  if (any(is.na(v))) {
    v <- bru_fill_missing(canopy, where, v)
  }
  v
}

#### rest camps ####
restcamp <- gcov$restcamp
restcamp<-scale(restcamp)
# Note: this method is usually not needed; the automatic invocation of
# `eval_spatial()` method by the component input evaluator is usually
# sufficient.
f.restcamp <- function(where) {
  # Extract the values
  v <- eval_spatial(restcamp, where, layer = "restcamp")
  # Fill in missing values; this example would work for
  # SpatialPixelsDataFrame data
  if (any(is.na(v))) {
    v <- bru_fill_missing(restcamp, where, v)
  }
  v
}
#### impala ####
impala<-gcov$impala
 
#### Specify Models ####
cmpnull<- geometry ~   mySmooth(geometry, model = matern) + Intercept(1)
cmp1<- geometry ~ ndvi(f.ndvi(.data.), model = "linear") + 
  mySmooth(geometry, model=matern) + Intercept(1)
cmp2<- geometry ~ vegetation(gcov$vegetation, model = "factor_contrast") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp3<- geometry ~ canopy(f.canopy(.data.), model = "linear") + 
  mySmooth(geometry, model=matern) + Intercept(1)
cmp4<- geometry ~ waterways(waterways, model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp5<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  mySmooth(geometry, model=matern) + Intercept(1)

cmp6<- geometry ~ vegetation(gcov$vegetation, model = "factor_contrast") +
  ndvi(f.ndvi(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)
veg3_mask<- gcov$vegetation == 3
ndvi_veg3 <- mask(ndvi, veg3_mask, maskvalue=FALSE)
veg5_mask<- gcov$vegetation == 5
ndvi_veg5 <- mask(ndvi, veg5_mask, maskvalue=FALSE)
cmp6b<- geometry ~ vegetation(gcov$vegetation, model = "factor_full") +
  ndvi_veg3(ndvi_veg3, model = "linear") +
  ndvi_veg5(ndvi_veg5, model = "linear") +
  mySmooth(geometry, model=matern) - Intercept(1)


cmp7<- geometry ~ canopy(f.canopy(.data.), model = "linear") +
  vegetation(gcov$vegetation, model = "factor_contrast") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp8<- geometry ~ waterways(f.waterways(.data.), model = "linear") +
  vegetation(gcov$vegetation, model = "factor_contrast") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp9<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  vegetation(gcov$vegetation, model = "factor_contrast") +
  mySmooth(geometry, model=matern) + Intercept(1)


cmp10<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  ndvi(f.ndvi(.data.), model = "linear") +
  canopy(f.canopy(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  vegetation(gcov$vegetation, model = "factor_contrast") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp11<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  ndvi(f.ndvi(.data.), model = "linear") +
  canopy(f.canopy(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

cmp12<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

cmp13<- geometry ~ impala(impala, model = "linear") + 
  mySmooth(geometry, model=matern) + Intercept(1)

canopy_x_impala <- canopy * impala
cmp14<- geometry ~ canopy(canopy, model = "linear") +
  impala(impala, model = "linear") +
  canopy_x_impala(canopy_x_impala, model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

impala_x_ndvi <- impala * ndvi
cmp15<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  ndvi(f.ndvi(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

waterways_x_ndvi <- waterways * ndvi
cmp16<- geometry ~ waterways_x_ndvi(waterways_x_ndvi, model = "linear") + 
  ndvi(f.ndvi(.data.), model = "linear") +
  canopy(f.canopy(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

waterways_x_restcamp <- waterways * restcamp
cmp17<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  waterways(f.waterways(.data.), model = "linear") +
  waterways_x_restcamp(waterways_x_restcamp, model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

waterways_x_canopy <- waterways * canopy
cmp18<- geometry ~ canopy(f.canopy(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  waterways_x_canopy(waterways_x_canopy, model = "linear") +
  vegetation(gcov$vegetation, model = "factor_contrast") +
  mySmooth(geometry, model=matern) + Intercept(1)


waterways_x_canopy <- waterways * canopy
cmp30<- geometry ~ canopy(f.canopy(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  waterways_x_canopy(waterways_x_canopy, model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

cmp31<- geometry ~ canopy(f.canopy(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

cmp32<- geometry ~ waterways_x_ndvi(waterways_x_ndvi, model = "linear") + 
  ndvi(f.ndvi(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp33<- geometry ~ 
  ndvi(f.ndvi(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp4<- geometry ~ waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

cmp34<- geometry ~ waterways(f.waterways(.data.), model = "linear") +
  impala(impala, model = "linear") +
  ndvi(f.ndvi(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp35<- geometry ~ waterways(waterways, model = "linear") +
  impala(impala, model = "linear") +
  ndvi(ndvi, model = "linear") +
  impala_x_ndvi(impala_x_ndvi, model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp40<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  ndvi(f.ndvi(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  vegetation(gcov$vegetation, model = "factor_contrast") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp41<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  ndvi(f.ndvi(.data.), model = "linear") +
  vegetation(gcov$vegetation, model = "factor_contrast") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp42<- geometry ~ restcamp(f.restcamp(.data.), model = "linear") + 
  ndvi(f.ndvi(.data.), model = "linear") +
  impala(impala, model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp43<- geometry ~ restcamp(restcamp, model = "linear") + 
  ndvi(ndvi, model = "linear") +
  impala(impala, model = "linear") +
  impala_x_ndvi(impala_x_ndvi, model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

cmp50<- geometry ~ ndvi(f.ndvi(.data.), model = "linear") +
  waterways(f.waterways(.data.), model = "linear") +
  waterways_x_ndvi(waterways_x_ndvi, model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)
cmp51<- geometry ~ canopy(f.canopy(.data.), model = "linear") + 
  waterways(f.waterways(.data.), model = "linear") +
  mySmooth(geometry, model=matern) + Intercept(1)

#### Fit models: lions ####
fitnull<- lgcp(cmpnull, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit1<- lgcp(cmp1, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit2<- lgcp(cmp2, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit3<- lgcp(cmp3, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit4<- lgcp(cmp4, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit5<- lgcp(cmp5, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))

fit10<-lgcp(cmp10, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit11<-lgcp(cmp11, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit13<-lgcp(cmp13, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit15<-lgcp(cmp15, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit16<-lgcp(cmp16, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
fit51<-lgcp(cmp51, cleanedlion2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))

deltaIC(fit3, fit4, fit51, fit10, fit11, fit16, criterion = c("DIC"))

allrange <- range(fit4.lp$all$mean)
allogrange<- range(fit4.lp$allog$mean)

csc <- scale_fill_viridis_c(option = "D", limits = allrange)
csclog<- scale_fill_viridis_c(option = "D", limits = allogrange)
ggplot() +
  gg(fit4.lp$all,geom="tile") +
  theme(legend.position = "bottom") +
  geom_sf(data=kruger, alpha = 0) +
  csc +
  theme_bw()

#### Posterior uncertainty ####
Lambda <- predict(
  fit4,
  fm_int(mesh, roadsampler),
  ~ sum(weight * exp(mySmooth + waterways + Intercept))
)
Lambda

total <- predict(
  fit4,
  fm_int(mesh, roadsampler),
  ~ data.frame(
    N = 5000:6000,
    dpois(5000:6000,
          lambda = sum(weight * exp(mySmooth + waterways + Intercept))
    )
  )
)
inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x = 5000:6000, y = total$mean))
inla.emarginal(identity, marginal = list(x = 5000:6000, y = total$mean))
total$plugin_estimate <- dpois(5000:6000, lambda = Lambda$mean)

ggplot(total) +
  geom_line(aes(x = N, y = mean, colour = "Posterior")) +
  geom_line(aes(x = N, y = plugin_estimate, colour = "Plugin"))

 
#### CARNIVORES ####
#### Extract Leopard Data####
#Extract only leopard sightings (for now)
leopard <- cleanedsf %>%
  filter(str_detect(animal, regex("leopard", ignore_case = TRUE)))
#Extract only 2024 sightings
leopard2024 <- leopard %>%
  filter(year== 2024)
#### Only take 1 point/day ####
#First step is to find which  grid the points are in
# Add date column
leopard2024$Date <- as.Date(leopard2024$tingTime)
# Join each leopard sighting to a grid cell
leopard2024km <- st_transform(leopard2024, "+proj=utm +zone=35 +south +units=km +datum=WGS84")
leopardgrid2024 <- st_join(leopard2024km, triangles, left = FALSE)
#Spatial and Temporal Bias
effortpergridday<-leopardgrid2024 %>%
  group_by(Date,grid_id) %>%
  summarise(usersingrid=n_distinct(username),.groups="drop")
totaleffortperday<-kruger2024km %>% ##not sure if I should put for kruger in general or just leopard sightings
  group_by(Date) %>%
  summarise(totalusers=n_distinct(username),.groups="drop")
effort<-left_join(effortpergridday,
                  st_drop_geometry(totaleffortperday), 
                  by="Date") %>%
  mutate(releffort=usersingrid/totalusers)
effortdaily <- leopardgrid2024 %>%
  st_join(effort, by = c("grid_id"="grid_id", "Date"="Date"))

View(effortdaily)
# Keep one sighting per grid cell per day
cleanedleopard2024km <- effortdaily %>%
  group_by(grid_id.x, Date.x) %>%
  slice(1) %>%
  ungroup()

#### leopard models ####

lfit8<-lgcp(cmp8, cleanedleopard2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
lfit10<-lgcp(cmp10, cleanedleopard2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
lfit11<-lgcp(cmp11, cleanedleopard2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
lfit14<-lgcp(cmp14, cleanedleopard2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
lfit18<-lgcp(cmp18, cleanedleopard2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
deltaIC(lfit8, lfit10, lfit11, lfit14, lfit18, criterion = c("DIC"))

#Predicting Surface Leopard
lfit14.lp <- predict(
  lfit14,
  pred.df,
  ~ list(
    all= exp(mySmooth + impala + canopy + canopy_x_impala + Intercept),
    allog= mySmooth + impala + canopy + canopy_x_impala + Intercept
  )
)
lallrange <- range(lfit14.lp$all$mean)
lallogrange<- range(lfit14.lp$allog$mean)
lcsc <- scale_fill_viridis_c(option = "D", limits = lallrange)
lcsclog<- scale_fill_viridis_c(option = "D", limits = lallogrange)
ggplot() +
  gg(lfit14.lp$all,geom="tile") +
  theme(legend.position = "bottom") +
  geom_sf(data=kruger, alpha = 0) +
  lcsc

#### Posterior Uncertainty for Leopard ####
lLambda <- predict(
  lfit14,
  fm_int(mesh, roadsampler),
  ~ sum(weight * exp(mySmooth + impala + canopy + canopy_x_impala + Intercept))
)
lLambda

ltotal <- predict(
  lfit14,
  fm_int(mesh, roadsampler),
  ~ data.frame(
    N = 2700:3700,
    dpois(2700:3700,
          lambda = sum(weight * exp(mySmooth + impala + canopy + canopy_x_impala + Intercept))
    )
  )
)
inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x = 2700:3700, y = ltotal$mean))
inla.emarginal(identity, marginal = list(x = 2700:3700, y = ltotal$mean))
ltotal$plugin_estimate <- dpois(2700:3700, lambda = lLambda$mean)

ggplot(ltotal) +
  geom_line(aes(x = N, y = mean, colour = "Posterior")) +
  geom_line(aes(x = N, y = plugin_estimate, colour = "Plugin"))

#### Extract Cheetah Data####
#Extract only cheetah sightings (for now)
cheetah <- cleanedsf %>%
  filter(str_detect(animal, regex("cheetah", ignore_case = TRUE)))
#Extract only 2024 sightings
cheetah2024 <- cheetah %>%
  filter(year== 2024)
#### Only take 1 point/day ####
#First step is to find which  grid the points are in
# Add date column
cheetah2024$Date <- as.Date(cheetah2024$tingTime)
# Join each cheetah sighting to a grid cell
cheetah2024km <- st_transform(cheetah2024, "+proj=utm +zone=35 +south +units=km +datum=WGS84")
cheetahgrid2024 <- st_join(cheetah2024km, triangles, left = FALSE)
#Spatial and Temporal Bias
effortpergridday<-cheetahgrid2024 %>%
  group_by(Date,grid_id) %>%
  summarise(usersingrid=n_distinct(username),.groups="drop")
totaleffortperday<-kruger2024km %>% ##not sure if I should put for kruger in general or just cheetah sightings
  group_by(Date) %>%
  summarise(totalusers=n_distinct(username),.groups="drop")
effort<-left_join(effortpergridday,
                  st_drop_geometry(totaleffortperday), 
                  by="Date") %>%
  mutate(releffort=usersingrid/totalusers)
effortdaily <- cheetahgrid2024 %>%
  st_join(effort, by = c("grid_id"="grid_id", "Date"="Date"))

View(effortdaily)
# Keep one sighting per grid cell per day
cleanedcheetah2024km <- effortdaily %>%
  group_by(grid_id.x, Date.x) %>%
  slice(1) %>%
  ungroup()

#### cheetah models ####

cfit4<- lgcp(cmp4, cleanedcheetah2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
cfit10<-lgcp(cmp10, cleanedcheetah2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
cfit11<-lgcp(cmp11, cleanedcheetah2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
cfit32<-lgcp(cmp32, cleanedcheetah2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
cfit33<-lgcp(cmp33, cleanedcheetah2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
cfit34<-lgcp(cmp34, cleanedcheetah2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
cfit35<-lgcp(cmp35, cleanedcheetah2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
deltaIC(cfit4, cfit10, cfit11, cfit32, cfit33, cfit34, cfit35, criterion = c("DIC"))

#Predicting Surface Cheetah
cfit35.lp <- predict(
  cfit35,
  pred.df,
  ~ list(
    all=exp(mySmooth + waterways + impala + ndvi + impala_x_ndvi + Intercept),
    allog=mySmooth + waterways + impala + ndvi + impala_x_ndvi + Intercept
  )
)
callrange <- range(cfit35.lp$all$mean)
callogrange<- range(cfit35.lp$allog$mean)
ccsc <- scale_fill_viridis_c(option = "D", limits = callrange)
ccsclog<- scale_fill_viridis_c(option = "D", limits = callogrange)
ggplot() +
  gg(cfit35.lp$all,geom="tile") +
  theme(legend.position = "bottom") +
  geom_sf(data=kruger, alpha = 0) +
  ccsc
#### Posterior Uncertainty for Cheetah ####
cLambda <- predict(
  cfit35,
  fm_int(mesh, roadsampler),
  ~ sum(weight * exp(mySmooth + waterways + impala + ndvi + impala_x_ndvi + Intercept))
)
cLambda

ctotal <- predict(
  cfit35,
  fm_int(mesh, roadsampler),
  ~ data.frame(
    N = 750:1000,
    dpois(750:1000,
          lambda = sum(weight * exp(mySmooth + waterways + impala + ndvi + impala_x_ndvi + Intercept))
    )
  )
)
inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x = 750:1000, y = ctotal$mean))
inla.emarginal(identity, marginal = list(x = 750:1000, y = ctotal$mean))
ctotal$plugin_estimate <- dpois(750:1000, lambda = cLambda$mean)

ggplot(ctotal) +
  geom_line(aes(x = N, y = mean, colour = "Posterior")) +
  geom_line(aes(x = N, y = plugin_estimate, colour = "Plugin"))
#### Extract Hyena Data####
#Extract only hyena sightings (for now)
hyena <- kruger2024 %>%
  filter(str_detect(animal, regex("hyena", ignore_case = TRUE)) &
           !str_detect(animal, regex("brown", ignore_case = TRUE)))
#Extract only 2024 sightings
hyena2024 <- hyena %>%
  filter(year== 2024)
#### Only take 1 point/day ####
#First step is to find which  grid the points are in
# Add date column
hyena2024$Date <- as.Date(hyena2024$tingTime)
# Join each hyena sighting to a grid cell
hyena2024km <- st_transform(hyena2024, "+proj=utm +zone=35 +south +units=km +datum=WGS84")
hyenagrid2024 <- st_join(hyena2024km, triangles, left = FALSE)
#Spatial and Temporal Bias
effortpergridday<-hyenagrid2024 %>%
  group_by(Date,grid_id) %>%
  summarise(usersingrid=n_distinct(username),.groups="drop")
totaleffortperday<-kruger2024km %>% ##not sure if I should put for kruger in general or just hyena sightings
  group_by(Date) %>%
  summarise(totalusers=n_distinct(username),.groups="drop")
effort<-left_join(effortpergridday,
                  st_drop_geometry(totaleffortperday), 
                  by="Date") %>%
  mutate(releffort=usersingrid/totalusers)
effortdaily <- hyenagrid2024 %>%
  st_join(effort, by = c("grid_id"="grid_id", "Date"="Date"))

View(effortdaily)
# Keep one sighting per grid cell per day
cleanedhyena2024km <- effortdaily %>%
  group_by(grid_id.x, Date.x) %>%
  slice(1) %>%
  ungroup()

#### hyena models ####

hfit10<-lgcp(cmp10, cleanedhyena2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
hfit11<-lgcp(cmp11, cleanedhyena2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
hfit40<-lgcp(cmp40, cleanedhyena2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
hfit41<-lgcp(cmp41, cleanedhyena2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
hfit42<-lgcp(cmp42, cleanedhyena2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
hfit34<-lgcp(cmp34, cleanedhyena2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
hfit12<-lgcp(cmp12, cleanedhyena2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
hfit43<-lgcp(cmp43, cleanedhyena2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))

deltaIC(hfit12,hfit10, hfit11,hfit40, hfit41, hfit42, hfit34, hfit43, criterion = c("DIC"))

#Predicting Surface hyena
hfit43.lp <- predict(
  hfit43,
  pred.df,
  ~ list(
    all=exp(mySmooth + restcamp + impala + ndvi + impala_x_ndvi + Intercept),
    allog=mySmooth + restcamp + impala + ndvi + impala_x_ndvi + Intercept
  )
)
hallrange <- range(hfit43.lp$all$mean)
hallogrange<- range(hfit43.lp$allog$mean)
hcsc <- scale_fill_viridis_c(option = "D", limits = hallrange)
hcsclog<- scale_fill_viridis_c(option = "D", limits = hallogrange)
ggplot() +
  gg(hfit43.lp$all,geom="tile") +
  theme(legend.position = "bottom") +
  geom_sf(data=kruger, alpha = 0) +
  hcsc +
  ggtitle("hyena (waterways restcamp)")

#### Posterior Uncertainty for hyena ####
hLambda <- predict(
  hfit43,
  fm_int(mesh, roadsampler),
  ~ sum(weight * exp(mySmooth + restcamp + impala + ndvi + impala_x_ndvi + Intercept))
)
hLambda

htotal <- predict(
  hfit43,
  fm_int(mesh, roadsampler),
  ~ data.frame(
    N = 2500:3700,
    dpois(2500:3700,
          lambda = sum(weight * exp(mySmooth + restcamp + impala + ndvi + impala_x_ndvi + Intercept))
    )
  )
)
inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x = 2500:3700, y = htotal$mean))
inla.emarginal(identity, marginal = list(x = 2500:3700, y = htotal$mean))
htotal$plugin_estimate <- dpois(2500:3700, lambda = hLambda$mean)

ggplot(htotal) +
  geom_line(aes(x = N, y = mean, colour = "Posterior")) +
  geom_line(aes(x = N, y = plugin_estimate, colour = "Plugin"))

 
#### HERBIVORE ####
#### Extract impala Data####
#Extract only impala sightings (for now)
impala <- cleanedsf %>%
  filter(str_detect(animal, regex("impala", ignore_case = TRUE)))
#Extract only 2024 sightings
impala2024 <- impala %>%
  filter(year== 2024)
#### Only take 1 point/day ####
#First step is to find which  grid the points are in
# Add date column
impala2024$Date <- as.Date(impala2024$tingTime)
# Join each impala sighting to a grid cell
impala2024km <- st_transform(impala2024, "+proj=utm +zone=35 +south +units=km +datum=WGS84")
impalagrid2024 <- st_join(impala2024km, triangles, left = FALSE)
#Spatial and Temporal Bias
effortpergridday<-impalagrid2024 %>%
  group_by(Date,grid_id) %>%
  summarise(usersingrid=n_distinct(username),.groups="drop")
totaleffortperday<-kruger2024km %>% ##not sure if I should put for kruger in general or just impala sightings
  group_by(Date) %>%
  summarise(totalusers=n_distinct(username),.groups="drop")
effort<-left_join(effortpergridday,
                  st_drop_geometry(totaleffortperday), 
                  by="Date") %>%
  mutate(releffort=usersingrid/totalusers)
effortdaily <- impalagrid2024 %>%
  st_join(effort, by = c("grid_id"="grid_id", "Date"="Date"))

View(effortdaily)
# Keep one sighting per grid cell per day
cleanedimpala2024km <- effortdaily %>%
  group_by(grid_id.x, Date.x) %>%
  slice(1) %>%
  ungroup()
#calculate effort
cleanedimpala2024km$weight <- 1/cleanedimpala2024km$releffort
cleanedimpala2024km$weight_scaled <- cleanedimpala2024km$weight / mean(cleanedimpala2024km$weight)

#### impala models ####

ifit10<-lgcp(cmp10, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit11<-lgcp(cmp11, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit1<- lgcp(cmp1, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit2<- lgcp(cmp2, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit3<- lgcp(cmp3, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit4<- lgcp(cmp4, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit5<- lgcp(cmp5, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit6<- lgcp(cmp6, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit30<- lgcp(cmp30, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit50<- lgcp(cmp50, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
ifit6b<- lgcp(cmp6b, cleanedimpala2024km, samplers = roadsampler, domain = list(geometry = fm_subdivide(mesh, 2)), options = list(
  control.compute = list(dic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)))
deltaIC(ifit1, ifit2, ifit3, ifit4, ifit5, ifit10,ifit11, ifit30, ifit50, ifit6, ifit6b, criterion = c("DIC"))

####Predicting surface Impala####

ifit6b.lp <- predict(
  ifit6b,
  pred.df,
  ~ list(
    all=exp(mySmooth + ndvi_veg3 + ndvi_veg5 + vegetation),
    allog=mySmooth + ndvi_veg3 + ndvi_veg5 + vegetation
  )
)
iallrange <- range(ifit6b.lp$all$mean,ifit6b.lp$ndvi5$mean,ifit6b.lp$ndvi3$mean,ifit6b.lp$vegetation$mean)
iallogrange<- range(ifit6b.lp$allog$mean)
icsc <- scale_fill_viridis_c(option = "D", limits = iallrange)
icsclog<- scale_fill_viridis_c(option = "D", limits = iallogrange)
ggplot() +
  gg(ifit6b.lp$all,geom="tile") +
  theme(legend.position = "bottom") +
  geom_sf(data=kruger, alpha = 0) +
  icsc +
  ggtitle("Predicted Impala Distribution")


#### Extract predicted impala map ####
coords <- st_coordinates(ifit6b.lp$allog)
pred_vals <- ifit6b.lp$allog$mean 
df_raster <- data.frame(x = coords[,1], y = coords[,2], value = pred_vals)
r <- rasterFromXYZ(df_raster)
crs(r) <- crs(ndvi)
plot(r)
res(r)
impala <- rast(r)
impala <- terra::resample(impala, ndvi, method = "bilinear")
# Preparing data for analysis
gcov$impala<- impala

f.impala <- function(where) {
  # Extract the values
  v <- eval_spatial(impala, where, layer = "impala")
  # Fill in missing values; this example would work for
  # SpatialPixelsDataFrame data
  if (any(is.na(v))) {
    v <- bru_fill_missing(impala, where, v)
  }
  v
}

#### Posterior Uncertainty for Impala ####
iLambda <- predict(
  ifit6b,
  fm_int(mesh, roadsampler),
  ~ sum(weight * exp(mySmooth + ndvi_veg3 + ndvi_veg5 + vegetation))
)
iLambda

itotal <- predict(
  ifit6b,
  fm_int(mesh, roadsampler),
  ~ data.frame(
    N = 150:350,
    dpois(150:350,
          lambda = sum(weight * exp(mySmooth + ndvi_veg3 + ndvi_veg5 + vegetation))
    )
  )
)
inla.qmarginal(c(0.025, 0.5, 0.975), marginal = list(x = 150:350, y = itotal$mean))
inla.emarginal(identity, marginal = list(x = 150:350, y = itotal$mean))
itotal$plugin_estimate <- dpois(150:350, lambda = iLambda$mean)

ggplot(itotal) +
  geom_line(aes(x = N, y = mean, colour = "Posterior")) +
  geom_line(aes(x = N, y = plugin_estimate, colour = "Plugin"))

 
  
#### Best Models ####
### Lion ###
#fit4
summary(fit4)
fit4.lp

csc <- scale_fill_viridis_c(option = "D", limits = allrange)
csclog<- scale_fill_viridis_c(option = "D", limits = allogrange)
ggplot() +
  gg(fit4.lp$all, geom = "tile") +
  theme(legend.position = "bottom") +
  geom_sf(data = kruger, alpha = 0) +
  csc +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()
  )

## Lion Occupancy
lionocurrence<-fit4.lp
lionocurrence$all$mean<-1-exp(-fit4.lp$all$mean)
ggplot() +
  gg(lionocurrence$all,geom="tile") +
  geom_sf(data=kruger, alpha = 0) +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()
  )

### Leopard ###
summary(lfit14)
lfit14.lp

allrange <- range(lfit14.lp$all$mean)
allogrange<- range(lfit14.lp$allog$mean)
csc <- scale_fill_viridis_c(option = "D", limits = allrange)
csclog<- scale_fill_viridis_c(option = "D", limits = allogrange)
ggplot() +
  gg(lfit14.lp$all, geom = "tile") +
  theme(legend.position = "bottom") +
  geom_sf(data = kruger, alpha = 0) +
  csc +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()
  )

## Leopard Occupancy
leopardocurrence<-lfit14.lp
leopardocurrence$all$mean<-1-exp(-lfit14.lp$all$mean)
ggplot() +
  gg(leopardocurrence$all,geom="tile") +
  geom_sf(data=kruger, alpha = 0) +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()
  )

### Cheetah ###
summary(cfit35)
cfit35.lp

allrange <- range(cfit35.lp$all$mean)
allogrange<- range(cfit35.lp$allog$mean)
csc <- scale_fill_viridis_c(option = "D", limits = allrange)
csclog<- scale_fill_viridis_c(option = "D", limits = allogrange)
ggplot() +
  gg(cfit35.lp$all, geom = "tile") +
  theme(legend.position = "bottom") +
  geom_sf(data = kruger, alpha = 0) +
  csc +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()  
  )

## Cheetah Occupancy
cheetahocurrence<-cfit35.lp
cheetahocurrence$all$mean<-1-exp(-cfit35.lp$all$mean)
ggplot() +
  gg(cheetahocurrence$all,geom="tile") +
  geom_sf(data=kruger, alpha = 0) +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()
  )

### Hyena ###
#hfit43
summary(hfit43)
hfit43.lp 

allrange <- range(hfit43.lp$all$mean)
allogrange<- range(hfit43.lp$allog$mean)
csc <- scale_fill_viridis_c(option = "D", limits = allrange)
csclog<- scale_fill_viridis_c(option = "D", limits = allogrange)
ggplot() +
  gg(hfit43.lp$all, geom = "tile") +
  theme(legend.position = "bottom") +
  geom_sf(data = kruger, alpha = 0) +
  csc +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()  
  )

hyenaocurrence<-hfit43.lp
hyenaocurrence$all$mean<-1-exp(-hfit43.lp$all$mean)
ggplot() +
  gg(hyenaocurrence$all,geom="tile") +
  geom_sf(data=kruger, alpha = 0) +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()
  )


### Impala ###
summary(ifit6b)
ifit6b.lp

allrange <- range(ifit6b.lp$all$mean)
allogrange<- range(ifit6b.lp$allog$mean)
csc <- scale_fill_viridis_c(option = "D", limits = allrange)
csclog<- scale_fill_viridis_c(option = "D", limits = allogrange)
ggplot() +
  gg(ifit6b.lp$all, geom = "tile") +
  theme(legend.position = "bottom") +
  geom_sf(data = kruger, alpha = 0) +
  csc +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()
  )

## Impala Occupancy
impalaocurrence<-ifit6b.lp
impalaocurrence$all$mean<-1-exp(-ifit6b.lp$all$mean)
ggplot() +
  gg(impalaocurrence$all,geom="tile") +
  geom_sf(data=kruger, alpha = 0) +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank()
  )

#### Combine all into one map ####
specieslist <- list(
  leopard = leopardocurrence$all,
  hyena = hyenaocurrence$all,
  lion = lionocurrence$all,
  cheetah = cheetahocurrence$all
  #impala = impalaocurrence$all
)

# Create a single data frame with coordinates and probabilities
combined_df <- bind_rows(
  lapply(names(specieslist), function(sp) {
    pred <- specieslist[[sp]]
    
    coords <- st_coordinates(pred)
    prob<- pred$mean
    
    data.frame(
      x = coords[, 1],
      y = coords[, 2],
      prob = prob,
      species = sp
    )
  })
)

# Plot combined map
ggplot() +
  geom_sf(data = krugerkm, fill = "midnightblue", color = "white", size = 1) +
  geom_tile(data = combined_df, aes(x = x, y = y, fill = species, alpha = prob), color = NA) +
  scale_fill_manual(values = c(
    leopard = "#FF00FF",  # you can pick bright colors per species here
    hyena = "#00FFFF",
    lion = "#FF4500",
    cheetah = "#FFFA3C"
  )) +
  scale_alpha(range = c(0, 1)) +  # min transparency almost invisible, max fully opaque
  coord_sf() +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    alpha = guide_legend(title = "Occurrence Probability"),
    fill = guide_legend(title = "Species")
  ) +
  labs(
    x = NULL, y = NULL
  )
#colours for reference
"#0B0B3B"
"midnightblue"
"#003153"
"#2B2D42"
"#264653"

#Facet map
ggplot(combined_df, aes(x = x, y = y, fill = prob)) +
  geom_tile() +
  facet_wrap(~ species,nrow=1) +
  scale_fill_viridis_c(option = "C") +
  coord_equal() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = NULL, y = NULL
  )






#### Predicting kill hotspots: hyenas and lions ####
# Create raster
#Extract lion map
coords <- st_coordinates(fit4.lp$all)
pred_vals <- fit4.lp$all$mean
df_raster <- data.frame(x = coords[,1], y = coords[,2], value = pred_vals)
li <- rasterFromXYZ(df_raster)
crs(li) <- crs(ndvi)
plot(li)
res(li)
lion <- rast(li)
lion <- terra::resample(lion, ndvi, method = "bilinear")

#Extract hyena map
hcoords <- st_coordinates(hfit43.lp$all)
hpred_vals <- hfit43.lp$all$mean
hdf_raster <- data.frame(x = hcoords[,1], y = hcoords[,2], value = hpred_vals)
h <- rasterFromXYZ(hdf_raster)
crs(h) <- crs(ndvi)
plot(h)
res(h)
hyena <- rast(h)
hyena <- terra::resample(hyena, ndvi, method = "bilinear")

# Normalise
max(lion)
max(hyena)
lion_rast <- lion / 0.549154639
hyena_rast <- hyena / 0.961865067

# Create overlap raster
lionxhyena<- lion_rast * hyena_rast

# Set threshold
threshold <- quantile(values(lionxhyena), probs = 0.90, na.rm = TRUE)
kill<- lionxhyena > threshold
killpoly<- st_as_stars(kill, values=TRUE)

#Check data for kill points
# Filter for any variation starting with "kill"
kills_data <- kruger2024 %>%
  filter(
    str_detect(action, regex("kill", ignore_case = TRUE)) |
      str_detect(tingDescription, regex("kill", ignore_case = TRUE))
  )
kills_data<-st_transform(kills_data, "+proj=utm +zone=35 +south +units=km +datum=WGS84")

# plot it
tmap_mode("plot")  # or "view" for interactive map

tm_shape(kill) +
  tm_raster(title = "Kill Hotspot", palette = "-inferno", alpha = 0.5) +
  tm_shape(kills_data) +
  tm_dots(col = "red", size = 0.1, title = "Reported Kills", legend.show=FALSE)

#Convert to polygons
killpoly<- as.polygons(kill)
killsf<- st_as_sf(killpoly)
killsfiltered<- killsf[killsf$value==1]

#Check for overlap
killinhotspot<- st_filter(kills_data,killsfiltered, .predicate = st_intersects)
obsprop<- nrow(killinhotspot)/nrow(kills_data)
plot(killsfiltered)

#now for random points generated
randompts<- st_sample(roadbuffer,size=nrow(kills_data), type = "random")
randompts<- st_as_sf(randompts)
randominhotspot<- st_filter(randompts,killsfiltered, .predicate = st_intersects)
exprop<- nrow(randominhotspot)/nrow(randompts)

#simulate replicates
set.seed(123)
repprops <- replicate(100, {
  randpts <- st_sample(roadbuffer, size = nrow(kills_data))
  randptssf <- st_as_sf(randpts)
  nrow(st_filter(randptssf, killsfiltered, .predicate = st_intersects)) / nrow(kills_data)
})

#randomisation test
p_value <- mean(repprops >= obsprop)

#Plot for visualization
hist(repprops, main = "Random Proportions of Points in Hotspot", xlab = "Proportion")
abline(v = obsprop, col = "red", lwd = 2)

effect_size <- obsprop-mean(repprops)

#CIs
ci_lower<- quantile(repprops,0.025)
ci_upper <- quantile(repprops,0.975)

#graph
xmin<-min(c(repprops,obsprop))
xmax<-max(c(repprops,obsprop))
hist(repprops, breaks=20, xlab="Proportion", col="lightgrey", xlim=c(0.92,1), ylim=c(0,19.5),main="")
abline(v=obsprop, col="firebrick",lwd = 2)
abline(v=ci_lower, col= "black", lty=2)
abline(v=ci_upper, col= "black", lty=2)


#### Calculate overlap between species (BA) ####
#Wrangle data and get the normalised prob
df_norm <- combined_df %>%
  group_by(species) %>%
  mutate(prob_norm = prob / sum(prob, na.rm = TRUE)) %>%
  ungroup()

df_wide <- df_norm %>%
  dplyr::select(x, y, species, prob_norm) %>%
  pivot_wider(names_from = species, values_from = prob_norm, values_fill = 0)


# Only include grids where road intersects
# average road width is 11m, we do 20m
roadbuffer <- st_buffer(roadsampler, dist = 1)
pts <- st_as_sf(df_wide, coords = c("x", "y"), crs = st_crs(roadbuffer))
mesh_on_roads <- st_filter(pts, roadbuffer, .predicate = st_intersects)
mesh_on_roads_df <- cbind(
  st_coordinates(mesh_on_roads),
  st_drop_geometry(mesh_on_roads)
)


#Lion & Hyena
BAlionxhyena2 <- sum(sqrt(mesh_on_roads_df$lion * mesh_on_roads_df$hyena), na.rm = TRUE)
#Lion & Cheetah
BAlionxcheetah2 <- sum(sqrt(mesh_on_roads_df$lion * mesh_on_roads_df$cheetah), na.rm = TRUE)
#Lion & Leopard
BAlionxleopard2 <- sum(sqrt(mesh_on_roads_df$lion * mesh_on_roads_df$leopard), na.rm = TRUE)
#Hyena & Cheetah
BAhyenaxcheetah2 <- sum(sqrt(mesh_on_roads_df$hyena * mesh_on_roads_df$cheetah), na.rm = TRUE)
#Hyena & Leopard
BAhyenaxleopard2 <- sum(sqrt(mesh_on_roads_df$hyena * mesh_on_roads_df$leopard), na.rm = TRUE)
#Cheetah & Leopard
BAcheetahxleopard2 <- sum(sqrt(mesh_on_roads_df$cheetah * mesh_on_roads_df$leopard), na.rm = TRUE)

# bootstrap to get CIs
calc_BA <- function(df, sp1, sp2) {
  sum(sqrt(df[[sp1]] * df[[sp2]]), na.rm = TRUE)
}

set.seed(123)  # for reproducibility
n_boot <- 1000  # number of bootstrap samples

# Prepare a data frame for bootstrap results
boot_results <- data.frame(
  BAlionxhyena = numeric(n_boot),
  BAlionxcheetah = numeric(n_boot),
  BAlionxleopard = numeric(n_boot),
  BAhyenaxcheetah = numeric(n_boot),
  BAhyenaxleopard = numeric(n_boot),
  BAcheetahxleopard = numeric(n_boot)
)

for (i in 1:n_boot) {
  # Sample rows with replacement from df_wide
  boot_sample <- mesh_on_roads_df %>% sample_n(size = nrow(mesh_on_roads_df), replace = TRUE)
  
  # Calculate BAs for this bootstrap sample
  boot_results$BAlionxhyena[i] <- calc_BA(boot_sample, "lion", "hyena")
  boot_results$BAlionxcheetah[i] <- calc_BA(boot_sample, "lion", "cheetah")
  boot_results$BAlionxleopard[i] <- calc_BA(boot_sample, "lion", "leopard")
  boot_results$BAhyenaxcheetah[i] <- calc_BA(boot_sample, "hyena", "cheetah")
  boot_results$BAhyenaxleopard[i] <- calc_BA(boot_sample, "hyena", "leopard")
  boot_results$BAcheetahxleopard[i] <- calc_BA(boot_sample, "cheetah", "leopard")
}

# Calculate 95% confidence intervals for each pairwise BA
ci_lower <- apply(boot_results, 2, quantile, probs = 0.025)
ci_upper <- apply(boot_results, 2, quantile, probs = 0.975)

# Combine summary stats
boot_summary <- data.frame(
  mean = colMeans(boot_results),
  lower_95 = ci_lower,
  upper_95 = ci_upper
)


# Convert to data frame
boot_df <- as.data.frame(boot_summary)
boot_df$comparison <- rownames(boot_df)

# Plot
ggplot(boot_df, aes(x = mean, y = reorder(comparison, mean))) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0.2, color = "firebrick") +
  labs(
    x = "Affinity",
    y = "Species Pair"
  ) +
  theme_bw()


#### Generate Figures ####
#### Numerical Integration Figure ####
# Define the normal distribution function (density)
f <- function(x) dnorm(x, mean = 0, sd = 1)  # standard normal distribution

# Define the range of x values
x_vals <- seq(-4, 4, length.out = 100)

# Compute the density values for the x_vals
y_vals <- f(x_vals)

# Plot the normal distribution curve
plot(x_vals, y_vals, type = "l", col = "black", lwd = 2, 
     xlab = "x", ylab = "Density", main = "")

# Rectangle Rule (Midpoint) Numerical Integration
rectangle_integral <- function(f, a, b, n) {
  # Compute step size
  h <- (b - a) / n
  # Midpoints of each subinterval
  x_mid <- seq(a + h/2, b - h/2, by = h)
  # Function values at the midpoints
  y_mid <- f(x_mid)
  # Sum up the areas of the rectangles
  integral <- h * sum(y_mid)
  return(integral)
}

# Define the range for integration
a <- -3  # left bound
b <- 3   # right bound
n <- 18   # number of subintervals

# Integrate using the Midpoint Rectangle Rule
integral_value <- rectangle_integral(f, a, b, n)
cat("Rectangle Rule (Midpoint) Approximation for Normal Distribution:", integral_value, "\n")

# Overlay the rectangles used by the Rectangle Rule
x_rect <- seq(a, b, length.out = n + 1)
y_rect <- f(x_rect)

# Draw the rectangles under the curve, now using the midpoints
for (i in 1:n) {
  # Midpoint of the current interval
  x_rect_mid <- (x_rect[i] + x_rect[i + 1]) / 2
  # Height of the rectangle at the midpoint
  y_rect_mid <- f(x_rect_mid)
  
  # Draw the rectangle with the correct height at the midpoint
  polygon(c(x_rect[i], x_rect[i], x_rect[i + 1], x_rect[i + 1]), 
          c(0, y_rect_mid, y_rect_mid, 0), col = rgb(0.2, 0.5, 1, 0.4), border = "steelblue")
}

#### Table for report (Naive intensity) ####

keywords <- c("lion", "elephant", "cheetah", "impala", 
              "giraffe", "buffalo", "hippo", "zebra", 
              "leopard", "wild dog","hyena","kudu")

top10_wildlife <- lapply(keywords, function(pattern) {
  data.frame(
    animal = pattern,
    n = nrow(filter(kruger2024, str_detect(animal, regex(pattern, ignore_case = TRUE))))
  )
}) %>%
  bind_rows() %>%
  slice_max(n, n = 10) %>%
  arrange(desc(n))

colnames(top10_wildlife)<- c("Animal","Count")

top10_wildlife %>%
  kbl(caption = "") %>%
  kable_classic(full_width = F, html_font = "Cambria")

#Now table for selected wildlife
totallion<- nrow(lion2024km)
totalleopard<- nrow(leopard2024km)
totalhyena<- nrow(hyena2024km)
totalcheetah<- nrow(cheetah2024km)
totalimpala<- nrow(impala2024km)
#Calculate naive intensity
totalroadlength<- sum(st_length(roadsampler))
nlion<- nrow(cleanedlion2024km)
nleopard<- nrow(cleanedleopard2024km)
nhyena<- nrow(cleanedhyena2024km)
ncheetah<- nrow(cleanedcheetah2024km)
nimpala<- nrow(cleanedimpala2024km)

intensitylion<- nlion/totalroadlength
intensityleopard<- nleopard/totalroadlength
intensityhyena<- nhyena/totalroadlength
intensitycheetah<- ncheetah/totalroadlength
intensityimpala<- nimpala/totalroadlength

selected_wildlife %>%
  kbl(caption = "") %>%
  kable_classic(full_width = F, html_font = "Cambria")


naive_table <- tibble::tibble(
  Animal = c("Lion", "Leopard", "Hyena","Cheetah","Impala"),
  `Total Sightings`= c(8541,3748,2881,1118,283),
  `Processed Sightings`= c(nlion, nleopard, nhyena, ncheetah, nimpala),
  `Naive intensity` = c(intensitylion, intensityleopard, intensityhyena, intensitycheetah, intensityimpala)
)

naive_table %>%
  kbl(caption = "") %>%
  kable_classic(full_width = F, html_font = "Cambria")

#### Model Evaluation (canned) ####
#K or L function
# For Lions: normal K function
cleanedlion_sp <- as_Spatial(cleanedlion2024km)
krugerkm_sp <- as_Spatial(krugerkm)
lp <- spatial.to.ppp(cleanedlion_sp, krugerkm_sp)
class(lp)
plot(lp, main="lions")

ke.lion <- Kest(lp)
plot(ke.lion)
env <- envelope(lp, Kest, nsim = 19, simulate = expression(rpoispp(lambda = intensity(lp))))
plot(env, main = "")

vrt<-fm_vertices(mesh)

# Model-based predictions on linear predictor scale
vint1 <- predict(fit16, vrt, formula = ~ exp(waterways + ndvi + waterways_x_ndvi + Intercept))
vint2 <- predict(fit16, vrt, formula = ~ exp(mySmooth + Intercept))
vint3 <- predict(fit16, vrt, formula = ~ exp(mySmooth + waterways + ndvi + waterways_x_ndvi + Intercept))
kruger_sp <- as(krugerkm, "Spatial")
#get road
roads_owin <- as.owin(road_buffer)

scaled_intensity1 <- vint1$mean * 20
scaled_intensity3 <- vint1$mean * 20
testsamples1 <- lapply(1:5, function(i) { spatial.to.ppp(sample.lgcp(mesh, log(scaled_intensity1), samplers = krugerkm), samplers = kruger_sp)})
samples1_filtered <- lapply(testsamples1, filter_to_roads)


testsamples3 <- lapply(1:5, function(i) { spatial.to.ppp(sample.lgcp(mesh, log(scaled_intensity3), samplers = krugerkm), samplers = kruger_sp)})

samples1 <- lapply(1:10, function(i) { spatial.to.ppp(sample.lgcp(mesh, log(vint1$mean), samplers = roadbuffer), samplers = kruger_sp)})
samples2 <- lapply(1:20, function(i) { spatial.to.ppp(sample.lgcp(mesh, log(vint2$mean), samplers = krugerkm), samplers = kruger_sp)})
samples3 <- lapply(1:100, function(i) { spatial.to.ppp(sample.lgcp(mesh, log(vint3$mean), samplers = krugerkm), samplers = kruger_sp)})


# Function to filter one ppp sample to points inside roads_owin
filter_to_roads <- function(ppp_sample) {
  inside_idx <- inside.owin(ppp_sample$x, ppp_sample$y, w = roads_owin)
  filtered_ppp <- ppp_sample[inside_idx]
  return(filtered_ppp)
}

# Apply filtering to each sample in samples1 list
samples1_filtered <- lapply(testsamples1, filter_to_roads)
samples3_filtered <- lapply(testsamples3, filter_to_roads)
samples2_filtered <- lapply(samples2, filter_to_roads)


plot(samples1_filtered[[2]], main = "")

env1K <- envelope(lp, Kest, simulate = samples1_filtered, nsim = 5)
plot(env1K)

env3K <- envelope(lp, Kest, simulate = samples3_filtered, nsim = 5)
plot(env3K)

env2K <- envelope(lp, Kest, simulate = samples2_filtered, nsim = 19)
plot(env2K)

env1L <- envelope(lp, Lest, simulate = samples1_filtered,nsim=5)
plot(env1L)
env3L <- envelope(lp, Lest, simulate = samples3_filtered,nsim=5)
plot(env3L)





#### Debugging section ####


#### Lion Graphs ####
w.pred <- predict(
  fit4,
  data.frame(waterways = seq(-1, 6, length.out = 1000)),
  formula = ~ exp(waterways + Intercept)
  # include = character(0) # Not needed from version 2.8.0
)
# recall: scaled = (value-mean)/sd
# to undo the process, original=scaled*sd + mean
w.pred$original<- w.pred$waterways* 10162.41 + 7265.31
ggplot(w.pred) +
  geom_line(aes(original, mean)) +
  geom_ribbon(
    aes(original,
        ymin = q0.025,
        ymax = q0.975
    ),
    alpha = 0.2
  ) +
  geom_ribbon(
    aes(original,
        ymin = mean - 1 * sd,
        ymax = mean + 1 * sd
    ),
    alpha = 0.2
  ) +
  labs(
    x = "Distance from Waterways (km)",
    y = expression(lambda)
  )

#### Leopard Graphs ####
# canopy range
canopy_seq <- seq(-0.5, 9, length.out = 200)

# impala values (low, medium, high)
impala_vals <- quantile(values(impala), probs = c(0.10, 0.5, 0.90), na.rm = TRUE)

# prediction grid
pred_df <- expand.grid(
  canopy = canopy_seq,
  impala = impala_vals
)
pred_df$canopy_x_impala <- pred_df$canopy * pred_df$impala

# Predict from the fitted model
leopard.pred <- predict(
  lfit14,
  pred_df,
  formula = ~ exp(Intercept + canopy + impala + canopy_x_impala)
)

# Add label for impala
leopard.pred$impala_cat <- cut(
  leopard.pred$impala,
  breaks = 3,
  labels = c("Low", "Medium", "High")
)
# recall: scaled = (value-mean)/sd
# to undo the process, original=scaled*sd + mean
m <- global(gcov$canopy, "mean", na.rm=TRUE)[[1]]
s <- global(gcov$canopy, "sd", na.rm=TRUE)[[1]]
leopard.pred$original<- leopard.pred$canopy*3.538917 + 1.127

# Plot with ribbon
ggplot(leopard.pred, aes(x = original, color = impala_cat, fill = impala_cat)) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), alpha = 0.2, color = NA) +
  geom_line(aes(y = mean), size = 1) +
  labs(
    x = "Canopy",
    y = expression(lambda),
    color = "Impala abundance",
    fill = "Impala abundance"
  ) +
  theme_bw()

# Leopard Heatmap
# Mean estimates
beta0 <- -0.155
beta_canopy <- -0.995
beta_impala <- 0.657
beta_int <- -0.185

# Range
canopy_seq <- seq(-0.3185027, 8.4412402, length.out = 200)
impala_seq <- seq(quantile(values(gcov$impala), .01, na.rm = TRUE),
                  quantile(values(gcov$impala), .99, na.rm = TRUE), length.out = 100)
grid <- expand.grid(canopy = canopy_seq, impala = impala_seq)

# plug-in linear predictor
grid$eta <- beta0 + beta_canopy * grid$canopy + beta_impala * grid$impala + beta_int * (grid$canopy * grid$impala)
grid$lambda <- exp(grid$eta)
grid$original <- grid$canopy*3.538917 + 1.127
ggplot(grid, aes(x = original, y = impala, fill = lambda)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log", name = expression(lambda)) +
  labs(x = "Canopy", y = "Impala (log intensity)") +
  theme_bw()

#### Hyena Graph ####
# restcamp graph
# Representative values for other covariates
mean_ndvi   <- 0
mean_impala <- -5.51

# Prediction grid
pred_df <- data.frame(
  restcamp = seq(-1.6, 4.5, length.out = 200),
  ndvi     = mean_ndvi,
  impala   = mean_impala,
  impala_x_ndvi = mean_impala * mean_ndvi
)

# Predict intensity
rh.pred <- predict(
  hfit43,
  pred_df,
  formula = ~ exp(restcamp + ndvi + impala + impala_x_ndvi + Intercept)
)

#
# recall: scaled = (value-mean)/sd
# to undo the process, original=scaled*sd + mean
rh.pred$original<- rh.pred$restcamp* 20791.9 + 31336.18
ggplot(rh.pred) +
  geom_line(aes(original, mean)) +
  geom_ribbon(
    aes(original,
        ymin = q0.025,
        ymax = q0.975
    ),
    alpha = 0.2
  ) +
  geom_ribbon(
    aes(original,
        ymin = mean - 1 * sd,
        ymax = mean + 1 * sd
    ),
    alpha = 0.2
  ) +
  labs(
    x = "Distance from Restcamps (km)",
    y = expression(lambda)
  )

# Impalaxndvi
# Representative impala values
impala_vals <- quantile(values(impala), probs = c(0.10, 0.5, 0.90), na.rm = TRUE)

# NDVI sequence
ndvi_seq <- seq(-2.892598, 4.759562, length.out = 200)

# Fixed restcamp
fixed_restcamp <- 0

# Build prediction data
pred_data <- expand.grid(
  ndvi = ndvi_seq,
  impala = impala_vals,
  restcamp = fixed_restcamp
)

# Compute the interaction
pred_data$impala_x_ndvi <- pred_data$ndvi * pred_data$impala

# Predict
pred <- predict(
  hfit43,
  pred_data,
  formula = ~ exp(restcamp + ndvi + impala + impala_x_ndvi + Intercept)
)

# Label impala levels for plotting
pred$impala_cat <- cut(pred$impala,
                       breaks = 3,
                       labels = c("Low", "Medium", "High"))
# recall: scaled = (value-mean)/sd
# to undo the process, original=scaled*sd + mean
pred$original<- pred$ndvi*0.09120805 + 0.167

# Plot
ggplot(pred, aes(x = original, color = impala_cat, fill = impala_cat)) +
  geom_line(aes(y = mean), size = 1) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), alpha = 0.2, color = NA) +
  labs(x = "NDVI",
       y = expression(lambda),
       color = "Impala abundance",
       fill = "Impala abundance") +
  theme_bw()

#### Cheetah Graph ####
# waterways graph
# Representative values for other covariates
mean_ndvi   <- 0
mean_impala <- -5.51

# Prediction grid
pred_df <- data.frame(
  waterways = seq(-1, 2, length.out = 200),
  ndvi     = mean_ndvi,
  impala   = mean_impala,
  impala_x_ndvi = mean_impala * mean_ndvi
)

# Predict intensity
cw.pred <- predict(
  cfit35,
  pred_df,
  formula = ~ exp(waterways + ndvi + impala + impala_x_ndvi + Intercept)
)

# recall: scaled = (value-mean)/sd
# to undo the process, original=scaled*sd + mean
cw.pred$original<- cw.pred$waterways* 10162.41 + 7265.31
# Plot
ggplot(cw.pred) +
  geom_line(aes(original, mean)) +
  geom_ribbon(
    aes(original,
        ymin = q0.025,
        ymax = q0.975
    ),
    alpha = 0.2
  ) +
  geom_ribbon(
    aes(original,
        ymin = mean - 1 * sd,
        ymax = mean + 1 * sd
    ),
    alpha = 0.2
  ) +
  labs(
    x = "Distance from Waterways (km)",
    y = expression(lambda)
  )

# Cheetah impalaxndvi graph
# Impalaxndvi
# impala values
impala_vals <- quantile(values(impala), probs = c(0.10, 0.5, 0.90), na.rm = TRUE)

# NDVI sequence
ndvi_seq <- seq(-2.892598, 4.759562, length.out = 200)

# Fixed restcamp
fixed_waterways <- 0

# prediction grid
cpred_data <- expand.grid(
  ndvi = ndvi_seq,
  impala = impala_vals,
  waterways = fixed_waterways
)

# interaction term
cpred_data$impala_x_ndvi <- cpred_data$ndvi * cpred_data$impala

# Predict
cpred <- predict(
  cfit35,
  cpred_data,
  formula = ~ exp(waterways + ndvi + impala + impala_x_ndvi + Intercept)
)

# Label impala levels for plotting
cpred$impala_cat <- cut(cpred$impala,
                       breaks = 3,
                       labels = c("Low", "Medium", "High"))

# recall: scaled = (value-mean)/sd
# to undo the process, original=scaled*sd + mean
cpred$original<- cpred$ndvi*0.09120805 + 0.167
# Plot
ggplot(cpred, aes(x = original, color = impala_cat, fill = impala_cat)) +
  geom_line(aes(y = mean), size = 1) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), alpha = 0.2, color = NA) +
  labs(x = "NDVI",
       y = expression(lambda),
       color = "Impala abundance",
       fill = "Impala abundance") +
  theme_bw()


### Impala Graph ####
# Get min and max values for the raster
range_ndvi3 <- global(ndvi_veg3, c("min", "max"), na.rm = TRUE)
range_ndvi5 <- global(ndvi_veg5, c("min", "max"), na.rm = TRUE)

# Define sequences
ndvi_veg3_seq <- seq(range_ndvi3$min, range_ndvi3$max, length.out = 200)
ndvi_veg5_seq <- seq(range_ndvi5$min, range_ndvi5$max, length.out = 200)

# Means of rasters
mean_ndvi_veg3 <- global(ndvi_veg3, "mean", na.rm = TRUE)[[1]]
mean_ndvi_veg5 <- global(ndvi_veg5, "mean", na.rm = TRUE)[[1]]


# Extract fixed effects
fixefs_mean  <- ifit6b$summary.fixed$mean
fixefs_lower <- ifit6b$summary.fixed$`0.025quant`
fixefs_upper <- ifit6b$summary.fixed$`0.975quant`

# Extract vegetation random effects
veg_re_mean  <- ifit6b$summary.random$vegetation$mean
veg_re_lower <- ifit6b$summary.random$vegetation$`0.025quant`
veg_re_upper <- ifit6b$summary.random$vegetation$`0.975quant`

# Function to compute linear predictor and credible intervals
compute_preds <- function(ndvi_veg3, ndvi_veg5, veg_index) {
  lin_pred_mean <- fixefs_mean[1] * ndvi_veg3 +
    fixefs_mean[2] * ndvi_veg5 +
    veg_re_mean[veg_index]
  
  lin_pred_lower <- fixefs_lower[1] * ndvi_veg3 +
    fixefs_lower[2] * ndvi_veg5 +
    veg_re_lower[veg_index]
  
  lin_pred_upper <- fixefs_upper[1] * ndvi_veg3 +
    fixefs_upper[2] * ndvi_veg5 +
    veg_re_upper[veg_index]
  
  data.frame(
    mean = lin_pred_mean,
    lower = lin_pred_lower,
    upper = lin_pred_upper
  )
}

# Prediction
# Lowveld (3)
df_veg3 <- data.frame(
  vegetation = factor("Lowveld", levels = c("Lowveld", "Mopane")),
  ndvi_veg3 = ndvi_veg3_seq,
  ndvi_veg5 = mean_ndvi_veg5
)
preds_veg3 <- compute_preds(df_veg3$ndvi_veg3, df_veg3$ndvi_veg5, veg_index = 1)
df_veg3 <- cbind(df_veg3, preds_veg3)
df_veg3$density <- exp(df_veg3$mean)
df_veg3$density_lower <- exp(df_veg3$lower)
df_veg3$density_upper <- exp(df_veg3$upper)
df_veg3$ndvi <- df_veg3$ndvi_veg3

# Mopane (5)
df_veg5 <- data.frame(
  vegetation = factor("Mopane", levels = c("Lowveld", "Mopane")),
  ndvi_veg3 = mean_ndvi_veg3,
  ndvi_veg5 = ndvi_veg5_seq
)
preds_veg5 <- compute_preds(df_veg5$ndvi_veg3, df_veg5$ndvi_veg5, veg_index = 2)
df_veg5 <- cbind(df_veg5, preds_veg5)
df_veg5$density <- exp(df_veg5$mean)
df_veg5$density_lower <- exp(df_veg5$lower)
df_veg5$density_upper <- exp(df_veg5$upper)
df_veg5$ndvi <- df_veg5$ndvi_veg5

# Plot
newdata <- rbind(df_veg3, df_veg5)
newdata$original<- newdata$ndvi*0.09120805 + 0.167
ggplot(newdata, aes(x = original, y = density, color = vegetation, fill = vegetation)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = density_lower, ymax = density_upper), alpha = 0.2, color = NA) +
  labs(
    x = "NDVI",
    y = expression(lambda),
    color = "Vegetation",
    fill = "Vegetation"
  ) +
  theme_bw()



