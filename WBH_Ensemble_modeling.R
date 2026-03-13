.....### Codes for White-bellied Heron habitat suitability mapping...................


.............### Install all the required packages....................................
install.packages("biomod2")
### all the packages below were installed prior to loading...........................


...........### Load the packages.....................................................


library(ggplot2)
library(raster)
library(sp)
library(rasterVis)
library(tidyterra)
library(randomForest)
library(xgboost)
library(maxnet)
library(biomod2)
library(terra)


..................### Load the data...............................................

ArdeaInsignis_data = read.csv("Ardea insignis_data.csv")

### Load and stack BioClimatic Variables
Bio_var = raster::stack(
  c(
    Bio_1 = "split_1.tif",
    Bio_3 = "split_3.tif",
    Bio_12 = "split_12.tif",
    Bio_15 = "split_15.tif",
    Bio_19 = "split_19.tif"
  )
)
Bio_var
plot(Bio_var)


................### format the data...............................................

wbh_data <- 
  BIOMOD_FormatingData(
    resp.var = ArdeaInsignis_data['Presence'],
    resp.xy = ArdeaInsignis_data[, c('Longitude', 'Latitude')],
    expl.var = Bio_var,
    resp.name = "Ardea. insignis",
    PA.nb.rep = 3,
    PA.nb.absences = 500,
    PA.strategy = 'random',
    filter.raster = TRUE
  )
wbh_data

plot(wbh_data)

.............## run the individual models..........................................

wbh_models <- BIOMOD_Modeling(bm.format = wbh_data,
                              modeling.id = 'AllModels',
                              models = c('RF', 'GLM', 'SRE', 'MAXNET', 
                                         'GBM', 'XGBOOST', 'ANN'),
                              CV.strategy = 'random',
                              CV.nb.rep = 3,
                              CV.perc = 0.8,
                              OPT.strategy = 'bigboss',
                              metric.eval = c('TSS','ROC'),
                              var.import = 2,
                              seed.val = 42)


wbh_models
summary(wbh_models)

length(get_built_models(wbh_models))
table(sub("_.*", "", get_built_models(wbh_models)))



........# Get evaluation scores & variables importance............................

get_evaluations(wbh_models)
modelval = get_evaluations(wbh_models)
write.csv(modelval, file = "modeleval.csv")

modelval_df <- as.data.frame(eval)



get_variables_importance(wbh_models)

var_imp = get_variables_importance(wbh_models)

..............# Represent evaluation scores........................................... 

bm_PlotEvalMean(bm.out = wbh_models, dataset = 'calibration')
bm_PlotEvalMean(bm.out = wbh_models, dataset = 'validation')
modelevaltable = bm_PlotEvalMean(bm.out = wbh_models, dataset = 'validation')
modelevaltable_df = as.data.frame(modelevaltable)
bm_PlotEvalBoxplot(bm.out = wbh_models, group.by = c('algo', 'run'))

...............# # Represent variables importance..........................................

bm_PlotVarImpBoxplot(bm.out = wbh_models, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = wbh_models, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = wbh_models, group.by = c('algo', 'expl.var', 'run'))

.........................# Represent response curves...................................

mods <- get_built_models(wbh_models, run = 'RUN1')
bm_PlotResponseCurves(bm.out = wbh_models, 
                      models.chosen = mods,
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = wbh_models, 
                      models.chosen = mods,
                      fixed.var = 'min')
mods <- get_built_models(wbh_models, full.name = 'Ardea..insignis_PA1_allRun_RF')
bm_PlotResponseCurves(bm.out = wbh_models, 
                      models.chosen = mods,
                      fixed.var = 'median',
                      do.bivariate = TRUE)


.................### Run the Ensemble Model.........................................

myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = wbh_models,
                                      models.chosen = 'all',
                                      em.by = 'algo',
                                      em.algo = c('EMmean', 'EMca'),
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      seed.val = 42)
myBiomodEM
summary(myBiomodEM)

.....................#Get evaluation scores & variables importance....................

get_evaluations(myBiomodEM)
Emodelval = get_evaluations(myBiomodEM)
write.csv(Emodelval, file = "Emodelval.csv")
get_variables_importance(myBiomodEM)
Emodelvar_imp = get_variables_importance(myBiomodEM)

..................# Represent evaluation scores................................

bm_PlotEvalMean(bm.out = myBiomodEM, dataset = 'calibration')

bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'algo'))

......................## Represent variables importance.................................

bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.PA'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.PA'))

.....................## Represent response curves..............................................

bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM),
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM, algo = 'EMmean'),
                      fixed.var = 'median',
                      do.bivariate = TRUE)


.............### Project Ensemble Models (Current Conditions)....................

Current_Ensemble_projection <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  proj.name = 'CurrentProjection',
  new.env = Bio_var,
  models.chosen = 'all',
  metric.binary = 'all',
  metric.filter = 'all',
  build.clamping.mask = TRUE
)

...............### Visualize the Current Projection................................

plot(Current_Ensemble_projection)




.................## Future Projection for the year 2050.................................................

........................### RCP4.5......................................................

Bio_2050_rcp_4.5 = raster::stack(
  c(
    Bio_1  = "Bio_1_SSP2-4.5.tif",
    Bio_3  = "Bio_3_SSP2-4.5.tif",
    Bio_12 = "Bio_12_SSP2-4.5.tif",
    Bio_15 = "Bio_15_SSP2-4.5.tif",
    Bio_19 = "Bio_19_SSP2-4.5.tif"
  )
)

Bio_2050_rcp_4.5
plot(Bio_2050_rcp_4.5)

FutureEM_projection_2050_rcp_4.5 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  proj.name = 'FutureProjection_4.5',
  new.env = Bio_2050_rcp_4.5,
  models.chosen = 'all',
  metric.binary = 'all',
  metric.filter = 'all',
  build.clamping.mask = TRUE
)

plot(FutureEM_projection_2050_rcp_4.5)


.....................#### RCP8.5 2050.......................................


Bio_2050_rcp_8.5 = raster::stack(
  c(
    Bio_1  = "Bio_1_SSP5-8.5.tif",
    Bio_3  = "Bio_3_SSP5-8.5.tif",
    Bio_12 = "Bio_12_SSP5-8.5.tif",
    Bio_15 = "Bio_15_SSP5-8.5.tif",
    Bio_19 = "Bio_19_SSP5-8.5.tif"
  )
)


Bio_2050_rcp_8.5

plot(Bio_2050_rcp_8.5)


FutureEM_projection_2050_rcp_8.5 <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  proj.name = 'FutureProjection_8.5',
  new.env = Bio_2050_rcp_8.5,
  models.chosen = 'all',
  metric.binary = 'all',
  metric.filter = 'all',
  build.clamping.mask = TRUE
)

plot(FutureEM_projection_2050_rcp_8.5) 


................## Species Range Shifts...................................................


wbh_bin_proj_current = raster::stack(
  c(
    ca = "proj_CurrentProjection_Ardea..insignis_ensemble_EMca_TSSbin.tif",
    em =  "proj_CurrentProjection_Ardea..insignis_ensemble_EMmean_TSSbin.tif"
  )
)
  


# ---- Future 2050 SSP2-4.5
wbh_bin_proj_2050_45 =  raster::stack(
  c(
    ca = "proj_FutureProjection_4.5_Ardea..insignis_ensemble_EMca_TSSbin.tif",
    em =  "proj_FutureProjection_4.5_Ardea..insignis_ensemble_EMmean_TSSbin.tif"
  )
)

# ---- Future 2050 SSP5-8.5
wbh_bin_proj_2050_85 = raster::stack(
  c(
    ca = "proj_FutureProjection_8.5_Ardea..insignis_ensemble_EMca_TSSbin.tif",
    em =  "proj_FutureProjection_4.5_Ardea..insignis_ensemble_EMmean_TSSbin.tif"
  )
)
  

..............### Aligning parameters.....................................................

crs(wbh_bin_proj_current)
res(wbh_bin_proj_current)
crs(wbh_bin_proj_current)
res(wbh_bin_proj_current)
# Dimensions (rows, cols, layers)
dim(wbh_bin_proj_current)
dim(wbh_bin_proj_2050_45)

# Extent (this is usually different)
extent(wbh_bin_proj_current)
extent(wbh_bin_proj_2050_45)

# Origin (MOST IMPORTANT & often ignored)
origin(wbh_bin_proj_current)
origin(wbh_bin_proj_2050_45)

wbh_bin_proj_2050_45_aligned <- resample(
  wbh_bin_proj_2050_45,
  wbh_bin_proj_current,
  method = "ngb"  
)

all.equal(dim(wbh_bin_proj_current), dim(wbh_bin_proj_2050_45_aligned))
all.equal(res(wbh_bin_proj_current), res(wbh_bin_proj_2050_45_aligned))
all.equal(extent(wbh_bin_proj_current), extent(wbh_bin_proj_2050_45_aligned))
all.equal(origin(wbh_bin_proj_current), origin(wbh_bin_proj_2050_45_aligned))


..............### Species Range Change for the year 2050 SSP2-4.5..........................

SRC_current_2050_45 <- BIOMOD_RangeSize(
  wbh_bin_proj_current,
  wbh_bin_proj_2050_45_aligned)
  
SRC_current_2050_45$Compt.By.Models
a_current_2050_45 = SRC_current_2050_45$Compt.By.Models
write.csv(a_current_2050_45, file = "a_current_2050_45.csv")

................### Species Range Change for the year 2050 SSP5-8.5.........................

wbh_bin_proj_2050_85_aligned <- resample(
  wbh_bin_proj_2050_85,
  wbh_bin_proj_current,
  method = "ngb"
)

SRC_current_2050_85 <- BIOMOD_RangeSize(
  wbh_bin_proj_current,
  wbh_bin_proj_2050_85_aligned
)

SRC_current_2050_85$Compt.By.Models
a_current_2050_85 = SRC_current_2050_85$Compt.By.Models
write.csv(a_current_2050_85, file = "a_current_2050_85.csv")



wbh_src_map <- c(
  SRC_current_2050_45$Diff.By.Pixel,
  SRC_current_2050_85$Diff.By.Pixel
)

names(wbh_src_map) <- c(
  "EMca | 2050 SSP2-4.5",
  "EMmean | 2050 SSP2-4.5",
  "EMca | 2050 SSP5-8.5",
  "EMmean | 2050 SSP5-8.5"
)

my.at <- seq(-2.5, 1.5, 1)

myColorkey <- list(
  at = my.at,
  labels = list(
    labels = c("lost", "pres", "abs", "gain"),
    at = my.at[-1] - 0.5
  )
)

plot(
  wbh_src_map,
  col = c("#f03b20", "#99d8c9", "#f0f0f0", "#2ca25f"),
  breaks = c(-2, -1, 0, 1, 2),
  main = names(wbh_src_map)
)

myColorkey <- list(
  at = my.at,
  labels = list(
    labels = c("Loss","Present", "Absent", "Gain"),
    at = my.at[-1] - 0.5
  )
)

labels = c("Loss","Present", "Absent", "Gain")
plot(
  wbh_src_map,
  breaks = c(-2, -1, 0, 1, 2),
  col = c("red", "yellow", "orange", "green"),
  main = names(wbh_src_map),
  legend = TRUE,
  plg = list(
    title = "Range change",
    text  = c("Loss", "Present", "Absent", "Gain")
  )
)


