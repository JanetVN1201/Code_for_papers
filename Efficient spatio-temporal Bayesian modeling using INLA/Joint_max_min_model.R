library(raster)
library(ggplot2)
library(terra)
library(tidyterra)
library(viridis)
library(colorspace)
library(geodata)
library(gganimate)
library(sp)
library(sf)

#Add districts to points data
mds <- gadm("ZAF", level = 2, path = "~/Env_book_ch/R" )
mds_poly <- sf::st_as_sf(mds)
dataSA <- readRDS("~/Env_book_ch/R/SA_clim_data_file_filtered.RDS")
crs0 <- st_crs(mds_poly)$proj4string

#Modeling using INLA 
library(INLA)

matern <-
  inla.spde2.pcmatern(mesh_fine,
                      prior.sigma = c(10, 0.5),
                      prior.range = c(18, 0.5)
  )

n.time = 36
iset_maxt <- inla.spde.make.index("imax", n.spde = matern$n.spde,
                             n.group = n.time)

iset_mint <- inla.spde.make.index("imin", n.spde = matern$n.spde,
                              n.group = n.time)

iset_common <- inla.spde.make.index("icommon", n.spde = matern$n.spde,
                                    n.group = n.time)
iset_common_min <- inla.spde.make.index("icommon_min", n.spde = matern$n.spde,
                                      n.group = n.time)

coo <- cbind(dataSA$lat, dataSA$lon)
Am <- inla.spde.make.A(mesh = mesh_fine, loc = coo, group = dataSA$time)

factor_matrix = model.matrix(dataSA$max_t ~ factor(dataSA$year))[,-1]

temp.stack.data_max <- inla.stack(data = list(y = cbind(as.vector(dataSA$max_t), NA)),
                                  A = list(1, Am, Am),
                                  effects = list(data.frame(Intercept_max = rep(1, nrow(dataSA)), 
                                                           year2022_max = factor_matrix[,1], year2023_max = factor_matrix[,2], btime_max = dataSA$time), 
                                                 i_max = iset_maxt, i_common = iset_common),
                                  tag = "temp.data_max")


 temp.stack.data_min <- inla.stack(data = list(y = cbind(NA, as.vector(dataSA$min_t))),
                                   A = list(1, Am, Am),
                                   effects = list(data.frame(Intercept_min = rep(1, nrow(dataSA)), 
                                                             year2022_min = factor_matrix[,1], year2023_min = factor_matrix[,2], btime_min = dataSA$time), 
                                                  i_min = iset_mint, i_common_min = iset_common_min),
                                                     tag = "temp.data_min")

temp.stack <- inla.stack(temp.stack.data_max, temp.stack.data_min)

 temp.f_co <- y ~ -1 + Intercept_max + Intercept_min + 
   year2022_max + year2023_max +
   year2022_min + year2023_min +
   f(imax, model = spde, group = imax.group, 
     control.group = list(model = "ar1")) +
    f(imin, model = spde, group = imin.group, 
      control.group = list(model = "ar1")) +
    f(icommon_min, copy = "imax", fixed = F, group = icommon_min.group)

  res_max_tco <- inla(temp.f_co,
                    data = inla.stack.data(temp.stack, spde = matern),
                    family = rep("gaussian",2),
                    control.predictor = list(A = inla.stack.A(temp.stack), compute = F),
                    control.inla = list(
                      int.strategy = "eb"
                    ), verbose = T)

