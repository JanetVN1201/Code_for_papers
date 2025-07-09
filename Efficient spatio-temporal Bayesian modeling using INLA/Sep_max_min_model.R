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

bnd1 <- inla.nonconvex.hull(matrix(c(unlist(coords_SA[,1:2])), ncol = 2), convex = 0.7)
mesh_fine <- inla.mesh.2d(max.edge = c(0.8,3), boundary = bnd1, offset = c(10,5))

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

coo <- cbind(dataSA$lat, dataSA$lon)
Am <- inla.spde.make.A(mesh = mesh_fine, loc = coo, group = dataSA$time)

factor_matrix = model.matrix(dataSA$max_t ~ factor(dataSA$year))[,-1]

temp.stack.data_maxonly <- inla.stack(data = list(y = dataSA$max_t),
                                  A = list(1, Am),
                                  effects = list(data.frame(Intercept_max = rep(1, nrow(dataSA)), 
                                                            year2022 = factor_matrix[,1], year2023 = factor_matrix[,2]), 
                                                 i_max = iset_maxt),
                                  tag = "temp.data_max")

temp.stack.data_minonly <- inla.stack(data = list(y = dataSA$min_t),
                                   A = list(1, Am),
                                   effects = list(data.frame(Intercept_min = rep(1, nrow(dataSA)), 
                                                             year2022 = factor_matrix[,1], year2023 = factor_matrix[,2], btime = dataSA$time), 
                                                  i_min = iset_mint),
                                                     tag = "temp.data_min")
 
 temp.f <- y ~ -1 + Intercept_max + year2022 + year2023 + #f(btime, model = "seasonal", season.length = 12) +
   f(imax, model = spde, group = imax.group,
     control.group = list(model = "ar1")) 
 
 temp.fmin <- y ~ -1 + Intercept_min + year2022 + year2023 + #f(btime, model = "seasonal", season.length = 12) +
   f(imin, model = spde, group = imin.group,
     control.group = list(model = "ar1")) 
   
res_max_t.pred <- inla(temp.f,
                   data = inla.stack.data(inla.stack(temp.stack.data_maxonly), spde = matern),
                   family = "gaussian",
                   control.predictor = list(A = inla.stack.A(inla.stack(temp.stack.data_maxonly)), compute = T),
                   control.inla = list(
                     int.strategy = "eb"), 
                   control.mode = list(theta = res_max_t$misc$theta.mode, restart = TRUE),
                     verbose = T,
                   control.compute = list(config = T))

res_min_t <- inla(temp.fmin,
                   data = inla.stack.data(inla.stack(temp.stack.data_minonly), spde = matern),
                   family = "gaussian",
                   control.predictor = list(A = inla.stack.A(temp.stack.data_minonly), compute = F),
                   control.inla = list(
                     int.strategy = "eb"
                   ), verbose = T)


