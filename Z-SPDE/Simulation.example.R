#SIMULATIONS
library(sp)
library(sf)
library(scales)
library(INLA)
library(INLAspacetime)
library(inlabru)
library(patchwork)
library(fmesher)
library(ggplot2)

##make polygons for first zonal configuration, two barriers with a canal

#geometry
smalldist <- 1.5
W <- c(3)
width <- c(W,W)
n <- c(100)
set.inla.seed <- 880

x.mid = 10
xlim.big  = c(0,20)
xlim.small = c(3,17)
ylim.big = c(0,20)
ylim.small = c(3,17)
#model parameters
range <- 2
sigma <- 1
max.edge.length <- 0.4
prior.range = c(1, 0.5)
prior.sigma = c(1, 0.1)
sigma.u = 1
sigma.epsilon = 0.2

#helpfull functions

#CORRELATION BETWEEN TWO POINTS
#make a square polygon
local.square.polygon <- function(xlim, ylim){
  # - output is a square
  xlim = range(xlim); ylim = range(ylim)
  corner1 = c(xlim[1], ylim[2])
  corner2 = c(xlim[2], ylim[1])
  poly = Polygon(rbind(corner1, c(corner1[1], corner2[2]), corner2, c(corner2[1], corner1[2]), corner1), hole = FALSE)
  return(SpatialPolygons(list(Polygons(list(poly), ID = runif(1)))))
}

#make a square polygon as sf, given xlim/ylim
local_square_sfc <- function(xlim, ylim) {
  xlim <- range(xlim); ylim <- range(ylim)
  xy <- rbind(
    c(xlim[1], ylim[2]),
    c(xlim[1], ylim[1]),
    c(xlim[2], ylim[1]),
    c(xlim[2], ylim[2]),
    c(xlim[1], ylim[2])  # close ring
  )
  st_sfc(st_polygon(list(xy)))
}

#find local correlation between points
localCorrel <- function(locs, mesh, Q) {
  nl <- nrow(locs)
  ii <- sapply(1:nl, function(i)
    which.min(rowSums(sweep(
      mesh$loc[, 1:ncol(locs)], 2, locs[i, ], "-")^2)))
  b <- matrix(0, nrow(Q), nl)
  for(i in 1:nl)
    b[ii[i], i] <- 1
  cc <- inla.qsolve(Q, b)
  s <- sqrt(diag(inla.qinv(Q)))
  for(i in 1:nl)
    cc[, i] <- cc[, i] / (s * s[ii[i]])
  return(drop(cc))
}


local.square.polygon_T = function(xlim, ylim){
  # - output is a square
  xlim = range(xlim); ylim = range(ylim)
  corner1 = c(xlim[1], ylim[2])
  corner2 = c(xlim[2], ylim[1])
  poly = Polygon(rbind(corner1, c(corner1[1], corner2[2]), corner2, c(corner2[1], corner1[2]), corner1), hole = TRUE)
  return(poly)
}

#example

barrier1 <- local.square.polygon(xlim=c(xlim.big[1]-1, x.mid-smalldist/2), 
                              ylim=x.mid+width[1]*c(-.5, .5))

barrier2 <- local.square.polygon(xlim=c(x.mid+smalldist/2, xlim.big[2]+1), 
                              ylim=x.mid+width[2]*c(-.5, .5))

barriers <- SpatialPolygons(c(barrier1@polygons, barrier2@polygons))

loc1 <- matrix(c(xlim.big[1],ylim.big[1], 
                 xlim.big[2],ylim.big[1], 
                 xlim.big[2],ylim.big[2], 
                 xlim.big[1],ylim.big[2]), 4, 2, byrow = T)
loc1.sp <- SpatialPolygons(list(Polygons(list(Polygon(loc1, FALSE)), '0')))

#your barriers (sf)
barrier1_ <- local_square_sfc(
  xlim = c(xlim.big[1] - 1, x.mid - smalldist/2),
  ylim = x.mid + width[1] * c(-.5, .5)
)

barrier2_ <- local_square_sfc(
  xlim = c(x.mid + smalldist/2, xlim.big[2] + 1),
  ylim = x.mid + width[2] * c(-.5, .5)
)

domain.xy_ <- rbind(loc1, loc1[1, ])   # close ring
domain_ <- st_sfc(st_polygon(list(domain.xy_)))

#multipolygon barriers
barriers_ <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier1_[[1]])[, 1:2])),
  st_polygon(list(st_coordinates(barrier2_[[1]])[, 1:2]))
)))

mesh <- fm_mesh_2d(
  loc.domain = domain.xy_, 
  max.edge = max.edge.length,
  offset = 1)

#multipolygon barriers
barriers_ <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier1_[[1]])[, 1:2])),
  st_polygon(list(st_coordinates(barrier2_[[1]])[, 1:2]))
)))

#mesh need all the barriers in one object
barriers_1 <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier1_[[1]])[, 1:2]))
)))

barriers_2 <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier2_[[1]])[, 1:2]))
)))

triBarrier_1and2 <- list(
  unlist(fm_contains(
    x = barriers_1, 
    y = mesh, 
    type = "centroid")),
  unlist(fm_contains(
    x = barriers_2, 
    y = mesh, 
    type = "centroid")))

bfem_1y2 <- mesh2fem.barrier(mesh, triBarrier_1and2)
ranges_true <- range * c(1, 0.01, 0.7)
Q <- inla.barrier.q(bfem_1y2, ranges = ranges_true, sigma = sigma)

#Data simulation
#true (simulated) field
#auxiliary plots first
gg0 <- ggplot() + 
  xlab("") + 
  ylab("") +
  theme_minimal() + 
  coord_fixed() 

gg.add <- list(
  scale_fill_viridis_c(
    option = "turbo",
    na.value = "transparent"
  )
)

#projection grids
bb <- matrix(c(3,17,3,17), nrow=2, ncol=2, byrow = T)

pgrid <- fm_evaluator(
  mesh,
  lattice = fm_evaluator_lattice(
    mesh,
    xlim = bb[1, ],
    ylim = bb[2, ],
    dims = c(300, 300)
  )
)
grid.df <- data.frame(
  x = rep(pgrid$x, times = length(pgrid$y)), 
  y = rep(pgrid$y, each = length(pgrid$x)))

#field
u <- inla.qsample(1, Q, seed = 1)[,1]
ugrid.sim <- fm_evaluate(pgrid, field = u)
grid.df$u <- as.vector(ugrid.sim)

gg0 + 
  geom_raster(
    data = grid.df, 
    aes(x = x, y = y, fill = u)) + 
   gg.add 

#Sampling and model the observations
poly1_h <- local.square.polygon_T(xlim=c(xlim.small[1], x.mid-smalldist/2), 
                                  ylim=x.mid+width[1]*c(-.5, .5))

poly2_h <- local.square.polygon_T(xlim=c(x.mid+smalldist/2, xlim.small[2]), 
                                  ylim=x.mid+width[2]*c(-.5, .5))

loc2 <- matrix(c(3,3, 17,3, 17,17, 3,17), 4, 2, byrow = T)

locp2 <- Polygon(loc2, hole = FALSE)

poly.water <- SpatialPolygons(list(Polygons(list(locp2, poly1_h, poly2_h), '0')))
poly.water_sf <- st_as_sf(poly.water)

set.seed(set.inla.seed)
loc.data <- spsample(x = poly.water, n = n, type = "random")
loc.data_sf <- st_as_sf(loc.data)
loc.data <- loc.data@coords

A.data <- inla.spde.make.A(mesh, loc.data)
u.data <- A.data %*% u

df <- data.frame(loc.data)
names(df) <- c('locx', 'locy')

df$y <- drop(sigma.u * u.data + sigma.epsilon * rnorm(nrow(df)))

#Model

stk <- inla.stack(data=list(y=df$y), 
                  A=list(A.data, 1),
                  effects=list(s=1:mesh$n, 
                               intercept=rep(1, nrow(df))), 
                  remove.unused = FALSE, 
                  tag='est')

zbm <- INLAspacetime::barrierModel.define(
  mesh = mesh, 
  barrier.triangles = triBarrier_1and2,
  prior.range = prior.range, 
  prior.sigma = prior.sigma, 
  range.fraction = c(0.01, 0.7))

formula <- y ~ 0 + intercept + f(s, model = zbm)

res.zbm <- inla(formula, data = inla.stack.data(stk),
                control.predictor = list(A = inla.stack.A(stk)),
                family = 'gaussian',
                control.family = list(hyper = list(
                  prec = list(prior = "pc.prec", fixed = FALSE, param = c(0.2,0.5)))),
                control.mode=list(restart=T, theta=c(3.2, 0.4, 1.6)))

grid.df$u.mean <- as.vector(
  fm_evaluate(
    pgrid,
    res.zbm$summary.random$s$mean))


gg0 + 
  geom_raster(
    data = grid.df, 
    aes(x = x, y = y, fill = u.mean)) +
  gg.add







