#CORRELATION BETWEEN TWO POINTS
library(sp)
library(sf)
library(scales)
library(INLA)
library(INLAspacetime)
library(inlabru)
library(patchwork)
library(fmesher)
library(ggplot2)
library(cowplot)
#make polygons for first zonal configuration, two barriers with a canal
smalldist = 1
width = c(2, 2)
#width = c(4, 4)
#max.edge.length = 0.2
max.edge.length = 0.4
n = 100
set.inla.seed = 602
x.mid <- 10
xlim.big <- c(0,20)
xlim.small <- c(3,17)
ylim.big <- c(0,20)
ylim.small <- c(3,17)

poly1 <- local.square.polygon(xlim=c(xlim.big[1]-1, x.mid-smalldist/2), 
                              ylim=x.mid+width[1]*c(-.5, .5))

poly2 <- local.square.polygon(xlim=c(x.mid+smalldist/2, xlim.big[2]+1), 
                              ylim=x.mid+width[2]*c(-.5, .5))

poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))

loc1 <- matrix(c(xlim.big[1],ylim.big[1], 
                 xlim.big[2],ylim.big[1], 
                 xlim.big[2],ylim.big[2], 
                 xlim.big[1],ylim.big[2]), 4, 2, byrow = T)

loc.out <- matrix(c(xlim.big[1],ylim.big[1], 
                    xlim.big[2],ylim.big[1], 
                    xlim.big[2],ylim.big[2], 
                    xlim.big[1],ylim.big[2],
                    xlim.big[1],ylim.big[1]), 5, 2, byrow = T)
loc.out.sp <- SpatialPolygons(list(Polygons(list(Polygon(loc.out, FALSE)), '0')))

loc.in <- matrix(c(xlim.small[1], ylim.small[1],
                   xlim.small[2], xlim.small[1],
                   xlim.small[2], ylim.small[2],
                   xlim.small[1], xlim.small[2],
                   xlim.small[1], ylim.small[1]), 5, 2, byrow = T)
loc.in.sp <- SpatialPolygons(list(Polygons(list(Polygon(loc.in, FALSE)), '0')))

plot(loc.out.sp)
plot(loc.in.sp)

loc.out.sf <- st_as_sf(loc.out.sp)
loc.in.sf  <- st_as_sf(loc.in.sp)
poly.out_sf <- st_difference(loc.out.sf, loc.in.sf)
plot(poly.out_sf, col = alpha("skyblue", 0.5))

seg <- inla.sp2segment(poly.original)

#your barriers (sf)
barrier1_ <- local_square_sfc(
  xlim = c(xlim.big[1] - 1, x.mid - smalldist/2),
  ylim = x.mid + width[1] * c(-.5, .5)
)

barrier2_ <- local_square_sfc(
  xlim = c(x.mid + smalldist/2, xlim.big[2] + 1),
  ylim = x.mid + width[2] * c(-.5, .5)
)

#your domain (sf) from loc1
loc1 <- matrix(c(
  xlim.big[1], ylim.big[1],
  xlim.big[2], ylim.big[1],
  xlim.big[2], ylim.big[2],
  xlim.big[1], ylim.big[2]
), 4, 2, byrow = TRUE)

domain.xy_ <- rbind(loc1, loc1[1, ])   # close ring
domain_ <- st_sfc(st_polygon(list(domain.xy_)))

#multipolygon barriers
barriers_ <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier1_[[1]])[, 1:2])),
  st_polygon(list(st_coordinates(barrier2_[[1]])[, 1:2]))
)))

#bounding box helper 
bb <- st_bbox(domain_)
bb_ <- matrix(c(0,20,0,20), nrow=2, ncol=2, byrow = T)
bb2 <- matrix(c(3,17,3,17), nrow=2, ncol=2, byrow=T)

mesh_ <- fm_mesh_2d(
  loc.domain = domain.xy_, 
  max.edge = 0.4,
  offset = 1)
mesh_$n

mesh <-  mesh_

sigma <- 1

#make coordinates list to choose reference points for illustration
width <- 3
y.mid <- 10
y.up <- y.mid + (width/2)
y.low <- y.mid - (width/2)
mesh <- mesh_

location <- matrix(c(c(15, 12.5, 10, 10, 10), rep(y.mid, 5)), ncol = 2)

max.edge.length = 0.4

#point 1, p1:
A <- list()
A[[1]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[1,1],y.up), nrow=1, ncol=2))
A[[2]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[1,1],(y.up + max.edge.length/2)), nrow=1, ncol=2))
A[[3]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[1,1],(y.up + max.edge.length)), nrow=1, ncol=2))

id.nodeA <- c()
id.nodeA[1] = which.max(A[[1]][1, ])
id.nodeA[2] = which.max(A[[2]][1, ])
id.nodeA[3] = which.max(A[[3]][1, ])

#just to check the coordinates
data.frame(x = c(mesh$loc[id.nodeA[1], 1], mesh$loc[id.nodeA[2], 1], mesh$loc[id.nodeA[3], 1]),
           y = c(mesh$loc[id.nodeA[1], 2], mesh$loc[id.nodeA[2], 2], mesh$loc[id.nodeA[3], 2]))

A.y.vector <- c(mesh$loc[id.nodeA[[1]],][2],
                mesh$loc[id.nodeA[[2]],][2],
                mesh$loc[id.nodeA[[3]],][2])

A.tmp <- A[[which(A.y.vector > y.up)[1]]]
id.node <- id.nodeA[which(A.y.vector > y.up)[1]]
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])
return.list.p1 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

#point 2, p2:
A <- list()
A[[1]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[2,1],y.up), nrow=1, ncol=2))
A[[2]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[2,1],(y.up + max.edge.length/2)), nrow=1, ncol=2))
A[[3]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[2,1],(y.up + max.edge.length)), nrow=1, ncol=2))

id.nodeA <- c()
id.nodeA[1] = which.max(A[[1]][1, ])
id.nodeA[2] = which.max(A[[2]][1, ])
id.nodeA[3] = which.max(A[[3]][1, ])

#just to check the coordinates and compare to the plot above
data.frame(x = c(mesh$loc[id.nodeA[1], 1], mesh$loc[id.nodeA[2], 1], mesh$loc[id.nodeA[3], 1]),
           y = c(mesh$loc[id.nodeA[1], 2], mesh$loc[id.nodeA[2], 2], mesh$loc[id.nodeA[3], 2]))

A.y.vector <- c(mesh$loc[id.nodeA[[1]],][2],
                mesh$loc[id.nodeA[[2]],][2],
                mesh$loc[id.nodeA[[3]],][2])

A.tmp <- A[[which(A.y.vector > y.up)[1]]]
id.node <- id.nodeA[which(A.y.vector > y.up)[1]]
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])
return.list.p2 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

#point 3, p3:
A <- list()
A[[1]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[3,1],y.up), nrow=1, ncol=2))
A[[2]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[3,1],(y.up + max.edge.length/2)), nrow=1, ncol=2))
A[[3]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[3,1],(y.up + max.edge.length)), nrow=1, ncol=2))

id.nodeA <- c()
id.nodeA[1] = which.max(A[[1]][1, ])
id.nodeA[2] = which.max(A[[2]][1, ])
id.nodeA[3] = which.max(A[[3]][1, ])

#just to check the coordinates and compare to the plot above
data.frame(x = c(mesh$loc[id.nodeA[1], 1], mesh$loc[id.nodeA[2], 1], mesh$loc[id.nodeA[3], 1]),
           y = c(mesh$loc[id.nodeA[1], 2], mesh$loc[id.nodeA[2], 2], mesh$loc[id.nodeA[3], 2]))

A.y.vector <- c(mesh$loc[id.nodeA[[1]],][2],
                mesh$loc[id.nodeA[[2]],][2],
                mesh$loc[id.nodeA[[3]],][2])

A.tmp <- A[[which(A.y.vector > y.up)[1]]]
id.node <- id.nodeA[which(A.y.vector > y.up)[1]]
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])
return.list.p3 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

#point 4, p4 (in the middle):
A.tmp <- inla.spde.make.A(mesh=mesh,
                          loc = matrix(c(location[4,1],location[4,2]), nrow=1, ncol=2))
id.node = which.max(A.tmp[1, ])
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])

return.list.p4 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

#point 5, p5 (in the middle down):
A <- list()
A[[1]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[5,1],y.low), nrow=1, ncol=2))
A[[2]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[5,1],(y.low - max.edge.length/2)), nrow=1, ncol=2))
A[[3]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[5,1],(y.low - max.edge.length)), nrow=1, ncol=2))

id.nodeA <- c()
id.nodeA[1] = which.max(A[[1]][1, ])
id.nodeA[2] = which.max(A[[2]][1, ])
id.nodeA[3] = which.max(A[[3]][1, ])


A.y.vector <- c(mesh$loc[id.nodeA[[1]],][2],
                mesh$loc[id.nodeA[[2]],][2],
                mesh$loc[id.nodeA[[3]],][2])

A.tmp <- A[[which(A.y.vector < y.low)[1]]]
id.node <- id.nodeA[which(A.y.vector < y.low)[1]]
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])
return.list.p5 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

#list of coordinates for all points
return.list.mesh1 <- list(return.list.p1, return.list.p2, return.list.p3, return.list.p4, return.list.p5)
id.coord.mesh1 <- as.matrix(rbind(return.list.mesh1[[1]]$id.coord, return.list.mesh1[[2]]$id.coord, return.list.mesh1[[3]]$id.coord, return.list.mesh1[[4]]$id.coord, return.list.mesh1[[5]]$id.coord))

#find correlation and make plots
locs_ <- id.coord.mesh1
#now Q, mcorrels gcorrels, ggcorrels for each point and range fraction
range <- 6
rf2 <- c(0.01, 0.2, 0.3, 0.5, 0.7, 0.8, 1)
rf1 <- c(0.01, 0.3)
#row1 <- range*c(1, rf1, rf2)
matrix <- list()
row <- list()

for (u in 1:length(rf1)) {
  row[[u]] <- list()
  for (f in 1:length(rf2)) {
    row[[u]][[f]] <- c(1, rf1[u], rf2[f])
  }
  matrix[[u]] <- matrix(range*unlist(row[[u]]), ncol = 3, byrow = T) 
}

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

triBarrier_1 <- unlist(fm_contains(
  x = barriers_1, 
  y = mesh_, 
  type = "centroid"))

triBarrier_2 <- unlist(fm_contains(
  x = barriers_2, 
  y = mesh_, 
  type = "centroid"))

triBarrier_1y2 <- list(
  unlist(fm_contains(
    x = barriers_1, 
    y = mesh_, 
    type = "centroid")),
  unlist(fm_contains(
    x = barriers_2, 
    y = mesh_, 
    type = "centroid")))

bfem_1y2 <- mesh2fem.barrier(mesh_, triBarrier_1y2)
Q <- list()

for (u in 1:length(rf1)) {
  Q[[u]] <- list()
  for (f in 1:length(rf2)) {
    Q[[u]][[f]] <- inla.barrier.q(bfem_1y2, ranges = matrix[[u]][f,], sigma = sigma)
  }
}

mcorrels <- list()
for (u in 1:length(rf1)) {
  mcorrels[[u]] <- list()
  for (f in 1:length(rf2)) {
    mcorrels[[u]][[f]] <- localCorrel(locs = locs_, mesh = mesh_, Q=Q[[u]][[f]])
  }
}

pgrid_ <- fm_evaluator(
  mesh_,
  lattice = fm_evaluator_lattice(
    mesh_,
    xlim = bb_[1, ],
    ylim = bb_[2, ],
    dims = c(600, 600)
  )
)

grid.df_ <- data.frame(
  x = rep(pgrid_$x, times = length(pgrid_$y)), 
  y = rep(pgrid_$y, each = length(pgrid_$x)))

gcorrels <- list()
for (u in 1:length(rf1)) {
  gcorrels[[u]] <- list()
  for (f in 1:length(rf2)) {
    gcorrels[[u]][[f]] <- as.matrix(fm_evaluate(
      pgrid_, field = mcorrels[[u]][[f]]
    ))
  }
}

ggcorrels <- list()
for (u in 1:length(rf1)) {
  ggcorrels[[u]] <- list()
  for (f in 1:length(rf2)) {
    ggcorrels[[u]][[f]] <- do.call(
      rbind, 
      lapply(1:5, function(l) #number of rows in locs
        data.frame(grid.df_, 
                   loc = paste(sprintf("%1.1f", locs_[l, 1]), 
                               sprintf("%1.1f", locs_[l, 2]), sep = ", "), 
                   correlation = gcorrels[[u]][[f]][, l])))
    
  }
}

gg0 <- ggplot() + 
  xlab("") + 
  ylab("") +
  theme_minimal() + 
  coord_fixed() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank() 
  )

gg.add <- list(
  scale_fill_viridis_c(
    option = "turbo",
    breaks = c(0.2, 0.4, 0.6, 0.8, 1),
    limits = c(0.2, 1),
    na.value = "transparent"
  )
)

#check first that data is alright with plots
gg0 + 
  geom_raster(
    data = ggcorrels[[2]][[4]][ggcorrels[[2]][[4]]$correlation>0.1, ], 
    mapping = aes(x = x, y = y, fill = correlation)) + 
  facet_wrap(~ loc) + 
  gg.add 

#plot with legend only
library(cowplot)

p_leg <- gg0 +
  geom_raster(data = ggcorrels[[1]][[1]], aes(x = x, y = y, fill = correlation)) +
  gg.add +
  theme_void()  # removes axes etc.

leg <- cowplot::get_legend(
  p_leg + theme(legend.position = "right")
)

#save legend as its own figure
if (!dir.exists("plots")) {
  dir.create("plots")
}
png("plots/legend_onlym1.png", width = 600, height = 800, res = 300) 
cowplot::ggdraw(leg)
dev.off()

rf1 <- c(0.01)
ggcorr2points <- list()
for (u in 1:length(rf1)) {
  ggcorr2points[[u]] <- list()
  for (f in 2:length(rf2)) {
    ggcorr2points[[u]][[f]] <- lapply(unique(ggcorrels[[u]][[f]]$loc), function(L) {
      gg0 +
        xlim(bb["xmin"], bb["xmax"]) +
        ylim(bb["ymin"], bb["ymax"]) +
        geom_sf(data = domain_, fill = "white", colour="black") +
        #        geom_sf(data = barriers_, fill = NA, colour="black") +
        geom_raster(
          data = subset(ggcorrels[[u]][[f]], loc == L & correlation > 0.2),
          aes(x = x, y = y, fill = correlation)
        ) +
        geom_sf(data = barriers_, fill = NA, colour="black") +
        geom_sf(data = barriers_1, fill = "white", colour="black") +
        geom_sf(data = poly.out_sf, fill = "white") +
        gg.add 
    })
  }
}


ggcorr2points[[1]][[1]] <- lapply(unique(ggcorrels[[1]][[1]]$loc), function(L) {
      gg0 +
        xlim(bb["xmin"], bb["xmax"]) +
        ylim(bb["ymin"], bb["ymax"]) +
        geom_sf(data = domain_, fill = "white", colour="black") +
        #        geom_sf(data = barriers_, fill = NA, colour="black") +
        geom_raster(
          data = subset(ggcorrels[[1]][[1]], loc == L & correlation > 0.2),
          aes(x = x, y = y, fill = correlation)
        ) +
        geom_sf(data = barriers_, fill = "white", colour="black") +
        geom_sf(data = barriers_1, fill = "white", colour="black") +
        geom_sf(data = poly.out_sf, fill = "white") +
        gg.add 
    })

for (u in 1:length(rf1)) {
  for (f in 1:length(rf2)) {
    for (m in 1:nrow(locs_)) {
      dir <- file.path("plots", paste0("2corr2points_m1point", m))
      if (!dir.exists(dir)) dir.create(dir) 
      fpp <- file.path(dir, (paste0(rf1[u], "_", rf2[f], ".png")))
      png(fpp)
      print( ggcorr2points[[u]][[f]][[m]] )
      dev.off()
    }
  }
}

#make polygons for second zonal configuration, two barriers with no canal
smalldist = 0
width = c(2, 2)
max.edge.length = 0.4
#n = 100
#set.inla.seed = 602
x.mid <- 10
xlim.big <- c(0,20)
xlim.small <- c(3,17)
ylim.big <- c(0,20)
ylim.small <- c(3,17)

poly1 <- local.square.polygon(xlim=c(xlim.big[1]-1, x.mid-smalldist/2), 
                              ylim=x.mid+width[1]*c(-.5, .5))

poly2 <- local.square.polygon(xlim=c(x.mid+smalldist/2, xlim.big[2]+1), 
                              ylim=x.mid+width[2]*c(-.5, .5))

poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))

loc1 <- matrix(c(xlim.big[1],ylim.big[1], 
                 xlim.big[2],ylim.big[1], 
                 xlim.big[2],ylim.big[2], 
                 xlim.big[1],ylim.big[2]), 4, 2, byrow = T)

loc.out <- matrix(c(xlim.big[1],ylim.big[1], 
                    xlim.big[2],ylim.big[1], 
                    xlim.big[2],ylim.big[2], 
                    xlim.big[1],ylim.big[2],
                    xlim.big[1],ylim.big[1]), 5, 2, byrow = T)
loc.out.sp <- SpatialPolygons(list(Polygons(list(Polygon(loc.out, FALSE)), '0')))

loc.in <- matrix(c(xlim.small[1], ylim.small[1],
                   xlim.small[2], xlim.small[1],
                   xlim.small[2], ylim.small[2],
                   xlim.small[1], xlim.small[2],
                   xlim.small[1], ylim.small[1]), 5, 2, byrow = T)
loc.in.sp <- SpatialPolygons(list(Polygons(list(Polygon(loc.in, FALSE)), '0')))

loc.out.sf <- st_as_sf(loc.out.sp)
loc.in.sf  <- st_as_sf(loc.in.sp)
poly.out_sf <- st_difference(loc.out.sf, loc.in.sf)

barrier1_ <- local_square_sfc(
  xlim = c(xlim.big[1] - 1, x.mid - smalldist/2),
  ylim = x.mid + width[1] * c(-.5, .5)
)

barrier2_ <- local_square_sfc(
  xlim = c(x.mid + smalldist/2, xlim.big[2] + 1),
  ylim = x.mid + width[2] * c(-.5, .5)
)

#your domain (sf) from loc1
loc1 <- matrix(c(
  xlim.big[1], ylim.big[1],
  xlim.big[2], ylim.big[1],
  xlim.big[2], ylim.big[2],
  xlim.big[1], ylim.big[2]
), 4, 2, byrow = TRUE)

domain.xy_ <- rbind(loc1, loc1[1, ])   # close ring
domain_ <- st_sfc(st_polygon(list(domain.xy_)))

#multipolygon barriers like the tutorial
barriers_ <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier1_[[1]])[, 1:2])),
  st_polygon(list(st_coordinates(barrier2_[[1]])[, 1:2]))
)))

gg0 +
  xlim(bb["xmin"], bb["xmax"]) +
  ylim(bb["ymin"], bb["ymax"]) +
  geom_sf(data = domain_, fill = "blue") +
  geom_sf(data = barriers_, fill = "gray")

mesh_ <- fm_mesh_2d(
  loc.domain = domain.xy_, 
  max.edge = 0.4,
  offset = 1)
mesh <- mesh_
mesh$n

width <- 2
y.mid <- 10
y.up <- y.mid + (width/2)
y.low <- y.mid - (width/2)

location <- matrix(c(c(15, 12.5, 10, 10), rep(y.mid, 4)), ncol = 2)
max.edge.length = 0.4

#p1:
A <- list()
A[[1]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[1,1],y.up), nrow=1, ncol=2))
A[[2]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[1,1],(y.up + max.edge.length/2)), nrow=1, ncol=2))
A[[3]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[1,1],(y.up + max.edge.length)), nrow=1, ncol=2))

id.nodeA <- c()
id.nodeA[1] = which.max(A[[1]][1, ])
id.nodeA[2] = which.max(A[[2]][1, ])
id.nodeA[3] = which.max(A[[3]][1, ])

A.y.vector <- c(mesh$loc[id.nodeA[[1]],][2],
                mesh$loc[id.nodeA[[2]],][2],
                mesh$loc[id.nodeA[[3]],][2])

A.tmp <- A[[which(A.y.vector > y.up)[1]]]
id.node <- id.nodeA[which(A.y.vector > y.up)[1]]
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])
return.list.p1 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

#p2: 
A <- list()
A[[1]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[2,1],y.up), nrow=1, ncol=2))
A[[2]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[2,1],(y.up + max.edge.length/2)), nrow=1, ncol=2))
A[[3]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[2,1],(y.up + max.edge.length)), nrow=1, ncol=2))

id.nodeA <- c()
id.nodeA[1] = which.max(A[[1]][1, ])
id.nodeA[2] = which.max(A[[2]][1, ])
id.nodeA[3] = which.max(A[[3]][1, ])

A.y.vector <- c(mesh$loc[id.nodeA[[1]],][2],
                mesh$loc[id.nodeA[[2]],][2],
                mesh$loc[id.nodeA[[3]],][2])

A.tmp <- A[[which(A.y.vector > y.up)[1]]]
id.node <- id.nodeA[which(A.y.vector > y.up)[1]]
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])
return.list.p2 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

#p3: 
A <- list()
A[[1]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[3,1],y.up), nrow=1, ncol=2))
A[[2]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[3,1],(y.up + max.edge.length/2)), nrow=1, ncol=2))
A[[3]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[3,1],(y.up + max.edge.length)), nrow=1, ncol=2))

id.nodeA <- c()
id.nodeA[1] = which.max(A[[1]][1, ])
id.nodeA[2] = which.max(A[[2]][1, ])
id.nodeA[3] = which.max(A[[3]][1, ])

A.y.vector <- c(mesh$loc[id.nodeA[[1]],][2],
                mesh$loc[id.nodeA[[2]],][2],
                mesh$loc[id.nodeA[[3]],][2])

A.tmp <- A[[which(A.y.vector > y.up)[1]]]
id.node <- id.nodeA[which(A.y.vector > y.up)[1]]
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])
return.list.p3 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

#p4:
A <- list()
A[[1]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[4,1],y.low), nrow=1, ncol=2))
A[[2]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[4,1],(y.low - max.edge.length/2)), nrow=1, ncol=2))
A[[3]] <- inla.spde.make.A(mesh=mesh,
                           loc = matrix(c(location[4,1],(y.low - max.edge.length)), nrow=1, ncol=2))

id.nodeA <- c()
id.nodeA[1] = which.max(A[[1]][1, ])
id.nodeA[2] = which.max(A[[2]][1, ])
id.nodeA[3] = which.max(A[[3]][1, ])

A.y.vector <- c(mesh$loc[id.nodeA[[1]],][2],
                mesh$loc[id.nodeA[[2]],][2],
                mesh$loc[id.nodeA[[3]],][2])

A.tmp <- A[[which(A.y.vector < y.low)[1]]]
id.node <- id.nodeA[which(A.y.vector < y.low)[1]]
id.coord <- c(mesh$loc[id.node, 1], mesh$loc[id.node, 2])
return.list.p4 <- list(A.tmp = A.tmp, id.coord = id.coord, id.node = id.node)

return.list.mesh2 <- list(return.list.p1, return.list.p2, return.list.p3, return.list.p4)
id.coord.mesh2 <- as.matrix(rbind(return.list.mesh2[[1]]$id.coord, return.list.mesh2[[2]]$id.coord, return.list.mesh2[[3]]$id.coord, return.list.mesh2[[4]]$id.coord))

locs_ <- id.coord.mesh2

pgrid_ <- fm_evaluator(
  mesh_,
  lattice = fm_evaluator_lattice(
    mesh_,
    xlim = bb_[1, ],
    ylim = bb_[2, ],
    dims = c(600, 600)
  )
)

grid.df_ <- data.frame(
  x = rep(pgrid_$x, times = length(pgrid_$y)), 
  y = rep(pgrid_$y, each = length(pgrid_$x)))


#multipolygon barriers
barriers_ <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier1_[[1]])[, 1:2])),
  st_polygon(list(st_coordinates(barrier2_[[1]])[, 1:2]))
)))

barriers_1 <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier1_[[1]])[, 1:2]))
)))

barriers_2 <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier2_[[1]])[, 1:2]))
)))

triBarrier_1 <- unlist(fm_contains(
  x = barriers_1, 
  y = mesh_, 
  type = "centroid"))

triBarrier_2 <- unlist(fm_contains(
  x = barriers_2, 
  y = mesh_, 
  type = "centroid"))

triBarrier_1y2 <- list(
  unlist(fm_contains(
    x = barriers_1, 
    y = mesh_, 
    type = "centroid")),
  unlist(fm_contains(
    x = barriers_2, 
    y = mesh_, 
    type = "centroid")))

bfem_1y2 <- mesh2fem.barrier(mesh_, triBarrier_1y2)

rf2 <- c(0.01, 0.2, 0.3, 0.5, 0.7, 0.8, 1)
rf1 <- c(0.01, 0.3)
matrix <- list()
row <- list()

for (u in 1:length(rf1)) {
  row[[u]] <- list()
  for (f in 1:length(rf2)) {
    row[[u]][[f]] <- c(1, rf1[u], rf2[f])
  }
  matrix[[u]] <- matrix(range*unlist(row[[u]]), ncol = 3, byrow = T) 
}

Q <- list()
#same bfem_1y2
for (u in 1:length(rf1)) {
  Q[[u]] <- list()
  for (f in 1:length(rf2)) {
    Q[[u]][[f]] <- inla.barrier.q(bfem_1y2, ranges = matrix[[u]][f,], sigma = sigma)
  }
}

#same locs_
mcorrels <- list()
for (u in 1:length(rf1)) {
  mcorrels[[u]] <- list()
  for (f in 1:length(rf2)) {
    mcorrels[[u]][[f]] <- localCorrel(locs = locs_, mesh = mesh_, Q=Q[[u]][[f]])
  }
}

#same pgrid_ 
gcorrels <- list()
for (u in 1:length(rf1)) {
  gcorrels[[u]] <- list()
  for (f in 1:length(rf2)) {
    gcorrels[[u]][[f]] <- as.matrix(fm_evaluate(
      pgrid_, field = mcorrels[[u]][[f]]
    ))
  }
}

#same grid.df_
ggcorrels <- list()
for (u in 1:length(rf1)) {
  ggcorrels[[u]] <- list()
  for (f in 1:length(rf2)) {
    ggcorrels[[u]][[f]] <- do.call(
      rbind, 
      lapply(1:nrow(locs_), function(l) #number of rows in locs_
        data.frame(grid.df_, 
                   loc = paste(sprintf("%1.1f", locs_[l, 1]), 
                               sprintf("%1.1f", locs_[l, 2]), sep = ", "), 
                   correlation = gcorrels[[u]][[f]][, l])))
    
  }
}

#check first that data is alright with plots
gg0 + 
  geom_raster(
    data = ggcorrels[[2]][[4]][ggcorrels[[2]][[4]]$correlation>0.1, ], 
    mapping = aes(x = x, y = y, fill = correlation)) + 
  facet_wrap(~ loc) + 
  gg.add ## look to the appendix for the code for this

p_leg <- gg0 +
  geom_raster(data = ggcorrels[[1]][[1]], aes(x = x, y = y, fill = correlation)) +
  gg.add +
  theme_void()  # removes axes etc.

leg <- cowplot::get_legend(
  p_leg + theme(legend.position = "right")
)

#save legend as its own figure
png("plots/legend_onlym2.png", width = 600, height = 800, res = 300) #change height for magick
cowplot::ggdraw(leg)
dev.off()

rf1 <- c(0.01)
ggcorr2points <- list()
for (u in 2:length(rf1)) {
  ggcorr2points[[u]] <- list()
  for (f in 1:length(rf2)) {
    ggcorr2points[[u]][[f]] <- lapply(unique(ggcorrels[[u]][[f]]$loc), function(L) {
      gg0 +
        xlim(bb["xmin"], bb["xmax"]) +
        ylim(bb["ymin"], bb["ymax"]) +
        geom_sf(data = domain_, fill = "white", colour="black") +
        #        geom_sf(data = barriers_, fill = NA, colour="black") +
        geom_raster(
          data = subset(ggcorrels[[u]][[f]], loc == L & correlation > 0.2),
          aes(x = x, y = y, fill = correlation)
        ) +
        geom_sf(data = barriers_, fill = NA, colour="black") +
        geom_sf(data = barriers_1, fill = "white", colour="black") +
        geom_sf(data = poly.out_sf, fill = "white") +
        gg.add 
    })
  }
}

ggcorr2points[[1]][[1]] <- lapply(unique(ggcorrels[[1]][[1]]$loc), function(L) {
  gg0 +
    xlim(bb["xmin"], bb["xmax"]) +
    ylim(bb["ymin"], bb["ymax"]) +
    geom_sf(data = domain_, fill = "white", colour="black") +
    #        geom_sf(data = barriers_, fill = NA, colour="black") +
    geom_raster(
      data = subset(ggcorrels[[1]][[1]], loc == L & correlation > 0.2),
      aes(x = x, y = y, fill = correlation)
    ) +
    geom_sf(data = barriers_, fill = "white", colour="black") +
    geom_sf(data = barriers_1, fill = "white", colour="black") +
    geom_sf(data = poly.out_sf, fill = "white") +
    gg.add 
})

for (u in 1:length(rf1)) {
  for (f in 1:length(rf2)) {
    for (m in 1:nrow(locs_)) {
      dir <- file.path("plots", paste0("2corr2points_m2point", m))
      if (!dir.exists(dir)) dir.create(dir) 
      fpp <- file.path(dir, (paste0(rf1[u], "_", rf2[f], ".png")))
      png(fpp)
      print( ggcorr2points[[u]][[f]][[m]] )
      dev.off()
    }
  }
}

#article figures
library(magick)
###0.01
r <- magick::image_read("plots/2corr2points_m1point1/0.01_0.01.png")
r <- 
  magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m1point2/0.01_0.01.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m1point3/0.01_0.01.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m1point4/0.01_0.01.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

i <- magick::image_read("plots/2corr2points_m1point5/0.01_0.01.png")
i <- magick::image_crop(i, geometry = "320x310+80+70")

a <- image_append(c(i,u,y,t), stack = F)
## 0.2
r <- magick::image_read("plots/2corr2points_m1point1/0.01_0.2.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m1point2/0.01_0.2.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m1point3/0.01_0.2.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m1point4/0.01_0.2.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

i <- magick::image_read("plots/2corr2points_m1point5/0.01_0.2.png")
i <- magick::image_crop(i, geometry = "320x310+80+70")

s <- image_append(c(i,u,y,t), stack = F)

## 0.01_0.3
r <- magick::image_read("plots/2corr2points_m1point1/0.01_0.3.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m1point2/0.01_0.3.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m1point3/0.01_0.3.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m1point4/0.01_0.3.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

i <- magick::image_read("plots/2corr2points_m1point5/0.01_0.3.png")
i <- magick::image_crop(i, geometry = "320x310+80+70")

d <- image_append(c(i,u,y,t), stack = F)

## 0.01_0.5
r <- magick::image_read("plots/2corr2points_m1point1/0.01_0.5.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m1point2/0.01_0.5.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m1point3/0.01_0.5.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m1point4/0.01_0.5.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

i <- magick::image_read("plots/2corr2points_m1point5/0.01_0.5.png")
i <- magick::image_crop(i, geometry = "320x310+80+70")

f <- image_append(c(i,u,y,t), stack = F)

## 0.01_0.7
r <- magick::image_read("plots/2corr2points_m1point1/0.01_0.7.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m1point2/0.01_0.7.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m1point3/0.01_0.7.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m1point4/0.01_0.7.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

i <- magick::image_read("plots/2corr2points_m1point5/0.01_0.7.png")
i <- magick::image_crop(i, geometry = "320x310+80+70")

g <- image_append(c(i,u,y,t), stack = F)

### 0.01_0.8

r <- magick::image_read("plots/2corr2points_m1point1/0.01_0.8.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m1point2/0.01_0.8.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m1point3/0.01_0.8.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m1point4/0.01_0.8.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

i <- magick::image_read("plots/2corr2points_m1point5/0.01_0.8.png")
i <- magick::image_crop(i, geometry = "320x310+80+70")

h <- image_append(c(i,u,y,t), stack = F)

##1
r <- magick::image_read("plots/2corr2points_m1point1/0.01_1.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m1point2/0.01_1.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m1point3/0.01_1.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m1point4/0.01_1.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

i <- magick::image_read("plots/2corr2points_m1point5/0.01_1.png")
i <- magick::image_crop(i, geometry = "320x310+80+70")

j <- image_append(c(i,u,y,t), stack = F)

aa <- image_append(c(j,h,g,f,d,s,a), stack = TRUE)

image_write(aa, path = "plots/2corr2points.m1_0.01.png", format = "png")

#for mesh 2
###0.01
r <- magick::image_read("plots/2corr2points_m2point1/0.01_0.01.png")
r <- 
  magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m2point2/0.01_0.01.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m2point3/0.01_0.01.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m2point4/0.01_0.01.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

a <- image_append(c(u,t), stack = F)
## 0.2
r <- magick::image_read("plots/2corr2points_m2point1/0.01_0.2.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m2point2/0.01_0.2.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m2point3/0.01_0.2.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m2point4/0.01_0.2.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

s <- image_append(c(u,t), stack = F)

## 0.01_0.3
r <- magick::image_read("plots/2corr2points_m2point1/0.01_0.3.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m2point2/0.01_0.3.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m2point3/0.01_0.3.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m2point4/0.01_0.3.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

d <- image_append(c(u,t), stack = F)

## 0.01_0.5
r <- magick::image_read("plots/2corr2points_m2point1/0.01_0.5.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m2point2/0.01_0.5.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m2point3/0.01_0.5.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m2point4/0.01_0.5.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

f <- image_append(c(u,t), stack = F)

## 0.01_0.7
r <- magick::image_read("plots/2corr2points_m2point1/0.01_0.7.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m2point2/0.01_0.7.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m2point3/0.01_0.7.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m2point4/0.01_0.7.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

g <- image_append(c(u,t), stack = F)

### 0.01_0.8

r <- magick::image_read("plots/2corr2points_m2point1/0.01_0.8.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m2point2/0.01_0.8.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m2point3/0.01_0.8.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m2point4/0.01_0.8.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

h <- image_append(c(u,t), stack = F)

##1
r <- magick::image_read("plots/2corr2points_m2point1/0.01_1.png")
r <- magick::image_crop(r, geometry = "320x310+80+70")

t <- magick::image_read("plots/2corr2points_m2point2/0.01_1.png")
t <- magick::image_crop(t, geometry = "320x310+80+70")

y <- magick::image_read("plots/2corr2points_m2point3/0.01_1.png")
y <- magick::image_crop(y, geometry = "320x310+80+70")

u <- magick::image_read("plots/2corr2points_m2point4/0.01_1.png")
u <- magick::image_crop(u, geometry = "320x310+80+70")

j <- image_append(c(u,t), stack = F)

ss <- image_append(c(j,h,g,f,d,s,a), stack = TRUE)

image_write(ss, path = "plots/2corr2points.m2_0.01.png", format = "png")

ss <- magick::image_read("plots/2corr2points.m2_0.01.png")

aa <- magick::image_read("plots/2corr2points.m1_0.01.png")

ssf <- ss |>
  image_border(color = "white", geometry = "35x35") |>  # frame thickness
  image_border(color = "black", geometry = "3x3") |>
  image_border(color = "white", geometry = "40x10")

aaf <- aa |>
  image_border(color = "white", geometry = "35x35") |>  # frame thickness
  image_border(color = "black", geometry = "3x3") |>
  image_border(color = "white", geometry = "40x10")

pp <- image_append(c(aaf,ssf), stack = F)

l <- magick::image_read("plots/legend_onlym1.png")
ll <- l |>
  image_scale("70%") |>
  image_scale("x700%") |>
  image_crop(geometry = "200x3360+00+1200") |>
  image_border(color = "white", geometry = "90x300") |>
  image_annotate(
    text = "1.0",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "+0-1320"
  ) |>
  image_annotate(
    text = "0.8",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-880"
  ) |>
  image_annotate(
    text = "0.6",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-440"
  ) |>
  image_annotate(
    text = "0.4",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-00"
  ) |>
  image_annotate(
    text = "0.2",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0+440"
  ) |>
  image_crop("330x3120+120+80")

img <- 
  image_append(c(pp,ll), stack = F)

image_write(img, path = "plots/2corr2points_0.01.new.png", format = "png")

b <- pp |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120+7
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120+320
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120+939
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120-303
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120-612
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120-921
"
  ) |> #################################
image_annotate(
  text = "0.01",
  size = 30,
  color = "black",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+1570+7
"
) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570+320
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570+939
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570-303
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570-612
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570-921
"
  ) |> #################################
image_annotate(
  text = "0.01",
  size = 30,
  color = "black",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+1895+7
"
) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1895+320
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1895+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1895+939
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1895-303
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1895-612
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1895-921
"
  ) |> #################################
image_annotate(
  text = "0.01",
  size = 30,
  color = "black",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+440+7
"
) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+440+320
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+440+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+440+939
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+440-303
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+440-612
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+440-921
"
  ) |> #################################
image_annotate(
  text = "0.01",
  size = 30,
  color = "black",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+760+7
"
) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+760+320
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+760+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+760+939
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+760-303
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+760-612
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+760-921
"
  ) |> #################################



image_annotate(
  text = "0.01",
  size = 30,
  color = "black",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+1080+7
"
) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1080+320
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1080+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1080+939
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1080-303
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1080-612
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1080-921
"
  ) #|> #################################

####
####
####

bb <- b |>
  image_annotate(
    text = "0.5",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340+7
"
  ) |>
  image_annotate(
    text = "0.3",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340+320
"
  ) |>
  image_annotate(
    text = "0.2",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+325+939
"
  ) |>
  image_annotate(
    text = "0.7",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340-303
"
  ) |>
  image_annotate(
    text = "0.8",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340-612
"
  ) |>
  image_annotate(
    text = "1",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+365-921
"
  ) 
####
####
####

bb <- bb |> 
  image_annotate(
    text = "0.5",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780+7
"
  ) |>
  image_annotate(
    text = "0.3",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780+320
"
  ) |>
  image_annotate(
    text = "0.2",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1763+939
"
  ) |>
  image_annotate(
    text = "0.7",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780-303
"
  ) |>
  image_annotate(
    text = "0.8",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780-612
"
  ) |>
  image_annotate(
    text = "1",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1803-921
"
  ) 

####
####
bb <- bb |> 
  image_annotate(
    text = "0.5",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+2143+7
"
  ) |>
  image_annotate(
    text = "0.3",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+2143+320
"
  ) |>
  image_annotate(
    text = "0.2",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+2143+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+2126+939
"
  ) |>
  image_annotate(
    text = "0.7",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+2143-303
"
  ) |>
  image_annotate(
    text = "0.8",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+2143-612
"
  ) |>
  image_annotate(
    text = "1",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+2156-921
"
  ) 


####
####
####

bb <- bb |> image_annotate(
  text = "0.5",
  size = 30,
  color = "red",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+662+7
"
) |>
  image_annotate(
    text = "0.3",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+662+320
"
  ) |>
  image_annotate(
    text = "0.2",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+662+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+650+939
"
  ) |>
  image_annotate(
    text = "0.7",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+662-303
"
  ) |>
  image_annotate(
    text = "0.8",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+662-612
"
  ) |>
  image_annotate(
    text = "1",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+685-921
"
  ) 
####
####

bb <- bb |> image_annotate(
  text = "0.5",
  size = 30,
  color = "red",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+983+7
"
) |>
  image_annotate(
    text = "0.3",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+983+320
"
  ) |>
  image_annotate(
    text = "0.2",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+983+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+968+939
"
  ) |>
  image_annotate(
    text = "0.7",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+983-303
"
  ) |>
  image_annotate(
    text = "0.8",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+983-612
"
  ) |>
  image_annotate(
    text = "1",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1003-921
"
  )#|> #################################



bb <- bb |> image_annotate(
  text = "0.5",
  size = 30,
  color = "red",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+1346+7
"
) |>
  image_annotate(
    text = "0.3",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1346+320
"
  ) |>
  image_annotate(
    text = "0.2",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1346+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1331+939
"
  ) |>
  image_annotate(
    text = "0.7",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1346-303
"
  ) |>
  image_annotate(
    text = "0.8",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1346-612
"
  ) |>
  image_annotate(
    text = "1",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1362-921
"
  ) #|> #################################

l <- magick::image_read("plots/legend_onlym1.png")
ll <- l |>
  image_scale("70%") |>
  image_scale("x700%") |>
  image_crop(geometry = "200x3360+00+1200") |>
  image_border(color = "white", geometry = "90x300") |>
  image_annotate(
    text = "1.0",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "+0-1320"
  ) |>
  image_annotate(
    text = "0.8",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-880"
  ) |>
  image_annotate(
    text = "0.6",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-440"
  ) |>
  image_annotate(
    text = "0.4",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-00"
  ) |>
  image_annotate(
    text = "0.2",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0+440"
  ) |>
  image_crop("330x3120+120+80")

img <- 
  image_append(c(bb,ll), stack = F)

image_write(img, path = "plots/2corr2points.bb_0.01leg.png", format = "png")
image_write(bb, path = "plots/2corr2points.bb_0.01.png", format = "png")

ss <- magick::image_read("plots/2corr2points.m2_0.01.png")

aa <- magick::image_read("plots/2corr2points.m1_0.01.png")

ssf <- ss |>
  image_border(color = "white", geometry = "35x35") |>  # frame thickness
  image_border(color = "black", geometry = "3x3") |>
  image_border(color = "white", geometry = "40x10")

aaf <- aa |>
  image_border(color = "white", geometry = "35x35") |>  # frame thickness
  image_border(color = "black", geometry = "3x3") |>
  image_border(color = "white", geometry = "40x10")

pp <- image_append(c(aaf,ssf), stack = F)

l <- magick::image_read("plots/legend_onlym1.png")
ll <- l |>
  image_scale("70%") |>
  image_scale("x700%") |>
  image_crop(geometry = "200x3360+00+1200") |>
  image_border(color = "white", geometry = "90x300") |>
  image_annotate(
    text = "1.0",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "+0-1320"
  ) |>
  image_annotate(
    text = "0.8",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-880"
  ) |>
  image_annotate(
    text = "0.6",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-440"
  ) |>
  image_annotate(
    text = "0.4",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-00"
  ) |>
  image_annotate(
    text = "0.2",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0+440"
  ) |>
  image_crop("330x3120+120+80")

img <- 
  image_append(c(pp,ll), stack = F)

image_write(img, path = "plots/2corr2points_0.01.new.png", format = "png")

b <- pp |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120+7
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120+320
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120+939
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120-303
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120-612
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+120-921
"
  ) |> #################################
image_annotate(
  text = "0.01",
  size = 30,
  color = "black",
  font = "Helvetica",
  #    strokewidth = 1.5,
  gravity = "west",      # top center
  location = "+1570+7
"
) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570+320
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570+939
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570-303
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570-612
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "black",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1570-921
"
  ) 

bb <- b |>
  image_annotate(
    text = "0.5",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340+7
"
  ) |>
  image_annotate(
    text = "0.3",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340+320
"
  ) |>
  image_annotate(
    text = "0.2",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+325+939
"
  ) |>
  image_annotate(
    text = "0.7",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340-303
"
  ) |>
  image_annotate(
    text = "0.8",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+340-612
"
  ) |>
  image_annotate(
    text = "1",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+365-921
"
  ) 
####
####
####

bb <- bb |> 
  image_annotate(
    text = "0.5",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780+7
"
  ) |>
  image_annotate(
    text = "0.3",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780+320
"
  ) |>
  image_annotate(
    text = "0.2",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780+629
"
  ) |>
  image_annotate(
    text = "0.01",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1763+939
"
  ) |>
  image_annotate(
    text = "0.7",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780-303
"
  ) |>
  image_annotate(
    text = "0.8",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1780-612
"
  ) |>
  image_annotate(
    text = "1",
    size = 30,
    color = "red",
    font = "Helvetica",
    #    strokewidth = 1.5,
    gravity = "west",      # top center
    location = "+1803-921
"
  ) 


l <- magick::image_read("plots/legend_onlym1.png")
ll <- l |>
  image_scale("70%") |>
  image_scale("x700%") |>
  image_crop(geometry = "200x3360+00+1200") |>
  image_border(color = "white", geometry = "90x300") |>
  image_annotate(
    text = "1.0",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "+0-1320"
  ) |>
  image_annotate(
    text = "0.8",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-880"
  ) |>
  image_annotate(
    text = "0.6",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-440"
  ) |>
  image_annotate(
    text = "0.4",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0-00"
  ) |>
  image_annotate(
    text = "0.2",
    size = 60,
    color = "black",
    font = "Helvetica",
    gravity = "east",      # top center
    location = "-0+440"
  ) |>
  image_crop("330x3120+120+80")

img <- 
  image_append(c(bb,ll), stack = F)

image_write(img, path = "plots/3corr2points.bb_0.01leg.png", format = "png")
image_write(bb, path = "plots/2corr2points.bb_0.01.png", format = "png")












