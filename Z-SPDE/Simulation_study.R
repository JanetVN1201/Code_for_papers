

local.plot.field = function(field, pal = plasma(50), ...){
  xlim = c(3, 17); ylim = xlim;
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, col = pal, ...)  
  # - Use image.plot to get nice colors and legend
}

#SECTION 4 RESULTS

################
################ TWO BARRIER CONFIGURATION
################

prior.range = c(1, 0.5); prior.sigma = c(1, 0.1)
zlim.sd <- c(0,1.5)
zlim.mean <- c(-2, 2)
xlim.marginal = c(0,6)
#geom 
SM <- c(1.5)
W <- c(3)
R <- 2
SS <- 880
N <- c(30)
#fr <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5,0.7,0.8,1)
fr <- c(0.01, 0.8, 1)

for (s in 1:length(SS)) {
  for (d in 1:length(SM)) {
    for (w in 1:length(W)) {
      for(r in 1:length(R)) {
        for (n in 1:length(N)) {
          
          
          
          
          smalldist = SM[d]
          width = c(W[w],W[w])
          range = R[r]
          max.edge.length = 0.4
          n = N[n]
          set.inla.seed = SS[s]
          x.mid = 10
          xlim.big  = c(0,20)
          xlim.small = c(3,17)
          ylim.big = c(0,20)
          ylim.small = c(3,17)
          
          prior.range = prior.range
          prior.sigma = c(1,0.1)
          sigma.u = 1
          sigma.epsilon = 0.2
          prior.range.st = prior.range
          prior.sigma.st = c(1,0.1)
          
          
          
          poly1 <- local.square.polygon(xlim=c(xlim.big[1]-1, x.mid-smalldist/2), 
                                        ylim=x.mid+width[1]*c(-.5, .5))
          
          poly2 <- local.square.polygon(xlim=c(x.mid+smalldist/2, xlim.big[2]+1), 
                                        ylim=x.mid+width[2]*c(-.5, .5))
          
          poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))
          
          loc1 <- matrix(c(xlim.big[1],ylim.big[1], 
                           xlim.big[2],ylim.big[1], 
                           xlim.big[2],ylim.big[2], 
                           xlim.big[1],ylim.big[2]), 4, 2, byrow = T)
          
          seg <- inla.sp2segment(poly.original)
          # - Transforms a SpatialPolygon to an "inla polygon"
          mesh <- inla.mesh.2d(loc=loc1, 
                               interior = seg, 
                               max.e = max.edge.length, 
                               offset=1)
          
          tl <- length(mesh$graph$tv[,1])
          # - the number of triangles in the mesh
          posTri <- matrix(0, tl, 2)
          
          for (t in 1:tl){
            temp = mesh$loc[mesh$graph$tv[t, ], ]
            posTri[t,] = colMeans(temp)[c(1,2)] 
          }
          
          posTri <- SpatialPoints(posTri)
          # - the positions of the triangle centers
          
          bar.original <- over(poly.original, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar.original <- unlist(bar.original)
          poly.bar.orginal <- inla.barrier.polygon(mesh, barrier.triangles = bar.original)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 1
          bar1 <- over(poly1, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar1 <- unlist(bar1)
          poly.bar1 <- inla.barrier.polygon(mesh, barrier.triangles = bar1)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 2
          bar2 <- over(poly2, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar2 <- unlist(bar2)
          poly.bar2 <- inla.barrier.polygon(mesh, barrier.triangles = bar2)
          
          mat <-  inla.barrier.fem.plus(mesh, list(bar1, bar2)) 
          fem <- mat 
          barrier.triangles <- list(bar1, bar2)
          
          # PLOTS
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
          
          poly.rect <- SpatialPolygons(list(Polygons(list(poly1_h, poly2_h), '0')))
          poly.rect_sf <- st_as_sf(poly.rect)
          poly.sq <- SpatialPolygons(list(Polygons(list(locp2), '0')))
          poly.sq_sf <- st_as_sf(poly.sq)
          
          sq_bbox <- st_bbox(poly.sq_sf) %>% st_as_sfc()
          poly1_bbox <- st_bbox(poly1) %>% st_as_sfc()
          poly2_bbox <- st_bbox(poly2) %>% st_as_sfc()
          
          
          
          tbm.frac <- list()
          trans <- list()
          
          
          for (p in 1:length(fr)) { #
            fp <- file.path("plots", paste0("sm", SM[d], "w", W[w], "r", R[r], "seed", SS[s], ".n30", "fr",fr[p], "/"))
            dir.create(fp, recursive = TRUE)
            gif.path = fp
            marginals.path = fp
            
            f <- fr[p]
            range.fraction <- c(0.01, f)
            tbm.frac[[p]] <- model.tbm.plus(mesh = mesh, fem = fem, 
                                            barrier.triangles = barrier.triangles,
                                            prior.range = prior.range, prior.sigma = prior.sigma,  
                                            range.fraction = range.fraction,
                                            range = range,
                                            set.inla.seed = set.inla.seed,
                                            loc.data = loc.data,
                                            sigma.u = 1, sigma.epsilon = 0.2,
                                            poly.original = poly.bar.orginal,
                                            prior.range.st = prior.range.st,
                                            prior.sigma.st = prior.sigma.st,
                                            return.list = TRUE)
            
            trans[[1]] <- tbm.frac[[p]]
            plot.model.tbm.plus.l.z(trans, 
                                    fr = fr[p],
                                    poly.bar.orginal = poly.bar.orginal,
                                    range = range,
                                    fps = 1,
                                    zlim.sd = zlim.sd,
                                    zlim.mean = zlim.mean,
                                    gif.path = fp)
            
            
            
            fpp <- file.path(paste0(marginals.path, "tbm", f, ".png"))
            png(fpp)
            #par(mfrow = c(3, 1))
            tmp = inla.tmarginal(function(x) exp(x), trans[[1]]$pos.tbm$res$marginals.hyperpar[[3]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = xlim.marginal)
            xvals = seq(0, 10, length.out=1000)
            lambda = 1.00; lines(xvals, exp(-lambda*xvals), lty='dashed')
            abline(v=range, col="blue")
            dev.off()
            
            
            
            fpp <- file.path(paste0(marginals.path, "st", f, ".png"))
            png(fpp)
            tmp = inla.tmarginal(function(x) (x), trans[[1]]$pos.st$res$marginals.hyperpar[[2]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = xlim.marginal)
            xvals = seq(0, 10, length.out=1000)
            abline(v=range, col="blue")
            dev.off()
            
            fpp <- file.path(paste0(marginals.path, "bm", f, ".png"))
            png(fpp)
            tmp = inla.tmarginal(function(x) exp(x), trans[[1]]$pos.bm$res$marginals.hyperpar[[3]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = xlim.marginal)
            xvals = seq(0, 10, length.out=1000)
            abline(v=range, col="blue")
            dev.off()
            
            
            
            
          }
          
          
        }
        
      }
    }
  }
  
}

#correct color of range fraction 0.01
new.u.tbm <- ifelse(trans[[1]]$pos.tbm$res$summary.random$s$mean > 2, 2, trans[[1]]$pos.tbm$res$summary.random$s$mean)
new.u.bm <- ifelse(trans[[1]]$pos.bm$res$summary.random$s$mean > 2, 2, trans[[1]]$pos.bm$res$summary.random$s$mean)
new.u.st <- ifelse(trans[[1]]$pos.st$res$summary.random$s$mean > 2, 2, trans[[1]]$pos.st$res$summary.random$s$mean)

fr <- 0.01
pal.pos.mean = turbo(50)

fp <- file.path("plots/color.corrected", paste0("sm", SM, "w", W, "r", R, "seed", SS, ".n30", "fr",fr, "/"))
dir.create(fp, recursive = TRUE)

## TBM
# POSTERIOR MEAN
fpp <- file.path(paste0(fp, "mean.tbm.png"))
png(fpp) 
local.plot.field(new.u.tbm, 
                 main="Spatial estimate for Transparent Barrier model",
                 sub = paste("range.fraction = c(", 
                             range.fraction[1], ", ", range.fraction[2], ")"),
                 axes = F, pal = pal.pos.mean,
                 zlim = zlim.mean)
plot(poly.bar.orginal, add=T)
dev.off()


## BM
# POSTERIOR MEAN
fpp <- file.path(paste0(fp, "mean.bm.png"))
png(fpp)
local.plot.field(new.u.bm, 
                 main="Spatial estimate for Transparent Barrier model",
                 sub = paste("range.fraction = c(", 
                             range.fraction[1], ", ", range.fraction[2], ")"),
                 axes = F, pal = pal.pos.mean,
                 zlim = zlim.mean)
plot(poly.bar.orginal, add=T)
dev.off()

## st
# POSTERIOR MEAN
fpp <- file.path(paste0(fp, "mean.st.png"))
png(fpp)
local.plot.field(new.u.st, 
                 main="Spatial estimate for Transparent Barrier model",
                 sub = paste("range.fraction = c(", 
                             range.fraction[1], ", ", range.fraction[2], ")"),
                 axes = F, pal = pal.pos.mean,
                 zlim = zlim.mean)
plot(poly.bar.orginal, add=T)
dev.off()

################
################ BARRIER NO CANAL
################


prior.range = c(1, 0.5); prior.sigma = c(1, 0.1)
#geom 
SM <- c(0)
W <- c(1)
R <- 4
SS <- 987
N <- c(50)
#fr <- c(0.01, 0.2, 0.3, 0.4, 0.5,0.7,0.8,1)
fr <- c(0.01,0.5,1)

for (s in 1:length(SS)) {
  for (d in 1:length(SM)) {
    for (w in 1:length(W)) {
      for(r in 1:length(R)) {
        for (n in 1:length(N)) {
          
          
          
          
          smalldist = SM[d]
          width = c(W[w],W[w])
          range = R[r]
          max.edge.length = 0.4
          n = N[n]
          set.inla.seed = SS[s]
          x.mid = 10
          xlim.big  = c(0,20)
          xlim.small = c(3,17)
          ylim.big = c(0,20)
          ylim.small = c(3,17)
          
          prior.range = prior.range
          prior.sigma = c(1,0.1)
          sigma.u = 1
          sigma.epsilon = 0.2
          prior.range.st = prior.range
          prior.sigma.st = c(1,0.1)
          
          
          
          poly1 <- local.square.polygon(xlim=c(xlim.big[1]-1, x.mid-smalldist/2), 
                                        ylim=x.mid+width[1]*c(-.5, .5))
          
          poly2 <- local.square.polygon(xlim=c(x.mid+smalldist/2, xlim.big[2]+1), 
                                        ylim=x.mid+width[2]*c(-.5, .5))
          
          poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))
          
          loc1 <- matrix(c(xlim.big[1],ylim.big[1], 
                           xlim.big[2],ylim.big[1], 
                           xlim.big[2],ylim.big[2], 
                           xlim.big[1],ylim.big[2]), 4, 2, byrow = T)
          
          seg <- inla.sp2segment(poly.original)
          # - Transforms a SpatialPolygon to an "inla polygon"
          mesh <- inla.mesh.2d(loc=loc1, 
                               interior = seg, 
                               max.e = max.edge.length, 
                               offset=1)
          
          tl <- length(mesh$graph$tv[,1])
          # - the number of triangles in the mesh
          posTri <- matrix(0, tl, 2)
          
          for (t in 1:tl){
            temp = mesh$loc[mesh$graph$tv[t, ], ]
            posTri[t,] = colMeans(temp)[c(1,2)] 
          }
          
          posTri <- SpatialPoints(posTri)
          # - the positions of the triangle centers
          
          bar.original <- over(poly.original, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar.original <- unlist(bar.original)
          poly.bar.orginal <- inla.barrier.polygon(mesh, barrier.triangles = bar.original)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 1
          bar1 <- over(poly1, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar1 <- unlist(bar1)
          poly.bar1 <- inla.barrier.polygon(mesh, barrier.triangles = bar1)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 2
          bar2 <- over(poly2, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar2 <- unlist(bar2)
          poly.bar2 <- inla.barrier.polygon(mesh, barrier.triangles = bar2)
          
          mat <-  inla.barrier.fem.plus(mesh, list(bar1, bar2)) 
          fem <- mat 
          barrier.triangles <- list(bar1, bar2)
          
          # PLOTS
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
          
          poly.rect <- SpatialPolygons(list(Polygons(list(poly1_h, poly2_h), '0')))
          poly.rect_sf <- st_as_sf(poly.rect)
          poly.sq <- SpatialPolygons(list(Polygons(list(locp2), '0')))
          poly.sq_sf <- st_as_sf(poly.sq)
          
          sq_bbox <- st_bbox(poly.sq_sf) %>% st_as_sfc()
          poly1_bbox <- st_bbox(poly1) %>% st_as_sfc()
          poly2_bbox <- st_bbox(poly2) %>% st_as_sfc()
          
          
          
          tbm.frac <- list()
          trans <- list()
          
          
          for (p in 1:length(fr)) { #
            fp <- file.path("plots", paste0("sm", SM[d], "w", W[w], "r", R[r], "seed", SS[s], ".n50", "fr", fr[p], "/"))
            dir.create(fp, recursive = TRUE)
            gif.path = fp
            marginals.path = fp
            
            f <- fr[p]
            range.fraction <- c(0.01, f)
            tbm.frac[[p]] <- model.tbm.plus(mesh = mesh, fem = fem, 
                                            barrier.triangles = barrier.triangles,
                                            prior.range = prior.range, prior.sigma = prior.sigma,  
                                            range.fraction = range.fraction,
                                            range = range,
                                            set.inla.seed = set.inla.seed,
                                            loc.data = loc.data,
                                            sigma.u = 1, sigma.epsilon = 0.2,
                                            poly.original = poly.bar.orginal,
                                            prior.range.st = prior.range.st,
                                            prior.sigma.st = prior.sigma.st,
                                            return.list = TRUE)
            
            trans[[1]] <- tbm.frac[[p]]
            plot.model.tbm.plus(trans, 
                                fr = fr[p],
                                poly.bar.orginal = poly.bar.orginal,
                                range = range,
                                fps = 1,
                                gif.path = fp)
            
            
            
            fpp <- file.path(paste0(marginals.path, "tbm", f, ".png"))
            png(fpp)
            #par(mfrow = c(3, 1))
            tmp = inla.tmarginal(function(x) exp(x), trans[[1]]$pos.tbm$res$marginals.hyperpar[[3]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = c(1,8))
            xvals = seq(0, 10, length.out=1000)
            
            abline(v=range, col="blue")
            dev.off()
            
            fpp <- file.path(paste0(marginals.path, "st", f, ".png"))
            png(fpp)
            tmp = inla.tmarginal(function(x) (x), trans[[1]]$pos.st$res$marginals.hyperpar[[2]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = c(1,8))
            xvals = seq(0, 10, length.out=1000)
            abline(v=range, col="blue")
            dev.off()
            
            fpp <- file.path(paste0(marginals.path, "bm", f, ".png"))
            png(fpp)
            tmp = inla.tmarginal(function(x) exp(x), trans[[1]]$pos.bm$res$marginals.hyperpar[[3]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = c(1,8))
            xvals = seq(0, 10, length.out=1000)
            
            abline(v=range, col="blue")
            dev.off()
            
            
            
            
          }
          
          
        }
        
      }
    }
  }
  
}




###########
###########
###########

#sm0w1r4seed987.n50
zlim.sd <- c(0,1.4)
zlim.mean <- c(-2, 2.822)


prior.range = c(1, 0.5); prior.sigma = c(1, 0.1)
#geom 
SM <- c(0)
W <- c(1)
R <- 4
SS <- 987
N <- c(50)
#fr <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5,0.7,0.8,1)

fr <- c(0.01,0.5,1)

for (s in 1:length(SS)) {
  for (d in 1:length(SM)) {
    for (w in 1:length(W)) {
      for(r in 1:length(R)) {
        for (n in 1:length(N)) {
          
          
          
          
          smalldist = SM[d]
          width = c(W[w],W[w])
          range = R[r]
          max.edge.length = 0.4
          n = N[n]
          set.inla.seed = SS[s]
          x.mid = 10
          xlim.big  = c(0,20)
          xlim.small = c(3,17)
          ylim.big = c(0,20)
          ylim.small = c(3,17)
          
          prior.range = prior.range
          prior.sigma = c(1,0.1)
          sigma.u = 1
          sigma.epsilon = 0.2
          prior.range.st = prior.range
          prior.sigma.st = c(1,0.1)
          
          
          
          poly1 <- local.square.polygon(xlim=c(xlim.big[1]-1, x.mid-smalldist/2), 
                                        ylim=x.mid+width[1]*c(-.5, .5))
          
          poly2 <- local.square.polygon(xlim=c(x.mid+smalldist/2, xlim.big[2]+1), 
                                        ylim=x.mid+width[2]*c(-.5, .5))
          
          poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))
          
          loc1 <- matrix(c(xlim.big[1],ylim.big[1], 
                           xlim.big[2],ylim.big[1], 
                           xlim.big[2],ylim.big[2], 
                           xlim.big[1],ylim.big[2]), 4, 2, byrow = T)
          
          seg <- inla.sp2segment(poly.original)
          # - Transforms a SpatialPolygon to an "inla polygon"
          mesh <- inla.mesh.2d(loc=loc1, 
                               interior = seg, 
                               max.e = max.edge.length, 
                               offset=1)
          
          tl <- length(mesh$graph$tv[,1])
          # - the number of triangles in the mesh
          posTri <- matrix(0, tl, 2)
          
          for (t in 1:tl){
            temp = mesh$loc[mesh$graph$tv[t, ], ]
            posTri[t,] = colMeans(temp)[c(1,2)] 
          }
          
          posTri <- SpatialPoints(posTri)
          # - the positions of the triangle centers
          
          bar.original <- over(poly.original, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar.original <- unlist(bar.original)
          poly.bar.orginal <- inla.barrier.polygon(mesh, barrier.triangles = bar.original)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 1
          bar1 <- over(poly1, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar1 <- unlist(bar1)
          poly.bar1 <- inla.barrier.polygon(mesh, barrier.triangles = bar1)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 2
          bar2 <- over(poly2, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar2 <- unlist(bar2)
          poly.bar2 <- inla.barrier.polygon(mesh, barrier.triangles = bar2)
          
          mat <-  inla.barrier.fem.plus(mesh, list(bar1, bar2)) 
          fem <- mat 
          barrier.triangles <- list(bar1, bar2)
          
          # PLOTS
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
          
          poly.rect <- SpatialPolygons(list(Polygons(list(poly1_h, poly2_h), '0')))
          poly.rect_sf <- st_as_sf(poly.rect)
          poly.sq <- SpatialPolygons(list(Polygons(list(locp2), '0')))
          poly.sq_sf <- st_as_sf(poly.sq)
          
          sq_bbox <- st_bbox(poly.sq_sf) %>% st_as_sfc()
          poly1_bbox <- st_bbox(poly1) %>% st_as_sfc()
          poly2_bbox <- st_bbox(poly2) %>% st_as_sfc()
          
          
          
          tbm.frac <- list()
          trans <- list()
          
          
          for (p in 1:length(fr)) { #
            fp <- file.path("plots", paste0("sm", SM[d], "w", W[w], "r", R[r], "seed", SS[s], ".n50", "fr", fr[p], "/"))
            dir.create(fp, recursive = TRUE)
            gif.path = fp
            marginals.path = fp
            
            f <- fr[p]
            range.fraction <- c(0.01, f)
            tbm.frac[[p]] <- model.tbm.plus(mesh = mesh, fem = fem, 
                                            barrier.triangles = barrier.triangles,
                                            prior.range = prior.range, prior.sigma = prior.sigma,  
                                            range.fraction = range.fraction,
                                            range = range,
                                            set.inla.seed = set.inla.seed,
                                            loc.data = loc.data,
                                            sigma.u = 1, sigma.epsilon = 0.2,
                                            poly.original = poly.bar.orginal,
                                            prior.range.st = prior.range.st,
                                            prior.sigma.st = prior.sigma.st,
                                            return.list = TRUE)
            
            trans[[1]] <- tbm.frac[[p]]
            plot.model.tbm.plus.l.z(trans, 
                                    fr = fr[p],
                                    poly.bar.orginal = poly.bar.orginal,
                                    range = range,
                                    fps = 1,
                                    zlim.sd = zlim.sd,
                                    zlim.mean = zlim.mean,
                                    gif.path = fp)
            
            
            
            fpp <- file.path(paste0(marginals.path, "tbm", f, ".png"))
            png(fpp)
            #par(mfrow = c(3, 1))
            tmp = inla.tmarginal(function(x) exp(x), trans[[1]]$pos.tbm$res$marginals.hyperpar[[3]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = xlim.marginal)
            xvals = seq(0, 10, length.out=1000)
            lambda = 1.00; lines(xvals, exp(-lambda*xvals), lty='dashed')
            abline(v=range, col="blue")
            dev.off()
            
            
            
            fpp <- file.path(paste0(marginals.path, "st", f, ".png"))
            png(fpp)
            tmp = inla.tmarginal(function(x) (x), trans[[1]]$pos.st$res$marginals.hyperpar[[2]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = xlim.marginal)
            xvals = seq(0, 10, length.out=1000)
            abline(v=range, col="blue")
            dev.off()
            
            fpp <- file.path(paste0(marginals.path, "bm", f, ".png"))
            png(fpp)
            tmp = inla.tmarginal(function(x) exp(x), trans[[1]]$pos.bm$res$marginals.hyperpar[[3]]) 
            plot(tmp, type = "l", xlab = "r", ylab = "Density", xlim = xlim.marginal)
            xvals = seq(0, 10, length.out=1000)
            abline(v=range, col="blue")
            dev.off()
            
            
            
            
          }
          
          
        }
        
      }
    }
  }
  
}







