## ============================================================================
## Barrier Model with Per-Barrier Ranges — ν = 1 and ν = 2
## ============================================================================
##
## Two barrier regions, each with its own range fraction.
## ranges = c(range_domain, range_barrier1, range_barrier2)
##
## We run FOUR model fits:
##   A) nu = 1, equal   range fractions  (range.fraction = c(0.01, 0.01))
##   B) nu = 1, unequal range fractions  (range.fraction = c(0.01, 0.7))
##   C) nu = 2, equal   range fractions  (range.fraction = c(0.01, 0.01))
##   D) nu = 2, unequal range fractions  (range.fraction = c(0.01, 0.7))
##
## Based on the minimal_sim_example.R geometry.
##
## FIX NOTE (rgeneric + sf crash):
##   If you get:  "Symbol not found: _sqlite3_enable_load_extension"
##   when running the rgeneric models, run in your terminal:
##     brew reinstall sqlite3 && brew reinstall proj && brew reinstall gdal
##   Then in R:  install.packages("sf", type = "source")
##   Or export DYLD_LIBRARY_PATH="$(brew --prefix sqlite3)/lib" before R.
## ============================================================================

library(sp)
library(sf)
library(scales)
library(INLA)
library(INLAspacetime)
library(inlabru)
library(patchwork)
library(fmesher)
library(ggplot2)
library(Matrix)

## ============================================================================
## PART 1: Domain geometry (from minimal_sim_example.R)
## ============================================================================

smalldist    <- 1.5
W            <- 3
width        <- c(W, W)
n_obs        <- 30
set.inla.seed <- 880

x.mid      <- 10
xlim.big   <- c(0, 20)
xlim.small <- c(3, 17)
ylim.big   <- c(0, 20)
ylim.small <- c(3, 17)

## Model parameters
range_domain <- 4
sigma_true   <- 1
sigma.u       <- 1
sigma.epsilon <- 0.2
max.edge.length <- 0.4

## Priors (shared by all models)
prior.range <- c(1, 0.5)
prior.sigma <- c(1, 0.1)

## Helper: square sf polygon
local_square_sfc <- function(xlim, ylim) {
  xlim <- range(xlim); ylim <- range(ylim)
  xy <- rbind(
    c(xlim[1], ylim[2]),
    c(xlim[1], ylim[1]),
    c(xlim[2], ylim[1]),
    c(xlim[2], ylim[2]),
    c(xlim[1], ylim[2])
  )
  st_sfc(st_polygon(list(xy)))
}

## Helper: sp square polygon
local.square.polygon <- function(xlim, ylim) {
  xlim <- range(xlim); ylim <- range(ylim)
  corner1 <- c(xlim[1], ylim[2])
  corner2 <- c(xlim[2], ylim[1])
  poly <- Polygon(rbind(corner1, c(corner1[1], corner2[2]),
                        corner2, c(corner2[1], corner1[2]), corner1),
                  hole = FALSE)
  SpatialPolygons(list(Polygons(list(poly), ID = runif(1))))
}

## Helper: sp square polygon (hole)
local.square.polygon_T <- function(xlim, ylim) {
  xlim <- range(xlim); ylim <- range(ylim)
  corner1 <- c(xlim[1], ylim[2])
  corner2 <- c(xlim[2], ylim[1])
  Polygon(rbind(corner1, c(corner1[1], corner2[2]),
                corner2, c(corner2[1], corner1[2]), corner1),
          hole = TRUE)
}

## Domain outline
loc1 <- matrix(c(
  xlim.big[1], ylim.big[1],
  xlim.big[2], ylim.big[1],
  xlim.big[2], ylim.big[2],
  xlim.big[1], ylim.big[2]
), 4, 2, byrow = TRUE)
domain.xy_ <- rbind(loc1, loc1[1, ])
domain_ <- st_sfc(st_polygon(list(domain.xy_)))

## Two separate barriers (sf)
barrier1_ <- local_square_sfc(
  xlim = c(xlim.big[1] - 1, x.mid - smalldist / 2),
  ylim = x.mid + width[1] * c(-0.5, 0.5))
barrier2_ <- local_square_sfc(
  xlim = c(x.mid + smalldist / 2, xlim.big[2] + 1),
  ylim = x.mid + width[2] * c(-0.5, 0.5))

## Individual barrier sfc objects (for fm_contains)
barriers_1 <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier1_[[1]])[, 1:2])))))
barriers_2 <- st_sfc(st_multipolygon(list(
  st_polygon(list(st_coordinates(barrier2_[[1]])[, 1:2])))))

## Mesh
mesh <- fm_mesh_2d(
  loc.domain = domain.xy_,
  max.edge   = max.edge.length,
  offset     = 1)
cat("Mesh nodes:", mesh$n, "\n")

## Identify barrier triangles — as a LIST (one vector per barrier)
triBarrier_1and2 <- list(
  unlist(fm_contains(x = barriers_1, y = mesh, type = "centroid")),
  unlist(fm_contains(x = barriers_2, y = mesh, type = "centroid"))
)
cat("Barrier 1 triangles:", length(triBarrier_1and2[[1]]), "\n")
cat("Barrier 2 triangles:", length(triBarrier_1and2[[2]]), "\n")

## ============================================================================
## PART 2: Build per-region FEM matrices
## ============================================================================

n_vert <- mesh$n
tv     <- mesh$graph$tv
coords <- mesh$loc[, 1:2]
n_tri  <- nrow(tv)

## Triangle areas
v1 <- coords[tv[, 1], ]; v2 <- coords[tv[, 2], ]; v3 <- coords[tv[, 3], ]
areas <- abs((v2[, 1] - v1[, 1]) * (v3[, 2] - v1[, 2]) -
               (v3[, 1] - v1[, 1]) * (v2[, 2] - v1[, 2])) / 2

## Classify triangles
allTri    <- seq_len(n_tri)
triB1     <- triBarrier_1and2[[1]]
triB2     <- triBarrier_1and2[[2]]
triDomain <- setdiff(allTri, union(triB1, triB2))

## Build lumped mass diagonals per region
Cd_diag  <- numeric(n_vert)
Cb1_diag <- numeric(n_vert)
Cb2_diag <- numeric(n_vert)

for (t in triDomain) {
  contrib <- areas[t] / 3
  for (j in 1:3) Cd_diag[tv[t, j]] <- Cd_diag[tv[t, j]] + contrib
}
for (t in triB1) {
  contrib <- areas[t] / 3
  for (j in 1:3) Cb1_diag[tv[t, j]] <- Cb1_diag[tv[t, j]] + contrib
}
for (t in triB2) {
  contrib <- areas[t] / 3
  for (j in 1:3) Cb2_diag[tv[t, j]] <- Cb2_diag[tv[t, j]] + contrib
}

Cd  <- Diagonal(n_vert, Cd_diag)
Cb1 <- Diagonal(n_vert, Cb1_diag)
Cb2 <- Diagonal(n_vert, Cb2_diag)
C_full_diag <- Cd_diag + Cb1_diag + Cb2_diag
C_full      <- Diagonal(n_vert, C_full_diag)

## Global stiffness matrix
fem_global <- fm_fem(mesh, order = 2)
G <- fem_global$g1

## Sanity check
cat("Max diff vs fm_fem c0:", max(abs(C_full_diag - fem_global$c0)), "\n")

## ============================================================================
## PART 3: Precision matrix builders (3 ranges)
## ============================================================================

#' Build barrier Q for nu = 1 (alpha = 2), 3 ranges
#'
#' Q = tau^2 * L * C^{-1} * L
#' where L = kappa_d^2 * C_d + kappa_b1^2 * C_b1 + kappa_b2^2 * C_b2 + G
#' and   tau^2 = 1 / (4*pi * kappa_d^2 * sigma^2)
build_barrier_Q_nu1 <- function(Cd, Cb1, Cb2, G, C_diag, ranges, sigma) {
  nu <- 1
  kappa_d  <- sqrt(8 * nu) / ranges[1]
  kappa_b1 <- sqrt(8 * nu) / ranges[2]
  kappa_b2 <- sqrt(8 * nu) / ranges[3]

  tau2 <- 1 / (4 * pi * kappa_d^2 * sigma^2)

  L    <- kappa_d^2 * Cd + kappa_b1^2 * Cb1 + kappa_b2^2 * Cb2 + G
  Cinv <- Diagonal(length(C_diag), 1 / C_diag)

  Q <- tau2 * L %*% Cinv %*% L
  Q <- (Q + t(Q)) / 2
  return(Q)
}

#' Build barrier Q for nu = 2 (alpha = 3), 3 ranges
#'
#' Q = tau^2 * L * C^{-1} * L * C^{-1} * L
#' where tau^2 = 1 / (8*pi * kappa_d^4 * sigma^2)
build_barrier_Q_nu2 <- function(Cd, Cb1, Cb2, G, C_diag, ranges, sigma) {
  nu <- 2
  kappa_d  <- sqrt(8 * nu) / ranges[1]
  kappa_b1 <- sqrt(8 * nu) / ranges[2]
  kappa_b2 <- sqrt(8 * nu) / ranges[3]

  tau2 <- 1 / (8 * pi * kappa_d^4 * sigma^2)

  L    <- kappa_d^2 * Cd + kappa_b1^2 * Cb1 + kappa_b2^2 * Cb2 + G
  Cinv <- Diagonal(length(C_diag), 1 / C_diag)

  A <- L %*% Cinv
  Q <- tau2 * A %*% A %*% L
  Q <- (Q + t(Q)) / 2
  return(Q)
}

## ============================================================================
## PART 4: Define the two range-fraction scenarios
## ============================================================================

## Scenario "eq":  equal fractions in both barriers
## Scenario "df":  different fractions per barrier
frac_eq <- c(0.01, 0.01)     # equal
frac_df <- c(0.01, 0.7)    # different

ranges_eq <- range_domain * c(1, frac_eq)
ranges_df <- range_domain * c(1, frac_df)

cat("\n--- True ranges (equal fractions) ---\n")
cat("Domain:", ranges_eq[1],
    " Barrier1:", ranges_eq[2],
    " Barrier2:", ranges_eq[3], "\n")
cat("\n--- True ranges (different fractions) ---\n")
cat("Domain:", ranges_df[1],
    " Barrier1:", ranges_df[2],
    " Barrier2:", ranges_df[3], "\n")

## ============================================================================
## PART 5: Build precision matrices — 4 variants
## ============================================================================

## Reference Q from INLA's built-in for nu=1 (for validation)
bfem_1y2 <- mesh2fem.barrier(mesh, triBarrier_1and2)
Q_nu1_ref_eq <- inla.barrier.q(bfem_1y2, ranges = ranges_eq, sigma = sigma_true)
Q_nu1_ref_df <- inla.barrier.q(bfem_1y2, ranges = ranges_df, sigma = sigma_true)

## Custom builders
Q_nu1_eq <- build_barrier_Q_nu1(Cd, Cb1, Cb2, G, C_full_diag, ranges_eq, sigma_true)
Q_nu1_df <- build_barrier_Q_nu1(Cd, Cb1, Cb2, G, C_full_diag, ranges_df, sigma_true)
Q_nu2_eq <- build_barrier_Q_nu2(Cd, Cb1, Cb2, G, C_full_diag, ranges_eq, sigma_true)
Q_nu2_df <- build_barrier_Q_nu2(Cd, Cb1, Cb2, G, C_full_diag, ranges_df, sigma_true)

cat("\n--- Dimension / nnz ---\n")
cat("Q_nu1_eq:", dim(Q_nu1_eq), " nnz:", nnzero(Q_nu1_eq), "\n")
cat("Q_nu1_df:", dim(Q_nu1_df), " nnz:", nnzero(Q_nu1_df), "\n")
cat("Q_nu2_eq:", dim(Q_nu2_eq), " nnz:", nnzero(Q_nu2_eq), "\n")
cat("Q_nu2_df:", dim(Q_nu2_df), " nnz:", nnzero(Q_nu2_df), "\n")

## Sanity: our nu=1 vs INLA reference
cat("\nMax |Q_nu1_eq - ref|:", max(abs(Q_nu1_eq - Q_nu1_ref_eq)), "\n")
cat("Max |Q_nu1_df - ref|:", max(abs(Q_nu1_df - Q_nu1_ref_df)), "\n")

## ============================================================================
## PART 6: Simulate fields (4 variants)
## ============================================================================

set.seed(1)
u_nu1_eq <- inla.qsample(1, Q_nu1_eq, seed = 1)[, 1]
u_nu1_df <- inla.qsample(1, Q_nu1_df, seed = 2)[, 1]
u_nu2_eq <- inla.qsample(1, Q_nu2_eq, seed = 3)[, 1]
u_nu2_df <- inla.qsample(1, Q_nu2_df, seed = 4)[, 1]

## ============================================================================
## PART 7: Visualise simulated fields
## ============================================================================

bb <- matrix(c(3, 17, 3, 17), nrow = 2, ncol = 2, byrow = TRUE)

pgrid <- fm_evaluator(
  mesh,
  lattice = fm_evaluator_lattice(
    mesh, xlim = bb[1, ], ylim = bb[2, ], dims = c(300, 300)))

grid.df <- data.frame(
  x = rep(pgrid$x, times = length(pgrid$y)),
  y = rep(pgrid$y, each  = length(pgrid$x)))

gg0 <- ggplot() + xlab("") + ylab("") + theme_minimal() + coord_fixed()
gg.add <- list(scale_fill_viridis_c(option = "turbo", na.value = "transparent"))

add.barriers <- list(
  geom_sf(data = barriers_1, fill = NA, colour = "black", linewidth = 0.5),
  geom_sf(data = barriers_2, fill = NA, colour = "black", linewidth = 0.5)
)

grid.df$u_nu1_eq <- as.vector(fm_evaluate(pgrid, field = u_nu1_eq))
grid.df$u_nu1_df <- as.vector(fm_evaluate(pgrid, field = u_nu1_df))
grid.df$u_nu2_eq <- as.vector(fm_evaluate(pgrid, field = u_nu2_eq))
grid.df$u_nu2_df <- as.vector(fm_evaluate(pgrid, field = u_nu2_df))

p1 <- gg0 + geom_raster(data = grid.df, aes(x, y, fill = u_nu1_eq)) +
  gg.add + add.barriers +
  ggtitle(expression(paste(nu, "=1, equal fracs (0.01, 0.01)")))
p2 <- gg0 + geom_raster(data = grid.df, aes(x, y, fill = u_nu1_df)) +
  gg.add + add.barriers +
  ggtitle(expression(paste(nu, "=1, diff fracs (0.01, 0.7)")))
p3 <- gg0 + geom_raster(data = grid.df, aes(x, y, fill = u_nu2_eq)) +
  gg.add + add.barriers +
  ggtitle(expression(paste(nu, "=2, equal fracs (0.01, 0.01)")))
p4 <- gg0 + geom_raster(data = grid.df, aes(x, y, fill = u_nu2_df)) +
  gg.add + add.barriers +
  ggtitle(expression(paste(nu, "=2, diff fracs (0.01, 0.7)")))

print((p1 + p2) / (p3 + p4))

## ============================================================================
## PART 8: Correlation comparison
## ============================================================================

localCorrel <- function(locs, mesh, Q) {
  nl <- nrow(locs)
  ii <- sapply(1:nl, function(i)
    which.min(rowSums(sweep(
      mesh$loc[, 1:ncol(locs)], 2, locs[i, ], "-")^2)))
  b <- matrix(0, nrow(Q), nl)
  for (i in 1:nl) b[ii[i], i] <- 1
  cc <- inla.qsolve(Q, b)
  s  <- sqrt(diag(inla.qinv(Q)))
  for (i in 1:nl) cc[, i] <- cc[, i] / (s * s[ii[i]])
  return(drop(cc))
}

## Reference points: one in each sub-domain + one in the canal
locs_ref <- cbind(
  c(5, 15, 10),
  c(8,  8, 10))

cat("\nComputing correlations (4 models)...\n")
corr_nu1_eq <- localCorrel(locs_ref, mesh, Q_nu1_eq)
corr_nu1_df <- localCorrel(locs_ref, mesh, Q_nu1_df)
corr_nu2_eq <- localCorrel(locs_ref, mesh, Q_nu2_eq)
corr_nu2_df <- localCorrel(locs_ref, mesh, Q_nu2_df)

gcorr_nu1_eq <- as.matrix(fm_evaluate(pgrid, field = corr_nu1_eq))
gcorr_nu1_df <- as.matrix(fm_evaluate(pgrid, field = corr_nu1_df))
gcorr_nu2_eq <- as.matrix(fm_evaluate(pgrid, field = corr_nu2_eq))
gcorr_nu2_df <- as.matrix(fm_evaluate(pgrid, field = corr_nu2_df))

make_corr_df <- function(gcorr, grid.df, locs, label) {
  do.call(rbind, lapply(1:nrow(locs), function(l)
    data.frame(grid.df[, c("x", "y")],
               loc = paste0("(", locs[l, 1], ", ", locs[l, 2], ")"),
               correlation = gcorr[, l],
               model = label)))
}

corr_all <- rbind(
  make_corr_df(gcorr_nu1_eq, grid.df, locs_ref, "nu=1 eq(0.01,0.01)"),
  make_corr_df(gcorr_nu1_df, grid.df, locs_ref, "nu=1 df(0.01,0.7)"),
  make_corr_df(gcorr_nu2_eq, grid.df, locs_ref, "nu=2 eq(0.01,0.01)"),
  make_corr_df(gcorr_nu2_df, grid.df, locs_ref, "nu=2 df(0.01,0.7)"))

p_corr <- gg0 +
  geom_raster(
    data = corr_all[corr_all$correlation > 0.1, ],
    mapping = aes(x = x, y = y, fill = correlation)) +
  facet_grid(model ~ loc) +
  scale_fill_viridis_c(option = "magma", na.value = "transparent") +
  gg.add +
  add.barriers +
  ggtitle("Correlation from reference points")
print(p_corr)

