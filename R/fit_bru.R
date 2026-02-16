library(INLA)
library(inlabru)
library(fmesher)
library(ggplot2)
library(patchwork)

source("./R/convert_units.R")
source("./R/load_data.R")
source("./R/plot_data.R")

# Load data
filename = './data/noaa-goes16/ABI-L2-SSTF/2023/297/12/OR_ABI-L2-SSTF-M6_G16_s20232971200209_e20232971259517_c20232971304395.nc'
limits = c(-110,-90,5,20)

data=load_data(filename,'DQF')
plot_data(data)
data_cropped = crop_data(data,limits)
plot_data(data_cropped)
dqf = data_cropped

data=load_data(filename,'SST')
data=data-273
plot_data(data)
data_cropped = crop_data(data,limits)
plot_data(data_cropped)
data_good = data_cropped*(dqf == 0)
plot_data(data_good)
data_bad = data_cropped*(dqf != 0)
plot_data(data_bad)
data_cropped = data_good
data_good = data_cropped
data_good[dqf != 0] = NA
plot_data(data_good)

# Pick some samples from SST data
xy_obs <- spatSample(data_good, 10000, method = 'random', as.points = TRUE)
xy_obs <- xy_obs[!is.na(values(xy_obs)),]

xy_obs1 <- extract(data_cropped, xy_obs)[,2]
xy_obs2 <- extract(data_good, xy_obs)[,2] # same as points are not on NA points

plot_data(data_good)
points(xy_obs, pch = 20, col = 'red')

# Prepare inputs for INLA
xy <- geom(xy_obs)[ , c('x','y')]
data_inla <- data.frame(sst = xy_obs1,
                        x = xy[ , 1],
                        y = xy[ , 2])
sst=xy_obs1

sst <- sf::st_as_sf(
  data.frame(easting = xy[,1], northing = xy[,2]),
  coords = c("easting", "northing")
)
sst$observed = xy_obs1

# mesh
mesh <- fm_mesh_2d(
  loc = xy,
  max.edge = c(5,0.5),
  cutoff = 0.1
)
meshpic <- ggplot() +
  geom_fm(data = mesh) +
  geom_point(aes(x = xy[,1], y = xy[,2]), size = 0.1)

matern <-
  inla.spde2.pcmatern(mesh,
                      prior.sigma = c(1, 0.01),
                      prior.range = c(1, 0.01)
  )

# model
formula <- observed ~ field(geometry, model = matern) + Intercept(1)

fit <- bru(formula, sst, family = "gaussian")
summary(fit)

# predict
pix <- fm_pixels(mesh, dims = c(100, 100))
pred <- predict(
  fit,
  pix,
  ~ field + Intercept
)

samp <- generate(
  fit,
  pix,
  ~ field + Intercept,
  n.samples = 4
)


pl_posterior_mean <- ggplot() +
  gg(pred, geom = "tile")+
ggtitle("Posterior mean")

pred$sample <- samp[, 1]
pl_posterior_mean1 <- ggplot() +
  gg(pred, aes(fill = sample), geom = "tile") +
  scale_color_gradientn(colors = terrain.colors(10), limits = c(10, 40))

pred$sample <- samp[, 2]
pl_posterior_mean2 <- ggplot() +
  gg(pred, aes(fill = sample), geom = "tile") +
  scale_color_gradientn(colors = terrain.colors(10), limits = c(10, 40))

pred$sample <- samp[, 3]
pl_posterior_mean3 <- ggplot() +
    gg(pred, aes(fill = sample), geom = "tile") +
  scale_color_gradientn(colors = terrain.colors(10), limits = c(10, 40))

pred$sample <- samp[, 4]
pl_posterior_mean4 <- ggplot() +
    gg(pred, aes(fill = sample), geom = "tile") +
  scale_color_viridis_c()

csc <- scale_color_gradientn(colors = terrain.colors(10), limits = c(10, 40))

(meshpic  | pl_posterior_mean2  ) / (pl_posterior_mean3  | pl_posterior_mean4  )

(pl_posterior_mean1 + csc | pl_posterior_mean2 + csc ) / (pl_posterior_mean3 + csc  | pl_posterior_mean4 + csc )
