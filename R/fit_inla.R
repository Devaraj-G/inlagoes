library(INLA)
library(fmesher)
library(terra)

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
# Create mesh
mesh <- fm_mesh_2d_inla(
  loc = xy,
  max.edge = c(2,0.75),
  cutoff = 0.1
)
plot(mesh); points(xy, pch = 19, col = 'red')

# SPDE model/stack

spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
A <- inla.spde.make.A(mesh, loc = xy)

fit_stack <- inla.stack(
  data = list(y = xy_obs1),
  A = list(A,matrix(1, nrow = nrow(xy), ncol = 1)),
  effects = list(spatial = 1:spde$n.spde, data.frame(intercept = 1))
)

formula <- y ~0 + f(spatial, model = spde)

# fit INLA

fit_inla <- inla(
  formula,
  data = inla.stack.data(fit_stack),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(fit_stack), compute = TRUE),
  control.compute=list(config=TRUE)
)

bru_components <- ~0 + f(spatial, model = spde)
bru_formula <- count ~ .

fit_bru <- bru(
bru_components,
bru_obs(
  formula = bru_formula,
family = 'gaussian',
data = xy_obs1,
),
options = list(
  control.compute = list(waic = TRUE, cpo = FALSE),
  control.inla = list(
    int.strategy = "eb"
  )
)
)


fit_inla$summary.fixed
fit_inla$summary.hyperpar

proj = inla.mesh.projector(mesh, xlim=c(-113,-89), ylim=c(4.9, 20.5), dim=c(100,100))
msurface = inla.mesh.project(proj, fit_inla$summary.random$spatial[,"mean"])
ssurface = inla.mesh.project(proj, fit_inla$summary.random$spatial[,"sd"])

out.p = c();
out.p$x = seq(-109, -90, length.out=100)
out.p$y = seq(5.5, 20, length.out=100)
out.p$z = msurface
out.s <- out.p
out.s$z = ssurface

#x <- as.matrix(cbind(out.p$x,out.p$y))
#z <- as.matrix(out.p$z)
ggmatplot()

ggplot() +
  gg(msurface) +
  scale_fill_viridis_c()
  geom_point(aes(mesh$loc[,1],mesh$loc[,2]))
ggplot() +
  gg(ssurface) +
  scale_fill_viridis_c()


sst_pred_mean <- fit_inla$summary.fitted.values$mean
sst_pred <- rast(xy_obs)
values(sst_pred) <- sst_pred_mean

plot(sst_pred, main = "predicted")
points(xy, pch=19, col="black")

# create grid for prediction
xy_pr <- spatSample(data_cropped, 10000, method = 'random', as.points = TRUE)

grid <- crds(xy_pr)
A_pred <- inla.spde.make.A(mesh, loc = as.matrix(grid))

stack_pred <- inla.stack(
  data = list(y = NA),
  A = list(A_pred, matrix(1, nrow = nrow(xy_pr), ncol = 1)),
  effects = list(spatial = 1:spde$n.spde,
                 data.frame(intercept = 1))
)

# predict (most probably redundant, info already in fit_inla)
pred_inla <- inla(
  formula,
  data = inla.stack.data(stack_pred),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(stack_pred), compute = TRUE)
)

sst_pred_mean <- pred_inla$summary.fitted.values$mean
sst_pred <- rast(xy_pr)
values(sst_pred) <- sst_pred_mean

plot(sst_pred, main = "predicted")
points(xy, pch=19, col="black")

# samples
samples <- inla.posterior.sample(5, fit_inla)

# Example: extract posterior mean field for each sample
mean_samples <- sapply(samples, function(s) mean(s$latent))
hist(mean_samples, col="skyblue", main="mu")

sst_sd <- rast(xy_pr)
sst_sd[] <- pred_inla$summary.fitted.values$sd
plot(sst_sd, main = "sd")

# Plot

field_samples <- sapply(samples, function(s) {
  s$latent
})
dim(field_samples)

pred_xy <- xyFromCell(data_cropped, 1:ncell(data_cropped))
A_pred <- inla.spde.make.A(mesh, loc = grid)

# Compute predictions for each sample
pred_matrix <- A_pred %*% field_samples

par(mfrow=c(1,3))
for (i in 1:3) {
  r_samp <- data_cropped
  values(r_samp) <- pred_matrix[, i]
  plot(r_samp, main=paste("Posterior sample", i))
}
