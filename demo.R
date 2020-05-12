library(radish2)
library(sp)
xx <- cbind(-5:5/10, (-5:5/10)^2)
coords <- cbind(0.5, 0:10/11+1/20)
uno  <- Lattice1D(xx, sp::SpatialPoints(as.matrix(coords)), radish_conductance_surface)
tres <- Lattice1D(xx, sp::SpatialPoints(as.matrix(coords)), radish_conductance_surface3)

E <- radish_distance(loglinear_conductance, uno, matrix(0, 1, 2), TRUE)$covariance[,,1]

S <- wishart_simulate_distance(1, E, 30)

# does original algorithm work
uno_fit <- radish_algorithm(loglinear_conductance, leastsquares, uno, S, objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE, validate=TRUE)
max(abs(uno_fit$gradient - uno_fit$num_gradient))
max(abs(uno_fit$hessian - uno_fit$num_hessian))
max(abs(uno_fit$partial_X - uno_fit$num_partial_X))
max(abs(uno_fit$partial_S[lower.tri(S)] - uno_fit$num_partial_S[lower.tri(S)]))

# what about via spam
dos_fit <- radish_algorithm2(loglinear_conductance, leastsquares, dos, S, objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE)

# what about via Matrix
tres_fit <- radish_algorithm3(loglinear_conductance, leastsquares, tres, S, objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE)

# speed test on more realistic data
library(radish2)
library(spam64)
library(sp)
NN <- 200^2
xx <- cbind(-NN:NN/NN, (-NN:NN/NN)^2)
rr <- raster::stack(raster::raster(as.matrix(xx[,1],200,200)), raster::raster(as.matrix(xx[,2],200,200)))
set.seed(1); coords <- cbind(runif(10),runif(10))
uno  <- radish_conductance_surface(rr, sp::SpatialPoints(as.matrix(coords)))
dos  <- radish_conductance_surface2(rr, sp::SpatialPoints(as.matrix(coords)))
tres  <- radish_conductance_surface3(rr, sp::SpatialPoints(as.matrix(coords)))

E <- radish_distance(loglinear_conductance, uno, matrix(0, 1, 2), TRUE)$covariance[,,1]
S <- wishart_simulate_distance(1, E, 30)

system.time(uno_fit <- radish_algorithm(loglinear_conductance, leastsquares, uno, S, objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE))
system.time(dos_fit <- radish_algorithm2(loglinear_conductance, leastsquares, dos, S, objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE)) #spam not worth it
system.time(tres_fit <- radish_algorithm3(loglinear_conductance, mlpe, tres, S, objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE))
