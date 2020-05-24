library(rgl)
library("car")
library("plot3D")
library(plot3Drgl)

phi <- function(x, y, z) {
  phi1 <- 3*x^2+y^2-1
  phi2 <- x^2+(y-0.5)^2-0.5
  phi3 <- 3*x^2+y^2-z-1
  return(max(c(phi1, phi2, phi3)))
}
bounds.ax = -0.75
bounds.bx = 1
bounds.ay = -1
bounds.by = 1
bounds.az = -1
bounds.bz = 0

grid.x <- seq(bounds.ax, bounds.bx, length = 40)
grid.y <- seq(bounds.ay, bounds.by, length = 30)
grid.z <- seq(bounds.az, bounds.bz, length = 10)

omega.x = c()
omega.y = c()
omega.z = c()
omega.color = c()

for (x in grid.x){
  for (y in grid.y){
    for(z in grid.z){
      if(phi(x,y,z) <= 0){
        omega.x = append(omega.x, x)
        omega.y = append(omega.y, y)
        omega.z = append(omega.z, z)
        omega.color = append(omega.color, phi(x, y, z))
      }
    }
  }
}

# also scatter3Drgl could be used
scatter3D(x= omega.x,y=omega.y,z=omega.z,  pch = 16, colvar=omega.color)
plotrgl(xlim = c(bounds.ax, bounds.bx),
        ylim = c(bounds.ay, bounds.by),
        zlim = c(bounds.az, bounds.bz),
        xlab = '', ylab = '', zlab = '',
        type = 'n')

# Add planes
planes3d(-4.5,-0.41421356, -1., 2.730, col = 'red', alpha = 0.6)
planes3d(-2.19280904, -0.41421356, -1., 1.4435, col = 'orange', alpha = 0.6)
planes3d(-0.97903939, -0.41421356, -1., 1.122, col = 'blue', alpha = 0.6)
planes3d(-1.58592422,  0.65671751, -1., 1.3174, col = 'green', alpha = 0.6)

