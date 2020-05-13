library(rgl)
library("plot3D")
library(plot3Drgl)

phi <- function(x, y, z) {
  phi1 <- 3*x^2+y^2-1
  phi2 <- x^2+(y-0.5)^2-0.5
  phi3 <- 3*x^2+y^2-z-1
  return(max(c(phi1, phi2, phi3)))
}
grid.x <- seq(-sqrt(3), sqrt(3), length = 20)
grid.y <- seq((1-sqrt(2))/2, 1, length = 10)
grid.z <- seq(-1, 0, length = 10)

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
        omega.color = append(omega.color, -phi(x, y, z))
      }
    }
  }
}

scatter3D(x= omega.x,y=omega.y,z=omega.z,  pch = 16, colvar=omega.color)
plotrgl(xlim = c(-sqrt(3), sqrt(3)),
        ylim = c((1-sqrt(2))/2, 1),
        zlim = c(-1, 0),
        xlab = '', ylab = '', zlab = '',
        type = 'n')

# Add planes
planes3d(0, -1, 0, -1, col = 'red', alpha = 0.6)
planes3d(0, 1, 0, 0, col = 'orange', alpha = 0.6)
planes3d(0, -2, 0, 0.25, col = 'blue', alpha = 0.6)
