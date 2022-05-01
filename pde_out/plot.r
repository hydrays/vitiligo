require(lattice)
require(rasterVis)
require(latticeExtra)
library(viridisLite)
coul <- viridis(100)
require("RColorBrewer")
##color Scheme I
chemokine.cols <- colorRampPalette(c('white', '#F8766D'))(n=100)
fb.col <- "red"
tcell.col <- "#1F78B4"

## Color Scheme II
##chemokine.cols <- colorRampPalette(c('white', '#1F78B4'))(n=100)
##fb.col <- "#0000FF"
##tcell.col <- "#F8766D"

## ## Color Scheme III
## chemokine.cols <- colorRampPalette(c('white', '#FFFF00'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"

## ## Color Scheme IV
## chemokine.cols <- colorRampPalette(c('white', '#7CAE00'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"

## ## Color Scheme V
## chemokine.cols <- colorRampPalette(c('white', '#F7931E'))(n=100)
## fb.col <- "#1F78B4"
## tcell.col <- "#F8766D"


##L <- 1024
L  <- 100
##L  <- 256
#L  <- 128

for( i in seq(0, 2000, by=1) )
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("phi_", padded_i, ".png", sep=''), height=600, width=1200)

    A <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
    p1  <- levelplot(z1)
    
    B <- matrix(unlist(read.csv(paste('eta', padded_i, '.dat', sep=''), header=F)), nrow=L)
    z2 = raster(B, xmn=0, xmx=L, ymn=0, ymx=L)
    p2 <- levelplot(z2)
    
    C <- matrix(unlist(read.csv(paste('psi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    z3 = raster(C, xmn=0, xmx=L, ymn=0, ymx=L)
    p3 <- levelplot(z3)
    
    print(p1, position=c(0, 0.0, 0.3, 1), more=TRUE)
    print(p2, position=c(0.3, 0.0, 0.6, 1), more=TRUE)
    print(p3, position=c(0.6, 0.0, 0.9, 1))    
    dev.off()
}

