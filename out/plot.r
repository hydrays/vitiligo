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
L  <- 512
##L  <- 256
#L  <- 128
colors <- c("#A7A7A7",
            "dodgerblue",
            "firebrick",
            "forestgreen",
            "gold")
myAt <- seq(5) - 1.5

for( i in seq(0, 2000, by=10) )
{
    cat(i, '\n')
    padded_i <- sprintf("%05d", i)
    png(paste("config_", padded_i, ".png", sep=''), height=600, width=1200)
    ##pdf(paste("config_", padded_i, ".pdf", sep=''), height=6, width=12)
    A <- matrix(unlist(read.csv(paste('c', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ##A[A==2] <- sample(c(0,2), length(A[A==2]), replace = TRUE, prob=c(0.9, 0.1))
    A[A==2] <- 0
    z1 = raster(A, xmn=0, xmx=L, ymn=0, ymx=L)
    znull = z1
    znull[] = 0
    p1  <- levelplot(znull, at=myAt, col.regions = c('white', '#C9CACA', '#F8766D', '#F8766D'), margin=FALSE,
                     colorkey=FALSE)
    coords <- xyFromCell(z1, which(z1[]==1))
    pp1 <- xyplot(coords[,2]~coords[,1], lwd=1, pch=8, cex=0.5, col="#C9CACA")
    coords <- xyFromCell(z1, which(z1[]==3))
    pp3 <- xyplot(coords[,2]~coords[,1], lwd=1, cex=0.25, col=tcell.col)    
    p1 <- p1 + as.layer(pp1) + as.layer(pp3)
    ##p1 <- as.layer(pp1) + p1
    
    ##B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    B <- matrix(unlist(read.csv(paste('phi', padded_i, '.dat', sep=''), header=F)), nrow=L)
    ##coords <- xyFromCell(z, which(z[]==3))
    z2 = raster(B, xmn=0, xmx=L, ymn=0, ymx=L)
    z2 = z2/5
    z2[1] = 1
    zlog2 <- log(z2)    
    p2 <- levelplot(zlog2, col.regions = chemokine.cols, margin=FALSE)
    C <- matrix(unlist(read.csv(paste('lambda', padded_i, '.dat', sep=''), header=F)), nrow=L)
    z3 = raster(C, xmn=0, xmx=L, ymn=0, ymx=L)    
    coords <- xyFromCell(z3, which(z3[]>1e-2))
    p3 <- xyplot(coords[,2]~coords[,1], cex=0.00001, col=fb.col)
    p4 <- p2 + as.layer(p3)
    
    print(p1, position=c(0, 0.0, .5, 1), more=TRUE)
    print(p4, position=c(0.5, 0.0, 1, 1))
    dev.off()
}

