# R Functions for computing the topographical variables for the 
# Tai Po Kau Forest Dynamics Plot, Hong Kong SAR

# By Jinlong Zhang 
# Tue Jan 31 14:29:58 2023


### 将规则点的数据， 转换为矩阵
### 必须是非常规则的点，x，y两个坐标， 都是等距的规则点， 该函数可以将有坐标的规则的点，排成矩阵
as.matrix.extract <- function(x, xcol = 1, ycol = 2, atrcol = 3){
   xcol.dat <- as.numeric(x[, xcol]) ## 提取横坐标列
   ycol.dat <- as.numeric(x[, ycol]) ## 提取纵坐标
   x <- x[order(ycol.dat, xcol.dat), ]  ### 按照行， 列， 排序
   res <- x[, atrcol]  ## 提取排序后的数据列
   if(is.factor(res)){
       res <- as.character(res)
   }
   dim(res) <- c(length(unique(xcol.dat)), length(unique(ycol.dat)))  ### 定义坐标维度
   tres <- t(res)  ## 行列转置
   ### image(tres)
   return(res)
}

################
## Calculating the slopes for each grid in a forest dynamic plot
## Provide the matrix of altitudes
## Provide the cell_size, the altitudes and cell_size must be in the same unit. 
## Algorithms could be found at ## 
## http://webhelp.esri.com/arcgiSDEsktop/9.3/body.cfm?tocVisable=1&ID=-1&TopicName=How%20Slope%20works
## Value: a matrix representing the degrees of each grid. 

get_slope <- function(z.mat, cell_size = 5){
     ### Extent the z.mat matrix using the value from the border, so that the slope for the out most cells could be calculated. 
     z.mat.colbind <- cbind(z.mat[,1], z.mat, z.mat[,ncol(z.mat)])
	 z.mat.rbind <- rbind(c(z.mat[1,1], z.mat[1,], z.mat[1,ncol(z.mat)]),z.mat.colbind,c(z.mat[nrow(z.mat),1], z.mat[nrow(z.mat),], z.mat[nrow(z.mat), ncol(z.mat)]) )
     
	 z.mat.ext <- z.mat.rbind
	 nrz <- nrow(z.mat.ext)
	 ncz <- ncol(z.mat.ext)
	 slope.mat <- rep(NA, nrz*ncz)
	 dim(slope.mat) <- c(nrz, ncz)  ## Create a Null matrix used in loop
	 for(m in 2:(nrz-1)){
	    for(n in 2:(ncz-1)){
		    a = z.mat.ext[m-1,n-1]
		    b = z.mat.ext[m-1,n]
		    c = z.mat.ext[m-1,n+1]
		    d = z.mat.ext[m,n-1]
		    e = z.mat.ext[m,n]
		    f = z.mat.ext[m,n+1]
		    g = z.mat.ext[m+1,n-1]
		    h = z.mat.ext[m+1,n]
		    i = z.mat.ext[m+1,n+1]
		    
		    dz_dx = ((c + 2*f + i) - (a + 2*d + g)) / (8 * cell_size)
		    dz_dy = ((g + 2*h + i) - (a + 2*b + c)) / (8 * cell_size)
			slope.mat[m,n] <- atan (sqrt ( dz_dx^2 + dz_dy^2 ))*(180/pi)  
		}
	}
	 return(slope.mat[c(-1, -nrow(slope.mat)),c(-1, -ncol(slope.mat))])
}
 


## test.dat.slope <- c(50, 30, 8, 45, 30, 10, 50, 30, 10)
## dim(test.dat.slope) <- c(3, 3)
## get_slope(test.dat.slope)

################################################################### 
## Calculating the aspect of each cell within a forest dynamic plot
## Adopted from ESRI http://webhelp.esri.com/arcgisdesktop/9.3/body.cfm?tocVisable=1&ID=-1&TopicName=how%20aspect%20works 
## the z.mat is a numerical matrix representing the altitude for each cell
## the output is a numerical matrix representing the degree of aspect. Starting from the north, 
## This function assums that The direction of the column, from greater to smaller columns points to the north. 

get_aspect <- function(z.mat){
     z.mat.colbind <- cbind(z.mat[,1], z.mat, z.mat[,ncol(z.mat)])
	 z.mat.rbind <- rbind(c(z.mat[1,1], z.mat[1,], z.mat[1,ncol(z.mat)]),z.mat.colbind,c(z.mat[nrow(z.mat),1], z.mat[nrow(z.mat),], z.mat[nrow(z.mat), ncol(z.mat)]) )
	 z.mat.ext <- z.mat.rbind
	 nrz <- nrow(z.mat.ext)
	 ncz <- ncol(z.mat.ext)
	 aspect.mat <- rep(NA, nrz*ncz)
	 dim(aspect.mat) <- c(nrz, ncz)
	 for(m in 2:(nrz-1)){
	    for(n in 2:(ncz-1)){
		    a = z.mat.ext[m-1,n-1]
		    b = z.mat.ext[m-1,n]
		    c = z.mat.ext[m-1,n+1]
		    d = z.mat.ext[m,n-1]
		    e = z.mat.ext[m,n]
		    f = z.mat.ext[m,n+1]
		    g = z.mat.ext[m+1,n-1]
		    h = z.mat.ext[m+1,n]
		    i = z.mat.ext[m+1,n+1]
            dz_dx = ((c + 2*f + i) - (a + 2*d + g))/8
            dz_dy = ((g + 2*h + i) - (a + 2*b + c))/8
            aspect.mat[m, n] = (180/pi) * atan2(dz_dy, -dz_dx)
        }
    }
        format.aspect <- function(aspect){
           cell = aspect
           for(i in 1:length(aspect)){
           if (aspect[i] < 0) {        
              cell[i] = 90.0 - aspect[i]        
              }        
           else if (aspect[i] > 90.0) {        
              cell[i] = 360.0 - aspect[i] + 90.0        
              }        
           else {        
              cell[i] = 90.0 - aspect[i]         
              }        
           }
           return(cell)        
        } 
    aspect.mat <- aspect.mat[c(-1, -nrow(aspect.mat)),c(-1, -ncol(aspect.mat))]
    return(format.aspect(aspect.mat))
}

###  test.dat.aspect <- c(101, 101, 101, 92, 92, 91, 85, 85, 84)
###  dim(test.dat.aspect) <- c(3, 3)
###  get_aspect(test.dat.aspect)



######################## Calculating  Convexity ################### 

get_convexity <- function(z.mat){
     z.mat.colbind <- cbind(rep(NA, nrow(z.mat)), z.mat, rep(NA, nrow(z.mat)))
	 z.mat.rbind <- rbind(rep(NA, ncol(z.mat)+2), z.mat.colbind, rep(NA, ncol(z.mat)+2) )
	 z.mat.ext <- z.mat.rbind
	 nrz <- nrow(z.mat.ext)
	 ncz <- ncol(z.mat.ext)
	 convexity.mat <- rep(NA, nrz*ncz)
	 dim(convexity.mat) <- c(nrz, ncz)
	 for(m in 2:(nrz-1)){
	    for(n in 2:(ncz-1)){
		    a = z.mat.ext[m-1,n-1]
		    b = z.mat.ext[m-1,n]
		    c = z.mat.ext[m-1,n+1]
		    d = z.mat.ext[m,n-1]
		    e = z.mat.ext[m,n]
		    f = z.mat.ext[m,n+1]
		    g = z.mat.ext[m+1,n-1]
		    h = z.mat.ext[m+1,n]
		    i = z.mat.ext[m+1,n+1]
			#### If na.....
			all.neighbour <- c(a, b, c, d, f, g, h, i)
            mean.neighbour <- mean(na.omit(all.neighbour))
            convexity.mat[m,n] <- (e - mean.neighbour)
        }
    }
    convexity.mat <- convexity.mat[c(-1, -nrow(convexity.mat)),c(-1, -ncol(convexity.mat))]
    return(convexity.mat)
}
