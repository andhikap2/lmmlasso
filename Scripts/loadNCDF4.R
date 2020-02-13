#Reading netCDF4 in raster r

# load package
 library(sp)
 library(raster)
 library(ncdf4)

 # read ncdf file
 nc<-nc_open("tn_ens_mean_0.1deg_reg_v19.0eHOM.nc")
)

 # extract variable name, size and dimension
 v <- nc$var[[1]]
 size <- v$varsize #the three dimensions (705, 465, 25171)
 dims <- v$ndims #3
 nt <- size[dims]              # length of time dimension
 lat <- nc$dim$latitude$vals   # latitude position
 lon <- nc$dim$longitude$vals  # longitude position

 # read sst variable
 r<-list()

 for (i in 1:nt) {
   start <- rep(1,dims)     # begin with start=(1,1,...,1)
   start[dims] <- i             # change to start=(1,1,...,i) to read    timestep i
   count <- size                # begin with count=(nx,ny,...,nt), reads entire var
   count[dims] <- 1             # change to count=(nx,ny,...,1) to read 1 tstep

   dt<-ncvar_get(nc, varid = 'tn', start = start, count = count)

   # convert to raster
   r[i]<-raster(dt)
 }

 # create layer stack with time dimension
 r<-stack(r)

 # transpose the raster to have correct orientation
 rt<-t(r)
 rt<-flip(rt,2) #flip to get the right shape

 extent(rt)<-extent(c(range(lon), range(lat)))

 # plot the result
 spplot(rt)
library(maps)
map('world',add=TRUE)



#writing multiRaster stack to tif
library(raster)
a=stack("tn_ens_mean_0.1deg_reg_v19.0eHOM.nc.1")
writeRaster(a,filename="tn_multilayer.tif") #takes forever