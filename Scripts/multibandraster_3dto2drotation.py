#How do I properly flatten a 3d array into 2d :)
import os
import rasterio
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")
tx=rasterio.open('tx_gdal.tif')
arraytx=tx.read() #numpy array of values with dimensions (time, lat, long)
day=20729 #days since Jan 1 1950
day2=day+203
maxtemp=arraytx[day:day2,:,:] #(203, 465, 705)

#need the dataset to be 465 x 705 x 203

a=maxtemp.transpose((1,2,0)) #465 x 705 x 203

#reshape intop a 2d array

a=a.reshape(-1,203) # 327825 x 203
a=a.astype('float64')

#need to make this shape back into 203 465 705

a=a.transpose() # 203 x 327825
a_reshaped=a.reshape(203,465,-1) #back to the original shape



os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/E-OBSv19.0HOM")


naip_data_ras = arraytx
naip_meta = tx.profile

naip_transform = naip_meta["transform"]
naip_crs = naip_meta["crs"]


naip_meta['count'] = 203
naip_meta['dtype'] = "float64"


with rasterio.open('test.tif', 'w', **naip_meta) as dst: #version 2: if any missing data, set to -9999
    dst.write(a_reshaped)
