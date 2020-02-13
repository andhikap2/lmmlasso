gdal_translate -of GTiff tx_ens_mean_0.1deg_reg_v19.0eHOM.nc tx_gdal.tif





gdal_translate -of GTiff tn_ens_mean_0.1deg_reg_v19.0eHOM.nc tn_gdal.tif





# K=2000
# for K in {2000..2018}
# do
# 	gdal_translate -of GTiff 'daymet_v3_tmin_'"$K"'_na.nc4' 'daymet_v3_tmin_'"$K"'_na.tif'
# done


gdal_translate -of GTiff daymet_v3_tmin_2018_na.nc4 daymet_v3_tmin_2018_na.tif


for K in {1979..2019}
do
	gdal_translate -of GTiff 'tmin.'"$K"'.nc' 'tmin'"$K"'.tif'
done