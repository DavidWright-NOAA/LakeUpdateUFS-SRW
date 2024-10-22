#!/bin/tcsh

# FVCOM output files
# if a lake forecast did not complete and there is no result for a lake,
# set 'None' or whatever filename that does not exist. NetCDF files
# are concatenated hourly output files. NCO used to concatenate.
#
set output_erie="2022111518z-112017z_leofs.nc"    #"202200512-00806_eri.nc"
set output_mhs="2022111518z-112017z_lmhofs.nc"    #"202200512-00806_mih.nc"
set output_ont="2022111518z-112017z_loofs.nc"     #"202200512-00806_ont.nc"
set output_sup="2022111518z-112017z_lsofs.nc"     #"202200512-00806_sup.nc"

# names of input files to the Python script
set infile_erie="erie_siglay0_tmp.nc"
set infile_mhs="lmhofs_siglay0_tmp.nc"
set infile_sup="sup_siglay0_tmp.nc"
set infile_ont="ont_siglay0_tmp.nc"

# subset the first layer,  hourly
echo ''
#echo 'subsetting surface layer only, one time slide from FVCOM output files....'
echo 'subsetting surface layer only, every hour from FVCOM output files....'
if (-e $output_erie) then
ncks -O -d time,,,1 -d siglay,0 -d siglev,0 $output_erie $infile_erie
set infile_time=$infile_erie
else
echo ''
echo 'No file found for Erie. Remap without Erie results...'
endif
if (-e $output_mhs) then
ncks -O -d time,,,1 -d siglay,0 -d siglev,0 $output_mhs $infile_mhs
set infile_time=$infile_mhs
else 
echo ''
echo 'No file found for Michigan-Huron. Remap without Mich-Huron results...'
endif
if (-e $output_sup) then
ncks -O -d time,,,1 -d siglay,0 -d siglev,0 $output_sup $infile_sup
set infile_time=$infile_sup
else
echo ''
echo 'No file found for Superior. Remap without Superior results...'
endif
if (-e $output_ont) then
ncks -O -d time,,,1 -d siglay,0 -d siglev,0 $output_ont $infile_ont
set infile_time=$infile_ont
else
echo ''
echo 'No file found for Ontario. Remap without Ontario results...'
endif

if (! -e $output_erie && ! -e  $output_mhs && ! -e  $output_sup && ! -e $output_ont ) then
echo ''
echo "None of output files is available. No remapping is done. "
echo ''
echo "Script ends here."
echo ''
exit
endif
 

## extract HRRR GL mask variables
echo 'Copying a file for GL mask on RRFS grid ...'
cp fv3mask_C3338.nc ./tmp.nc

# get time info
echo ''
echo 'getting time information from a FVCOM output file ...'
ncks -O -h -v time,Times $infile_time  time.nc
ncrename -O -h -d time,Time time.nc # rename 
ncks -h -A time.nc tmp.nc # append
rm time.nc

## add a variable container for tsfc
echo ''
echo 'adding variables (blank for now) to output file ...'
#ncap2  -O -h -s "twsfc=array(0.0f,0.0f,/$Time,$lat,$lon/)" tmp.nc out_fv3grid.nc
#ncap2  -O -h -s "twsfc=glmask" tmp.nc out_fv3grid.nc
ncap2 -O -h -s 'twsfc[$Time,$lat,$lon]=glmask' tmp.nc out_fv3grid.nc
ncatted -O -h -a long_name,twsfc,o,c,water_surface_temperature out_fv3grid.nc
ncatted -O -h -a units,twsfc,o,c,degC out_fv3grid.nc
ncatted -O -h -a description,twsfc,o,c,"water surface temperature" out_fv3grid.nc
#ncap2  -O -h -s "tisfc=array(0.0f,0.0f,/$Time,$lat,$lon/)" out_fv3grid.nc out_fv3grid.nc
#ncap2  -O -h -s "tisfc=glmask" out_fv3grid.nc out_fv3grid.nc
ncap2 -O -h -s 'tisfc[$Time,$lat,$lon]=glmask' out_fv3grid.nc out_fv3grid.nc
ncatted -O -h -a long_name,twsfc,o,c,water_surface_temperature out_fv3grid.nc
ncatted -O -h -a long_name,tisfc,o,c,ice_surface_temperature out_fv3grid.nc
ncatted -O -h -a units,tisfc,o,c,degC out_fv3grid.nc
ncatted -O -h -a description,tisfc,o,c,"ice surface temperature" out_fv3grid.nc
ncap2 -O -h -s 'aice[$Time,$lat,$lon]=glmask' out_fv3grid.nc out_fv3grid.nc
ncatted -O -h -a long_name,aice,o,c,ice_concentration out_fv3grid.nc
ncatted -O -h -a units,aice,o,c,- out_fv3grid.nc
ncatted -O -h -a description,aice,o,c,"ice fraction [0-1]" out_fv3grid.nc
ncap2 -O -h -s 'vice[$Time,$lat,$lon]=glmask' out_fv3grid.nc out_fv3grid.nc
ncatted -O -h -a long_name,vice,o,c,ice_thickness out_fv3grid.nc
ncatted -O -h -a units,vice,o,c,m out_fv3grid.nc
ncatted -O -h -a description,vice,o,c,"mean ice thickness [m]" out_fv3grid.nc
ncap2 -O -h -s 'tsfc[$Time,$lat,$lon]=glmask' out_fv3grid.nc out_fv3grid.nc
ncatted -O -h -a long_name,tsfc,o,c,lake_skin_temperature out_fv3grid.nc
ncatted -O -h -a units,tsfc,o,c,degC out_fv3grid.nc
ncatted -O -h -a description,tsfc,o,c,"skin temperature of lake and ice surfaces (weighted-mean based on ice concentration)" out_fv3grid.nc
# sensible heat flux
ncap2 -O -h -s 'hs[$Time,$lat,$lon]=glmask' out_fv3grid.nc out_fv3grid.nc
ncatted -O -h -a long_name,hs,o,c,sensible_heat_flux out_fv3grid.nc
ncatted -O -h -a units,hs,o,c,W/m2 out_fv3grid.nc
ncatted -O -h -a description,hs,o,c,"sensible heat flux at the water surface (not the ice surface), downward positive" out_fv3grid.nc
#latent heat flux
ncap2 -O -h -s 'hl[$Time,$lat,$lon]=glmask' out_fv3grid.nc out_fv3grid.nc
ncatted -O -h -a long_name,hl,o,c,latent_heat_flux out_fv3grid.nc
ncatted -O -h -a units,hl,o,c,W/m2 out_fv3grid.nc
ncatted -O -h -a description,hl,o,c,"latent heat flux at the water surface (not the ice surface), downward positive" out_fv3grid.nc


# edit mask attributes
ncatted -O -h -a description,glmask,o,c,"Great Lakes mask (1 if overwater 0 otherwise" out_fv3grid.nc

# run interpolation script
echo ''
echo 'now running the Python script for remapping ...'
python ./interp_skintemp_5lakes_ice_fv3grid_ver6.py $infile_erie $infile_mhs $infile_sup $infile_ont 

# remove time dimension from mask, lon&lat info.
#ncwa -O -h -v LAKEMASK -a Time  out_fv3grid.nc maskinfo.nc
# extract skin temp info
ncks -O -h -v geolon,geolat,glmask,twsfc,tisfc,aice,vice,tsfc,hs,hl,time,Times  out_fv3grid.nc tsfc_C3338_grid.nc
ncap2 -O -s 'geolon=double(geolon);geolat=double(geolat);time=double(time)' tsfc_C3338_grid.nc tsfc_C3338_grid.nc

echo 'skin temp info extracted.'

# rename Time dimension (for Ferret use)
#ncrename -O -h -d Time,time tsfc_rrfs_grid.nc
#echo 'renamed'

# remove temporary files
echo ''
echo 'removing temporary files .....'
rm -f $infile_erie $infile_mhs $infile_sup $infile_ont 
rm -f tmp.nc out_fv3grid.nc
echo ''
echo 'done.'



