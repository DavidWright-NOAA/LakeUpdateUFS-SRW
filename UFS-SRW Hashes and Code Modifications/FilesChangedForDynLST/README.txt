Places to move files to before compiling:
src/ufs-weather-model/FV3/atmos_cubed_sphere/model/fv_regional_bc.F90
#Routine to read in new lower_bc variables

src/ufs-weather-model/FV3/atmos_cubed_sphere/model/fv_control.F90
#Controls the reading in of new namelist section/variables

src/ufs-weather-model/FV3/atmos_cubed_sphere/model/fv_arrays.F90
#Creates new variables in the Atm object to read new lake conditions into

src/ufs-weather-model/FV3/atmos_cubed_sphere/driver/fvGFS/atmosphere.F90
#Calls update lowbc and moves variable to SfcProp object



#Namelist section:
&lowbc_nml
    lowbc_base = 'fvcom_f'
    lowbc_dir = '/PATH/TO/LAKE/DATA/lake_stage/'
    lowbc_int = 1
    lowbcupdate = .true.
/

#lowbc_base: Base name of the files UFS will look for to update lower boundary conditions.
#    HARD CODED 'fvcom_f' as the file name, so changes to this variable currently won't do anything.
#lowbc_dir: File to look for FVCOM files in.
#    THIS CURRENTLY DOES NOT WORK. Files should be linked into the YYYYMMDDHH/INPUT dir or the UFS-SRW run.
#lowbc_int: Interval (hours) between lowbc update calls
#lowbcupdate: T/F, T if to turn on update to the lower boundary conditions.
