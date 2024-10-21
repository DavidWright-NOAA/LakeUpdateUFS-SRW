# -*- coding: utf-8 -*-

import numpy as np
import pygrib
import netCDF4 as nc
from datetime import datetime as dt
from datetime import timedelta as delta

### Global constants
filterValue = 0.0

fontsize = 8 #Font Size
xpos = 0.08 #Panel Label location
ypos = 0.15
MH_Divide = -84.7540+360. #Dividing longitude taken from Fujisaki-Manome's 2020 paper to divide Michigan and Huron

Round_Value = 2
#Blank Lists to append to
mi_ctl_lh = []
mi_sta_lh = []
mi_dyn_lh = []

hu_ctl_lh = []
hu_sta_lh = []
hu_dyn_lh = []

su_ctl_lh = []
su_sta_lh = []
su_dyn_lh = []

er_ctl_lh = []
er_sta_lh = []
er_dyn_lh = []

on_ctl_lh = []
on_sta_lh = []
on_dyn_lh = []

plot_time = []
### File Paths and File Loads
base_fp = ''
mask_fp = 'CreateMask/GreatLakesRegion_Mask_UFS-SRW-OUTPUTGrid_v1.nc'
ctl_fp = 'Feb2021CaseStudy/UFS-SRW/'
sta_fp = ctl_fp


fp = open('Feb2021_lh_spatialMean.txt','w')

startTime = dt(year=2021,month=2,day=5,hour=0,minute=00)

### Read in Masks from external file and divide MIH into Michigan and Huron
masks = nc.Dataset(mask_fp)
sup_mask = masks['supmask'][:,:]
mih_mask = masks['mhsmask'][:,:]
ont_mask = masks['ontmask'][:,:]
eri_mask = masks['erimask'][:,:]
lat_mask = masks['geolat'][:,:]
lon_mask = masks['geolon'][:,:]

mic_mask = np.where((mih_mask == 1) & (lon_mask >= MH_Divide),0,mih_mask)
hur_mask = np.where((mih_mask == 1) & (lon_mask < MH_Divide),0,mih_mask)

fp.write('Forecast Time   \tMIC-C\tMIC-S\tMIC-D\tHUR-C\tHUR-S\tHUR-D\tSUP-C\tSUP-S\tSUP-D\tERI-C\tERI-S\tERI-D\tONT-C\tONT-S\tONT-D\n')
for fhr in np.arange(0,55,1):
    #plot_time.append(fhr)
    deltaT= delta(hours=int(fhr))
    curTime = startTime + deltaT
    plot_time.append(curTime.strftime("%Y/%m/%d %H:%M"))
    
    ctl_ds = pygrib.open(base_fp + ctl_fp + f'2021020500_CTRL/rrfs.t00z.prslevf0{fhr:02d}.tm00.grib2')
    
    sta_ds = pygrib.open(base_fp + sta_fp + f'2021020500_STAT/rrfs.t00z.prslevf0{fhr:02d}.tm00.grib2')
    
    dyn_ds = pygrib.open(base_fp + sta_fp + f'2021020500_DYN/rrfs.t00z.prslevf0{fhr:02d}.tm00.grib2')
    
    ###Get variables and take difference
    ctl_lh = ctl_ds.select(name='Latent heat net flux', typeOfLevel='surface', forecastTime=fhr)[0].values
    sta_lh = sta_ds.select(name='Latent heat net flux', typeOfLevel='surface', forecastTime=fhr)[0].values
    dyn_lh = dyn_ds.select(name='Latent heat net flux', typeOfLevel='surface', forecastTime=fhr)[0].values

    
    mi_ctl_lh.append(np.ma.mean(np.ma.array(ctl_lh, mask=mic_mask == 0)))
    mi_sta_lh.append(np.ma.mean(np.ma.array(sta_lh, mask=mic_mask == 0)))
    mi_dyn_lh.append(np.ma.mean(np.ma.array(dyn_lh, mask=mic_mask == 0)))

    hu_ctl_lh.append(np.ma.mean(np.ma.array(ctl_lh, mask=hur_mask == 0)))
    hu_sta_lh.append(np.ma.mean(np.ma.array(sta_lh, mask=hur_mask == 0)))
    hu_dyn_lh.append(np.ma.mean(np.ma.array(dyn_lh, mask=hur_mask == 0)))
    
    su_ctl_lh.append(np.ma.mean(np.ma.array(ctl_lh, mask=sup_mask == 0)))
    su_sta_lh.append(np.ma.mean(np.ma.array(sta_lh, mask=sup_mask == 0)))
    su_dyn_lh.append(np.ma.mean(np.ma.array(dyn_lh, mask=sup_mask == 0)))
    
    er_ctl_lh.append(np.ma.mean(np.ma.array(ctl_lh, mask=eri_mask == 0)))
    er_sta_lh.append(np.ma.mean(np.ma.array(sta_lh, mask=eri_mask == 0)))
    er_dyn_lh.append(np.ma.mean(np.ma.array(dyn_lh, mask=eri_mask == 0)))
    
    on_ctl_lh.append(np.ma.mean(np.ma.array(ctl_lh, mask=ont_mask == 0)))
    on_sta_lh.append(np.ma.mean(np.ma.array(sta_lh, mask=ont_mask == 0)))
    on_dyn_lh.append(np.ma.mean(np.ma.array(dyn_lh, mask=ont_mask == 0)))

    fp.write(str(plot_time[-1]) +'\t' + str(round(mi_ctl_lh[-1],Round_Value)) + '\t' + str(round(mi_sta_lh[-1],Round_Value)) + '\t' + str(round(mi_dyn_lh[-1],Round_Value)) + '\t'
             +  str(round(hu_ctl_lh[-1],Round_Value)) + '\t' + str(round(hu_sta_lh[-1],Round_Value)) + '\t' + str(round(hu_dyn_lh[-1],Round_Value)) +'\t'
             +  str(round(su_ctl_lh[-1],Round_Value)) + '\t' + str(round(su_sta_lh[-1],Round_Value)) + '\t' + str(round(su_dyn_lh[-1],Round_Value)) +'\t'
             +  str(round(er_ctl_lh[-1],Round_Value)) + '\t' + str(round(er_sta_lh[-1],Round_Value)) + '\t' + str(round(er_dyn_lh[-1],Round_Value)) +'\t'
             +  str(round(on_ctl_lh[-1],Round_Value)) + '\t' + str(round(on_sta_lh[-1],Round_Value)) + '\t' + str(round(on_dyn_lh[-1],Round_Value)) + '\n')

fp.close()

