# -*- coding: utf-8 -*-

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pygrib
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import netCDF4 as nc
import pygrib

### Draw maps using Cartopy Features
def plotMap(ax):
    ax.coastlines(resolution='10m', linewidth=1.0)
    #ax.add_feature(cfeature.LAKES.with_scale('10m'),facecolor='none', edgecolor='black', linewidth=0.75)
    ax.add_feature(cfeature.STATES.with_scale('10m'),facecolor='white', edgecolor='black', linewidth=1.0)
    return ax


### Global constants
trans = ccrs.PlateCarree() #map transformation/projection
fontsize = 8 #Font Size
xpos = 0.02 #Panel Label location
ypos = 0.07

### Zoom Region for Map
lonMin_Z = -92.5 #Great Lakes: -92.5 #Lake Huron: -85.1 #Lake Erie:-83.8 #Lake Michigan: -88.0# # #-92.5#-94.0#-82.083####
lonMax_Z = -75.75 #Great Lakes: -75.75 #Lake Huron: -79.1 #Lake Erie:-75.0 #Lake Michigan: -84.0# # #-84.0# -74.0#-77.695####
latMin_Z = 41.0 #Great Lakes: 41.0 #Lake Huron: 42.8 #Lake Erie: 41.0 #Lake Michigan: 41.5# # #45.2#40.0#41.012####
latMax_Z = 49.5 #Great Lakes: 49.5 #Lake Huron: 46.5 #Lake Erie: 44.8 #Lake Michigan: 44.5# # #49.3#50.0#43.469####
 
img_extent = [lonMin_Z, lonMax_Z, latMin_Z, latMax_Z] 

### Color map and range
cmap = plt.cm.jet
cLevels = [-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]#[0,.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0] #np.arange(0,13,.2)
norm = BoundaryNorm(cLevels, ncolors=cmap.N, clip=True)

cmap_diff = plt.cm.bwr
cLevels_diff = [-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5,3.0] #np.arange(-2.0,2.2,.2)
norm_diff = BoundaryNorm(cLevels_diff, ncolors=cmap_diff.N, clip=True)

### File Paths and File Loads
base_fp = ''
dec17_gs_Ids = nc.Dataset(base_fp + 'Dec2017CaseStudy/GLSEA3/2017_358_glsea_sst.nc')
dec17_gs_Fds = nc.Dataset(base_fp + 'Dec2017CaseStudy/GLSEA3/2017_360_glsea_sst.nc')
feb21_gs_Ids = nc.Dataset(base_fp + 'Feb2021CaseStudy/GLSEA3/2021_036_glsea_sst.nc')
feb21_gs_Fds = nc.Dataset(base_fp + 'Feb2021CaseStudy/GLSEA3/2021_038_glsea_sst.nc')
jan22_gs_Ids = nc.Dataset(base_fp + 'Jan2022CaseStudy/GLSEA3/2022_005_glsea_sst.nc')
jan22_gs_Fds = nc.Dataset(base_fp + 'Jan2022CaseStudy/GLSEA3/2022_008_glsea_sst.nc')
nov22_gs_Ids = nc.Dataset(base_fp + 'Nov2022CaseStudy/GLSEA3/2022_320_glsea_sst.nc')
nov22_gs_Fds = nc.Dataset(base_fp + 'Nov2022CaseStudy/GLSEA3/2022_323_glsea_sst.nc')

#Dec 2017 Simulation
dec17_srw_Cds = nc.Dataset(base_fp + 'Dec2017CaseStudy/UFS-SRW/2017122412_CTRL/phyf000.nc')
dec17_srw_Ids = nc.Dataset(base_fp + 'Dec2017CaseStudy/UFS-SRW/2017122412_DYN/phyf000.nc')
dec17_srw_Fds = nc.Dataset(base_fp + 'Dec2017CaseStudy/UFS-SRW/2017122412_DYN/phyf042.nc')

#Feb 2021 Simulation
feb21_srw_Cds = pygrib.open(base_fp + 'Feb2021CaseStudy/UFS-SRW/2021020500_CTRL/rrfs.t00z.prslevf000.tm00.grib2')
feb21_srw_Ids = pygrib.open(base_fp + 'Feb2021CaseStudy/UFS-SRW/2021020500_DYN/rrfs.t00z.prslevf000.tm00.grib2')
feb21_srw_Fds = pygrib.open(base_fp + 'Feb2021CaseStudy/UFS-SRW/2021020500_DYN/rrfs.t00z.prslevf054.tm00.grib2')

#Jan 2022 Simulation
jan22_srw_Cds = nc.Dataset(base_fp + 'Jan2022CaseStudy/UFS-SRW/22010512_CTRL/phyf000.nc')
jan22_srw_Ids = nc.Dataset(base_fp + 'Jan2022CaseStudy/UFS-SRW/22010512_DYN/phyf000.nc')
jan22_srw_Fds = nc.Dataset(base_fp + 'Jan2022CaseStudy/UFS-SRW/22010512_DYN/phyf066.nc')

#Nov 2022 Simulation
nov22_srw_Cds = nc.Dataset(base_fp + 'Nov2022CaseStudy/UFS-SRW/22111612_CTRL/phyf000.nc')
nov22_srw_Ids = nc.Dataset(base_fp + 'Nov2022CaseStudy/UFS-SRW/22111612_DYN/phyf000.nc')
nov22_srw_Fds = nc.Dataset(base_fp + 'Nov2022CaseStudy/UFS-SRW/22111612_DYN/phyf066.nc')


#fp = open('grib2Dump.txt','w')
#for n in feb21_srw_Cds:
#    fp.write(str(n)+'\n')
#fp.close()    
    

###Get Variables
#Dec 2017
dec17_gs_I = dec17_gs_Ids['sst'][0] #(lat, lon)
dec17_gs_F = dec17_gs_Fds['sst'][0]
dec17_srw_C= dec17_srw_Cds['tmpsfc'][0] - 273.15
dec17_srw_I= dec17_srw_Ids['tmpsfc'][0] - 273.15
dec17_srw_F= dec17_srw_Fds['tmpsfc'][0] - 273.15

dec17_srw_ice_C= dec17_srw_Cds['icec'][0]
dec17_srw_ice_I= dec17_srw_Ids['icec'][0]
dec17_srw_ice_F= dec17_srw_Fds['icec'][0]

#Feb 2021
feb21_gs_I = feb21_gs_Ids['sst'][0]
feb21_gs_F = feb21_gs_Fds['sst'][0]
feb21_srw_C= feb21_srw_Cds.select(name='Sea surface temperature')[0].values - 273.15
feb21_srw_I= feb21_srw_Ids.select(name='Sea surface temperature')[0].values - 273.15
feb21_srw_F= feb21_srw_Fds.select(name='Sea surface temperature')[0].values - 273.15

feb21_srw_ice_C = feb21_srw_Cds.select(name='Sea ice area fraction')[0].values
feb21_srw_ice_I = feb21_srw_Ids.select(name='Sea ice area fraction')[0].values
feb21_srw_ice_F = feb21_srw_Fds.select(name='Sea ice area fraction')[0].values

#Jan 2022
jan22_gs_I = jan22_gs_Ids['sst'][0]
jan22_gs_F = jan22_gs_Fds['sst'][0]
jan22_srw_C= jan22_srw_Cds['tmpsfc'][0] - 273.15
jan22_srw_I= jan22_srw_Ids['tmpsfc'][0] - 273.15
jan22_srw_F= jan22_srw_Fds['tmpsfc'][0] - 273.15

jan22_srw_ice_C = jan22_srw_Cds['icec'][0]
jan22_srw_ice_I = jan22_srw_Ids['icec'][0]
jan22_srw_ice_F = jan22_srw_Fds['icec'][0]

#Nov 2022
nov22_gs_I = nov22_gs_Ids['sst'][0]
nov22_gs_F = nov22_gs_Fds['sst'][0]
nov22_srw_C= nov22_srw_Cds['tmpsfc'][0] - 273.15
nov22_srw_I= nov22_srw_Ids['tmpsfc'][0] - 273.15
nov22_srw_F= nov22_srw_Fds['tmpsfc'][0] - 273.15



#Lat/Lon variables
gs_lat = dec17_gs_Ids['lat'][:]
gs_lon = dec17_gs_Ids['lon'][:]

srw_lat= nov22_srw_Ids['lat'][:,:]
srw_lon= nov22_srw_Ids['lon'][:,:]

grb_lat,grb_lon = feb21_srw_Cds[1].latlons()

#Change over time
dec17_gs_chg = dec17_gs_F - dec17_gs_I
dec17_srw_chg= dec17_srw_F- dec17_srw_I
dec17_ice_chg = dec17_srw_ice_F - dec17_srw_ice_I

feb21_gs_chg = feb21_gs_F - feb21_gs_I
feb21_srw_chg= feb21_srw_F - feb21_srw_I

jan22_gs_chg = jan22_gs_F - jan22_gs_I
jan22_srw_chg= jan22_srw_F- jan22_srw_I

nov22_gs_chg = nov22_gs_F - nov22_gs_I
nov22_srw_chg= nov22_srw_F- nov22_srw_I

### Plotting
plt.rcParams['hatch.linewidth'] = .5 #Changes hatching width

fig, axarr = plt.subplots(nrows=4,ncols=5, figsize=(8.0,4.0), subplot_kw={'projection': trans})

axlist = axarr.flatten()

## Row 1 (Dec 2017)
outplt = axlist[0].pcolormesh(srw_lon, srw_lat, dec17_srw_C, cmap=cmap, norm=norm, transform=trans) #, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[0].contourf(srw_lon, srw_lat, dec17_srw_ice_C, levels=[0.01,1.0], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[0])
axlist[0].set_extent(img_extent, trans)
axlist[0].set_title('Control', fontdict={'fontsize':fontsize,'fontweight':'bold'})
axlist[0].text(-.1,.5, "Dec. 2017", horizontalalignment='center',verticalalignment='center',transform = axlist[0].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[0].text(xpos,ypos, "a", transform = axlist[0].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[1].pcolormesh(srw_lon, srw_lat, dec17_srw_I, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[1].contourf(srw_lon, srw_lat, dec17_srw_ice_I, levels=[0.01,1.0], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[1])
axlist[1].set_extent(img_extent, trans)
axlist[1].set_title('Static', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[1].text(xpos,ypos, "b", transform = axlist[1].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[2].pcolormesh(gs_lon, gs_lat, dec17_gs_I, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[2].contourf(gs_lon, gs_lat, dec17_gs_I, levels=[0.0,0.3], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[2])
axlist[2].set_extent(img_extent, trans)
axlist[2].set_title('GLSEA3 Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[2].text(xpos,ypos, "c", transform = axlist[2].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[3].pcolormesh(srw_lon, srw_lat, dec17_srw_chg, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[3].contourf(srw_lon,srw_lat,dec17_srw_ice_F, levels=[0.01,1.0], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[3])
axlist[3].set_extent(img_extent, trans)
axlist[3].set_title('Dynamic Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[3].text(xpos,ypos, "d", transform = axlist[3].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[4].pcolormesh(gs_lon, gs_lat, dec17_gs_chg, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[4].contourf(gs_lon, gs_lat, dec17_gs_F, levels=[0.0,0.3], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[4])
axlist[4].set_extent(img_extent, trans)
axlist[4].set_title('GLSEA3 Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[4].text(xpos,ypos, "e", transform = axlist[4].transAxes, weight='bold',size=fontsize, backgroundcolor="none")


## Row 2 (Feb 2021)
outplt = axlist[5].pcolormesh(grb_lon, grb_lat, feb21_srw_C, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[5].contourf(grb_lon, grb_lat, feb21_srw_ice_C, levels=[0.01,1.0,], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[5])
axlist[5].set_extent(img_extent, trans)
axlist[5].text(-.1,.5, "Feb. 2021", horizontalalignment='center',verticalalignment='center',transform = axlist[5].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')
#axlist[7].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[5].text(xpos,ypos, "f", transform = axlist[5].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[6].pcolormesh(grb_lon, grb_lat, feb21_srw_I, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[6].contourf(grb_lon, grb_lat, feb21_srw_ice_I, levels=[0.01,1.0,], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
#for i, collection in enumerate(outplt.collections):
#    collection.set_edgecolor('pink')
plotMap(axlist[6])
axlist[6].set_extent(img_extent, trans)
#axlist[7].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[6].text(xpos,ypos, "g", transform = axlist[6].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[7].pcolormesh(gs_lon, gs_lat, feb21_gs_I, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[7].contourf(gs_lon, gs_lat, feb21_gs_I, levels=[0.0, 0.3], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[7])
axlist[7].set_extent(img_extent, trans)
#axlist[7].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[7].text(xpos,ypos, "h", transform = axlist[7].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[8].pcolormesh(grb_lon, grb_lat, feb21_srw_chg, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[8].contourf(grb_lon, grb_lat, feb21_srw_ice_F, levels=[0.01,1.0], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
#for i, collection in enumerate(outplt.collections):
#    collection.set_edgecolor('white')
plotMap(axlist[8])
axlist[8].set_extent(img_extent, trans)
#axlist[7].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[8].text(xpos,ypos, "i", transform = axlist[8].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[9].pcolormesh(gs_lon, gs_lat, feb21_gs_chg, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[9].contourf(gs_lon, gs_lat, feb21_gs_F, levels=[0.0, 0.3], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[9])
axlist[9].set_extent(img_extent, trans)
#axlist[9].set_title('GLSEA Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[9].text(xpos,ypos, "j", transform = axlist[9].transAxes, weight='bold',size=fontsize, backgroundcolor="none")


## Row 3 (Jan 2022)
outplt = axlist[10].pcolormesh(srw_lon, srw_lat, jan22_srw_C, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[10].contourf(srw_lon, srw_lat, jan22_srw_ice_C, levels=[0.01,1.0,], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[10])
axlist[10].set_extent(img_extent, trans)
axlist[10].text(-.1,.5, "Jan. 2022", horizontalalignment='center',verticalalignment='center',transform = axlist[10].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[10].text(xpos,ypos, "k", transform = axlist[10].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[11].pcolormesh(srw_lon, srw_lat, jan22_srw_I, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[11].contourf(srw_lon, srw_lat, jan22_srw_ice_I, levels=[0.01,1.0,], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[11])
axlist[11].set_extent(img_extent, trans)
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[11].text(xpos,ypos, "l", transform = axlist[11].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[12].pcolormesh(gs_lon, gs_lat, jan22_gs_I, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[12].contourf(gs_lon, gs_lat, jan22_gs_I, levels=[0.0,0.3], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[12])
axlist[12].set_extent(img_extent, trans)
#axlist[12].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[12].text(xpos,ypos, "m", transform = axlist[12].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[13].pcolormesh(srw_lon, srw_lat, jan22_srw_chg, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[13].contourf(srw_lon, srw_lat, jan22_srw_ice_F, levels=[0.01,1.0,], hatches=['xxxxxx','xxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[13])
axlist[13].set_extent(img_extent, trans)
#axlist[14].set_title('GLSEA Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[13].text(xpos,ypos, "n", transform = axlist[13].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[14].pcolormesh(gs_lon, gs_lat, jan22_gs_chg, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt = axlist[14].contourf(gs_lon, gs_lat, jan22_gs_F, levels=[0.0,0.3], hatches=['xxxxxx','xxxxxxxx'], colors='lightgrey', transform=trans)
plotMap(axlist[14])
axlist[14].set_extent(img_extent, trans)
#axlist[14].set_title('GLSEA Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[14].text(xpos,ypos, "o", transform = axlist[14].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

## Row 4 (Nov 2022)
outplt = axlist[15].pcolormesh(srw_lon, srw_lat, nov22_srw_C, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[15])
axlist[15].set_extent(img_extent, trans)
axlist[15].text(-.1,.5, "Nov. 2022", horizontalalignment='center',verticalalignment='center',transform = axlist[15].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')
#axlist[17].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[15].text(xpos,ypos, "p", transform = axlist[15].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[16].pcolormesh(srw_lon, srw_lat, nov22_srw_I, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[16])
axlist[16].set_extent(img_extent, trans)
#axlist[17].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[16].text(xpos,ypos, "q", transform = axlist[16].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt2 = axlist[17].pcolormesh(gs_lon, gs_lat, nov22_gs_I, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[17])
axlist[17].set_extent(img_extent, trans)
#axlist[17].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[17].text(xpos,ypos, "r", transform = axlist[17].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[18].pcolormesh(srw_lon, srw_lat, nov22_srw_chg, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[18])
axlist[18].set_extent(img_extent, trans)
#axlist[17].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[18].text(xpos,ypos, "s", transform = axlist[18].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[19].pcolormesh(gs_lon, gs_lat, nov22_gs_chg, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[19])
axlist[19].set_extent(img_extent, trans)
#axlist[19].set_title('GLSEA Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[19].text(xpos,ypos, "t", transform = axlist[19].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

#Colorbars for figures

#ax1_Test = make_axes_locatable(axlist[1])
#cax1 = ax1_Test.append_axes("bottom",size="7%", pad="2%")

cax1 = axlist[15].inset_axes([0.15,-.3,2.75,.2]) #[X-location, Y-Location, Width, Height]  https://matplotlib.org/stable/gallery/subplots_axes_and_figures/colorbar_placement.html
cax2 = axlist[18].inset_axes([0.15,-.3,1.75,.2])

cb1 = fig.colorbar(outplt2, ax=axlist[15:18], cax=cax1, orientation='horizontal', ticks=cLevels[1::2]) #, location = 'bottom')#, shrink=.8)
cb1.set_label(label=r'${Deg. C}$', size=6)
cb1.ax.tick_params(labelsize=6)
cb2 = fig.colorbar(outplt, ax=axlist[18:20], cax=cax2, orientation='horizontal', ticks=[-3.0,-2.0,-1.0,-.5,.5,1.0,2.0,3.0])#, location = 'bottom', shrink=.8)
cb2.set_label(label=r'${Deg. C}$', size=6)
cb2.ax.tick_params(labelsize=6)
### Save Figure
plt.subplots_adjust(hspace = 0.0, wspace=0.0) #left=0.05,bottom=.1,top=.95,right=0.99,   Lake Erie: (left=0.013,right=0.987,top=0.826,bottom=0.219,hspace=0.053,wspace=0.030)#Lake Michigan: (left=0.013,right=0.987,top=0.956,bottom=0.049,hspace=0.098,wspace=0.030)#(left=0.025,right=0.987,top=0.750,bottom=0.100,hspace=0.225,wspace=0.081) #.034 #hspace = 0.05 for Erie/Ontario
#fig.tight_layout()
plt.savefig(base_fp +'Publication/Figures/LST_GLSEA_Comparison_20panels_WithIce.png', dpi=1800)

###Close all datasets
dec17_gs_Ids.close()
dec17_gs_Fds.close()
feb21_gs_Ids.close() 
feb21_gs_Fds.close() 
jan22_gs_Ids.close() 
jan22_gs_Fds.close() 
nov22_gs_Ids.close() 
nov22_gs_Fds.close() 

#Dec 2017 Simulation
dec17_srw_Cds.close() 
dec17_srw_Ids.close()
dec17_srw_Fds.close() 

#Feb 2021 Simulation
feb21_srw_Cds.close()
feb21_srw_Ids.close()
feb21_srw_Fds.close()

#Jan 2022 Simulation
jan22_srw_Cds.close()
jan22_srw_Ids.close()
jan22_srw_Fds.close()

#Nov 2022 Simulation
nov22_srw_Cds.close()
nov22_srw_Ids.close()
nov22_srw_Fds.close()