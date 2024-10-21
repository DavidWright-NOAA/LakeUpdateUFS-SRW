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
    ax.add_feature(cfeature.STATES.with_scale('10m'),facecolor='none', edgecolor='black', linewidth=1.0)
    return ax


### Global constants
filterValue = 0.0
trans = ccrs.PlateCarree() #map transformation/projection
fontsize = 8 #Font Size
xpos = 0.03 #Panel Label location
ypos = 0.10

### Zoom Region for Map
lonMin_Z = -90.8 #Lake Huron: -85.1 #Lake Erie:-83.8 #Lake Michigan: -88.0# # #-92.5#-94.0#-82.083####
lonMax_Z = -83.5 #Lake Huron: -79.1 #Lake Erie:-75.0 #Lake Michigan: -84.0# # #-84.0# -74.0#-77.695####
latMin_Z = 44.9 #Lake Huron: 42.8 #Lake Erie: 41.0 #Lake Michigan: 41.5# # #45.2#40.0#41.012####
latMax_Z = 48.0 #Lake Huron: 46.5 #Lake Erie: 44.8 #Lake Michigan: 44.5# # #49.3#50.0#43.469####
 
img_extent = [lonMin_Z, lonMax_Z, latMin_Z, latMax_Z] 

### Color map and range
cmap = plt.cm.jet
cLevels = np.arange(0.0,42.0,2)
norm = BoundaryNorm(cLevels, ncolors=cmap.N, clip=True)

cmap_diff = plt.cm.bwr
cLevels_diff = [-8.,-6.,-4.,-2.,-1.,1.,2.,4.,6.,8.] #np.arange(-2.0,2.2,.2)
norm_diff = BoundaryNorm(cLevels_diff, ncolors=cmap_diff.N, clip=True)

### File Paths and File Loads
base_fp = ''

swe_ctl_Ids = pygrib.open(base_fp + 'UFS-SRW/2021020500_CTRL/rrfs.t00z.prslevf006.tm00.grib2')
swe_ctl_Fds = pygrib.open(base_fp + 'UFS-SRW/2021020500_CTRL/rrfs.t00z.prslevf054.tm00.grib2')

swe_sta_Ids = pygrib.open(base_fp + 'UFS-SRW/2021020500_STAT/rrfs.t00z.prslevf006.tm00.grib2')
swe_sta_Fds = pygrib.open(base_fp + 'UFS-SRW/2021020500_STAT/rrfs.t00z.prslevf054.tm00.grib2')

swe_dyn_Ids = pygrib.open(base_fp + 'UFS-SRW/2021020500_DYN/rrfs.t00z.prslevf006.tm00.grib2')
swe_dyn_Fds = pygrib.open(base_fp + 'UFS-SRW/2021020500_DYN/rrfs.t00z.prslevf054.tm00.grib2')

snodas_ds = nc.Dataset(base_fp + 'SNODAS/SNODAS_2021_0207-0205.nc')

###Get variables and take difference

swe_ctlI = swe_ctl_Ids.select(name='Water equivalent of accumulated snow depth (deprecated)')[0].values
swe_ctlF = swe_ctl_Fds.select(name='Water equivalent of accumulated snow depth (deprecated)')[0].values

swe_staI = swe_sta_Ids.select(name='Water equivalent of accumulated snow depth (deprecated)')[0].values
swe_staF = swe_sta_Fds.select(name='Water equivalent of accumulated snow depth (deprecated)')[0].values

swe_dynI = swe_dyn_Ids.select(name='Water equivalent of accumulated snow depth (deprecated)')[0].values
swe_dynF = swe_dyn_Fds.select(name='Water equivalent of accumulated snow depth (deprecated)')[0].values

ctl_acc  = swe_ctlF - swe_ctlI
ctl_acc  = np.ma.array(ctl_acc, mask=ctl_acc<=filterValue)

sta_acc  = swe_staF - swe_staI
sta_acc  = np.ma.array(sta_acc, mask=sta_acc<=filterValue)

dyn_acc  = swe_dynF - swe_dynI
dyn_acc  = np.ma.array(dyn_acc, mask=dyn_acc<=filterValue)

sta_ctl_diff = sta_acc - ctl_acc
dyn_sta_diff = dyn_acc - sta_acc

srw_lat,srw_lon = swe_ctl_Ids[1].latlons()


snodas_acc = snodas_ds['Band1'][:]
snodas_lat = snodas_ds['lat']
snodas_lon = snodas_ds['lon']
snodas_acc = np.ma.array(snodas_acc, mask=snodas_acc<=filterValue)

### Plotting
fig, axarr = plt.subplots(nrows=2,ncols=3, figsize=(8.0,3.5), subplot_kw={'projection': trans})
axlist = axarr.flatten()

#Panel 0
outplt = axlist[0].pcolormesh(srw_lon, srw_lat, ctl_acc, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[0])
axlist[0].set_extent(img_extent, trans)
axlist[0].set_title('Control', fontdict={'fontsize':fontsize,'fontweight':'bold'})
text = axlist[0].text(xpos,ypos, "a", transform = axlist[0].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
text.set_bbox(dict(facecolor='lightgrey',alpha=.75, linewidth=0))
#Panel 1
outplt = axlist[1].pcolormesh(srw_lon, srw_lat, sta_acc, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[1])
axlist[1].set_extent(img_extent, trans)
axlist[1].set_title('Static', fontdict={'fontsize':fontsize,'fontweight':'bold'})
text = axlist[1].text(xpos,ypos, "b", transform = axlist[1].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
text.set_bbox(dict(facecolor='lightgrey',alpha=.75, linewidth=0))
#Panel 2
outplt = axlist[2].pcolormesh(srw_lon, srw_lat, dyn_acc, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[2])
axlist[2].set_extent(img_extent, trans)
axlist[2].set_title('Dynamic', fontdict={'fontsize':fontsize,'fontweight':'bold'})
text = axlist[2].text(xpos,ypos, "c", transform = axlist[2].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
text.set_bbox(dict(facecolor='lightgrey',alpha=.75, linewidth=0))
#Panel 3
outplt = axlist[3].pcolormesh(snodas_lon, snodas_lat, snodas_acc, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[3])
axlist[3].set_extent(img_extent, trans)
axlist[3].set_title('SNODAS', fontdict={'fontsize':fontsize,'fontweight':'bold'})
text = axlist[3].text(xpos,ypos, "d", transform = axlist[3].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
text.set_bbox(dict(facecolor='lightgrey',alpha=.75, linewidth=0))
#Panel 4
outplt2 = axlist[4].pcolormesh(srw_lon, srw_lat, sta_ctl_diff, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[4])
axlist[4].set_extent(img_extent, trans)
axlist[4].set_title('Static - Control', fontdict={'fontsize':fontsize,'fontweight':'bold'})
text = axlist[4].text(xpos,ypos, "e", transform = axlist[4].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
text.set_bbox(dict(facecolor='lightgrey',alpha=.75, linewidth=0))
#Panel 5
outplt2 = axlist[5].pcolormesh(srw_lon, srw_lat, dyn_sta_diff, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[5])
axlist[5].set_extent(img_extent, trans)
axlist[5].set_title('Dynamic - Static', fontdict={'fontsize':fontsize,'fontweight':'bold'})
text = axlist[5].text(xpos,ypos, "f", transform = axlist[5].transAxes, weight='bold',size=fontsize, backgroundcolor="white")
text.set_bbox(dict(facecolor='lightgrey',alpha=.75, linewidth=0))

cax1 = axlist[2].inset_axes([1.05,0.0,.05,1.0]) #[X-location, Y-Location, Width, Height]  https://matplotlib.org/stable/gallery/subplots_axes_and_figures/colorbar_placement.html
cax2 = axlist[5].inset_axes([1.05,0.0,.05,1.0])

cb1 = fig.colorbar(outplt, ax=axlist[2], cax=cax1, orientation='vertical')#, ticks=cLevels[::2]) #, location = 'bottom')#, shrink=.8)
cb1.set_label(label=r'${mm}$', size=8)
cb1.ax.tick_params(labelsize=8)

cb2 = fig.colorbar(outplt2, ax=axlist[5], cax=cax2, orientation='vertical', ticks=[-8.,-6.,-4.,-2,-1.,1.,2.,4.,6.,8.])#, ticks=cLevels[::2]) #, location = 'bottom')#, shrink=.8)
cb2.set_label(label=r'${mm}$', size=8)
cb2.ax.tick_params(labelsize=8)


plt.subplots_adjust(hspace = -0.1, wspace=0.0) #left=0.05,bottom=.1,top=.95,right=0.99,
#plt.show()
plt.savefig('Publication/Figures/SWE_Feb2021_6panels_ZOOM_LakeSuperior.png', dpi=1800, bbox_inches='tight')

###Close Files
swe_ctl_Ids.close() 
swe_ctl_Fds.close()

swe_sta_Ids.close()
swe_sta_Fds.close()

swe_dyn_Ids.close()
swe_dyn_Fds.close()

snodas_ds.close()