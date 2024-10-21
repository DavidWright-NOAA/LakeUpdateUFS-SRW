# -*- coding: utf-8 -*-

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pygrib
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import netCDF4 as nc


### Draw maps using Cartopy Features
def plotMap(ax):
    ax.coastlines(resolution='10m', linewidth=1.0)
    #ax.add_feature(cfeature.LAKES.with_scale('10m'),facecolor='none', edgecolor='black', linewidth=0.75)
    ax.add_feature(cfeature.STATES.with_scale('10m'),facecolor='white', edgecolor='black', linewidth=1.0)
    return ax


### Global constants
trans = ccrs.PlateCarree() #map transformation/projection
fontsize = 8 #Font Size
fontsize_title = 10 #Font Size for Titles
xpos = 0.02 #Panel Label location
ypos = 0.07

### Zoom Region for Map
lonMin_Z = -92.5 #Lake Huron: -85.1 #Lake Erie:-83.8 #Lake Michigan: -88.0# # #-92.5#-94.0#-82.083####
lonMax_Z = -75.75 #Lake Huron: -79.1 #Lake Erie:-75.0 #Lake Michigan: -84.0# # #-84.0# -74.0#-77.695####
latMin_Z = 41.0 #Lake Huron: 42.8 #Lake Erie: 41.0 #Lake Michigan: 41.5# # #45.2#40.0#41.012####
latMax_Z = 49.5 #Lake Huron: 46.5 #Lake Erie: 44.8 #Lake Michigan: 44.5# # #49.3#50.0#43.469####
 
img_extent = [lonMin_Z, lonMax_Z, latMin_Z, latMax_Z] 

### Color map and range
cmap = plt.cm.jet
cLevels = [0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400] #np.arange(0,13,.2)
norm = BoundaryNorm(cLevels, ncolors=cmap.N, clip=True)

cmap_diff = plt.cm.bwr
cLevels_diff = [-100,-75,-50,-25,-10,10,25,50,75,100]#[-50,-40,-30,-20,-10,10,20,30,40,50] #np.arange(-2.0,2.2,.2)
norm_diff = BoundaryNorm(cLevels_diff, ncolors=cmap_diff.N, clip=True)

### File Paths and File Loads
base_fp = ''

#Dec 2017 Simulation
dec17_ctl = nc.Dataset(base_fp + 'Dec2017CaseStudy/UFS-SRW/2017122412_CTRL/phyf038.nc')
dec17_sta = nc.Dataset(base_fp + 'Dec2017CaseStudy/UFS-SRW/2017122412_STAT/phyf038.nc')
dec17_dyn = nc.Dataset(base_fp + 'Dec2017CaseStudy/UFS-SRW/2017122412_DYN/phyf038.nc')

#Feb 2021 Simulation
feb21_ctl = pygrib.open(base_fp + 'Feb2021CaseStudy/UFS-SRW/2021020500_CTRL/rrfs.t00z.prslevf033.tm00.grib2')
feb21_sta = pygrib.open(base_fp + 'Feb2021CaseStudy/UFS-SRW/2021020500_STAT/rrfs.t00z.prslevf033.tm00.grib2')
feb21_dyn = pygrib.open(base_fp + 'Feb2021CaseStudy/UFS-SRW/2021020500_DYN/rrfs.t00z.prslevf033.tm00.grib2')

#Jan 2022 Simulation
jan22_ctl = nc.Dataset(base_fp + 'Jan2022CaseStudy/UFS-SRW/22010512_CTRL/phyf019.nc')
jan22_sta = nc.Dataset(base_fp + 'Jan2022CaseStudy/UFS-SRW/22010512_STAT/phyf019.nc')
jan22_dyn = nc.Dataset(base_fp + 'Jan2022CaseStudy/UFS-SRW/22010512_DYN/phyf019.nc')

#Nov 2022 Simulation
nov22_ctl = nc.Dataset(base_fp + 'Nov2022CaseStudy/UFS-SRW/22111612_CTRL/phyf047.nc')
nov22_sta = nc.Dataset(base_fp + 'Nov2022CaseStudy/UFS-SRW/22111612_STATIC/phyf047.nc')
nov22_dyn = nc.Dataset(base_fp + 'Nov2022CaseStudy/UFS-SRW/22111612_DYN/phyf047.nc')
  

###Get Variables
#Dec 2017
dec17_ctl_lh= dec17_ctl['lhtfl'][0]
dec17_sta_lh= dec17_sta['lhtfl'][0]
dec17_dyn_lh= dec17_dyn['lhtfl'][0]

#Feb 2021
feb21_ctl_lh= feb21_ctl.select(name='Latent heat net flux', typeOfLevel='surface', forecastTime=33)[0].values
feb21_sta_lh= feb21_sta.select(name='Latent heat net flux', typeOfLevel='surface', forecastTime=33)[0].values
feb21_dyn_lh= feb21_dyn.select(name='Latent heat net flux', typeOfLevel='surface', forecastTime=33)[0].values

#Jan 2022
jan22_ctl_lh= jan22_ctl['lhtfl'][0]
jan22_sta_lh= jan22_sta['lhtfl'][0]
jan22_dyn_lh= jan22_dyn['lhtfl'][0]

#Nov 2022
nov22_ctl_lh= nov22_ctl['lhtfl'][0]
nov22_sta_lh= nov22_sta['lhtfl'][0]
nov22_dyn_lh= nov22_dyn['lhtfl'][0]


#Lat/Lon variables
srw_lat= nov22_ctl['lat']
srw_lon= nov22_ctl['lon']

grb_lat,grb_lon = feb21_ctl[1].latlons()

#Change over time
dec17_sta_ctl = dec17_sta_lh - dec17_ctl_lh
dec17_dyn_sta = dec17_dyn_lh - dec17_sta_lh

feb21_sta_ctl = feb21_sta_lh - feb21_ctl_lh
feb21_dyn_sta = feb21_dyn_lh - feb21_sta_lh

jan22_sta_ctl = jan22_sta_lh - jan22_ctl_lh
jan22_dyn_sta = jan22_dyn_lh - jan22_sta_lh

nov22_sta_ctl = nov22_sta_lh - nov22_ctl_lh
nov22_dyn_sta = nov22_dyn_lh - nov22_sta_lh

### Plotting
fig, axarr = plt.subplots(nrows=5,ncols=4, figsize=(8.5,6.0), subplot_kw={'projection': trans})

axlist = axarr.flatten()

## Row 1 (Control)
outplt2 = axlist[0].pcolormesh(srw_lon, srw_lat, dec17_ctl_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[0])
axlist[0].set_extent(img_extent, trans)
axlist[0].set_title("Dec. 2017", fontdict={'fontsize':fontsize_title,'fontweight':'bold'},)
#axlist[0].set_title('Control', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#axlist[0].text(-.1,.5, "Dec. 2017", horizontalalignment='center',verticalalignment='center',transform = axlist[0].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[0].text(xpos,ypos, "a", transform = axlist[0].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[1].pcolormesh(grb_lon, grb_lat, feb21_ctl_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[1])
axlist[1].set_extent(img_extent, trans)
#axlist[1].text(-.1,.5, "Feb. 2021", horizontalalignment='center',verticalalignment='center',transform = axlist[5].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')
#axlist[7].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[1].text(xpos,ypos, "b", transform = axlist[1].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
axlist[1].set_title("Feb. 2021", fontdict={'fontsize':fontsize_title,'fontweight':'bold'},)

outplt = axlist[2].pcolormesh(srw_lon, srw_lat, jan22_ctl_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[2])
axlist[2].set_extent(img_extent, trans)
#axlist[2].text(-.1,.5, "Jan. 2022", horizontalalignment='center',verticalalignment='center',transform = axlist[10].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[2].text(xpos,ypos, "c", transform = axlist[2].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
axlist[2].set_title("Jan. 2022", fontdict={'fontsize':fontsize_title,'fontweight':'bold'},)

outplt = axlist[3].pcolormesh(srw_lon, srw_lat, nov22_ctl_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[3])
axlist[3].set_extent(img_extent, trans)
#axlist[3].text(-.1,.5, "Nov. 2022", horizontalalignment='center',verticalalignment='center',transform = axlist[15].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')
#axlist[17].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[3].text(xpos,ypos, "d", transform = axlist[3].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
axlist[3].text(1.1,.5,'Control', horizontalalignment='center',verticalalignment='center',transform = axlist[3].transAxes, weight='bold',size=fontsize_title, backgroundcolor="none", rotation='vertical')
axlist[3].set_title("Nov. 2022", fontdict={'fontsize':fontsize_title,'fontweight':'bold'},)

## Row 2 (Static)
outplt = axlist[4].pcolormesh(srw_lon, srw_lat, dec17_sta_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[4])
axlist[4].set_extent(img_extent, trans)
#axlist[4].set_title('Static', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[4].text(xpos,ypos, "e", transform = axlist[4].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[5].pcolormesh(grb_lon, grb_lat, feb21_sta_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[5])
axlist[5].set_extent(img_extent, trans)
#axlist[7].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[5].text(xpos,ypos, "f", transform = axlist[5].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[6].pcolormesh(srw_lon, srw_lat, jan22_sta_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[6])
axlist[6].set_extent(img_extent, trans)
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[6].text(xpos,ypos, "g", transform = axlist[6].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[7].pcolormesh(srw_lon, srw_lat, nov22_sta_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[7])
axlist[7].set_extent(img_extent, trans)
#axlist[17].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[7].text(xpos,ypos, "h", transform = axlist[7].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
axlist[7].text(1.1,.5,'Static', horizontalalignment='center',verticalalignment='center',transform = axlist[7].transAxes, weight='bold',size=fontsize_title, backgroundcolor="none", rotation='vertical')


## Row 3 (Dynamic)
outplt = axlist[8].pcolormesh(srw_lon, srw_lat, dec17_dyn_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[8])
axlist[8].set_extent(img_extent, trans)
#axlist[8].set_title('Dynamic', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[8].text(xpos,ypos, "i", transform = axlist[8].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[9].pcolormesh(grb_lon, grb_lat, feb21_dyn_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[9])
axlist[9].set_extent(img_extent, trans)
#axlist[7].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[9].text(xpos,ypos, "j", transform = axlist[9].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[10].pcolormesh(srw_lon, srw_lat, jan22_dyn_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[10])
axlist[10].set_extent(img_extent, trans)
#axlist[12].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[10].text(xpos,ypos, "k", transform = axlist[10].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[11].pcolormesh(srw_lon, srw_lat, nov22_dyn_lh, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[11])
axlist[11].set_extent(img_extent, trans)
#axlist[12].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[11].text(xpos,ypos, "l", transform = axlist[11].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
axlist[11].text(1.1,.5,'Dynamic', horizontalalignment='center',verticalalignment='center',transform = axlist[11].transAxes, weight='bold',size=fontsize_title, backgroundcolor="none", rotation='vertical')


## Row 4 (Static - Control)
outplt = axlist[12].pcolormesh(srw_lon, srw_lat, dec17_sta_ctl, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[12])
axlist[12].set_extent(img_extent, trans)
#axlist[12].set_title('Static - Control', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[12].text(xpos,ypos, "m", transform = axlist[12].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[13].pcolormesh(grb_lon, grb_lat, feb21_sta_ctl, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[13])
axlist[13].set_extent(img_extent, trans)
#axlist[7].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[13].text(xpos,ypos, "n", transform = axlist[13].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[14].pcolormesh(srw_lon, srw_lat, jan22_sta_ctl, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[14])
axlist[14].set_extent(img_extent, trans)
#axlist[14].set_title('GLSEA Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[14].text(xpos,ypos, "o", transform = axlist[14].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[15].pcolormesh(srw_lon, srw_lat, nov22_sta_ctl, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[15])
axlist[15].set_extent(img_extent, trans)
#axlist[17].set_title('GLSEA Initial', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[15].text(xpos,ypos, "p", transform = axlist[15].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
axlist[15].text(1.1,.5,'Stat - Ctrl', horizontalalignment='center',verticalalignment='center',transform = axlist[15].transAxes, weight='bold',size=fontsize_title, backgroundcolor="none", rotation='vertical')


## Row 5 (Dynamic - Static)
outplt = axlist[16].pcolormesh(srw_lon, srw_lat, dec17_dyn_sta, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[16])
axlist[16].set_extent(img_extent, trans)
#axlist[16].set_title('Dynamic - Static', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[16].text(xpos,ypos, "q", transform = axlist[16].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[17].pcolormesh(grb_lon, grb_lat, feb21_dyn_sta, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[17])
axlist[17].set_extent(img_extent, trans)
#axlist[9].set_title('GLSEA Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[17].text(xpos,ypos, "r", transform = axlist[17].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[18].pcolormesh(srw_lon, srw_lat, jan22_dyn_sta, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[18])
axlist[18].set_extent(img_extent, trans)
#axlist[14].set_title('GLSEA Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[18].text(xpos,ypos, "s", transform = axlist[18].transAxes, weight='bold',size=fontsize, backgroundcolor="none")

outplt = axlist[19].pcolormesh(srw_lon, srw_lat, nov22_dyn_sta, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
plotMap(axlist[19])
axlist[19].set_extent(img_extent, trans)
#axlist[19].set_title('GLSEA Change', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#cb = plt.colorbar(outplt, ax=axlist[1], orientation='horizontal', pad=0.01)# ,shrink=.75 No shrink for erie-Ontario
#cb.ax.tick_params(labelsize=6)
#cb.set_label(label='',size=6)
axlist[19].text(xpos,ypos, "t", transform = axlist[19].transAxes, weight='bold',size=fontsize, backgroundcolor="none")
axlist[19].text(1.1,.5,'Dyn - Stat', horizontalalignment='center',verticalalignment='center',transform = axlist[19].transAxes, weight='bold',size=fontsize_title, backgroundcolor="none", rotation='vertical')
#Colorbars for figures

#ax1_Test = make_axes_locatable(axlist[1])
#cax1 = ax1_Test.append_axes("bottom",size="7%", pad="2%")

cax1 = axlist[16].inset_axes([0.15,-.3,1.75,.2]) #[X-location, Y-Location, Width, Height]  https://matplotlib.org/stable/gallery/subplots_axes_and_figures/colorbar_placement.html
cax2 = axlist[18].inset_axes([0.15,-.3,1.75,.2])

cb1 = fig.colorbar(outplt2, ax=axlist[16:18], cax=cax1, orientation='horizontal', ticks=cLevels[::2]) #, location = 'bottom')#, shrink=.8)
cb1.set_label(label=r'${W/m^2}$', size=8)
cb1.ax.tick_params(labelsize=8)
cb2 = fig.colorbar(outplt, ax=axlist[18:20], cax=cax2, orientation='horizontal', ticks=cLevels_diff)#, location = 'bottom', shrink=.8)
cb2.set_label(label=r'${W/m^2}$', size=8)
cb2.ax.tick_params(labelsize=8)
### Save Figure
plt.subplots_adjust(hspace = -0.05, wspace=0.05, left=0.01,right=0.93,top=.93) #left=0.05,bottom=.1,top=.95,right=0.99,   Lake Erie: (left=0.013,right=0.987,top=0.826,bottom=0.219,hspace=0.053,wspace=0.030)#Lake Michigan: (left=0.013,right=0.987,top=0.956,bottom=0.049,hspace=0.098,wspace=0.030)#(left=0.025,right=0.987,top=0.750,bottom=0.100,hspace=0.225,wspace=0.081) #.034 #hspace = 0.05 for Erie/Ontario

#fig.tight_layout()
plt.savefig(base_fp +'Publication/Figures/LH_Spatial_Comparison_20panels_Transpose.png', dpi=1800)

###Close all datasets
dec17_ctl.close()
dec17_sta.close()
dec17_dyn.close()
feb21_ctl.close() 
feb21_sta.close()
feb21_dyn.close()
jan22_ctl.close() 
jan22_sta.close()
jan22_dyn.close()
nov22_ctl.close() 
nov22_sta.close()
nov22_dyn.close()
