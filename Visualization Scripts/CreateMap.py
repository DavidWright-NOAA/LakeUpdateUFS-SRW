# -*- coding: utf-8 -*-


import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pygrib
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import netCDF4 as nc
import pygrib
import os


### Draw maps using Cartopy Features
def plotMap(ax):
    ax.coastlines(resolution='10m', linewidth=1.0)
    ax.add_feature(cfeature.STATES.with_scale('10m'),facecolor=[0.9375,0.9375,0.859375], edgecolor='black', linewidth=1.0)
    ax.add_feature(cfeature.LAKES.with_scale('10m'),facecolor=[0.59375,0.71484375,0.8828125], edgecolor='black', linewidth=0.75)
    return ax


### Global constants
base_fp = ''
trans = ccrs.LambertConformal(central_longitude=-84.365, central_latitude = 45.5, standard_parallels=(45.5,45.5)) #ccrs.PlateCarree() #map transformation/projection
pc = ccrs.PlateCarree()
fontsize = 8 #Font Size
xpos = 0.08 #Panel Label location
ypos = 0.15

#Read in Data to get domain bounds
ctl_ds = nc.Dataset(base_fp + 'sfc_data.tile7.halo0.nc')

srw_lat= ctl_ds['geolat'][:,:]
srw_lon= ctl_ds['geolon'][:,:] - 360.0

### Zoom Region for Map
lonMin_Z = srw_lon.min() + 1.0 #-93.8 #Lake Huron: -85.1 #Lake Erie:-83.8 #Lake Michigan: -88.0# # #-92.5#-94.0#-82.083####
lonMax_Z = srw_lon.max() - 1.0 #-74.9 #Lake Huron: -79.1 #Lake Erie:-75.0 #Lake Michigan: -84.0# # #-84.0# -74.0#-77.695####
latMin_Z = srw_lat.min() - 0.5 #39.2 #Lake Huron: 42.8 #Lake Erie: 41.0 #Lake Michigan: 41.5# # #45.2#40.0#41.012####
latMax_Z = srw_lat.max() + 0.1 #51.0 #Lake Huron: 46.5 #Lake Erie: 44.8 #Lake Michigan: 44.5# # #49.3#50.0#43.469####
 
img_extent = [lonMin_Z, lonMax_Z, latMin_Z, latMax_Z] 


polygon_points_latlon = [[srw_lon[0,0],srw_lat[0,0]],[srw_lon[0,-1],srw_lat[0,-1]], [srw_lon[-1,-1],srw_lat[-1,-1]], [srw_lon[-1,0],srw_lat[-1,0]], [srw_lon[0,0],srw_lat[0,0]]]
#polygon = Polygon([])

### Plotting
fig, axlist = plt.subplots(figsize=(8.0,4.0), subplot_kw={'projection': trans})


polygon_points_ll = trans.transform_point(*polygon_points_latlon[0], src_crs=pc)
polygon_points_ul = trans.transform_point(*polygon_points_latlon[1], src_crs=pc)
polygon_points_ur = trans.transform_point(*polygon_points_latlon[2], src_crs=pc)
polygon_points_lr = trans.transform_point(*polygon_points_latlon[3], src_crs=pc)
polygon = Polygon([polygon_points_ll, polygon_points_ul, polygon_points_ur, polygon_points_lr, polygon_points_ll], facecolor='none', edgecolor='red')
plotMap(axlist)

axlist.set_extent(img_extent, pc)

###Set Lake Labels
alpha = .6
#Lake Superior
text_point_LS = trans.transform_point(-90.0,47.5,src_crs=pc)
axlist.text(*text_point_LS, "Lake Superior", color='black', weight='bold', fontsize = 7, bbox=dict(boxstyle='Round, pad=0.1',facecolor='lightgrey',alpha=alpha, linewidth=0))
#Lake Huron
text_point_LH = trans.transform_point(-83.9,45.8,src_crs=pc)
axlist.text(*text_point_LH, "Lake Huron", color='black', weight='bold', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=-55., bbox=dict(boxstyle='Round, pad=0.1',facecolor='lightgrey',alpha=alpha, linewidth=0))
#Lake Michigan
text_point_LM = trans.transform_point(-87.01,42.46,src_crs=pc)
axlist.text(*text_point_LM, "Lake Michigan", color='black', weight='bold', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=78., bbox=dict(boxstyle='Round, pad=0.1',facecolor='lightgrey',alpha=alpha, linewidth=0))
#Lake Erie
text_point_LE = trans.transform_point(-82.2,41.58,src_crs=pc)
axlist.text(*text_point_LE, "Lake Erie", color='black', weight='bold', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=25., bbox=dict(boxstyle='Round, pad=0.1',facecolor='lightgrey',alpha=alpha, linewidth=0))
#Lake Ontario
text_point_LO = trans.transform_point(-79.7,43.35,src_crs=pc)
axlist.text(*text_point_LO, "Lake Ontario", color='black', weight='bold', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=10., bbox=dict(boxstyle='Round, pad=0.1',facecolor='lightgrey',alpha=alpha, linewidth=0))

###Add City/Obs points
#BUFN6
text_point_BUFN6 = trans.transform_point(-78.890,42.878,src_crs=pc)
axlist.plot(*text_point_BUFN6,'g^',markersize=5,alpha=.7,markeredgewidth=0.0)
axlist.text(text_point_BUFN6[0] + 10000, text_point_BUFN6[1] + 10000, "BUFN6", fontsize = 7, ha='left', va='bottom', color='green')
#MKGM4
text_point_MKGM4 = trans.transform_point(-86.339,43.227,src_crs=pc)
axlist.plot(*text_point_MKGM4,'g^',markersize=5,alpha=.7,markeredgewidth=0.0)
axlist.text(text_point_MKGM4[0] + 10000, text_point_MKGM4[1] + 10000, "MKGM4", fontsize = 7, ha='left', va='bottom', color='green')
#TWCO1
text_point_TWCO1 = trans.transform_point(-83.259,41.699,src_crs=pc)
axlist.plot(*text_point_TWCO1,'g^',markersize=5,alpha=.7,markeredgewidth=0.0)
axlist.text(text_point_TWCO1[0] - 10000, text_point_TWCO1[1] + 10000, "TWCO1", fontsize = 7, ha='right', va='bottom', color='green')

###Add State Names
#Michigan
text_point_MI = trans.transform_point(-85.8,42.6,src_crs=pc)
axlist.text(*text_point_MI, "Michigan", color='black', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=00.)
#Ohio
text_point_OH = trans.transform_point(-83.2,40.8,src_crs=pc)
axlist.text(*text_point_OH, "Ohio", color='black', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=00.)
#Pennsylvania
text_point_PA = trans.transform_point(-79.7,41.1,src_crs=pc)
axlist.text(*text_point_PA, "Pennsylvania", color='black', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=00.)
#New York
text_point_NY = trans.transform_point(-78.65,42.15,src_crs=pc)
axlist.text(*text_point_NY, "New\nYork", color='black', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=00.)
#Indiana
text_point_IN = trans.transform_point(-87.1,40.8,src_crs=pc)
axlist.text(*text_point_IN, "Indiana", color='black', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=00.)
#Illinois
text_point_IL = trans.transform_point(-89.7,41.4,src_crs=pc)
axlist.text(*text_point_IL, "Illinois", color='black', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=00.)
#Wisconsin
text_point_WI = trans.transform_point(-91.0,44.5,src_crs=pc)
axlist.text(*text_point_WI, "Wisconsin", color='black', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=00.)
#Minnesota
text_point_MN = trans.transform_point(-94.0,47.4,src_crs=pc)
axlist.text(*text_point_MN, "Minnesota", color='black', fontsize = 7, ha='left', va='bottom', rotation_mode='anchor', rotation=00.)



#Add Polygon for domain
axlist.add_patch(polygon)

### Save Figure

plt.savefig(base_fp +'Publication/Figures/Map_With_Labels_States.png', dpi=1800, bbox_inches='tight')

