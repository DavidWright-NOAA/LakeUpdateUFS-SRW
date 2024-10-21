# -*- coding: utf-8 -*-

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import csv
import string
from datetime import datetime as dt

### Global constants
filterValue = 0.0

fontsize = 8 #Font Size
xpos = 5 #Panel Label location
#ypos = -16
yadd = 10 #Offset from the lowest y-position
markersize = 1.2 #.4
linewidth = .8 #.2


MH_Divide = -84.7540+360. #Dividing longitude taken from Ayumi's 2020 paper to divide Michigan and Huron


mi_sta_ctl = []
mi_dyn_sta = []
hu_sta_ctl = []
hu_dyn_sta = []
su_sta_ctl = []
su_dyn_sta = []
er_sta_ctl = []
er_dyn_sta = []
on_sta_ctl = []
on_dyn_sta = []

### File Paths and File Loads
for file in ['Dec2017_sh_spatialMean.txt','Feb2021_sh_spatialMean.txt','Jan2022_sh_spatialMean.txt','Nov2022_sh_spatialMean.txt']:
    with open(file,'r', newline = '\n') as fp:
        all_lines = csv.reader(fp, delimiter='\t')
        next(all_lines) #skip the header row
        #Blank Lists to append to
        mi_ctl_sh = []
        mi_sta_sh = []
        mi_dyn_sh = []

        hu_ctl_sh = []
        hu_sta_sh = []
        hu_dyn_sh = []

        su_ctl_sh = []
        su_sta_sh = []
        su_dyn_sh = []

        er_ctl_sh = []
        er_sta_sh = []
        er_dyn_sh = []

        on_ctl_sh = []
        on_sta_sh = []
        on_dyn_sh = []

        plot_time = []
        for x in all_lines:
            plot_time.append(dt.strptime(x[0],'%Y/%m/%d %H:%M'))
            mi_ctl_sh.append(float(x[1]))
            mi_sta_sh.append(float(x[2]))
            mi_dyn_sh.append(float(x[3]))
    
            hu_ctl_sh.append(float(x[4]))
            hu_sta_sh.append(float(x[5]))
            hu_dyn_sh.append(float(x[6]))
    
            su_ctl_sh.append(float(x[7]))
            su_sta_sh.append(float(x[8]))
            su_dyn_sh.append(float(x[9]))
    
            er_ctl_sh.append(float(x[10]))
            er_sta_sh.append(float(x[11]))
            er_dyn_sh.append(float(x[12]))
    
            on_ctl_sh.append(float(x[13]))
            on_sta_sh.append(float(x[14]))
            on_dyn_sh.append(float(x[15]))            
    
    
    ###Take differences between values, convert lists to numpy array for easy subtraction
    mi_ctl_sh = np.array(mi_ctl_sh)
    mi_sta_sh = np.array(mi_sta_sh)
    mi_dyn_sh = np.array(mi_dyn_sh)
    
    hu_ctl_sh = np.array(hu_ctl_sh)
    hu_sta_sh = np.array(hu_sta_sh)
    hu_dyn_sh = np.array(hu_dyn_sh)
    
    su_ctl_sh = np.array(su_ctl_sh)
    su_sta_sh = np.array(su_sta_sh)
    su_dyn_sh = np.array(su_dyn_sh)
    
    er_ctl_sh = np.array(er_ctl_sh)
    er_sta_sh = np.array(er_sta_sh)
    er_dyn_sh = np.array(er_dyn_sh)
    
    on_ctl_sh = np.array(on_ctl_sh)
    on_sta_sh = np.array(on_sta_sh)
    on_dyn_sh = np.array(on_dyn_sh)
    
    mi_sta_ctl.append(mi_sta_sh - mi_ctl_sh)  
    mi_dyn_sta.append(mi_dyn_sh - mi_sta_sh)
    
    hu_sta_ctl.append(hu_sta_sh - hu_ctl_sh)
    hu_dyn_sta.append(hu_dyn_sh - hu_sta_sh)
    
    su_sta_ctl.append(su_sta_sh - su_ctl_sh)
    su_dyn_sta.append(su_dyn_sh - su_sta_sh)
    
    er_sta_ctl.append(er_sta_sh - er_ctl_sh)
    er_dyn_sta.append(er_dyn_sh - er_sta_sh)
    
    on_sta_ctl.append(on_sta_sh - on_ctl_sh)
    on_dyn_sta.append(on_dyn_sh - on_sta_sh)
    
    #Clear the lists to use with next file
    del(mi_ctl_sh)
    del(mi_sta_sh)
    del(mi_dyn_sh)
    del(hu_ctl_sh)
    del(hu_sta_sh)
    del(hu_dyn_sh)
    del(su_ctl_sh)
    del(su_sta_sh)
    del(su_dyn_sh)
    del(er_ctl_sh)
    del(er_sta_sh)
    del(er_dyn_sh)
    del(on_ctl_sh)
    del(on_sta_sh)
    del(on_dyn_sh)

    plot_time.clear()



### Plotting
fig, axarr = plt.subplots(nrows=5,ncols=4, figsize=(6.5,6.0), sharey='row', sharex=True) #6.5, 6.0 ..... 8.5, 4.0
fig.tight_layout()
axlist = axarr.flatten()

x_time_dec = np.arange(0,su_sta_ctl[0].size,1)
x_time_feb = np.arange(0,su_sta_ctl[1].size,1)
x_time_jan = np.arange(0,su_sta_ctl[2].size,1)
x_time_nov = np.arange(0,su_sta_ctl[3].size,1)
#Lake Superior

outplt0 = axlist[0].plot(x_time_dec, su_sta_ctl[0],'#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, su_dyn_sta, '#b367ff')#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt0 = axlist[0].plot(x_time_dec, su_dyn_sta[0],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt1 = axlist[1].plot(x_time_feb, su_sta_ctl[1],'#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, su_dyn_sta, '#b367ff')#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt1 = axlist[1].plot(x_time_feb, su_dyn_sta[1],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt2 = axlist[2].plot(x_time_jan, su_sta_ctl[2],'#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, su_dyn_sta, '#b367ff')#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt2 = axlist[2].plot(x_time_jan, su_dyn_sta[2],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt3 = axlist[3].plot(x_time_nov, su_sta_ctl[3],'#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, su_dyn_sta, '#b367ff')#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
outplt3 = axlist[3].plot(x_time_nov, su_dyn_sta[3],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)
#axlist[0].set_ylabel('Dec. 2017')
axlist[3].text(1.1,.5, "Superior", horizontalalignment='center',verticalalignment='center',transform = axlist[3].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical') #-.38,.5

#Lake Michigan
outplt4 = axlist[4].plot(x_time_dec, mi_sta_ctl[0], '#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, mi_dyn_sta, '#b367ff')
outplt4 = axlist[4].plot(x_time_dec, mi_dyn_sta[0],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt5 = axlist[5].plot(x_time_feb, mi_sta_ctl[1], '#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, mi_dyn_sta, '#b367ff')
outplt5 = axlist[5].plot(x_time_feb, mi_dyn_sta[1],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt6 = axlist[6].plot(x_time_jan, mi_sta_ctl[2], '#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, mi_dyn_sta, '#b367ff')
outplt6 = axlist[6].plot(x_time_jan, mi_dyn_sta[2],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt7 = axlist[7].plot(x_time_nov, mi_sta_ctl[3], '#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, mi_dyn_sta, '#b367ff')
outplt7 = axlist[7].plot(x_time_nov, mi_dyn_sta[3],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

x_time = np.arange(0,su_sta_ctl[1].size,1)
axlist[7].text(1.1,.5, "Michigan", horizontalalignment='center',verticalalignment='center',transform = axlist[7].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')


#Lake Huron
outplt8 = axlist[8].plot(x_time_dec, hu_sta_ctl[0], '#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, hu_dyn_sta, '#b367ff')
outplt8 = axlist[8].plot(x_time_dec, hu_dyn_sta[0],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt9 = axlist[9].plot(x_time_feb, hu_sta_ctl[1], '#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, hu_dyn_sta, '#b367ff')
outplt9 = axlist[9].plot(x_time_feb, hu_dyn_sta[1],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt10 = axlist[10].plot(x_time_jan, hu_sta_ctl[2], '#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, hu_dyn_sta, '#b367ff')
outplt10 = axlist[10].plot(x_time_jan, hu_dyn_sta[2],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt11 = axlist[11].plot(x_time_nov, hu_sta_ctl[3], '#7fffbe', marker='s',markersize=markersize, linewidth=linewidth)#, x_time, hu_dyn_sta, '#b367ff')
outplt11 = axlist[11].plot(x_time_nov, hu_dyn_sta[3],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

axlist[11].text(1.1,.5, "Huron", horizontalalignment='center',verticalalignment='center',transform = axlist[11].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')


#Lake Erie
outplt12 = axlist[12].plot(x_time_dec, er_sta_ctl[0], '#7fffbe',marker='s',markersize=markersize, linewidth=linewidth)#, x_time, er_dyn_sta, '#b367ff')
outplt12 = axlist[12].plot(x_time_dec, er_dyn_sta[0],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt13 = axlist[13].plot(x_time_feb, er_sta_ctl[1], '#7fffbe',marker='s',markersize=markersize, linewidth=linewidth)#, x_time, er_dyn_sta, '#b367ff')
outplt13 = axlist[13].plot(x_time_feb, er_dyn_sta[1],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt14 = axlist[14].plot(x_time_jan, er_sta_ctl[2], '#7fffbe',marker='s',markersize=markersize, linewidth=linewidth)#, x_time, er_dyn_sta, '#b367ff')
outplt14 = axlist[14].plot(x_time_jan, er_dyn_sta[2],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt15 = axlist[15].plot(x_time_nov, er_sta_ctl[3], '#7fffbe',marker='s',markersize=markersize, linewidth=linewidth)#, x_time, er_dyn_sta, '#b367ff')
outplt15 = axlist[15].plot(x_time_nov, er_dyn_sta[3],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

axlist[15].text(1.1,.5, "Erie", horizontalalignment='center',verticalalignment='center',transform = axlist[15].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')

#Lake Ontario
outplt16 = axlist[16].plot(x_time_dec, on_sta_ctl[0], '#7fffbe',marker='s',markersize=markersize, linewidth=linewidth)#, x_time, on_dyn_sta, '#b367ff')
outplt16 = axlist[16].plot(x_time_dec, on_dyn_sta[0],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt17 = axlist[17].plot(x_time_feb, on_sta_ctl[1], '#7fffbe',marker='s',markersize=markersize, linewidth=linewidth)#, x_time, on_dyn_sta, '#b367ff')
outplt17 = axlist[17].plot(x_time_feb, on_dyn_sta[1],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt18 = axlist[18].plot(x_time_jan, on_sta_ctl[2], '#7fffbe',marker='s',markersize=markersize, linewidth=linewidth)#, x_time, on_dyn_sta, '#b367ff')
outplt18 = axlist[18].plot(x_time_jan, on_dyn_sta[2],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

outplt19 = axlist[19].plot(x_time_nov, on_sta_ctl[3], '#7fffbe',marker='s',markersize=markersize, linewidth=linewidth)#, x_time, on_dyn_sta, '#b367ff')
outplt19 = axlist[19].plot(x_time_nov, on_dyn_sta[3],'#b367ff', marker='^',markersize=markersize, linewidth=linewidth)

axlist[19].text(1.1,.5, "Ontario", horizontalalignment='center',verticalalignment='center',transform = axlist[19].transAxes, weight='bold',size=fontsize, backgroundcolor="white", rotation='vertical')

### Set the same parameters for each subplot, including x and y limits, tick marks, and set subplot labels
axlist[0].set_title('Dec. 2017',weight='bold',size=fontsize,backgroundcolor='none')
axlist[1].set_title('Feb. 2021',weight='bold',size=fontsize,backgroundcolor='none')
axlist[2].set_title('Jan. 2022',weight='bold',size=fontsize,backgroundcolor='none')
axlist[3].set_title('Nov. 2022',weight='bold',size=fontsize,backgroundcolor='none')
#axlist[4].set_title('Lake Ontario',weight='bold',size=fontsize,backgroundcolor='none')

cnt = 0
letters = string.ascii_lowercase
for y in axlist[0:20]:
    y.set_ylim([-60,60])
    ypos = y.get_ylim()[0] + yadd #Y position of letter label
    y.set_xlim([0,66])
    y.set_xticks([])
    y.set_yticks([-60,-40,-20,0,20,40,60])
    y.tick_params(left=True,right=True)
    y.tick_params(axis='y', labelsize=6)
    y.tick_params(axis='x', labelsize=6)
    minor_locator = AutoMinorLocator(2)
    y.yaxis.set_minor_locator(minor_locator)
    y.tick_params(axis='y',which='minor',right=True)
    y.text(xpos,ypos,letters[cnt],weight='bold',size=fontsize, backgroundcolor="none", va='center',ha='center')#"white")
    y.grid(linewidth=0.2)
    cnt = cnt+1
'''
for y in axlist[5:10]:
    y.set_ylim([-22,20])
    ypos = y.get_ylim()[0] + yadd #Y position of letter label
    y.set_xlim([0,66])
    y.set_xticks([])
    y.set_yticks([-20,0,20,40])
    minor_locator = AutoMinorLocator(2)
    y.yaxis.set_minor_locator(minor_locator)
    y.text(xpos,ypos,letters[cnt],weight='bold',size=fontsize, backgroundcolor="white")
    cnt = cnt+1

for y in axlist[10:15]:
    y.set_ylim([-22,50])
    ypos = y.get_ylim()[0] + yadd #Y position of letter label
    y.set_xlim([0,66])
    y.set_xticks([])
    y.set_yticks([-20,0,20,40])
    minor_locator = AutoMinorLocator(2)
    y.yaxis.set_minor_locator(minor_locator)
    y.text(xpos,ypos,letters[cnt],weight='bold',size=fontsize, backgroundcolor="white")
    cnt = cnt+1

for y in axlist[15:20]:
    y.set_ylim([-60,60])
    ypos = y.get_ylim()[0] + yadd*4 #Y position of letter label
    y.set_xlim([0,66])
    y.set_xticks([])
    y.set_yticks([-40,-20,0,20,40])
    minor_locator = AutoMinorLocator(2)
    y.yaxis.set_minor_locator(minor_locator)
    y.text(xpos,ypos,letters[cnt],weight='bold',size=fontsize, backgroundcolor="white")
    cnt = cnt+1
'''
axlist[15].set_xticks([10,20,30,40,50,60])
axlist[16].set_xticks([10,20,30,40,50,60])
axlist[17].set_xticks([10,20,30,40,50,60])
axlist[18].set_xticks([10,20,30,40,50,60])
axlist[19].set_xticks([10,20,30,40,50,60])

fig.text(0.02,.535, r"${W/m^2}$", horizontalalignment='center',verticalalignment='center', size=fontsize, backgroundcolor="none", rotation='vertical') #0.057, .58
fig.text(0.5,.06, r"${Forecast}$ ${Hour}$", horizontalalignment='center',verticalalignment='center', size=fontsize, backgroundcolor="none", rotation='horizontal')

leg = fig.legend(['Static - Control', 'Dynamic - Static'], ncol=2, loc='center', bbox_to_anchor=(.5,.03), markerscale=8.0)
for line in leg.get_lines():
    line.set_linewidth(1.5)

#axlist[0].set_title('Control', fontdict={'fontsize':fontsize,'fontweight':'bold'})
#axlist[0].text(xpos,ypos, "A", transform = axlist[0].transAxes, weight='bold',size=fontsize, backgroundcolor="white")
'''
#Panel 1
outplt = axlist[1].pcolormesh(srw_lon, srw_lat, sta_acc, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
axlist[1].set_extent(img_extent, trans)
axlist[1].set_title('Static', fontdict={'fontsize':fontsize,'fontweight':'bold'})
axlist[1].text(xpos,ypos, "B", transform = axlist[1].transAxes, weight='bold',size=fontsize, backgroundcolor="white")

#Panel 2
outplt = axlist[2].pcolormesh(srw_lon, srw_lat, dyn_acc, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
axlist[2].set_extent(img_extent, trans)
axlist[2].set_title('Dynamic', fontdict={'fontsize':fontsize,'fontweight':'bold'})
axlist[2].text(xpos,ypos, "C", transform = axlist[2].transAxes, weight='bold',size=fontsize, backgroundcolor="white")

#Panel 3
outplt = axlist[3].pcolormesh(snodas_lon, snodas_lat, snodas_acc, cmap=cmap, norm=norm, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
axlist[3].set_extent(img_extent, trans)
axlist[3].set_title('SNODAS', fontdict={'fontsize':fontsize,'fontweight':'bold'})
axlist[3].text(xpos,ypos, "D", transform = axlist[3].transAxes, weight='bold',size=fontsize, backgroundcolor="white")

#Panel 4
outplt2 = axlist[4].pcolormesh(srw_lon, srw_lat, sta_ctl_diff, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
axlist[4].set_extent(img_extent, trans)
axlist[4].set_title('Static - Control', fontdict={'fontsize':fontsize,'fontweight':'bold'})
axlist[4].text(xpos,ypos, "E", transform = axlist[4].transAxes, weight='bold',size=fontsize, backgroundcolor="white")

#Panel 5
outplt2 = axlist[5].pcolormesh(srw_lon, srw_lat, dyn_sta_diff, cmap=cmap_diff, norm=norm_diff, transform=trans)#, clip_box=(img_extent[0],img_extent[2],img_extent[1],img_extent[3]))
axlist[5].set_extent(img_extent, trans)
axlist[5].set_title('Dynamic - Static', fontdict={'fontsize':fontsize,'fontweight':'bold'})
axlist[5].text(xpos,ypos, "F", transform = axlist[5].transAxes, weight='bold',size=fontsize, backgroundcolor="white")


#plt.subplots_adjust(hspace = -.3, wspace=0.0, top=.99) #left=0.05,bottom=.1,top=.95,right=0.99,
#plt.show()
'''
plt.subplots_adjust(hspace = 0.15, wspace=0.0, bottom=.10, left=0.08) #bottom = .19, left = .1
plt.savefig('Publication/Figures/SH_Line_20Panel_ROTATED.png', dpi=1800)

plt.show()