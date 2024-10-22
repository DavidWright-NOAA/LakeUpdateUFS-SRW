import sys
import numpy as np
from netCDF4 import Dataset
import math
import scipy.interpolate as interp
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

brdr=0.05
missing_value = -99.9999

# get input file
argvs = sys.argv
ncfile_fvcom_eri  = argvs[1]
ncfile_fvcom_mhs  = argvs[2]
ncfile_fvcom_sup  = argvs[3]
ncfile_fvcom_ont  = argvs[4]
print ncfile_fvcom_eri
print ncfile_fvcom_mhs
print ncfile_fvcom_sup
print ncfile_fvcom_ont


# read "containar" file
ncfile='./out_fv3grid.nc'
fh = Dataset(ncfile, mode='a')

lakemask = fh.variables['glmask'][:]
erimask = fh.variables['erimask'][:]
ontmask = fh.variables['ontmask'][:]
supmask = fh.variables['supmask'][:]
mhsmask = fh.variables['mhsmask'][:]
twsfc      = fh.variables['twsfc'][:]#.reshape((1,399,500))
tisfc      = fh.variables['tisfc'][:]#.reshape((1,399,500))
aicec      = fh.variables['aice'][:]#.reshape((1,399,500))
vicec      = fh.variables['vice'][:]#.reshape((1,399,500))
tsfc      = fh.variables['tsfc'][:]#.reshape((1,399,500))
hsc      = fh.variables['hs'][:]#.reshape((1,399,500))
hlc      = fh.variables['hl'][:]#.reshape((1,399,500))
lon      = fh.variables['geolon'][:]
lat      = fh.variables['geolat'][:]
time      = fh.variables['time'][:]
#lon=lon0[0,:,:]#+360.
#lat=lat0[0,:,:]
#lakemask=lakemask0[0,:,:]
#erimask=erimask0[0,:,:]
#mhsmask=mhsmask0[0,:,:]
#supmask=supmask0[0,:,:]
#ontmask=ontmask0[0,:,:]

# read FVCOM outputs
try: # check if erie outputs exist
	fh_fvcom_eri     = Dataset(ncfile_fvcom_eri,mode='r')
	lon_fvcom_eri    = fh_fvcom_eri.variables['lon'][:]
	lat_fvcom_eri    = fh_fvcom_eri.variables['lat'][:]
	temp_fvcom_eri   = fh_fvcom_eri.variables['temp'][:]
	hs_fvcom_eri   = fh_fvcom_eri.variables['sensible_heat_flux'][:]
	hl_fvcom_eri   = fh_fvcom_eri.variables['latent_heat_flux'][:]
	eri_exist = True
except:
	lon_fvcom_eri    = np.array([missing_value])
        lat_fvcom_eri    = np.array([missing_value])
        temp_fvcom_eri   = np.ones((len(time),1,1))*missing_value
        tice_fvcom_eri   = np.ones((len(time),1))*missing_value
        aice_fvcom_eri   = np.ones((len(time),1))*missing_value
        vice_fvcom_eri   = np.ones((len(time),1))*missing_value
        hs_fvcom_eri   = np.ones((len(time),1))*missing_value
        hl_fvcom_eri   = np.ones((len(time),1))*missing_value
	eri_exist = False
if ( eri_exist ): # check if ice outputs exist
	try: 
		tice_fvcom_eri   = fh_fvcom_eri.variables['tsfc'][:]
		aice_fvcom_eri   = fh_fvcom_eri.variables['aice'][:]
		vice_fvcom_eri   = fh_fvcom_eri.variables['vice'][:]
	except:
		tice_fvcom_eri   = np.zeros((len(time),len(lon_fvcom_eri)))
        	aice_fvcom_eri   = np.zeros((len(time),len(lon_fvcom_eri)))
        	vice_fvcom_eri   = np.zeros((len(time),len(lon_fvcom_eri)))


try: # check if mich-huron outputs exist
	fh_fvcom_mhs     = Dataset(ncfile_fvcom_mhs,mode='r')
	lon_fvcom_mhs    = fh_fvcom_mhs.variables['lon'][:]
	lat_fvcom_mhs    = fh_fvcom_mhs.variables['lat'][:]
	temp_fvcom_mhs   = fh_fvcom_mhs.variables['temp'][:]
#	tice_fvcom_mhs   = fh_fvcom_mhs.variables['tsfc'][:]
#	aice_fvcom_mhs   = fh_fvcom_mhs.variables['aice'][:]
        hs_fvcom_mhs   = fh_fvcom_mhs.variables['sensible_heat_flux'][:]
        hl_fvcom_mhs   = fh_fvcom_mhs.variables['latent_heat_flux'][:]
	mhs_exist = True
except:
        lon_fvcom_mhs    = np.array([missing_value])
        lat_fvcom_mhs    = np.array([missing_value])
        temp_fvcom_mhs   = np.ones((len(time),1,1))*missing_value
        tice_fvcom_mhs   = np.ones((len(time),1))*missing_value
        aice_fvcom_mhs   = np.ones((len(time),1))*missing_value
        vice_fvcom_mhs   = np.ones((len(time),1))*missing_value
        hs_fvcom_mhs   = np.ones((len(time),1))*missing_value
        hl_fvcom_mhs   = np.ones((len(time),1))*missing_value
	mhs_exist = False 
if ( mhs_exist ): # check if ice outputs exist
        try: 
                tice_fvcom_mhs   = fh_fvcom_mhs.variables['tsfc'][:]
                aice_fvcom_mhs   = fh_fvcom_mhs.variables['aice'][:]
                vice_fvcom_mhs   = fh_fvcom_mhs.variables['vice'][:]
        except:
                tice_fvcom_mhs   = np.zeros((len(time),len(lon_fvcom_mhs)))
                aice_fvcom_mhs   = np.zeros((len(time),len(lon_fvcom_mhs)))
                vice_fvcom_mhs   = np.zeros((len(time),len(lon_fvcom_mhs)))

try: # check if sup outputs exist
	fh_fvcom_sup     = Dataset(ncfile_fvcom_sup,mode='r')
	lon_fvcom_sup    = fh_fvcom_sup.variables['lon'][:]
	lat_fvcom_sup    = fh_fvcom_sup.variables['lat'][:]
	temp_fvcom_sup   = fh_fvcom_sup.variables['temp'][:]
#	tice_fvcom_sup   = fh_fvcom_sup.variables['tsfc'][:]
#	aice_fvcom_sup   = fh_fvcom_sup.variables['aice'][:]
        hs_fvcom_sup   = fh_fvcom_sup.variables['sensible_heat_flux'][:]
        hl_fvcom_sup   = fh_fvcom_sup.variables['latent_heat_flux'][:]
	sup_exist = True
except:
        lon_fvcom_sup    = np.array([missing_value])
        lat_fvcom_sup    = np.array([missing_value])
        temp_fvcom_sup   = np.ones((len(time),1,1))*missing_value
        tice_fvcom_sup   = np.ones((len(time),1))*missing_value
        aice_fvcom_sup   = np.ones((len(time),1))*missing_value
        vice_fvcom_sup   = np.ones((len(time),1))*missing_value
        hs_fvcom_sup   = np.ones((len(time),1))*missing_value
        hl_fvcom_sup   = np.ones((len(time),1))*missing_value
	sup_exist = False
if ( sup_exist ): # then check if ice outputs exist
        try: 
                tice_fvcom_sup   = fh_fvcom_sup.variables['tsfc'][:]
                aice_fvcom_sup   = fh_fvcom_sup.variables['aice'][:]
                vice_fvcom_sup   = fh_fvcom_sup.variables['vice'][:]
        except:
                tice_fvcom_sup   = np.zeros((len(time),len(lon_fvcom_sup)))
                aice_fvcom_sup   = np.zeros((len(time),len(lon_fvcom_sup)))
                vice_fvcom_sup   = np.zeros((len(time),len(lon_fvcom_sup)))


try: # check if ont outputs exist
	fh_fvcom_ont     = Dataset(ncfile_fvcom_ont,mode='r')
	lon_fvcom_ont    = fh_fvcom_ont.variables['lon'][:]
	lat_fvcom_ont    = fh_fvcom_ont.variables['lat'][:]
	temp_fvcom_ont   = fh_fvcom_ont.variables['temp'][:]
#	tice_fvcom_ont   = fh_fvcom_ont.variables['tsfc'][:]
#	aice_fvcom_ont   = fh_fvcom_ont.variables['aice'][:]
        hs_fvcom_ont   = fh_fvcom_ont.variables['sensible_heat_flux'][:]
        hl_fvcom_ont   = fh_fvcom_ont.variables['latent_heat_flux'][:]
	ont_exist = True
except:
        lon_fvcom_ont    = np.array([missing_value])
        lat_fvcom_ont    = np.array([missing_value])
        temp_fvcom_ont   = np.ones((len(time),1,1))*missing_value
        tice_fvcom_ont   = np.ones((len(time),1))*missing_value
        aice_fvcom_ont   = np.ones((len(time),1))*missing_value
        vice_fvcom_ont   = np.ones((len(time),1))*missing_value
        hs_fvcom_ont   = np.ones((len(time),1))*missing_value
        hl_fvcom_ont   = np.ones((len(time),1))*missing_value
	ont_exist = False
if ( ont_exist ): # check if ice outputs exist
        try: 
                tice_fvcom_ont   = fh_fvcom_ont.variables['tsfc'][:]
                aice_fvcom_ont   = fh_fvcom_ont.variables['aice'][:]
                vice_fvcom_ont   = fh_fvcom_ont.variables['vice'][:]
        except:
                tice_fvcom_ont   = np.zeros((len(time),len(lon_fvcom_ont)))
                aice_fvcom_ont   = np.zeros((len(time),len(lon_fvcom_ont)))
                vice_fvcom_ont   = np.zeros((len(time),len(lon_fvcom_ont)))

print eri_exist, mhs_exist, sup_exist, ont_exist
print 'array shapes'
print temp_fvcom_eri.shape
print temp_fvcom_sup.shape
print temp_fvcom_ont.shape
print temp_fvcom_mhs.shape

# concatenate 
lon_fvcom0 = np.concatenate([lon_fvcom_eri,lon_fvcom_mhs,lon_fvcom_sup,lon_fvcom_ont])
lat_fvcom0 = np.concatenate([lat_fvcom_eri,lat_fvcom_mhs,lat_fvcom_sup,lat_fvcom_ont])
temp_fvcom0 = np.concatenate((temp_fvcom_eri,temp_fvcom_mhs,temp_fvcom_sup,temp_fvcom_ont),axis=2)
tice_fvcom0 = np.concatenate((tice_fvcom_eri,tice_fvcom_mhs,tice_fvcom_sup,tice_fvcom_ont),axis=1)
aice_fvcom0 = np.concatenate((aice_fvcom_eri,aice_fvcom_mhs,aice_fvcom_sup,aice_fvcom_ont),axis=1)
vice_fvcom0 = np.concatenate((vice_fvcom_eri,vice_fvcom_mhs,vice_fvcom_sup,vice_fvcom_ont),axis=1)
hs_fvcom0 = np.concatenate((hs_fvcom_eri,hs_fvcom_mhs,hs_fvcom_sup,hs_fvcom_ont),axis=1)
hl_fvcom0 = np.concatenate((hl_fvcom_eri,hl_fvcom_mhs,hl_fvcom_sup,hl_fvcom_ont),axis=1)
lon_fvcom = np.ma.masked_where( lon_fvcom0 == missing_value, lon_fvcom0 )
lat_fvcom = np.ma.masked_where( lat_fvcom0 == missing_value, lat_fvcom0 )
temp_fvcom = np.ma.masked_where( temp_fvcom0 == missing_value, temp_fvcom0 )
tice_fvcom = np.ma.masked_where( tice_fvcom0 == missing_value, tice_fvcom0 )
aice_fvcom = np.ma.masked_where( aice_fvcom0 == missing_value, aice_fvcom0 )
vice_fvcom = np.ma.masked_where( vice_fvcom0 == missing_value, vice_fvcom0 )
hs_fvcom = np.ma.masked_where( hs_fvcom0 == missing_value, hs_fvcom0 )
hl_fvcom = np.ma.masked_where( hl_fvcom0 == missing_value, hl_fvcom0 )


#print temp_fvcom.shape

for nn in range(len(time)):
#for nn in range(1):

	print nn,str(time[nn])

	twsfc_fvcom = np.array(temp_fvcom[nn,0,:])
	# RRFS grid GL subset
	#lon_wrfsubset = lon
	#lat_wrfsubset = lat
	#twsfc_wrf0 = np.array(0. * twsfc[nn,:,:]) # make it all zero. overwritten by interpolated values later. 
	#twsfc_wrf = np.array(twsfc_wrf0) # GL subset
#!set region/j=634:910/i=996:1440
        lon_wrfsubset = lon#[648:925,1005:1450]
        lat_wrfsubset = lat#[648:925,1005:1450]
        twsfc_wrf0 = np.array(0. * twsfc[nn,:,:]) # make it all zero. overwritten by interpolated values later. 
        #twsfc_wrf = np.array(twsfc_wrf0[648:925,1005:1450]) # GL subset
        twsfc_wrf = np.array(twsfc_wrf0) # GL subset


	#print lon.shape,twsfc_wrf0.shape, twsfc_wrf.shape,twsfc_fvcom.shape

        tisfc_fvcom = np.array(tice_fvcom[nn,:])
	# HRRR grid GL subset
        tisfc_wrf0 = np.array(0. * tisfc[nn,:,:]) # make it all zero. overwritten by interpolated values later.
        tisfc_wrf = np.array(tisfc_wrf0) # GL subset

        aicec_fvcom = np.array(aice_fvcom[nn,:])
        vicec_fvcom = np.array(vice_fvcom[nn,:])
	# HRRR grid GL subset
        aicec_wrf0 = np.array(0. * aicec[nn,:,:]) # make it all zero. overwritten by interpolated values later.
        aicec_wrf  = np.array(aicec_wrf0) # GL subset
        vicec_wrf0 = np.array(0. * vicec[nn,:,:]) # make it all zero. overwritten by interpolated values later.
        vicec_wrf  = np.array(vicec_wrf0) # GL subset


        hsc_fvcom = np.array(hs_fvcom[nn,:])
        # HRRR grid GL subset
        hsc_wrf0 = np.array(0. * hsc[nn,:,:]) # make it all zero. overwritten by interpolated values later.
        hsc_wrf  = np.array(hsc_wrf0) # GL subset

        hlc_fvcom = np.array(hl_fvcom[nn,:])
        # HRRR grid GL subset
        hlc_wrf0 = np.array(0. * hlc[nn,:,:]) # make it all zero. overwritten by interpolated values later.
        hlc_wrf  = np.array(hlc_wrf0) # GL subset

	# make a map
	map = Basemap(projection='merc',llcrnrlon=-93+360.,llcrnrlat=39.,urcrnrlon=-73.+360.,urcrnrlat=50.,resolution='i')
	map.drawcoastlines()
	map.drawcountries()
	parallels = np.arange(30,50,5.) # make latitude lines ever 5 degrees from 30N-50N
	meridians = np.arange(-95,-70,5.) # make longitude lines every 5 degrees from 95W to 70W
	map.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
	map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

	xwrf,ywrf=map(lon_wrfsubset,lat_wrfsubset)
	x_fvcom,y_fvcom=map(lon_fvcom,lat_fvcom)
	twsfc_wrf=interp.griddata((x_fvcom,y_fvcom),twsfc_fvcom,(xwrf,ywrf),method='nearest')#'linear')
	tisfc_wrf=interp.griddata((x_fvcom,y_fvcom),tisfc_fvcom,(xwrf,ywrf),method='nearest')#'linear')
	aicec_wrf=interp.griddata((x_fvcom,y_fvcom),aicec_fvcom,(xwrf,ywrf),method='nearest')#'linear')
	vicec_wrf=interp.griddata((x_fvcom,y_fvcom),vicec_fvcom,(xwrf,ywrf),method='nearest')#'linear')
	hsc_wrf=interp.griddata((x_fvcom,y_fvcom),hsc_fvcom,(xwrf,ywrf),method='nearest')#'linear')
	hlc_wrf=interp.griddata((x_fvcom,y_fvcom),hlc_fvcom,(xwrf,ywrf),method='nearest')#'linear')

        #print twsfc.max(), hsc_fvcom.max()

	#twsfc_wrf0=twsfc_wrf
	#tisfc_wrf0=tisfc_wrf
	#aicec_wrf0=aicec_wrf

        #twsfc_wrf0[648:925,1005:1450]=twsfc_wrf
        #tisfc_wrf0[648:925,1005:1450]=tisfc_wrf
        #aicec_wrf0[648:925,1005:1450]=aicec_wrf
        twsfc_wrf0=twsfc_wrf
        tisfc_wrf0=tisfc_wrf
        aicec_wrf0=aicec_wrf
        vicec_wrf0=vicec_wrf
        hsc_wrf0=hsc_wrf
        hlc_wrf0=hlc_wrf

        	
	#print x_fvcom.shape,twsfc_fvcom.shape,xwrf.shape
	#print twsfc_wrf.shape, twsfc_wrf0.shape, lakemask.shape
	twsfc_wrf0[lakemask==0.] = missing_value
	tisfc_wrf0[lakemask==0.] = missing_value
	aicec_wrf0[lakemask==0.] = missing_value
	vicec_wrf0[lakemask==0.] = missing_value
	hsc_wrf0[lakemask==0.] = missing_value
	hlc_wrf0[lakemask==0.] = missing_value
	if ( not eri_exist): 
	#	print 'erie_mask'
	        twsfc_wrf0[erimask==1.] =  missing_value
        	tisfc_wrf0[erimask==1.] =  missing_value
        	aicec_wrf0[erimask==1.] =  missing_value
        	vicec_wrf0[erimask==1.] =  missing_value
        	hsc_wrf0[erimask==1.] =  missing_value
        	hlc_wrf0[erimask==1.] =  missing_value
        if ( not mhs_exist):
	#	print 'mhs_mask'
                twsfc_wrf0[mhsmask==1.] =  missing_value
                tisfc_wrf0[mhsmask==1.] =  missing_value
                aicec_wrf0[mhsmask==1.] =  missing_value
                vicec_wrf0[mhsmask==1.] =  missing_value
                hsc_wrf0[mhsmask==1.] =  missing_value
                hlc_wrf0[mhsmask==1.] =  missing_value
        if ( not sup_exist):
	#	print 'sup_mask'
                twsfc_wrf0[supmask==1.] =  missing_value
                tisfc_wrf0[supmask==1.] =  missing_value
                aicec_wrf0[supmask==1.] =  missing_value
                vicec_wrf0[supmask==1.] =  missing_value
                hsc_wrf0[supmask==1.] =  missing_value
                hlc_wrf0[supmask==1.] =  missing_value
        if ( not ont_exist):
	#	print 'ont_mask'
                twsfc_wrf0[ontmask==1.] =  missing_value
                tisfc_wrf0[ontmask==1.] =  missing_value
                aicec_wrf0[ontmask==1.] =  missing_value
                vicec_wrf0[ontmask==1.] =  missing_value
                hsc_wrf0[ontmask==1.] =  missing_value
                hlc_wrf0[ontmask==1.] =  missing_value

 
	#map.pcolor(xwrf,ywrf,tsfc_wrf,cmap=plt.cm.jet)
	#plt.savefig("mask.png")
	twsfc[nn,:,:] = twsfc_wrf0
	tisfc[nn,:,:] = tisfc_wrf0
	aicec[nn,:,:] = aicec_wrf0
	vicec[nn,:,:] = vicec_wrf0
	hsc[nn,:,:] = hsc_wrf0
	hlc[nn,:,:] = hlc_wrf0
	tsfc[nn,:,:]=twsfc_wrf0*(1.0-aicec_wrf0)+tisfc_wrf0*aicec_wrf0


fh.variables['twsfc'][:]=twsfc
fh.variables['tisfc'][:]=tisfc
fh.variables['aice'][:]=aicec
fh.variables['vice'][:]=vicec
fh.variables['tsfc'][:]=tsfc
fh.variables['hs'][:]=hsc
fh.variables['hl'][:]=hlc

fh.close()

if ( eri_exist ):
	fh_fvcom_eri.close()
if ( mhs_exist ):
	fh_fvcom_mhs.close()
if ( sup_exist ):
	fh_fvcom_sup.close()
if ( ont_exist ):
	fh_fvcom_ont.close()

print "Python script completed. "


