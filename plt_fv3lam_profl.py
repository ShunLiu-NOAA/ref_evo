#!/usr/bin/env python3
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
import argparse
import glob
import os
import pandas as pd

def plot_world_map(lon, lat, lont, latt, plotpath):
    # plot generic world map
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
#   ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
#   ax.set_extent([-138, -56.5, 17.5, 60.0],crs=ccrs.PlateCarree())
    ax.set_extent([-100.25, -100.0, 40, 40.25],crs=ccrs.PlateCarree())
    cmap = 'viridis'
    cbarlabel = 'grid'
#       cmap = 'bwr'
    cs = ax.scatter(lon, lat,s=35,marker="o",c='r')
    cs = ax.scatter(lont, latt,s=35,marker="s",c='b')
#   cs = ax.pcolormesh(lons, lats, data,vmin=vmin,vmax=vmax,cmap=cmap)
#   cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04)
#   cb.set_label(cbarlabel, fontsize=12)

    plttitle = 'JEDI FV3 grid in 0.25x0.25 box by %s' % (os.environ['LOGNAME'])
    plt.title(plttitle)
    plt.savefig(plotpath,bbox_inches='tight',dpi=100)
    plt.close('all')

def read_var(geopath):
    tmpdata = nc.Dataset(geopath,'r')
    tmplatt = tmpdata.variables['grid_latt'][:]
    tmplat = tmpdata.variables['grid_lat'][:]
    tmpdata.close()

    arrayshapet = tmplatt.shape
    lontout = np.empty(arrayshapet)
    lattout = np.empty(arrayshapet)

    arrayshape = tmplat.shape
    lonout = np.empty(arrayshape)
    latout = np.empty(arrayshape)

    geonc = nc.Dataset(geopath)
    lat = geonc.variables['grid_lat'][:]
    lon = geonc.variables['grid_lon'][:]
    latt = geonc.variables['grid_latt'][:]
    lont = geonc.variables['grid_lont'][:]
    geonc.close()

    latout[:,:] = lat
    lonout[:,:] = lon
    lattout[:,:] = latt
    lontout[:,:] = lont
    return lonout, latout, lontout, lattout

def find_loc(lon,lat,stalon,stalat):
    locx=0
    locy=0
    arrayshape=lon.shape
    print(arrayshape[0])
    nx=arrayshape[0]
    ny=arrayshape[1]
    tmp=np.empty(arrayshape)
    tmp=0.0
         
    tmplon=np.empty(arrayshape)
    tmplon=lon-stalon
    tmplat=np.empty(arrayshape)
    tmplat=lat-stalat
    tmp=np.empty(arrayshape)
    tmp=tmplon*tmplon+tmplat*tmplat

    a=np.where(tmp==np.min(tmp))
    print(a)
    locx=a[0]
    locy=a[1]
    print(lon[locx,locy])
    print(lat[locx,locy])
    return locx, locy
    
def gen_figure(geopath):
    # read the files to get the 2D array to plot
    lon, lat, lont, latt = read_var(geopath)
    plotpath ='fv3grid.png'
    plot_world_map(lon, lat, lont, latt, plotpath)
  
def readfield(corefilepath,tracerfilepath,geopath):

    tmpdata = nc.Dataset(corefilepath,'r')
    u = tmpdata.variables['u'][:]
    v = tmpdata.variables['v'][:]
    W = tmpdata.variables['W'][:]
    T = tmpdata.variables['T'][:]
    DZ = tmpdata.variables['DZ'][:]
    delp = tmpdata.variables['delp'][:]
    tmpdata.close()

    tmpdata = nc.Dataset(tracerfilepath,'r')
    sphum = tmpdata.variables['sphum'][:]
    liq_wat = tmpdata.variables['liq_wat'][:]
    ice_wat = tmpdata.variables['ice_wat'][:]
    rainwat = tmpdata.variables['rainwat'][:]
    snowwat = tmpdata.variables['snowwat'][:]
    graupel = tmpdata.variables['graupel'][:]
    ice_nc = tmpdata.variables['ice_nc'][:]
    rain_nc = tmpdata.variables['rain_nc'][:]
    sgs_tke = tmpdata.variables['sgs_tke'][:]
    tmpdata.close()

    lon, lat, lont, latt = read_var(geopath)

    sta_name='bwi'
    stalon=283.31
    stalat=39.18
    ixloc,iyloc=find_loc(lon,lat,stalon,stalat)
    ix=ixloc[0]
    iy=iyloc[0]
    writeprofl(corefilepath,sta_name,ix,iy,\
               u,v,W,T,DZ,delp,\
               sphum,liq_wat,ice_wat,rainwat,snowwat,graupel,ice_nc,rain_nc,sgs_tke)

    sta_name='iad'
    stalon=282.55
    stalat=38.95
    ixloc,iyloc=find_loc(lon,lat,stalon,stalat)
    ix=ixloc[0]
    iy=iyloc[0]
    writeprofl(corefilepath,sta_name,ix,iy,\
               u,v,W,T,DZ,delp,\
               sphum,liq_wat,ice_wat,rainwat,snowwat,graupel,ice_nc,rain_nc,sgs_tke)

    sta_name='dca'
    stalon=292.96
    stalat=38.85
    ixloc,iyloc=find_loc(lon,lat,stalon,stalat)
    ix=ixloc[0]
    iy=iyloc[0]
    writeprofl(corefilepath,sta_name,ix,iy,\
               u,v,W,T,DZ,delp,\
               sphum,liq_wat,ice_wat,rainwat,snowwat,graupel,ice_nc,rain_nc,sgs_tke)


def writeprofl(corefilepath,sta_name,ix,iy,u,v,W,T,DZ,delp,sphum,liq_wat,ice_wat,rainwat,snowwat,graupel,ice_nc,rain_nc,sgs_tke):
    uprof=pd.Series(u[0,:,ix,iy])
    vprof=pd.Series(v[0,:,ix,iy])
    Wprof=pd.Series(W[0,:,ix,iy])
    Tprof=pd.Series(T[0,:,ix,iy])
    DZprof=pd.Series(DZ[0,:,ix,iy])
    delpprof=pd.Series(delp[0,:,ix,iy])
    sphumprof=pd.Series(sphum[0,:,ix,iy])
    liq_watprof=pd.Series(liq_wat[0,:,ix,iy])
    ice_watprof=pd.Series(ice_wat[0,:,ix,iy])
    rainwatprof=pd.Series(rainwat[0,:,ix,iy])
    snowwatprof=pd.Series(snowwat[0,:,ix,iy])
    graupelprof=pd.Series(graupel[0,:,ix,iy])
    ice_ncprof=pd.Series(ice_nc[0,:,ix,iy])
    rain_ncprof=pd.Series(rain_nc[0,:,ix,iy])
    sgs_tkeprof=pd.Series(sgs_tke[0,:,ix,iy])

    df=pd.DataFrame(
       {
         "u": uprof,
         "v": vprof,
         "W": Wprof,
         "T": Tprof,
         "DZ": DZprof,
         "delp": delpprof,
         "sphum": sphumprof,
         "liq_wat": liq_watprof,
         "ice_wat": ice_watprof,
         "rainwat": rainwatprof,
         "snowwat": snowwatprof,
         "graupel": graupelprof,
         "ice_nc": ice_ncprof,
         "rain_nc": rain_ncprof,
         "sgs_tke": sgs_tkeprof,
       }
    )

    df.to_csv(sta_name + '.csv')
#   print(df.columns[1])

    plt_profl(corefilepath,df,sta_name)

def plt_profl(corefilepath,df,sta_name):

    fileinfo=corefilepath.split('/')
    anafile=fileinfo[-2]
    cyc=fileinfo[-3]
    idate=fileinfo[-4]

#   print(df)
    y=np.linspace(0, 59, 60)
    fig = plt.figure(figsize=(18,9))

    ncol=8
    ax0 = fig.add_subplot(2, ncol, 1)  #row, column, index
    ax0.scatter(np.flip(df.u,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[0])

    ax1 = fig.add_subplot(2, ncol, 2)
    ax1.scatter(np.flip(df.v,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[1])

    ax2 = fig.add_subplot(2, ncol, 3)
    ax2.scatter(np.flip(df.W,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[2])

    ax3 = fig.add_subplot(2, ncol, 4)
    ax3.scatter(np.flip(df['T'],0),y,s=10,marker="o",c='b')
    plt.title(df.columns[3])

    ax4 = fig.add_subplot(2, ncol, 5)
    ax4.scatter(np.flip(df['DZ'],0),y,s=10,marker="o",c='b')
    plt.title(df.columns[4])

    ax5 = fig.add_subplot(2, ncol, 6)
    ax5.scatter(np.flip(df['delp'],0),y,s=10,marker="o",c='b')
    plt.title(df.columns[5])

    ax6 = fig.add_subplot(2, ncol, 7)
    ax6.scatter(np.flip(df['sphum']*1000,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[6])

    ax7 = fig.add_subplot(2, ncol, 8)
    ax7.scatter(np.flip(df['liq_wat']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[7])

    ax8 = fig.add_subplot(2, ncol, 9)
    ax8.scatter(np.flip(df['ice_wat']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[8])

    ax9 = fig.add_subplot(2, ncol, 10)
    ax9.scatter(np.flip(df['rainwat']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[9])

    ax10 = fig.add_subplot(2, ncol, 11)
    ax10.scatter(np.flip(df['snowwat']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[10])

    ax11 = fig.add_subplot(2, ncol, 12)
    ax11.scatter(np.flip(df['graupel']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[11])

    ax12 = fig.add_subplot(2, ncol, 13)
    ax12.scatter(np.flip(df['ice_nc'],0),y,s=10,marker="o",c='b')
    plt.title(df.columns[12])

    ax13 = fig.add_subplot(2, ncol, 14)
    ax13.scatter(np.flip(df['rain_nc'],0),y,s=10,marker="o",c='b')
    plt.title(df.columns[13])

    ax14 = fig.add_subplot(2, ncol, 15)
    ax14.scatter(np.flip(df['sgs_tke'],0),y,s=10,marker="o",c='b')
    plt.title(df.columns[14])

    output_infor=sta_name + '_' + idate + "_cyc" + cyc + '_' + anafile
    fig.text(0.5, 0.04, output_infor, ha='center')
    plt.savefig(sta_name+'_profl.png',bbox_inches='tight',dpi=100)

    print(df.columns[0])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-g', '--geoin', help="path to prefix of input files with geolat/geolon", required=True)
    ap.add_argument('-c', '--core', help="path to prefix of input files with core", required=True)
    ap.add_argument('-t', '--tracer', help="path to prefix of input files with tracer", required=True)
    MyArgs = ap.parse_args()
#   gen_figure(MyArgs.geoin)
#   gen_location(MyArgs.geoin)
    readfield(MyArgs.core,MyArgs.tracer,MyArgs.geoin)
