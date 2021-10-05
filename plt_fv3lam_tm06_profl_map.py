#!/usr/bin/env python3
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mpl_toolkits
mpl_toolkits.__path__.append('/gpfs/dell2/emc/modeling/noscrub/gwv/py/lib/python/basemap-1.2.1-py3.6-linux-x86_64.egg/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap, maskoceans
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
import argparse
import glob
import os
import pandas as pd

def plot_world_map(lon, lat, lont, latt, data, plotpath):
    # plot generic world map
    fig = plt.figure(figsize=(12,12))
    #ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=240))
    ax = fig.add_subplot(1, 1, 1)

    llcrnrlon = -161.0
    llcrnrlat = 8.5
    urcrnrlon = -22.75
    urcrnrlat = 40.25
    lat_0 = 45.0
    lon_0 = -120.0
    lat_ts = 30.0
    xscale=0.08
    yscale=0.12
    m = Basemap(ax=ax,projection='stere',lat_ts=lat_ts,lat_0=lat_0,lon_0=lon_0,\
                  llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,\
                  urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
                  rsphere=6371229,resolution='l')
    m.fillcontinents(color='LightGrey',zorder=0)
    m.drawcoastlines(linewidth=0.75)
    m.drawstates(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
#   cmap = 'bwr'
#   cbarlabel = 'grid'
#   cs = ax.scatter(lon, lat,s=35,marker="o",c='r')
#   cs = ax.scatter(lont, latt,s=35,marker="s",c='b')

    vmin=np.min(data)
    vmax=np.max(data)
    x,y = m(lon,lat)
    tmp2m_1 = (data - 273.15)*1.8 + 32.0

    units = '\xb0''F'
    clevs = np.linspace(-46,134,61)
    clevsdif = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]
    cm = cmap_t2m()
    difcolors = ['blue','#1874CD','dodgerblue','deepskyblue','turquoise','white','white','#EEEE00','#EEC900','darkorange','orangered','red']
    cmdif = matplotlib.colors.ListedColormap(difcolors)
    norm = matplotlib.colors.BoundaryNorm(clevs, cm.N)
    normdif = matplotlib.colors.BoundaryNorm(clevsdif, cmdif.N)

    cs = m.pcolormesh(x, y, tmp2m_1,cmap=cm,norm=norm)
    cb = m.colorbar(cs, location='bottom',pad=0.05,ticks=[-40,-28,-16,-4,8,20,32,44,56,68,80,92,104,116,128],extend='both')

    plttitle = 'FV3 T at level 10 by %s' % (os.environ['LOGNAME'])
    plt.title(plttitle)
    plt.savefig(plotpath,bbox_inches='tight',dpi=100)
    plt.close('all')

def cmap_t2m():
 # Create colormap for 2-m temperature
 # Modified version of the ncl_t2m colormap from Jacob's ncepy code
    r=np.array([255,128,0,  70, 51, 0,  255,0, 0,  51, 255,255,255,255,255,171,128,128,36,162,255])
    g=np.array([0,  0,  0,  70, 102,162,255,92,128,185,255,214,153,102,0,  0,  0,  68, 36,162,255])
    b=np.array([255,128,128,255,255,255,255,0, 0,  102,0,  112,0,  0,  0,  56, 0,  68, 36,162,255])
    xsize=np.arange(np.size(r))
    r = r/255.
    g = g/255.
    b = b/255.
    red = []
    green = []
    blue = []
    for i in range(len(xsize)):
        xNorm=np.float(i)/(np.float(np.size(r))-1.0)
        red.append([xNorm,r[i],r[i]])
        green.append([xNorm,g[i],g[i]])
        blue.append([xNorm,b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    cmap_t2m_coltbl = matplotlib.colors.LinearSegmentedColormap('CMAP_T2M_COLTBL',colorDict)
    return cmap_t2m_coltbl

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

def read_var_tm06(geopath):
    tmpdata = nc.Dataset(geopath,'r')
    tmplatt = tmpdata.variables['geolat_s'][:]
    tmplat = tmpdata.variables['geolat'][:]
    tmpdata.close()

    arrayshapet = tmplatt.shape
    lontout = np.empty(arrayshapet)
    lattout = np.empty(arrayshapet)

    arrayshape = tmplat.shape
    lonout = np.empty(arrayshape)
    latout = np.empty(arrayshape)

    geonc = nc.Dataset(geopath)
    lat = geonc.variables['geolat'][:]
    lon = geonc.variables['geolon'][:]
    latt = geonc.variables['geolat_s'][:]
    lont = geonc.variables['geolon_s'][:]
    geonc.close()

    print("lat and lon range::")
    print(np.max(lat),np.min(lat))
    print(np.max(lon),np.min(lon))

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
    plot_world_map(lon, lat, lont, latt, data, plotpath)
  
def readfield(corefilepath,tracerfilepath,geopath):

    tmpdata = nc.Dataset(corefilepath,'r')
    u = tmpdata.variables['u_s'][:]
    v = tmpdata.variables['v_s'][:]
    W = tmpdata.variables['w'][:]
    T = tmpdata.variables['t'][:]
    delp = tmpdata.variables['delp'][:]

    sphum = tmpdata.variables['sphum'][:]
    liq_wat = tmpdata.variables['liq_wat'][:]
    ice_wat = tmpdata.variables['ice_wat'][:]
    rainwat = tmpdata.variables['rainwat'][:]
    snowwat = tmpdata.variables['snowwat'][:]
    graupel = tmpdata.variables['graupel'][:]
#   ice_nc = tmpdata.variables['ice_nc'][:]
#   rain_nc = tmpdata.variables['rain_nc'][:]
#   sgs_tke = tmpdata.variables['sgs_tke'][:]
#   tmpdata.close()

#   lon, lat, lont, latt = read_var(geopath)
    lon, lat, lont, latt = read_var_tm06(geopath)

    plotpath ='fv3grid.png'
    arrayshape=lon.shape
    data=np.empty(arrayshape)
    data=T[65,:,:]
    plot_world_map(lon, lat, lont, latt, data, plotpath)
    exit()

    sta_name='bwi'
    stalon=283.31
    stalat=39.18
    ixloc,iyloc=find_loc(lon,lat,stalon,stalat)
    ix=ixloc[0]
    iy=iyloc[0]
    writeprofl(corefilepath,sta_name,ix,iy,\
               u,v,W,T,delp,\
               sphum,liq_wat,ice_wat,rainwat,snowwat,graupel)

#   sta_name='iad'
#   stalon=282.55
#   stalat=38.95
#   ixloc,iyloc=find_loc(lon,lat,stalon,stalat)
#   ix=ixloc[0]
#   iy=iyloc[0]
#   writeprofl(corefilepath,sta_name,ix,iy,\
#              u,v,W,T,DZ,delp,\
#              sphum,liq_wat,ice_wat,rainwat,snowwat,graupel,ice_nc,rain_nc,sgs_tke)

#   sta_name='dca'
#   stalon=292.96
#   stalat=38.85
#   ixloc,iyloc=find_loc(lon,lat,stalon,stalat)
#   ix=ixloc[0]
#   iy=iyloc[0]
#   writeprofl(corefilepath,sta_name,ix,iy,\
#              u,v,W,T,DZ,delp,\
#              sphum,liq_wat,ice_wat,rainwat,snowwat,graupel,ice_nc,rain_nc,sgs_tke)


def writeprofl(corefilepath,sta_name,ix,iy,u,v,W,T,delp,sphum,liq_wat,ice_wat,rainwat,snowwat,graupel):
    uprof=pd.Series(u[:,ix,iy])
    vprof=pd.Series(v[:,ix,iy])
    Wprof=pd.Series(W[:,ix,iy])
    Tprof=pd.Series(T[:,ix,iy])
    delpprof=pd.Series(delp[:,ix,iy])
    sphumprof=pd.Series(sphum[:,ix,iy])
    liq_watprof=pd.Series(liq_wat[:,ix,iy])
    ice_watprof=pd.Series(ice_wat[:,ix,iy])
    rainwatprof=pd.Series(rainwat[:,ix,iy])
    snowwatprof=pd.Series(snowwat[:,ix,iy])
    graupelprof=pd.Series(graupel[:,ix,iy])

    df=pd.DataFrame(
       {
         "u": uprof,
         "v": vprof,
         "W": Wprof,
         "T": Tprof,
         "delp": delpprof,
         "sphum": sphumprof,
         "liq_wat": liq_watprof,
         "ice_wat": ice_watprof,
         "rainwat": rainwatprof,
         "snowwat": snowwatprof,
         "graupel": graupelprof,
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
    y=np.linspace(0, 65, 66)
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
    ax4.scatter(np.flip(df['T'],0),y,s=10,marker="o",c='b')
    plt.title(df.columns[3])

    ax5 = fig.add_subplot(2, ncol, 6)
    ax5.scatter(np.flip(df['delp'],0),y,s=10,marker="o",c='b')
    plt.title(df.columns[4])

    ax6 = fig.add_subplot(2, ncol, 7)
    ax6.scatter(np.flip(df['sphum']*1000,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[5])

    ax7 = fig.add_subplot(2, ncol, 8)
    ax7.scatter(np.flip(df['liq_wat']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[6])

    ax8 = fig.add_subplot(2, ncol, 9)
    ax8.scatter(np.flip(df['ice_wat']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[7])

    ax9 = fig.add_subplot(2, ncol, 10)
    ax9.scatter(np.flip(df['rainwat']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[8])

    ax10 = fig.add_subplot(2, ncol, 11)
    ax10.scatter(np.flip(df['snowwat']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[9])

    ax11 = fig.add_subplot(2, ncol, 12)
    ax11.scatter(np.flip(df['graupel']*1000.0,0),y,s=10,marker="o",c='b')
    plt.title(df.columns[10])

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
