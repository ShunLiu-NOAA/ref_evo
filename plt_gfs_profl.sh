
#itile=$1
# conus domain
itile=5
python plt_gfs_profl.py -g /work/noaa/da/Cory.R.Martin/noscrub/UFO_eval/global-workflow/fix/fix_fv3_gmted2010/C768/C768_oro_data.tile${itile}.nc \
   -c /work/noaa/da/Cory.R.Martin/noscrub/UFO_eval/GFS/para/c768_v15ics_v16tag/2020121500/gdas.20201214/18/atmos/RESTART/20201215.030000.fv_core.res.tile6.nc \
   -t /work/noaa/da/Cory.R.Martin/noscrub/UFO_eval/GFS/para/c768_v15ics_v16tag/2020121500/gdas.20201214/18/atmos/RESTART/20201215.030000.fv_tracer.res.tile6.nc
