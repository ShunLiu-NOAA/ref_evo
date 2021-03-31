#!/bin/bash -l
#set -x

echo $SHELL
echo `date`
source ~/bin/loadp.sh

execdir=/gpfs/dell2/emc/modeling/save/Shun.Liu/code/ref_evo
thisdate=`$NDATE -24 |cut -c1-8`
rundir=/gpfs/dell2/emc/modeling/noscrub/Shun.Liu/test/ref_evo/test_$thisdate
mkdir -p $rundir
cd $rundir

allcyc='00 12'
alltmmark='tm06'
allanl='guess anl'
for cyc in $allcyc
do
  for ianl in $allanl
  do
  for tmmark in $alltmmark
  do
  idatadir=/gpfs/dell5/ptmp/emc.campara/fv3lamda/fv3lamda.$thisdate/$cyc
  echo $idatadir/$ianl.$tmmark
  python $execdir/plt_fv3lam_tm06_profl.py -g $idatadir/guess.$tmmark/gfs_data.tile7.nc \
                          -c $idatadir/$ianl.$tmmark/gfs_data.tile7.nc \
                          -t $idatadir/$ianl.$tmmark/gfs_data.tile7.nc
  mv bwi.csv bwi_${thisdate}_${cyc}_${ianl}_$tmmark.csv
  mv iad.csv iad_${thisdate}_${cyc}_${ianl}_$tmmark.csv
  mv dca.csv dca_${thisdate}_${cyc}_${ianl}_$tmmark.csv
  mv bwi_profl.png bwi_profl_${thisdate}_${cyc}_${ianl}_$tmmark.png
  mv iad_profl.png iad_profl_${thisdate}_${cyc}_${ianl}_$tmmark.png
  mv dca_profl.png dca_profl_${thisdate}_${cyc}_${ianl}_$tmmark.png
  exit
  done
  done
done
