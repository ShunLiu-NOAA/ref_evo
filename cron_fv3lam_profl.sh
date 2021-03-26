#!/bin/bash -l
#set -x

echo $SHELL
echo `date`
source ~/bin/loadp.sh

execdir=/gpfs/dell2/emc/modeling/save/Shun.Liu/code/ref_evo
thisdate=`$NDATE |cut -c1-8`
rundir=/gpfs/dell2/emc/modeling/noscrub/Shun.Liu/test/ref_evo/$thisdate
mkdir -p $rundir
cd $rundir

allcyc='00 12'
alltmmark='tm00 tm01 tm02 tm03 tm04 tm05'
allanl='guess anl'
for cyc in $allcyc
do
  for ianl in $allanl
  do
  for tmmark in $alltmmark
  do
  idatadir=/gpfs/dell5/ptmp/emc.campara/fv3lamda/fv3lamda.$thisdate/$cyc
  echo $idatadir/$ianl.$tmmark
  python $execdir/plt_fv3lam_profl.py -g $idatadir/guess.$tmmark/grid_spec.nc \
                          -c $idatadir/$ianl.$tmmark/fv_core.res.tile1.nc \
                          -t $idatadir/$ianl.$tmmark/fv_tracer.res.tile1.nc
  mv bwi.csv bwi_${thisdate}_${cyc}_${ianl}_$tmmark.csv
  mv iad.csv iad_${thisdate}_${cyc}_${ianl}_$tmmark.csv
  mv dca.csv dca_${thisdate}_${cyc}_${ianl}_$tmmark.csv
  #exit
  done
  done
done
