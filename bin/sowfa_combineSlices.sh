#!/bin/bash
set -e
if [ -z "$1" ]; then
    echo 'specify SAMPLENAME, for the following data structure:'
    echo '  ./*/var1_${SAMPLENAME}.vtk'
    echo '  ./*/var2_${SAMPLENAME}.vtk'
    echo '      ...'
    echo '  ./*/varN_${SAMPLENAME}.vtk'
    exit
fi
name="$1"
postprocess='vtkCombine.py'

i=0
for d in `foamOutputDirs.py`; do
    filelist=`ls -m $d/*_${name}.vtk`
    if [ -z "$filelist" ]; then continue; fi
    filelist=`echo $filelist | sed 's/,//g'`
    echo "Processing $filelist"
    $postprocess $filelist
    newfile="${name}_${i}.vtp"
    mv combined.vtp $newfile
    echo "  wrote $newfile" 
    i=$((i+1))
done
