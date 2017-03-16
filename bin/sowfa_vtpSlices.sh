#!/bin/bash
set -e
if [ -z "$1" ]; then
    echo 'specify SAMPLENAME, for the following data structure:'
    echo '  ./*/${SAMPLENAME}.vtk'
    exit
fi
name="$1"
postprocess='vtk2vtp.py'

name=${name%.vtk}

i=0
for d in `foamOutputDirs.py`; do
    fname="$d/$name.vtk"
    if [ ! -f "$fname" ]; then continue; fi
    $postprocess $fname
    newfile="${name}_${i}.vtp"
    mv $d/$name.vtp $newfile
    echo "  wrote $newfile" 
    i=$((i+1))
done
