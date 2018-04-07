#!/bin/bash
#
# Run from inside postProcessing/inflowPlanes directory
# Assumes 30 sampled planes, ID'd from 0-29
#
set -e
converter="$HOME/inflow/utilities/ensight_planes_to_hawc.py"
mkdir -p inflowPlanes
for i in `seq -f '%02g' 0 29`; do
    time python $converter "inflowPlane_${i}_U" 8.0 &> log.inflowPlane.$i
    mv u.bin inflowPlanes/u_${i}.bin
    mv v.bin inflowPlanes/v_${i}.bin
    mv w.bin inflowPlanes/w_${i}.bin
    mv InflowWind_from_SOWFA.dat InflowWind_${i}.dat
    sed -i "s/u.bin/u_${i}.bin/" InflowWind_${i}.dat
    sed -i "s/v.bin/v_${i}.bin/" InflowWind_${i}.dat
    sed -i "s/w.bin/w_${i}.bin/" InflowWind_${i}.dat
done
