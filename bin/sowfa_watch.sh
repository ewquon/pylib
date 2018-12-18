#!/bin/bash
#logfile=`ls log.*Solver* | tail -n 1`
solver=`grep application system/controlDict | awk '{print $2}'`
solver=${solver%;*}
logfile=`ls -t log.*$solver | head -n 1`
echo "Following $logfile"
echo 'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv'
tail -f $logfile \
    | grep --line-buffered -e '^Time =' -e '^Courant' -e '-Local Cell' -e 'maximum:' -e '-Boundary Flux' -e 'total - flux'  -e 'No Iterations 0' -e 'Re-reading object' \
    | sed -u 's/^Courant/   -Courant/g' \
    | sed 's/\/ area.*$//g'
