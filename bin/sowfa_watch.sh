#!/bin/bash
logfile=`ls log.*Solver | tail -n 1`
echo "Following $logfile"
echo 'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv'
tail -f $logfile | grep -e '^Time =' -e '^Courant' -e '-Local Cell' -e 'maximum:' -e '-Boundary Flux' -e 'total - flux'
