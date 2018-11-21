#!/bin/bash
logfile=`ls log.*Solver | tail -n 1`
tail -f $logfile | grep -e '^Time =' -e '^Courant' -e '-Local Cell' -e 'maximum:' -e '-Boundary Flux' -e 'total - flux'
