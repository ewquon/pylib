#!/usr/local/bin/python
import sys
import turbsim_field

if len(sys.argv) <= 1:
    print 'USAGE:',sys.argv[0],'input_prefix [output_prefix]'
    print 'Default output_prefix is "vtk/input_prefix" such that files vtk/input_prefix_*.vtk will be written out'
    sys.exit()
prefix = sys.argv[1]
if len(sys.argv) > 2: output = sys.argv[2]
else: output = None

field = turbsim_field.turbsim_field(prefix,verbose=True)
field.writeVTKSeries(prefix=output)
