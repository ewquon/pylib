#!/usr/local/bin/python
import sys, os
import turbsim_field

if len(sys.argv) <= 1:
    print 'USAGE:',sys.argv[0],'input_prefix [step] [output_dir]'
    print '  Default output step is 1, i.e. every time step will be output'
    print '  Default output_dir is "vtk" such that files vtk/{input_prefix}_*.vtk will be written out'
    sys.exit()
prefix = sys.argv[1]
if len(sys.argv) > 2: output_step = int(sys.argv[2])
else: output_step = 1
if len(sys.argv) > 3: output_dir = sys.argv[3]
else: output_dir = 'vtk'

if not os.path.isdir(output_dir): 
    print 'Creating output path:',output_dir
    os.makedirs(output_dir)
output_prefix = output_dir + os.sep + prefix

field = turbsim_field.turbsim_field(prefix,verbose=True)
field.writeVTKSeries(prefix=output_prefix,step=output_step)
