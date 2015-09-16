#!/usr/local/bin/python
import sys
import shlex #for smart splitting to preserve enquoted text
from pylab import *

plotSymbols = False
savePlot = False

# stop reading file if change in any of the data values is > THRESHOLD
#THRESHOLD = 9e9
#THRESHOLD = 1e3
#THRESHOLD = 10

def plot(fname,plotlist):

    data = []
    with open(fname) as f:

        #line = f.readline().split()
        line = shlex.split(f.readline())
        print line
        Nvar = len(line) - 1
        try:
            float(line[0]) # first line is data
            varNames = ['x'] + ['y'+str(i) for i in range(1,Nvar+1)]
        except ValueError: # first line is header
            varNames = line
            line = f.readline().split()
        print Nvar,'dependent variables in file :',varNames

        # save first line of data
        for val in line:
            data.append([float(val)])

        # read rest of file
        prevValues = []
        for line in f:
            values = [float(s) for s in line.split()]
    #        if len(prevValues) > 0:
    #            delta = [abs(v2/v1) for v1,v2 in zip(prevValues,values)]
    #            #if any(delta > THRESHOLD):
    #            if any([r > THRESHOLD for r in delta]):
    #                print 'WARNING: solution blow up detected! (TRESHOLD='+str(THRESHOLD)+')'
    #                print '  previous values',prevValues
    #                print '  current values ',values
    #                break
            # save data
            for col,val in zip(data,values):
                col.append(float(val))
            prevValues = values

    if Nvar==0: #single column of data
        Nvar = 1
        N = len(data[0])
        data = [range(1,N+1), data[0]]
        varNames += ['y']

    # allow specification of specific columns to plot
#    if len(sys.argv[2:]) > 0: 
#        plotlist = [int(var) for var in sys.argv[2:]]
#    else: plotlist = range(1,Nvar+1)
    if len(plotlist)==0: plotlist = range(1,Nvar+1)
    print 'Plotting variables:',[varNames[v] for v in plotlist]

    legendStrs = []
    for var in plotlist:
        if plotSymbols: plot(data[0],data[var],'s-')
        else: plot(data[0],data[var])
        legendStrs.append(varNames[var])
    xlabel(varNames[0])
    if Nvar>1: legend(legendStrs,loc='best')

    if savePlot: 
        fname = '.'.join(sys.argv[1].split('.')[:-1])
        for var in plotlist: fname += '_' + str(var)
        savefig(fname+'.png')
    show()

################################################################################
if __name__ == "__main__":

    if len(sys.argv) <= 1:
        print 'USAGE: ', sys.argv[0],'data_file [col1 col2 ...]'
        print '    where data_file has data in columns with an optional one-line header.'
        print '    Data are delimited by white-space. The first column is assumed to be'
        print '    the independent variable, and an optional list of dependent data'
        print '    columns may be specified for plotting.'
        sys.exit()

    plotlist = [int(var) for var in sys.argv[2:]]
    plot( sys.argv[1], plotlist )
