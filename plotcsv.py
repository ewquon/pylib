#!/usr/bin/python
import sys
import matplotlib.pyplot as plt

# stop reading file if change in any of the data values is > THRESHOLD
#THRESHOLD = 9e9
#THRESHOLD = 1e3
#THRESHOLD = 10

curFig = None
curAxes = None

def new():
    global curFig, curAxes
    curFig = plt.figure()
    curAxes = curFig.add_subplot(111)
    plt.draw()

def plot(fname, plotlist=[], legendnames=[], \
         xname='x', \
         yname='', \
         title='', \
         style='-', \
         savefig='', \
         verbose=False):
    global curFig,curAxes
    if curAxes==None or not plt.fignum_exists(curFig.number): new()

    data = []
    with open(fname) as f:

        # check for single header line
        line = f.readline().strip().split(',')
        if verbose: print line
        Nvar = len(line) - 1
        try:
            float(line[0]) # first line is data
            varNames = [xname] + ['y'+str(i) for i in range(1,Nvar+1)]
        except ValueError: # first line is header
            varNames = line
            line = f.readline().split(',') # get next line
        if verbose: print Nvar,'dependent variables in file :',varNames

        # save first line of data
        for val in line:
            data.append([float(val)])

        # read rest of file
        prevValues = []
        for line in f:
            values = [float(s) for s in line.split(',')]
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
    try:
        if len(plotlist)==0: plotlist = range(1,Nvar+1)
    except TypeError:
        plotlist = [plotlist]
    if verbose: print 'Plotting variables:',[varNames[var] for var in plotlist]

    if len(legendnames)==0:
        legendnames = [ varNames[var] for var in plotlist ]

    for var,lstr in zip(plotlist,legendnames):
        curAxes.plot(data[0],data[var],style,label=lstr)
    curAxes.legend(loc='best')
    curAxes.set_xlabel(varNames[0])
    curAxes.set_ylabel(yname)
    curAxes.set_title(title)

    if not savefig=='':
        plt.savefig(savefig)

    plt.gcf().canvas.draw()


def title(tstr):
    global curFig, curAxes
    curAxes.set_title(tstr)
    curFig.canvas.draw()

################################################################################
if __name__ == "__main__":

    if len(sys.argv) <= 1:
        print 'USAGE: ', sys.argv[0],'data_file [col1 col2 ...] [-save imageName.png]'
        print '    where data_file has data in columns with an optional one-line header.'
        print '    Data are delimited by white-space. The first column is assumed to be'
        print '    the independent variable, and an optional list of dependent data'
        print '    columns may be specified for plotting.'
        sys.exit()

    savename=''
    for i in range(1,len(sys.argv)):
        if sys.argv[i]=='-save':
            sys.argv.pop(i)
            savename = sys.argv.pop(i)
            #print 'saving file to',savename
            break

    #plotlist = [int(var) for var in sys.argv[2:]]
    plotlist = []
    varnames = []
    #style = ''
    for var in sys.argv[2:]:
        try:
            plotlist.append( int(var) )
        except ValueError:
            #style += var
            varnames.append(var)
    #if style=='': style = '-'
    if len(varnames) < len(plotlist): varnames = []

    #plot( sys.argv[1], plotlist, style=style, savefig=savename )
    plot( sys.argv[1], plotlist, legendnames=varnames, savefig=savename )
    #plt.show(block=True)
    if not savename: plt.show()
