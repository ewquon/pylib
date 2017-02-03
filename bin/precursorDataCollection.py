#!/usr/bin/env python
#
# automated postprocessing script to extract input/output quantities for any precursor
# - will attempt to scrape 'setUp' or 'initialConditions' file for input information
#
import sys,os
import SOWFA.postProcessing.averaging
import SOWFA.include
import refs
import datetime
verbose = False
DBname = 'precursorTrainingData.csv'

if len(sys.argv) > 1:
    DBpath = sys.argv[1]
else:
    DBpath = os.environ['HOME'] + os.sep + DBname
if 'VERBOSE' in os.environ:
    verbose = os.environ['VERBOSE']

now = datetime.datetime.now()
timestamp = str(now)
path = os.getcwd()

inputVars = [
        'timestamp','path',
        'zRef','URef','dir',
        'TRef',
        'momSourceType','tempSourceType',
        'rotationPeriod','latitude',
        'heatingRate','z0',
        'betaM','gammaM','betaH','gammaH','alphaH',
        'domainHeight','zInversion','inversionWidth',
        ]
outputVars = [
        'simTime','dt',
        'shearCoeff',
        'veerAngle',
        'Tsurf',
        'dTdz_surf',
        'dTdz_inv',
        'dTdz_upper',
        'TI',
        'Ri',
        ]

refFile = 'refs.yaml' # optional definitions
defaults = { #{{{ 
        # wind setup
        'U0Mag':            8.0,    # desired planar-average wind
        'windHeight':       90.0,
        'dir':              270.0,
        # ABL properties
        'momentumSourceType':       'n/a',
        'temperatureSourceType':    'n/a',
        'rotationPeriod':           1.0e30,
        'latitude':                 0.0,
        # fixedHeatingRate boundary parameters
        'T0':               300.0,  # initial/surface temperature
        'heatingRate':      0.0,    # surface heating rate
        'z0':               0.2,    # surface roughness
        'betaM':            15.0,   # M-O scaling unstable constant
        'gammaM':           4.7,    # M-O scaling stable constant
        'betaH':            9.0,    # M-O scaling unstable constant
        'gammaH':           4.7,    # M-O scaling stable constant
        'alphaH':           0.74,   # M-O scaling constant
        # initialization
        'zMax':             1000.0,
        'zInversion':       600.0,
        'inversionWidth':   1e-5,
        'TBottom':          300.0,
        'TTop':             300.0,
        'TGradUpper':       0.003,
        # turbine properties (for estimating output)
        'D':                126.0,
        # misc
        'g':                9.81,
        }#}}}

# extract input parameters from SOWFA include files
params = None
if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
    params,comments = SOWFA.include.read( sys.argv[1], verbose=verbose )
elif os.path.isfile('setUp'):
    params,comments = SOWFA.include.read( 'setUp', verbose=verbose )
elif os.path.isfile('0.original/include/initialConditions'):
    params,comments = SOWFA.include.read( '0.original/include/initialConditions', verbose=verbose )

if not params is None:
    print 'Checking scraped parameters' 
    for p in params.keys():
        if verbose:
            if p in defaults.keys():
                print '  overwriting default',p,'=',defaults[p],'with',params[p]
            else:
                print '  setting',p,'=',params[p]
        defaults[p] = params[p]

# get reference parameters from file (optional)
ref = refs.read( defaults, refFile=refFile, verbose=verbose ) # allows values to be overridden 
ref.readABLProperties(verbose=verbose)
zRef            = ref.windHeight#{{{
URef            = ref.U0Mag
dir             = ref.dir
TRef            = ref.T0
momSourceType   = ref.momentumSourceType
tempSourceType  = ref.temperatureSourceType
rotationPeriod  = ref.rotationPeriod
latitude        = ref.latitude
heatingRate     = ref.heatingRate
z0              = ref.z0
betaM           = ref.betaM
gammaM          = ref.gammaM
betaH           = ref.betaH
gammaH          = ref.gammaH
alphaH          = ref.alphaH
domainHeight    = ref.zMax
zInversion      = ref.zInversion
inversionWidth  = ref.inversionWidth
TBottom         = ref.TBottom
TTop            = ref.TTop
TGradUpper      = ref.TGradUpper#}}}
if tempSourceType=='computed': heatingRate = 'n/a'

#
# read averaging data
#
data = SOWFA.postProcessing.averaging.read( 'postProcessing/averaging' )
simTime = data.t[-1]
dt = data.dt[-1]

Hhub = zRef
D = ref.D
g = ref.g

#
# calculate all output quantities of interest (at final step)
#
# - temperature gradients (inversion, upper)
dTdz_inv, dTdz_upper = data.calcTGradients(zi=zInversion)

# - shear, veer from mean wind
shearCoeff = data.calcShear( heights=[Hhub/4,Hhub/2,Hhub], zref=Hhub, Uref=URef, verbose=verbose )
veerAngle = data.calcVeer( hmax=Hhub+D/2, verbose=verbose )

# - turbulence level
TIx,TIy,TIz,TIdir,TIxyz,TKE = data.calcTI( heights=Hhub, SFS=True, verbose=verbose )
TI = TIdir

# - stability
Ri = data.calcRichardsonNumber( g=g, zref=Hhub, D=D, verbose=verbose )
Tsurf = data.Tsurf # from calcRichardsonNumber
dTdz_surf = data.dTdz_surf # from calcRichardsonNumber

# - model the initial 

#
# screen output
#
def echo(var):#{{{
    val = globals()[var]
    if isinstance(val,str):
        print '{:25s}: "{:s}"'.format(var,val)   
    else:
        print '{:25s}: {:g}'.format(var,val) #}}}
print 80*'='
print 'INPUTS'
print '------'
for v in inputVars: echo(v)
print ''
print 'OUTPUTS (at end of precursor)'
print '-----------------------------'
for v in outputVars: echo(v)
print ''

#
# update training database (csv file)
#
def writeDB(f,var):#{{{
    val = globals()[var]
    if isinstance(val,str):
        f.write('"{:s}",'.format(val))
    else:
        f.write('{:g},'.format(val)) #}}}
DBexists = os.path.isfile(DBpath) and os.path.getsize(DBpath) > 0
with open(DBpath,'a') as f:
    if not DBexists: # write header
        f.write( ','.join(inputVars+outputVars) + '\n' )
    for v in inputVars: writeDB(f,v)
    for v in outputVars: writeDB(f,v)
    f.write('\n')
print 'Appended training data to',DBpath
