#!/usr/bin/env python
from __future__ import print_function
import sys
import os
print('Running Python',sys.version)

import numpy as np
import matplotlib.pyplot as plt

from SOWFA.postProcessing.sample_set import uniform, SampleCollection

for dpath in sys.argv[1:]:

    linesTransverse = uniform(dpath) # every 1000 m in x, spanning y, at hub height
    print(linesTransverse)

    # get sampling locations
    xlist = []
    for name in linesTransverse.sampleNames:
        # sample name: hline_x1000_U
        if not name.endswith('_U'): continue
        locstr = name.split('_')[1] # e.g. "x1000"
        xval = float(locstr[1:])
        xlist.append(xval)
    xsamples = np.array(sorted(xlist))

    # set up the collection of data samples
    transverseData = SampleCollection(xsamples,linesTransverse,'hline_x{:d}')
    transverseData.sample_all(sample_list=['U'])

