#!/usr/bin/env python
import numpy as np
import NWTC.datatools.SOWFA.timeVaryingMappedBC as bc
points = bc.read_boundary_points('points')
np.savez_compressed('points.npz', points=points)
