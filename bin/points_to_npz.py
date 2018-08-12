#!/usr/bin/env python
import numpy as np
import NWTC.datatools.SOWFA.timeVaryingMappedBC as bc
y, z, is_structured = bc.read_boundary_points('points')
np.savez_compressed('points.npz', y=y, z=z)
