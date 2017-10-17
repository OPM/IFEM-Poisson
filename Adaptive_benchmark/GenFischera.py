from splipy import *
from splipy.IO import *
import splipy.curve_factory as cf
import splipy.surface_factory as sf
import splipy.volume_factory as vf
import numpy as np

v = Volume() - [1,1,1];

with G2("Fischera.g2") as f:
    for i in range(2):
        for j in range(2):
            for k in range(2):
                if i==j==k==1: # skip first octant
                    continue
                f.write(v + [i,j,k])
