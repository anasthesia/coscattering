import numpy as np
from acropolis.models import DecayModel
from acropolis.scans import ScanParameter, BufferedScanner

res = BufferedScanner( DecayModel,
                       mphi = ScanParameter(0,3,200),
                       tau = 10.**7,
                       temp0 = 10.,
                       n0a = ScanParameter(-14,-3,200,fast=True),
                       bree = 1.,
                       braa = 0.).perform_scan()

f=open("tst_ee.dat","w")
np.savetxt(f,res)
f.close()

