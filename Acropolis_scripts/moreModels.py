# math
from math import pi, exp, log, log10, sqrt, ceil
import cmath as cm
# numpy
import numpy as np
# scipy
from scipy.linalg import expm
# abc
from abc import ABC, abstractmethod

# input
from acropolis.input import InputInterface, locate_sm_file
# nucl
from acropolis.nucl import NuclearReactor, MatrixGenerator
# params
from acropolis.params import zeta3
from acropolis.params import hbar, c_si, me2, alpha, tau_t, me
mZ = np.longdouble(92000)
from acropolis.params import Emin
from acropolis.params import NY
# pprint
from acropolis.pprint import print_info, print_warning
# abstract model
from acropolis.models import AbstractModel

# call
from subprocess import call
# time
from time import time
# interpolation
from scipy.interpolate import griddata

class ourDecayModel(AbstractModel):

    def __init__(self, MD1, MD2, MZp, ttheta, epsilon, aX):
        # Initialize the Input_Interface
        self._sII   = InputInterface( locate_sm_file() )

        # The masses of our particles
        self._sMD1 = MD1            # in MeV
        self._sMD2 = MD2            # in MeV
        self._sMZp = MZp            # in MeV
        # The couplings
        self._sTtheta  = ttheta
        self._sEps = epsilon
        self._sAlphaX = aX
        # The injection energy, we take the maximal energy possible
        self._sE0   = (self._sMD2-self._sMD1)/2  # in MeV

        # The number density of the mediator
        # (relative to photons) ...
        # self._sN0a  = n0a
        # ... at T = temp0 ...
        # self._sT0   = temp0           # in MeV
        # ... corresponding to t = t(temp0)
        # self._st0   = self._sII.time(self._sT0)

        # Call the super constructor
        super(ourDecayModel, self).__init__(self._sE0, self._sII)

        #We only calculate the number density if the lifetime is above 1e3 s
        if self._tau() < 1e3:
            print_warning("The lifetime of the decaying particle is less than 1e3 s, so no photodisintegrationeffect is expected. It's number density is not calculated.")
        else:
            # Calculate number density
            self.calc_number_density()


    # DEPENDENT QUANTITIES ##############################################################

    def calc_number_density(self):
        #Write all parameters to file
        print("\x1B[1;32mINFO   \x1B[0m: Calculating the number denisty")
        start_time = time()
        f = open("param.m","w")
        f.write("m1 = \"%s\"\n"%(self._sMD1))
        f.write("m2 = \"%s\"\n"%(self._sMD2))
        f.write("mZp = \"%s\"\n"%(self._sMZp))
        f.write("tant = \"%s\"\n"%(self._sTtheta))
        f.write("eps = \"%s\"\n"%(self._sEps))
        f.write("ax = \"%s\"\n"%(self._sAlphaX))
        f.write("xmin = 1\n")
        f.write("xmax = %s\n"%(ceil(self._sMD1/min(self._sTrg))))
        #f.write("xmax = 20")
        f.close()

        #Use mathematica script to calculate the number density
        call(['wolframscript', '-script', 'number_density_script.m'])

        #Load the number density from a file
        data = np.loadtxt("number_density.dat")
        self.xRange = data.T[0]
        self.numDens = data.T[1]
        
        end_time = time()
        print_info("Finished after " + str( int( (end_time - start_time)*10 )/10 ) + "s.")

    def _number_density_new(self,T):
        x = self._sMD1/T

        # Logarithmic interpolation of the number density
        res = 10**griddata(self.xRange, np.log10(self.numDens), ([x]), method='linear')[0]
        return res

    def _tau(self): 
        Gamma = (5*self._sAlphaX*self._sEps**2*self._sTtheta**2*(-((-self._sMD1**2 + self._sMD2**2)*(1000*mZ**2 + 80*(25*self._sMD1**2 - 75*self._sMD1*self._sMD2 + 25*self._sMD2**2 - 16*mZ**2 - 50*self._sMZp**2) - (101*mZ**2*(-self._sMD1**4 - self._sMD2**4 + 6*self._sMD1*self._sMD2*mZ**2 - self._sMD2**2*mZ**2 + 2*mZ**4 + self._sMD1**2*(2*self._sMD2**2 - mZ**2)))/(mZ**2 - self._sMZp**2)**2 + ((961*mZ**4 - 1860*mZ**2*self._sMZp**2 + 1000*self._sMZp**4)*(self._sMD1**4 + self._sMD2**4 - 6*self._sMD1*self._sMD2*self._sMZp**2 + self._sMD2**2*self._sMZp**2 - 2*self._sMZp**4 + self._sMD1**2*(-2*self._sMD2**2 + self._sMZp**2)))/(-(mZ**2*self._sMZp) + self._sMZp**3)**2)) - 2*(3000*self._sMD1**3*self._sMD2 - 241*mZ**4 - 280*mZ**2*self._sMZp**2 - 3000*self._sMZp**4 + 60*self._sMD1*self._sMD2*(50*self._sMD2**2 - 7*mZ**2 - 100*self._sMZp**2) + 30*self._sMD1**2*(7*mZ**2 + 100*self._sMZp**2) + 30*self._sMD2**2*(7*mZ**2 + 100*self._sMZp**2))*cm.log(-2*self._sMD1**2) + 2*(3000*self._sMD1**3*self._sMD2 - 241*mZ**4 - 280*mZ**2*self._sMZp**2 - 3000*self._sMZp**4 + 60*self._sMD1*self._sMD2*(50*self._sMD2**2 - 7*mZ**2 - 100*self._sMZp**2) + 30*self._sMD1**2*(7*mZ**2 + 100*self._sMZp**2) + 30*self._sMD2**2*(7*mZ**2 + 100*self._sMZp**2))*cm.log(-2*self._sMD1*self._sMD2) + (2*mZ**2*(-self._sMD1**2 + 2*self._sMD1*self._sMD2 - self._sMD2**2 + mZ**2)*(241*mZ**8 - 443*mZ**6*self._sMZp**2 - self._sMD1**6*(31*mZ**2 + 70*self._sMZp**2) - 2*self._sMD1**5*self._sMD2*(31*mZ**2 + 70*self._sMZp**2) + self._sMD1**4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2) - self._sMD2**6*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(-210*mZ**6 + 513*mZ**4*self._sMZp**2) + self._sMD1**3*self._sMD2*(-179*mZ**4 + 583*mZ**2*self._sMZp**2 + 4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2)) - self._sMD1*self._sMD2*(62*mZ**6 + 140*mZ**4*self._sMZp**2 + 2*self._sMD2**4*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(179*mZ**4 - 583*mZ**2*self._sMZp**2)) + self._sMD1**2*(-210*mZ**6 + 513*mZ**4*self._sMZp**2 + self._sMD2**4*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(-358*mZ**4 + 1166*mZ**2*self._sMZp**2)))*cm.log(-mZ**2))/(cm.sqrt(self._sMD1**4 + (self._sMD2**2 - mZ**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + mZ**2))*(mZ**2 - self._sMZp**2)**3) - (2*mZ**2*(-self._sMD1**2 + 2*self._sMD1*self._sMD2 - self._sMD2**2 + mZ**2)*(241*mZ**8 - 443*mZ**6*self._sMZp**2 - self._sMD1**6*(31*mZ**2 + 70*self._sMZp**2) - 2*self._sMD1**5*self._sMD2*(31*mZ**2 + 70*self._sMZp**2) + self._sMD1**4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2) - self._sMD2**6*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(-210*mZ**6 + 513*mZ**4*self._sMZp**2) + self._sMD1**3*self._sMD2*(-179*mZ**4 + 583*mZ**2*self._sMZp**2 + 4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2)) - self._sMD1*self._sMD2*(62*mZ**6 + 140*mZ**4*self._sMZp**2 + 2*self._sMD2**4*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(179*mZ**4 - 583*mZ**2*self._sMZp**2)) + self._sMD1**2*(-210*mZ**6 + 513*mZ**4*self._sMZp**2 + self._sMD2**4*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(-358*mZ**4 + 1166*mZ**2*self._sMZp**2)))*cm.log((self._sMD1 - self._sMD2)**2 - mZ**2))/(cm.sqrt(self._sMD1**4 + (self._sMD2**2 - mZ**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + mZ**2))*(mZ**2 - self._sMZp**2)**3) + (2*mZ**2*(-self._sMD1**2 + 2*self._sMD1*self._sMD2 - self._sMD2**2 + mZ**2)*(241*mZ**8 - 443*mZ**6*self._sMZp**2 - self._sMD1**6*(31*mZ**2 + 70*self._sMZp**2) - 2*self._sMD1**5*self._sMD2*(31*mZ**2 + 70*self._sMZp**2) + self._sMD1**4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2) - self._sMD2**6*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(-210*mZ**6 + 513*mZ**4*self._sMZp**2) + self._sMD1**3*self._sMD2*(-179*mZ**4 + 583*mZ**2*self._sMZp**2 + 4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2)) - self._sMD1*self._sMD2*(62*mZ**6 + 140*mZ**4*self._sMZp**2 + 2*self._sMD2**4*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(179*mZ**4 - 583*mZ**2*self._sMZp**2)) + self._sMD1**2*(-210*mZ**6 + 513*mZ**4*self._sMZp**2 + self._sMD2**4*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(-358*mZ**4 + 1166*mZ**2*self._sMZp**2)))*cm.log(2*self._sMD1*self._sMD2*(self._sMD1**2 - 2*self._sMD1*self._sMD2 + self._sMD2**2 - mZ**2)))/(cm.sqrt(self._sMD1**4 + (self._sMD2**2 - mZ**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + mZ**2))*(mZ**2 - self._sMZp**2)**3) - (2*mZ**2*(-self._sMD1**2 + 2*self._sMD1*self._sMD2 - self._sMD2**2 + mZ**2)*(241*mZ**8 - 443*mZ**6*self._sMZp**2 - self._sMD1**6*(31*mZ**2 + 70*self._sMZp**2) - 2*self._sMD1**5*self._sMD2*(31*mZ**2 + 70*self._sMZp**2) + self._sMD1**4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2) - self._sMD2**6*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(-210*mZ**6 + 513*mZ**4*self._sMZp**2) + self._sMD1**3*self._sMD2*(-179*mZ**4 + 583*mZ**2*self._sMZp**2 + 4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2)) - self._sMD1*self._sMD2*(62*mZ**6 + 140*mZ**4*self._sMZp**2 + 2*self._sMD2**4*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(179*mZ**4 - 583*mZ**2*self._sMZp**2)) + self._sMD1**2*(-210*mZ**6 + 513*mZ**4*self._sMZp**2 + self._sMD2**4*(31*mZ**2 + 70*self._sMZp**2) + self._sMD2**2*(-358*mZ**4 + 1166*mZ**2*self._sMZp**2)))*cm.log(self._sMD1**4 + self._sMD2**4 - self._sMD2**2*mZ**2 - self._sMD1**2*(2*self._sMD2**2 + mZ**2) + (-self._sMD1**2 + self._sMD2**2)*cm.sqrt(self._sMD1**4 + (self._sMD2**2 - mZ**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + mZ**2))))/(cm.sqrt(self._sMD1**4 + (self._sMD2**2 - mZ**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + mZ**2))*(mZ**2 - self._sMZp**2)**3) + (4*(self._sMD1**2 - 2*self._sMD1*self._sMD2 + self._sMD2**2 - self._sMZp**2)*(2*self._sMD1**5*self._sMD2*mZ**2*(31*mZ**2 + 70*self._sMZp**2) - self._sMD1**4*self._sMD2**2*mZ**2*(31*mZ**2 + 70*self._sMZp**2) + self._sMD1**6*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMD2**6*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMZp**4*(2883*mZ**6 - 8401*mZ**4*self._sMZp**2 + 8720*mZ**2*self._sMZp**4 - 3000*self._sMZp**6) + self._sMD2**2*(-2883*mZ**6*self._sMZp**2 + 8370*mZ**4*self._sMZp**4 - 8790*mZ**2*self._sMZp**6 + 3000*self._sMZp**8) - self._sMD1**3*self._sMD2*(2883*mZ**6 - 8339*mZ**4*self._sMZp**2 + 8860*mZ**2*self._sMZp**4 - 3000*self._sMZp**6 + 4*self._sMD2**2*(31*mZ**4 + 70*mZ**2*self._sMZp**2)) - self._sMD1**2*(2883*mZ**6*self._sMZp**2 - 8370*mZ**4*self._sMZp**4 + 8790*mZ**2*self._sMZp**6 - 3000*self._sMZp**8 + self._sMD2**4*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + 2*self._sMD2**2*(2883*mZ**6 - 8339*mZ**4*self._sMZp**2 + 8860*mZ**2*self._sMZp**4 - 3000*self._sMZp**6)) + self._sMD1*self._sMD2*(62*mZ**4*self._sMZp**4 + 140*mZ**2*self._sMZp**6 + 2*self._sMD2**4*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMD2**2*(-2883*mZ**6 + 8339*mZ**4*self._sMZp**2 - 8860*mZ**2*self._sMZp**4 + 3000*self._sMZp**6)))*cm.log(self._sMZp))/((-mZ**2 + self._sMZp**2)**3*cm.sqrt(self._sMD1**4 + (self._sMD2**2 - self._sMZp**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + self._sMZp**2))) + (2*(self._sMD1**2 - 2*self._sMD1*self._sMD2 + self._sMD2**2 - self._sMZp**2)*(2*self._sMD1**5*self._sMD2*mZ**2*(31*mZ**2 + 70*self._sMZp**2) - self._sMD1**4*self._sMD2**2*mZ**2*(31*mZ**2 + 70*self._sMZp**2) + self._sMD1**6*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMD2**6*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMZp**4*(2883*mZ**6 - 8401*mZ**4*self._sMZp**2 + 8720*mZ**2*self._sMZp**4 - 3000*self._sMZp**6) + self._sMD2**2*(-2883*mZ**6*self._sMZp**2 + 8370*mZ**4*self._sMZp**4 - 8790*mZ**2*self._sMZp**6 + 3000*self._sMZp**8) - self._sMD1**3*self._sMD2*(2883*mZ**6 - 8339*mZ**4*self._sMZp**2 + 8860*mZ**2*self._sMZp**4 - 3000*self._sMZp**6 + 4*self._sMD2**2*(31*mZ**4 + 70*mZ**2*self._sMZp**2)) - self._sMD1**2*(2883*mZ**6*self._sMZp**2 - 8370*mZ**4*self._sMZp**4 + 8790*mZ**2*self._sMZp**6 - 3000*self._sMZp**8 + self._sMD2**4*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + 2*self._sMD2**2*(2883*mZ**6 - 8339*mZ**4*self._sMZp**2 + 8860*mZ**2*self._sMZp**4 - 3000*self._sMZp**6)) + self._sMD1*self._sMD2*(62*mZ**4*self._sMZp**4 + 140*mZ**2*self._sMZp**6 + 2*self._sMD2**4*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMD2**2*(-2883*mZ**6 + 8339*mZ**4*self._sMZp**2 - 8860*mZ**2*self._sMZp**4 + 3000*self._sMZp**6)))*cm.log(2*self._sMD1*self._sMD2*(self._sMD1**2 - 2*self._sMD1*self._sMD2 + self._sMD2**2 - self._sMZp**2)))/((-mZ**2 + self._sMZp**2)**3*cm.sqrt(self._sMD1**4 + (self._sMD2**2 - self._sMZp**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + self._sMZp**2))) - (2*(self._sMD1**2 - 2*self._sMD1*self._sMD2 + self._sMD2**2 - self._sMZp**2)*(2*self._sMD1**5*self._sMD2*mZ**2*(31*mZ**2 + 70*self._sMZp**2) - self._sMD1**4*self._sMD2**2*mZ**2*(31*mZ**2 + 70*self._sMZp**2) + self._sMD1**6*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMD2**6*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMZp**4*(2883*mZ**6 - 8401*mZ**4*self._sMZp**2 + 8720*mZ**2*self._sMZp**4 - 3000*self._sMZp**6) + self._sMD2**2*(-2883*mZ**6*self._sMZp**2 + 8370*mZ**4*self._sMZp**4 - 8790*mZ**2*self._sMZp**6 + 3000*self._sMZp**8) - self._sMD1**3*self._sMD2*(2883*mZ**6 - 8339*mZ**4*self._sMZp**2 + 8860*mZ**2*self._sMZp**4 - 3000*self._sMZp**6 + 4*self._sMD2**2*(31*mZ**4 + 70*mZ**2*self._sMZp**2)) - self._sMD1**2*(2883*mZ**6*self._sMZp**2 - 8370*mZ**4*self._sMZp**4 + 8790*mZ**2*self._sMZp**6 - 3000*self._sMZp**8 + self._sMD2**4*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + 2*self._sMD2**2*(2883*mZ**6 - 8339*mZ**4*self._sMZp**2 + 8860*mZ**2*self._sMZp**4 - 3000*self._sMZp**6)) + self._sMD1*self._sMD2*(62*mZ**4*self._sMZp**4 + 140*mZ**2*self._sMZp**6 + 2*self._sMD2**4*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMD2**2*(-2883*mZ**6 + 8339*mZ**4*self._sMZp**2 - 8860*mZ**2*self._sMZp**4 + 3000*self._sMZp**6)))*cm.log(-(self._sMD1 - self._sMD2)**2 + self._sMZp**2))/((-mZ**2 + self._sMZp**2)**3*cm.sqrt(self._sMD1**4 + (self._sMD2**2 - self._sMZp**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + self._sMZp**2))) - (2*(self._sMD1**2 - 2*self._sMD1*self._sMD2 + self._sMD2**2 - self._sMZp**2)*(2*self._sMD1**5*self._sMD2*mZ**2*(31*mZ**2 + 70*self._sMZp**2) - self._sMD1**4*self._sMD2**2*mZ**2*(31*mZ**2 + 70*self._sMZp**2) + self._sMD1**6*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMD2**6*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMZp**4*(2883*mZ**6 - 8401*mZ**4*self._sMZp**2 + 8720*mZ**2*self._sMZp**4 - 3000*self._sMZp**6) + self._sMD2**2*(-2883*mZ**6*self._sMZp**2 + 8370*mZ**4*self._sMZp**4 - 8790*mZ**2*self._sMZp**6 + 3000*self._sMZp**8) - self._sMD1**3*self._sMD2*(2883*mZ**6 - 8339*mZ**4*self._sMZp**2 + 8860*mZ**2*self._sMZp**4 - 3000*self._sMZp**6 + 4*self._sMD2**2*(31*mZ**4 + 70*mZ**2*self._sMZp**2)) - self._sMD1**2*(2883*mZ**6*self._sMZp**2 - 8370*mZ**4*self._sMZp**4 + 8790*mZ**2*self._sMZp**6 - 3000*self._sMZp**8 + self._sMD2**4*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + 2*self._sMD2**2*(2883*mZ**6 - 8339*mZ**4*self._sMZp**2 + 8860*mZ**2*self._sMZp**4 - 3000*self._sMZp**6)) + self._sMD1*self._sMD2*(62*mZ**4*self._sMZp**4 + 140*mZ**2*self._sMZp**6 + 2*self._sMD2**4*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + self._sMD2**2*(-2883*mZ**6 + 8339*mZ**4*self._sMZp**2 - 8860*mZ**2*self._sMZp**4 + 3000*self._sMZp**6)))*cm.log(self._sMD1**4 + self._sMD2**4 - self._sMD2**2*self._sMZp**2 - self._sMD1**2*(2*self._sMD2**2 + self._sMZp**2) + (-self._sMD1**2 + self._sMD2**2)*cm.sqrt(self._sMD1**4 + (self._sMD2**2 - self._sMZp**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + self._sMZp**2))))/((-mZ**2 + self._sMZp**2)**3*cm.sqrt(self._sMD1**4 + (self._sMD2**2 - self._sMZp**2)**2 - 2*self._sMD1**2*(self._sMD2**2 + self._sMZp**2)))))/(7.374714e6*self._sMD2**3*pi)
        # print(Gamma)
        if abs(Gamma.imag) != 0:
            print_warning("The decay width has an imaginary part Im(Gamma)={} which has been ignored.".format(Gamma.imag))
        if Gamma.real < 0:
            print_warning("The decay width is negative, Gamma={}".format(Gamma.real))
        return 6.6e-22/Gamma.real #in s

    # ABSTRACT METHODS ##################################################################

    def _temperature_range(self):
        # The number of degrees-of-freedom to span
        mag = 2.
        # Calculate the approximate decay temperature
        Td = self._sII.temperature( self._tau() )
        # Calculate Tmin and Tmax from Td
        Td_ofm = log10(Td)
        # Here we choose -1.5 (+0.5) orders of magnitude
        # below (above) the approx. decay temperature,
        # since the main part happens after t = \tau
        Tmin = 10.**(Td_ofm - 3.*mag/4.)
        Tmax = 10.**(Td_ofm + 1.*mag/4.)

        return (Tmin, Tmax)

    # def _number_density(self, T):
    #     sf_ratio = self._sII.scale_factor(self._sT0)/self._sII.scale_factor(T)

    #     delta_t = self._sII.time(T) - self._st0
    #     n_gamma = (2.*zeta3)*(self._sT0**3.)/(pi**2.)

    #     return self._sN0a * n_gamma * sf_ratio**3. * exp( -delta_t/self._tau() )

    def _source_photon_0(self,T):
        return 0.

    def _source_photon_c(self, E, T):
        return 0.

    def _source_electron_0(self,T):
        return 0.

    def _source_electron_c(self, E, T):
        dGamma_dE = ((640*self._sAlphaX*self._sEps**2*pi**2*self._sTtheta**2*(-1000*(2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) - cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + mZ**2) + (101*mZ**4*(2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + 2*self._sMD1*self._sMD2*mZ**2 + 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**2 + mZ**4 + self._sMD1**2*(2*self._sMD2**2 - 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - mZ**2) - self._sMD2**2*(2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + mZ**2)))/((2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) - cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + mZ**2)*(mZ**2 - self._sMZp**2)**2) + ((961*mZ**4 - 1860*mZ**2*self._sMZp**2 + 1000*self._sMZp**4)*(2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + 2*self._sMD1*self._sMD2*self._sMZp**2 + 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*self._sMZp**2 + self._sMZp**4 + self._sMD1**2*(2*self._sMD2**2 - 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - self._sMZp**2) - self._sMD2**2*(2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + self._sMZp**2)))/((mZ**2 - self._sMZp**2)**2*(2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) - cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + self._sMZp**2)) + (mZ**2*(self._sMD1*self._sMD2*(78*mZ**4 - 482*mZ**2*self._sMZp**2) + self._sMD2**2*(124*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**2 - 39*mZ**4 + 280*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*self._sMZp**2 + 241*mZ**2*self._sMZp**2) - 2*(62*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2*mZ**2 - 39*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**4 - 70*mZ**6 + 140*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2*self._sMZp**2 + 241*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**2*self._sMZp**2 + 171*mZ**4*self._sMZp**2) + self._sMD1**2*(124*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**2 - 39*mZ**4 + 280*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*self._sMZp**2 + 241*mZ**2*self._sMZp**2 - 4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2)))*cm.log(2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) - cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + mZ**2))/(mZ**2 - self._sMZp**2)**3 + ((2*self._sMD1*self._sMD2*(961*mZ**6 - 2759*mZ**4*self._sMZp**2 + 3000*mZ**2*self._sMZp**4 - 1000*self._sMZp**6) + self._sMD1**2*(-961*mZ**6 + 2759*mZ**4*self._sMZp**2 - 3000*mZ**2*self._sMZp**4 + 1000*self._sMZp**6 + 4*self._sMD2**2*(31*mZ**4 + 70*mZ**2*self._sMZp**2) - 4*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*(31*mZ**4 + 70*mZ**2*self._sMZp**2)) - self._sMD2**2*(961*mZ**6 - 2759*mZ**4*self._sMZp**2 + 3000*mZ**2*self._sMZp**4 - 1000*self._sMZp**6 + 4*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*(31*mZ**4 + 70*mZ**2*self._sMZp**2)) + 2*(961*mZ**6*self._sMZp**2 - 2790*mZ**4*self._sMZp**4 + 2930*mZ**2*self._sMZp**6 - 1000*self._sMZp**8 + 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + (-2*E*self._sMD2 + self._sMD2**2 + me**2)*(961*mZ**6 - 2759*mZ**4*self._sMZp**2 + 3000*mZ**2*self._sMZp**4 - 1000*self._sMZp**6)))*cm.log(2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) - cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + self._sMZp**2))/(mZ**2 - self._sMZp**2)**3))/1.229119e6 - (640*self._sAlphaX*self._sEps**2*pi**2*self._sTtheta**2*(-1000*(2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) + cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + mZ**2) + (101*mZ**4*(2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + 2*self._sMD1*self._sMD2*mZ**2 + 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**2 + mZ**4 + self._sMD1**2*(2*self._sMD2**2 - 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - mZ**2) - self._sMD2**2*(2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + mZ**2)))/((2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) + cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + mZ**2)*(mZ**2 - self._sMZp**2)**2) + ((961*mZ**4 - 1860*mZ**2*self._sMZp**2 + 1000*self._sMZp**4)*(2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + 2*self._sMD1*self._sMD2*self._sMZp**2 + 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*self._sMZp**2 + self._sMZp**4 + self._sMD1**2*(2*self._sMD2**2 - 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - self._sMZp**2) - self._sMD2**2*(2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + self._sMZp**2)))/((mZ**2 - self._sMZp**2)**2*(2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) + cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + self._sMZp**2)) + (mZ**2*(self._sMD1*self._sMD2*(78*mZ**4 - 482*mZ**2*self._sMZp**2) + self._sMD2**2*(124*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**2 - 39*mZ**4 + 280*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*self._sMZp**2 + 241*mZ**2*self._sMZp**2) - 2*(62*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2*mZ**2 - 39*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**4 - 70*mZ**6 + 140*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2*self._sMZp**2 + 241*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**2*self._sMZp**2 + 171*mZ**4*self._sMZp**2) + self._sMD1**2*(124*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*mZ**2 - 39*mZ**4 + 280*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*self._sMZp**2 + 241*mZ**2*self._sMZp**2 - 4*self._sMD2**2*(31*mZ**2 + 70*self._sMZp**2)))*cm.log(2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) + cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + mZ**2))/(mZ**2 - self._sMZp**2)**3 + ((2*self._sMD1*self._sMD2*(961*mZ**6 - 2759*mZ**4*self._sMZp**2 + 3000*mZ**2*self._sMZp**4 - 1000*self._sMZp**6) + self._sMD1**2*(-961*mZ**6 + 2759*mZ**4*self._sMZp**2 - 3000*mZ**2*self._sMZp**4 + 1000*self._sMZp**6 + 4*self._sMD2**2*(31*mZ**4 + 70*mZ**2*self._sMZp**2) - 4*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*(31*mZ**4 + 70*mZ**2*self._sMZp**2)) - self._sMD2**2*(961*mZ**6 - 2759*mZ**4*self._sMZp**2 + 3000*mZ**2*self._sMZp**4 - 1000*self._sMZp**6 + 4*(-2*E*self._sMD2 + self._sMD2**2 + me**2)*(31*mZ**4 + 70*mZ**2*self._sMZp**2)) + 2*(961*mZ**6*self._sMZp**2 - 2790*mZ**4*self._sMZp**4 + 2930*mZ**2*self._sMZp**6 - 1000*self._sMZp**8 + 2*(-2*E*self._sMD2 + self._sMD2**2 + me**2)**2*(31*mZ**4 + 70*mZ**2*self._sMZp**2) + (-2*E*self._sMD2 + self._sMD2**2 + me**2)*(961*mZ**6 - 2759*mZ**4*self._sMZp**2 + 3000*mZ**2*self._sMZp**4 - 1000*self._sMZp**6)))*cm.log(2*me**2 - (self._sMD2**2*me**2 - me**4 + self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) + 2*me**2*(-2*E*self._sMD2 + self._sMD2**2 + me**2) - (-2*E*self._sMD2 + self._sMD2**2 + me**2)**2 + self._sMD1**2*(-2*E*self._sMD2 + 2*me**2) + cm.sqrt(self._sMD1**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD1**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2))*cm.sqrt(self._sMD2**4 + (2*E*self._sMD2 - self._sMD2**2)**2 - 2*self._sMD2**2*(-2*E*self._sMD2 + self._sMD2**2 + 2*me**2)))/(2.*(-2*E*self._sMD2 + self._sMD2**2 + me**2)) + self._sMZp**2))/(mZ**2 - self._sMZp**2)**3))/1.229119e6)/(64.*self._sMD2*pi**3)
        if abs(dGamma_dE.imag) != 0:
            print("WARNING: The decay width has an imaginary part Im(Gamma)={} which has been ignored.".format(dGamma_dE.imag))
        return self._number_density_new(T)*dGamma_dE.real