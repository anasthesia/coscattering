# math
from math import pi, exp, log, log10, sqrt, ceil
import cmath as cm
from mpmath import mp
mp.dps=100
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
mZ = 92000.
mE = 0.511
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
        self._sE0   = (self._sMD2**2-self._sMD1**2-2*mE*self._sMD1)/(2*self._sMD2)  # in MeV

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
        elif (MD2-MD1) < 2*0.511:
            print_warning("The decay of Chi2 to electrons is kinematically forbidden due to the small mass splitting. The number density is not calculated.")
        else:
            # Calculate number density
            self.calc_number_density()

        # The number density of the mediator taken from the decay model
        # (relative to photons) ...
        #self._sN0a  = n0a
        # ... at T = temp0 ...
        #self._sT0   = 10          # in MeV
        # ... corresponding to t = t(temp0)
        #self._st0   = self._sII.time(self._sT0)

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

    def _number_density(self, T):
        sf_ratio = self._sII.scale_factor(self._sT0)/self._sII.scale_factor(T)

        delta_t = self._sII.time(T) - self._st0
        n_gamma = (2.*zeta3)*(self._sT0**3.)/(pi**2.)

        return self._sN0a * n_gamma * sf_ratio**3. * exp( -delta_t/self._tau() )

    def _tau(self): 
        epsilon = mp.mpmathify(self._sEps)
        aX      = mp.mpmathify(self._sAlphaX)
        ttheta  = mp.mpmathify(self._sTtheta)
        MD1     = mp.mpmathify(self._sMD1)
        MD2     = mp.mpmathify(self._sMD2)
        MZp     = mp.mpmathify(self._sMZp)
        MZ      = mp.mpmathify(mZ)

        Gamma=(5*aX*epsilon**2*ttheta**2*((MD1**2 - MD2**2)*(1000*MZ**2 + 80*(25*MD1**2 - 75*MD1*MD2 + 25*MD2**2 - 16*MZ**2 - 50*MZp**2) - (101*MZ**2*(-MD1**4 - MD2**4 + 6*MD1*MD2*MZ**2 - MD2**2*MZ**2 + 2*MZ**4 + MD1**2*(2*MD2**2 - MZ**2)))/(MZ**2 - MZp**2)**2 + ((961*MZ**4 - 1860*MZ**2*MZp**2 + 1000*MZp**4)*(MD1**4 + MD2**4 - 6*MD1*MD2*MZp**2 + MD2**2*MZp**2 - 2*MZp**4 + MD1**2*(-2*MD2**2 + MZp**2)))/(-(MZ**2*MZp) + MZp**3)**2) - 2*(3000*MD1**3*MD2 - 241*MZ**4 - 280*MZ**2*MZp**2 - 3000*MZp**4 + 60*MD1*MD2*(50*MD2**2 - 7*MZ**2 - 100*MZp**2) + 30*MD1**2*(7*MZ**2 + 100*MZp**2) + 30*MD2**2*(7*MZ**2 + 100*MZp**2))*mp.log(MD1/MD2) - (2*MZ**2*(-MD1**2 + 2*MD1*MD2 - MD2**2 + MZ**2)*(241*MZ**8 - 443*MZ**6*MZp**2 - MD1**6*(31*MZ**2 + 70*MZp**2) - 2*MD1**5*MD2*(31*MZ**2 + 70*MZp**2) + MD1**4*MD2**2*(31*MZ**2 + 70*MZp**2) - MD2**6*(31*MZ**2 + 70*MZp**2) + MD2**2*(-210*MZ**6 + 513*MZ**4*MZp**2) + MD1**3*MD2*(-179*MZ**4 + 583*MZ**2*MZp**2 + 4*MD2**2*(31*MZ**2 + 70*MZp**2)) - MD1*MD2*(62*MZ**6 + 140*MZ**4*MZp**2 + 2*MD2**4*(31*MZ**2 + 70*MZp**2) + MD2**2*(179*MZ**4 - 583*MZ**2*MZp**2)) + MD1**2*(-210*MZ**6 + 513*MZ**4*MZp**2 + MD2**4*(31*MZ**2 + 70*MZp**2) + MD2**2*(-358*MZ**4 + 1166*MZ**2*MZp**2)))*mp.log(-(MD1**4 + MD2**2*(MD2**2 - MZ**2 + mp.sqrt(MD1**4 + (MD2**2 - MZ**2)**2 - 2*MD1**2*(MD2**2 + MZ**2))) - MD1**2*(2*MD2**2 + MZ**2 + mp.sqrt(MD1**4 + (MD2**2 - MZ**2)**2 - 2*MD1**2*(MD2**2 + MZ**2))))/(2.*MD1*MD2*MZ**2)))/(mp.sqrt(MD1**4 + (MD2**2 - MZ**2)**2 - 2*MD1**2*(MD2**2 + MZ**2))*(MZ**2 - MZp**2)**3) - (2*(MD1**2 - 2*MD1*MD2 + MD2**2 - MZp**2)*(2*MD1**5*MD2*MZ**2*(31*MZ**2 + 70*MZp**2) - MD1**4*MD2**2*MZ**2*(31*MZ**2 + 70*MZp**2) + MD1**6*(31*MZ**4 + 70*MZ**2*MZp**2) + MD2**6*(31*MZ**4 + 70*MZ**2*MZp**2) + MZp**4*(2883*MZ**6 - 8401*MZ**4*MZp**2 + 8720*MZ**2*MZp**4 - 3000*MZp**6) + MD2**2*(-2883*MZ**6*MZp**2 + 8370*MZ**4*MZp**4 - 8790*MZ**2*MZp**6 + 3000*MZp**8) - MD1**3*MD2*(2883*MZ**6 - 8339*MZ**4*MZp**2 + 8860*MZ**2*MZp**4 - 3000*MZp**6 + 4*MD2**2*(31*MZ**4 + 70*MZ**2*MZp**2)) - MD1**2*(2883*MZ**6*MZp**2 - 8370*MZ**4*MZp**4 + 8790*MZ**2*MZp**6 - 3000*MZp**8 + MD2**4*(31*MZ**4 + 70*MZ**2*MZp**2) + 2*MD2**2*(2883*MZ**6 - 8339*MZ**4*MZp**2 + 8860*MZ**2*MZp**4 - 3000*MZp**6)) + MD1*MD2*(62*MZ**4*MZp**4 + 140*MZ**2*MZp**6 + 2*MD2**4*(31*MZ**4 + 70*MZ**2*MZp**2) + MD2**2*(-2883*MZ**6 + 8339*MZ**4*MZp**2 - 8860*MZ**2*MZp**4 + 3000*MZp**6)))*mp.log(-(MD1**4 + MD2**2*(MD2**2 - MZp**2 + mp.sqrt(MD1**4 + (MD2**2 - MZp**2)**2 - 2*MD1**2*(MD2**2 + MZp**2))) - MD1**2*(2*MD2**2 + MZp**2 + mp.sqrt(MD1**4 + (MD2**2 - MZp**2)**2 - 2*MD1**2*(MD2**2 + MZp**2))))/(2.*MD1*MD2*MZp**2)))/((-MZ**2 + MZp**2)**3*mp.sqrt(MD1**4 + (MD2**2 - MZp**2)**2 - 2*MD1**2*(MD2**2 + MZp**2)))))/(7.374714e6*MD2**3*mp.pi)
        
        res = float(Gamma)

        if abs(res.imag) != 0:
            print_warning("The decay width has an imaginary part Im(Gamma)={} which has been ignored.".format(res.imag))
        if res.real < 0:
            print_warning("The decay width is negative, Gamma={}".format(Gamma.real))
        return 6.6e-22/res.real #in s

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

    def _source_photon_0(self,T):
        return 0.

    def _source_photon_c(self, E, T):
        return 0.

    def _source_electron_0(self,T):
        return 0.

    def _source_electron_c(self, E, T):

        if E <= mE:
            return 0
        elif E >= (self._sMD2**2-self._sMD1**2-2*mE*self._sMD1)/(2*self._sMD2):
            return 0
        else:
            epsilon = mp.mpmathify(self._sEps)
            aX      = mp.mpmathify(self._sAlphaX)
            ttheta  = mp.mpmathify(self._sTtheta)
            MD1     = mp.mpmathify(self._sMD1)
            MD2     = mp.mpmathify(self._sMD2)
            MZp     = mp.mpmathify(self._sMZp)
            MZ      = mp.mpmathify(mZ)
            ME      = mp.mpmathify(mE)
            Ee      = mp.mpmathify(E)

            dGamma_dE = (5*aX*epsilon**2*ttheta**2*((-2000*MD2*mp.sqrt((Ee - ME)*(Ee + ME))*mp.sqrt((MD1**2 + (2*Ee - MD2)*MD2)**2 - 4*MD1**2*ME**2))/(-2*Ee*MD2 + MD2**2 + ME**2) - (202*MD2*mp.sqrt((Ee - ME)*(Ee + ME))*(2*Ee*MD2 - MD2**2 - ME**2)*mp.sqrt((MD1**2 + (2*Ee - MD2)*MD2)**2 - 4*MD1**2*ME**2)*MZ**4*(2*(2*Ee*MD2 - ME**2)*(MD1**2 + 2*Ee*MD2 - MD2**2 - ME**2) + (-MD1**2 - 4*Ee*MD2 + 2*MD1*MD2 + MD2**2 + 2*ME**2)*MZ**2 + MZ**4))/((2*Ee**2*MD2**2 - MD1**2*ME**2 + 2*ME**4 + MD2*mp.sqrt((Ee - ME)*(Ee + ME))*mp.sqrt((MD1**2 + 2*Ee*MD2 - MD2**2)**2 - 4*MD1**2*ME**2) + MD2**2*MZ**2 + ME**2*(MD2**2 + MZ**2) - Ee*MD2*(-MD1**2 + MD2**2 + 4*ME**2 + 2*MZ**2))*(-2*Ee**2*MD2**2 + MD1**2*ME**2 - 2*ME**4 + MD2*mp.sqrt((Ee - ME)*(Ee + ME))*mp.sqrt((MD1**2 + 2*Ee*MD2 - MD2**2)**2 - 4*MD1**2*ME**2) - MD2**2*MZ**2 - ME**2*(MD2**2 + MZ**2) + Ee*MD2*(-MD1**2 + MD2**2 + 4*ME**2 + 2*MZ**2))*(MZ - MZp)**2*(MZ + MZp)**2) - (2*MD2*mp.sqrt((Ee - ME)*(Ee + ME))*(2*Ee*MD2 - MD2**2 - ME**2)*mp.sqrt((MD1**2 + (2*Ee - MD2)*MD2)**2 - 4*MD1**2*ME**2)*(2*(2*Ee*MD2 - ME**2)*(MD1**2 + 2*Ee*MD2 - MD2**2 - ME**2) + (-MD1**2 - 4*Ee*MD2 + 2*MD1*MD2 + MD2**2 + 2*ME**2)*MZp**2 + MZp**4)*(961*MZ**4 - 1860*MZ**2*MZp**2 + 1000*MZp**4))/((MZ - MZp)**2*(MZ + MZp)**2*(2*Ee**2*MD2**2 - MD1**2*ME**2 + 2*ME**4 + MD2*mp.sqrt((Ee - ME)*(Ee + ME))*mp.sqrt((MD1**2 + 2*Ee*MD2 - MD2**2)**2 - 4*MD1**2*ME**2) + MD2**2*MZp**2 + ME**2*(MD2**2 + MZp**2) - Ee*MD2*(-MD1**2 + MD2**2 + 4*ME**2 + 2*MZp**2))*(-2*Ee**2*MD2**2 + MD1**2*ME**2 - 2*ME**4 + MD2*mp.sqrt((Ee - ME)*(Ee + ME))*mp.sqrt((MD1**2 + 2*Ee*MD2 - MD2**2)**2 - 4*MD1**2*ME**2) - MD2**2*MZp**2 - ME**2*(MD2**2 + MZp**2) + Ee*MD2*(-MD1**2 + MD2**2 + 4*ME**2 + 2*MZp**2))) + (MZ**2*(-124*(2*Ee*MD2 - ME**2)*(MD1**2 + 2*Ee*MD2 - MD2**2 - ME**2)*MZ**2 - 39*(MD1**2 + 4*Ee*MD2 - 2*MD1*MD2 - MD2**2 - 2*ME**2)*MZ**4 + 140*MZ**6 - (280*(2*Ee*MD2 - ME**2)*(MD1**2 + 2*Ee*MD2 - MD2**2 - ME**2) - 241*(MD1**2 + 4*Ee*MD2 - 2*MD1*MD2 - MD2**2 - 2*ME**2)*MZ**2 + 342*MZ**4)*MZp**2)*mp.log(-((2*Ee**2*MD2**2 - MD1**2*ME**2 + 2*ME**4 + MD2*mp.sqrt(Ee**2 - ME**2)*mp.sqrt((MD1**2 + (2*Ee - MD2)*MD2)**2 - 4*MD1**2*ME**2) + MD2**2*MZ**2 + ME**2*(MD2**2 + MZ**2) - Ee*MD2*(-MD1**2 + MD2**2 + 4*ME**2 + 2*MZ**2))/(-2*Ee**2*MD2**2 + MD1**2*ME**2 - 2*ME**4 + MD2*mp.sqrt(Ee**2 - ME**2)*mp.sqrt((MD1**2 + (2*Ee - MD2)*MD2)**2 - 4*MD1**2*ME**2) - MD2**2*MZ**2 - ME**2*(MD2**2 + MZ**2) + Ee*MD2*(-MD1**2 + MD2**2 + 4*ME**2 + 2*MZ**2)))))/(MZ**2 - MZp**2)**3 + ((124*(2*Ee*MD2 - ME**2)*(MD1**2 + 2*Ee*MD2 - MD2**2 - ME**2)*MZ**4 - 961*(MD1**2 + 4*Ee*MD2 - 2*MD1*MD2 - MD2**2 - 2*ME**2)*MZ**6 + MZ**2*(280*(2*Ee*MD2 - ME**2)*(MD1**2 + 2*Ee*MD2 - MD2**2 - ME**2) + 2759*(MD1**2 + 4*Ee*MD2 - 2*MD1*MD2 - MD2**2 - 2*ME**2)*MZ**2 + 1922*MZ**4)*MZp**2 - 60*MZ**2*(50*(MD1**2 + 4*Ee*MD2 - 2*MD1*MD2 - MD2**2 - 2*ME**2) + 93*MZ**2)*MZp**4 + 20*(50*(MD1**2 + 4*Ee*MD2 - 2*MD1*MD2 - MD2**2 - 2*ME**2) + 293*MZ**2)*MZp**6 - 2000*MZp**8)*mp.log(-((2*Ee**2*MD2**2 - MD1**2*ME**2 + 2*ME**4 + MD2*mp.sqrt(Ee**2 - ME**2)*mp.sqrt((MD1**2 + (2*Ee - MD2)*MD2)**2 - 4*MD1**2*ME**2) + MD2**2*MZp**2 + ME**2*(MD2**2 + MZp**2) - Ee*MD2*(-MD1**2 + MD2**2 + 4*ME**2 + 2*MZp**2))/(-2*Ee**2*MD2**2 + MD1**2*ME**2 - 2*ME**4 + MD2*mp.sqrt(Ee**2 - ME**2)*mp.sqrt((MD1**2 + (2*Ee - MD2)*MD2)**2 - 4*MD1**2*ME**2) - MD2**2*MZp**2 - ME**2*(MD2**2 + MZp**2) + Ee*MD2*(-MD1**2 + MD2**2 + 4*ME**2 + 2*MZp**2)))))/(MZ**2 - MZp**2)**3))/(1.229119e6*MD2**2*mp.pi)

            res = float(dGamma_dE)
            if abs(res.imag) != 0:
                print("WARNING: The decay width has an imaginary part Im(Gamma)={} which has been ignored.".format(res.imag))
            return self._number_density_new(T)*res.real


#Model with monochromatic spectrum with energy at our peak energy in order to compare results
class ourDecayModelTest(AbstractModel):

    def __init__(self, MD1, MD2, MZp, ttheta, epsilon, aX, n0a):
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
        self._sE0   = 0.5*(self._sMD2**2-self._sMD1**2-2*mE*self._sMD1)/(2*self._sMD2)  # in MeV

        # The number density of the mediator
        # (relative to photons) ...
        # self._sN0a  = n0a
        # ... at T = temp0 ...
        # self._sT0   = temp0           # in MeV
        # ... corresponding to t = t(temp0)
        # self._st0   = self._sII.time(self._sT0)

        # Call the super constructor
        super(ourDecayModelTest, self).__init__(self._sE0, self._sII)

        # The number density of the mediator
        # (relative to photons) ...
        self._sN0a  = n0a
        # ... at T = temp0 ...
        self._sT0   = 10          # in MeV
        # ... corresponding to t = t(temp0)
        self._st0   = self._sII.time(self._sT0)

    # DEPENDENT QUANTITIES ##############################################################

    def _number_density(self, T):
        sf_ratio = self._sII.scale_factor(self._sT0)/self._sII.scale_factor(T)

        delta_t = self._sII.time(T) - self._st0
        n_gamma = (2.*zeta3)*(self._sT0**3.)/(pi**2.)

        return self._sN0a * n_gamma * sf_ratio**3. * exp( -delta_t/self._tau() )

    def _tau(self): 
        epsilon = mp.mpmathify(self._sEps)
        aX      = mp.mpmathify(self._sAlphaX)
        ttheta  = mp.mpmathify(self._sTtheta)
        MD1     = mp.mpmathify(self._sMD1)
        MD2     = mp.mpmathify(self._sMD2)
        MZp     = mp.mpmathify(self._sMZp)
        MZ      = mp.mpmathify(mZ)

        Gamma=(5*aX*epsilon**2*ttheta**2*((MD1**2 - MD2**2)*(1000*MZ**2 + 80*(25*MD1**2 - 75*MD1*MD2 + 25*MD2**2 - 16*MZ**2 - 50*MZp**2) - (101*MZ**2*(-MD1**4 - MD2**4 + 6*MD1*MD2*MZ**2 - MD2**2*MZ**2 + 2*MZ**4 + MD1**2*(2*MD2**2 - MZ**2)))/(MZ**2 - MZp**2)**2 + ((961*MZ**4 - 1860*MZ**2*MZp**2 + 1000*MZp**4)*(MD1**4 + MD2**4 - 6*MD1*MD2*MZp**2 + MD2**2*MZp**2 - 2*MZp**4 + MD1**2*(-2*MD2**2 + MZp**2)))/(-(MZ**2*MZp) + MZp**3)**2) - 2*(3000*MD1**3*MD2 - 241*MZ**4 - 280*MZ**2*MZp**2 - 3000*MZp**4 + 60*MD1*MD2*(50*MD2**2 - 7*MZ**2 - 100*MZp**2) + 30*MD1**2*(7*MZ**2 + 100*MZp**2) + 30*MD2**2*(7*MZ**2 + 100*MZp**2))*mp.log(MD1/MD2) - (2*MZ**2*(-MD1**2 + 2*MD1*MD2 - MD2**2 + MZ**2)*(241*MZ**8 - 443*MZ**6*MZp**2 - MD1**6*(31*MZ**2 + 70*MZp**2) - 2*MD1**5*MD2*(31*MZ**2 + 70*MZp**2) + MD1**4*MD2**2*(31*MZ**2 + 70*MZp**2) - MD2**6*(31*MZ**2 + 70*MZp**2) + MD2**2*(-210*MZ**6 + 513*MZ**4*MZp**2) + MD1**3*MD2*(-179*MZ**4 + 583*MZ**2*MZp**2 + 4*MD2**2*(31*MZ**2 + 70*MZp**2)) - MD1*MD2*(62*MZ**6 + 140*MZ**4*MZp**2 + 2*MD2**4*(31*MZ**2 + 70*MZp**2) + MD2**2*(179*MZ**4 - 583*MZ**2*MZp**2)) + MD1**2*(-210*MZ**6 + 513*MZ**4*MZp**2 + MD2**4*(31*MZ**2 + 70*MZp**2) + MD2**2*(-358*MZ**4 + 1166*MZ**2*MZp**2)))*mp.log(-(MD1**4 + MD2**2*(MD2**2 - MZ**2 + mp.sqrt(MD1**4 + (MD2**2 - MZ**2)**2 - 2*MD1**2*(MD2**2 + MZ**2))) - MD1**2*(2*MD2**2 + MZ**2 + mp.sqrt(MD1**4 + (MD2**2 - MZ**2)**2 - 2*MD1**2*(MD2**2 + MZ**2))))/(2.*MD1*MD2*MZ**2)))/(mp.sqrt(MD1**4 + (MD2**2 - MZ**2)**2 - 2*MD1**2*(MD2**2 + MZ**2))*(MZ**2 - MZp**2)**3) - (2*(MD1**2 - 2*MD1*MD2 + MD2**2 - MZp**2)*(2*MD1**5*MD2*MZ**2*(31*MZ**2 + 70*MZp**2) - MD1**4*MD2**2*MZ**2*(31*MZ**2 + 70*MZp**2) + MD1**6*(31*MZ**4 + 70*MZ**2*MZp**2) + MD2**6*(31*MZ**4 + 70*MZ**2*MZp**2) + MZp**4*(2883*MZ**6 - 8401*MZ**4*MZp**2 + 8720*MZ**2*MZp**4 - 3000*MZp**6) + MD2**2*(-2883*MZ**6*MZp**2 + 8370*MZ**4*MZp**4 - 8790*MZ**2*MZp**6 + 3000*MZp**8) - MD1**3*MD2*(2883*MZ**6 - 8339*MZ**4*MZp**2 + 8860*MZ**2*MZp**4 - 3000*MZp**6 + 4*MD2**2*(31*MZ**4 + 70*MZ**2*MZp**2)) - MD1**2*(2883*MZ**6*MZp**2 - 8370*MZ**4*MZp**4 + 8790*MZ**2*MZp**6 - 3000*MZp**8 + MD2**4*(31*MZ**4 + 70*MZ**2*MZp**2) + 2*MD2**2*(2883*MZ**6 - 8339*MZ**4*MZp**2 + 8860*MZ**2*MZp**4 - 3000*MZp**6)) + MD1*MD2*(62*MZ**4*MZp**4 + 140*MZ**2*MZp**6 + 2*MD2**4*(31*MZ**4 + 70*MZ**2*MZp**2) + MD2**2*(-2883*MZ**6 + 8339*MZ**4*MZp**2 - 8860*MZ**2*MZp**4 + 3000*MZp**6)))*mp.log(-(MD1**4 + MD2**2*(MD2**2 - MZp**2 + mp.sqrt(MD1**4 + (MD2**2 - MZp**2)**2 - 2*MD1**2*(MD2**2 + MZp**2))) - MD1**2*(2*MD2**2 + MZp**2 + mp.sqrt(MD1**4 + (MD2**2 - MZp**2)**2 - 2*MD1**2*(MD2**2 + MZp**2))))/(2.*MD1*MD2*MZp**2)))/((-MZ**2 + MZp**2)**3*mp.sqrt(MD1**4 + (MD2**2 - MZp**2)**2 - 2*MD1**2*(MD2**2 + MZp**2)))))/(7.374714e6*MD2**3*mp.pi)
        
        res = float(Gamma)

        if abs(res.imag) != 0:
            print_warning("The decay width has an imaginary part Im(Gamma)={} which has been ignored.".format(res.imag))
        if res.real < 0:
            print_warning("The decay width is negative, Gamma={}".format(Gamma.real))
        return 6.6e-22/res.real #in s

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

    def _source_photon_0(self,T):
        return 0.

    def _source_photon_c(self, E, T):
        return 0.

    def _source_electron_0(self, T):
        return self._number_density(T) * (hbar/self._tau())

    def _source_electron_c(self, E, T):
        return 0.