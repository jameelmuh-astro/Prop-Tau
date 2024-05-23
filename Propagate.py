#!/usr/local/bin/python
# -*- coding: utf-8 -*-
"""
File : Propagate.py 
Definite version for class Propagate
"""
import os
import sys

# Import other libraries for Prop_Tau
from Constants import *
from Particle import eParticleType,Particle
import Output
import random
import math
from enum import Enum

######### ENUMERATIONS #########
class eWhere(Enum):
    """
    Enumeration for defining the position of a particle relative to Earth's geometry and predefined altitude limits.
    """
    INEARTH = -2
    OUTOFCUTS = -1 
    INEARTH_AND_WITHINCUTS = 0
    OUTEARTH_AND_WITHINCUTS = 1
    LOST = 2
    ANDES = 3
    
class eProcType(Enum):
    """
    Enumeration for defining different types of processes for particle propagation.

    Attributes:
        Earth : Process type for Earth interaction.
        CC : Process type for charged current interaction.
        NC : Process type for neutral current interaction.
        Decay : Process type for particle decay.
    """
    Earth = 0
    CC = 1
    NC = 2
    Decay = 3
    
class eModelType(Enum):
    """
    Enumeration for defining different model types for particle propagation.
    """
    OnlyProton = 1
    Salomon = 2
    Gauss = 3
    Gauss2 = 4
    BreitWigner = 5
    BreitWigner2 = 6

# Global variables for tau decays
read = False
_tauola_fraction = []

ALTITUDE_MAX = 11  # Arbitrary value for ALTITUDE_MAX in km
ALTITUDE_MIN = -0.2  # Arbitrary value for ALTITUDE_MIN in km

class Propagate:
    def __init__(self, model=None, LossType=None, Stoch=None, Bdecay=None, PairProd=None):
        """
        Initialize the Propagation class attributes.

        Parameters:
            model (str): The model type for particle propagation.
            LossType (str): The loss type for particle propagation.
            Stoch (bool): A flag indicating whether to use stochastic interactions.
            Bdecay (bool): A flag indicating whether to include beta decay in the propagation.
            PairProd (bool): A flag indicating whether to include pair production in the propagation.

        Attributes:
            gModel (None): The model object for particle propagation.
            nFactor (int): The number of factors for particle propagation.
            aFactor (int): The number of factors for particle propagation.
            fLosses (int): A flag indicating whether to include losses in the propagation.
            epsIBL_min (float): The minimum value for the interaction length.
            epsIBL_max (float): The maximum value for the interaction length.
            fZ (list): A list of Z values for particle propagation.
            fA (list): A list of A values for particle propagation.
            k_com (float): The inelasticity parameter for particle propagation.
            Gam_com (float): The inelasticity parameter for particle propagation.
            z_com (float): The initial z coordinate for particle propagation.
            s_com (float): The initial s coordinate for particle propagation.
            dtdz_com (float): The differential cross section for particle propagation.
            _tot_length_in_earth (float): The total length of the propagation path in Earth.
            _running_length (float): The current length of the propagation path in Earth.
            _horizontal_angle (float): The horizontal angle of the propagation path.

        """
        # Initializing attributes specific to particle propagation
        self.gModel = None
        self.nFactor = 0
        self.aFactor = 0
        self.fLosses = 0
        self.epsIBL_min = 0
        self.epsIBL_max = 0
        self.fZ = []
        self.fA = []
        self.k_com = 0
        self.Gam_com = 0
        self.z_com = 0
        self.s_com = 0
        self.dtdz_com = 0
        self._tot_length_in_earth = 0
        self._running_length = 0
        self._horizontal_angle = 0
        
        # Initialize basic attributes for particle propagation
        self._E = 0
        self._x = 0
        self._y = 0
        self._z = 0
        self._l_decay = 0
        self._n_bang = 0
        
        # Initialize attributes specific to interactions
        if model is not None and LossType is not None and Stoch is not None and Bdecay is not None and PairProd is not None:
            self.fType = model
            self.fStoch = Stoch
            self.fBdecay = Bdecay
            self.fPairProd = PairProd
            self.fe_1 = 0.0
            self.fe_max = 0.0
            self.feth_1 = 0.0
            self.feth_2 = 0.0
            self.fe0_1 = 0.0
            self.fcsi_1 = 0.0
            self.fdelta_1 = 0.0
            self.fe0_2 = 0.0
            self.fdelta_2 = 0.0
            self.fzita = 0.0
            self.fnelements = 0
            # Lists for nuclear particle properties
            self.fNuclt = []                      # kinetic energy threshold
            self.fNuclh1, self.fNuclh2 = [], []   # height of potential barrier
            self.fNuclx1, self.fNuclx2 = [], []   # position parameter
            self.fNuclw1, self.fNuclw2 = [], []   # weighting factors
            self.fNuclc = []                      # constants for nuclear interactions
            # Lists for alpha particle properties
            self.fAlpht = []                      # kinetic energy threshold
            self.fAlphh1, self.fAlphh2 = [], []   # height of potential barrier
            self.fAlphx1, self.fAlphx2 = [], []   # position parameter
            self.fAlphw1, self.fAlphw2 = [], []   # weighting factors
            self.fAlphc = []                      # constants for alpha particle interactions
            # Magnetic field list
            self.fB = []
    
    

    # Class methods
    def GetNElements(self):
        """ Obtaining the number of elements
        """
        return self.fnelements

    def GetLossType(self):
        """ Obtaining the type of losses
        """
        return self.fLosses

    # Other methods for the class

    def set_x(self, value):
        """ Set the value of x """
        self._x = value
    
    def set_y(self, value):
        """ Set the value of y """
        self._y = value
        
    def set_z(self, value):
        """ Set the value of z """
        self._z = value
    
    def get_x(self):
        """ Get the value of x """
        return self._x
    
    def get_y(self):
        """ Get the value of y """
        return self._y
    
    def get_z(self):
        """ Get the value of z """
        return self._z
    
    def increment_running_length(self, v):
        """ Incrementing the traversal length
        """
        self._running_length += v

    def set_tot_length_in_earth(self, v):
        """ Defining the total length in Earth
        """
        self._tot_length_in_earth = v

    def set_horizontal_angle(self, v):
        """ Defining the horizontal angle
        """
        self._horizontal_angle = v

    def get_tot_length_in_earth(self):
        """ Obtaining the total lenth in Earth
        """
        return self._tot_length_in_earth

    def get_running_length(self):
        """ Obtaining running length
        """
        return self._running_length

    def get_horizontal_angle(self):
        """ Obtaining horizontal angle
        """
        return self._horizontal_angle

    def rand_(self):
        """ Generates a random number between 0 and 1
        """
        return random.random()

    def test_earth_geometry(self):
        """ 
        Tests the position of the particle relative to Earth's geometry and predefined altitude limits.
        """
        if self._running_length >= self._tot_length_in_earth + ALTITUDE_MAX / math.sin(self._horizontal_angle):
            return eWhere.OUTOFCUTS
        elif self._running_length >= self._tot_length_in_earth + ALTITUDE_MIN / math.sin(self._horizontal_angle) and self._running_length < self._tot_length_in_earth:  # ALTITUDE_MIN is defined negative
            return eWhere.INEARTH_AND_WITHINCUTS
        elif self._running_length >= self._tot_length_in_earth and self._running_length < self._tot_length_in_earth + ALTITUDE_MAX / math.sin(self._horizontal_angle):
            return eWhere.OUTEARTH_AND_WITHINCUTS
        else:
            return eWhere.INEARTH

    def earth_density(self, R):
        """ Calculation of Earth's density based on the Preliminary Reference Earth Model
        (from Dziewonski and Anderson 1981)
        """
        x = R / kREarth
        
        # Atmosphere density
        if R > kREarth:
            dens = 0.00224
            A = 0.20*16+0.80*14  # 1 g/mol for nucleons 
        
        # Inside the Earth
        if R < 1221.500:   # core
            dens = 13.0885 - 8.8381 * x * x
            A=0.90*52+0.10*58
        elif R < 3480.000:
            dens = 12.5815 - x * (1.2638 + x * (3.6426 + x * 5.5281))
            A=0.90*52+0.10*58
        elif R < 5701.000:   # mantle
            dens = 7.9565 - x * (6.4761 - x * (5.5283 - x * 3.0807))
            A=0.44*16+0.22*22+0.22*28+0.06*52+0.03*40+0.03*26
        elif R < 5771.000:
            dens = 5.3197 - 1.4836 * x
            A=0.44*16+0.22*22+0.22*28+0.06*52+0.03*40+0.03*26
        elif R < 5971.000:
            dens = 11.2494 - 8.0298 * x
            A=0.44*16+0.22*22+0.22*28+0.06*52+0.03*40+0.03*26
        elif R < 6151.000:
            dens = 7.1089 - 3.8045 * x
            A=0.44*16+0.22*22+0.22*28+0.06*52+0.03*40+0.03*26
        elif R < 6346.600:    # crust
            dens = 2.691 + 0.6924 * x
            A=0.46*16+0.28*28.+0.083*26.+0.056*52.+0.042*40.+0.025*22.+0.024*24.+0.02*28.
        elif R < 6356.000:
            dens = 2.9
            A=0.46*16+0.28*28.+0.083*26.+0.056*52.+0.042*40.+0.025*22.+0.024*24.+0.02*28.
        elif R <= kREarth:
            dens = 2.6
            A=0.46*16+0.28*28.+0.083*26.+0.056*52.+0.042*40.+0.025*22.+0.024*24.+0.02*28.

        return dens * kNAvogadro / A  # nucleons/cm^3
    
    def PropagateParticle(self, input_particle):
        # Displaying the type of particle currently being propagated
        print("... propagating {} ({}, Ecurr {})".format(input_particle.GetName(), input_particle.GetMassNum(), input_particle.GetEprod()))

        # Case where the particle is emitted at the Earth's surface
        if input_particle.GetBranch() == 0:
            # Resetting the particle's characteristics
            input_particle.SetMassNum(0)
            input_particle.SetCharge(0)
            input_particle.SetZprod(0)

            # Generating random coordinates at the Earth's surface
            F = random.random()
            r_min = kREarth * math.cos(alpha_max)
            r_max = kREarth * math.cos(alpha_min)
            K = 2 / (r_max ** 2 - r_min ** 2)
            r = math.sqrt(2 * F / K + r_min ** 2)
            z = -math.sqrt(kREarth ** 2 - r ** 2) 
            phi = 2 * math.pi * random.random()
            x = r * math.cos(phi)
            y = r * math.sin(phi)
            print("rnd = {}; r = {}; z = {} ; phi = {}; x = {}; y = {}".format(F, r, z, phi, x, y))

            # Defining new coordinates for the particle and the horizontal angle
            self.set_x(x)
            self.set_y(y)
            self.set_z(z)
            self.set_tot_length_in_earth(-2 * z)
            self.set_horizontal_angle(math.acos(r / kREarth))

            # Resetting the traversal lenght
            self._running_length = 0

        # Checking propagation conditions for certain special particles
        if input_particle.GetType() == eParticleType.Tau and ((self.test_earth_geometry() == eWhere.INEARTH_AND_WITHINCUTS or self.test_earth_geometry() == eWhere.OUTEARTH_AND_WITHINCUTS)):
            input_particle.SetType(eParticleType.Null)
        if input_particle.GetType() == eParticleType.Neutrino_tau and (self.test_earth_geometry() == eWhere.OUTOFCUTS):
            input_particle.SetType(eParticleType.Null)
        if input_particle.GetType() == eParticleType.Antineutrino_tau and (self.test_earth_geometry() == eWhere.OUTOFCUTS):
            input_particle.SetType(eParticleType.Null)
        if input_particle.GetType() == eParticleType.Neutrino_mu and (self.test_earth_geometry() == eWhere.OUTOFCUTS):
            input_particle.SetType(eParticleType.Null)
        if input_particle.GetType() == eParticleType.Antineutrino_mu and (self.test_earth_geometry() == eWhere.OUTOFCUTS):
            input_particle.SetType(eParticleType.Null)
        if input_particle.GetType() == eParticleType.Neutrino_e and (self.test_earth_geometry() == eWhere.OUTOFCUTS):
            input_particle.SetType(eParticleType.Null)
        if input_particle.GetType() == eParticleType.Antineutrino_e and (self.test_earth_geometry() == eWhere.OUTOFCUTS):
            input_particle.SetType(eParticleType.Null)
        

        # Initializing output list of particles
        output = []

        # Propagation based on the particle type
        if input_particle.GetType() == eParticleType.Pion:
            output = self.PropagatePion(input_particle)
        elif input_particle.GetType() == eParticleType.Tau:
            output = self.PropagateTau(input_particle)
        elif input_particle.GetType() == eParticleType.Muon:
            output = self.PropagateMuon(input_particle)
        elif input_particle.GetType() in [eParticleType.Neutrino_e, eParticleType.Antineutrino_e, eParticleType.Neutrino_mu, eParticleType.Antineutrino_mu]:
            output = self.PropagateNeutrino(input_particle)
            
        # Propagation of a Neutrino Tau   # TEST
        elif input_particle.GetType() in [eParticleType.Neutrino_tau, eParticleType.Antineutrino_tau]:
            output = self.PropagateNeutrino(input_particle)

        # Returns the list of propagated particles
        return output

    def PropagateNeutrino(self, input_particle):
        output = []

        # Converting of the energy of the particle in GeV
        E = input_particle.GetEprod() / 1e9
        A = input_particle.GetMassNum()
        Z = input_particle.GetCharge()
        z = input_particle.GetZint()

        # Updating characteristics of the input particle after the propagation process
        input_particle = self.GetProcess(input_particle, 0)
        E = input_particle.GetEprod()
        nBr = input_particle.GetBranch() + 1

        # Adding the generated particle in the list
        output.append(Particle(Energy=0, branch=nBr, type=input_particle.GetType(), charge=input_particle.GetCharge(), Eprod=E, Zprod=z)) # SOUCI

        # Updating characteristics of the input particle
        input_particle.SetEint(E)
        input_particle.SetZprod(self._running_length)

        # Returns list of generated particles during propagation
        return output

    def PropagateTau(self, input_particle):
        # Initializing output list
        output = []

        # Converting energy in GeV
        E = input_particle.GetEprod() / 1e9
        A = input_particle.GetMassNum()
        Z = input_particle.GetCharge()
        z = input_particle.GetZint()

        # Updating characteristics of the input particle after propagation
        input_particle = self.GetProcess(input_particle, 0)
        E = input_particle.GetEprod()
        print("Energy in prop tau =", E)  # Displaying energy after propagation

        # Branch of the particle
        nBr = input_particle.GetBranch() + 1

        # Adding generated particle to output list
        output.append(Particle(Energy=0, charge=input_particle.GetCharge(), branch=nBr, Eprod=E, Zprod=z, type=input_particle.GetType()))

        # Updating characteristics of input particle
        input_particle.SetEint(E)
        input_particle.SetZprod(self._running_length)

        # Returning output list of particles generated during propagation
        return output
    
    """def PropagateNucleus(self, input_particle):
        # Characteristics of input particle
        Acurr = input_particle.GetMassNum()
        Zcurr = input_particle.GetCharge()

        # Verifying if particle can decay via beta decay
        if Zcurr != GetBetaDecayStableZ(Acurr):
            if fBdecay:
                return BetaDecay(input_particle)
            else:
                input_particle.SetCharge(GetBetaDecayStableZ(Acurr))
                if input_particle.GetType() == eParticleType.Neutron:
                 input_particle.SetType(eParticleType.Proton)

        # Initializing output list
        output = []

        # Vérification du type de la particule et du mode de propagation
        if input_particle.GetType() == eParticleType.Proton and fStoch <= 0:
            DetermProton(input_particle)
            return output

        # Détermination de la branche de la particule
        nBr = input_particle.GetBranch() + 1

        # Détermination du processus de propagation
        proc = GetProcess(input_particle, electrons if fPairProd else 0)

        # Récupération des caractéristiques finales de la particule
        zfin = input_particle.GetZint()
        Gfin = input_particle.GetGint()

        # Initialisation des variables pour les protons et les neutrons
        Afin, Zfin, intType, protons, neutrons = 0, 0, 0, 0, 0

        # Sélection du processus de propagation
        if proc == eProc.Earth:
            return electrons
        elif proc == eAlphaProd:
            Afin = Acurr - 4
            Zfin = Zcurr - 2
            intType = 4
            protons, neutrons = 0, 0
        elif proc in [eNucleonProd, eHadron, eHadron_IR]:
            # Gestion des processus de production de nucléons et de hadrons
            # (certaines parties sont commentées car des fonctions spécifiques ne sont pas définies ici)
            pass
        elif proc == eDisi:
            # Gestion du processus de dissociation
            # (certaines parties sont commentées car des fonctions spécifiques ne sont pas définies ici)
            pass

        # Gestion de la production de pions
        if proc == eHadron or proc == eHadron_IR:
            # (certaines parties sont commentées car des fonctions spécifiques ne sont pas définies ici)
            pass

        # Gestion des cas particuliers
        if Afin == 1:
            if Zfin == 1:
                protons += 1
            else:
                neutrons += 1
            Afin = 0
        if not fBdecay:
            if Afin > 1:
                Zfin = GetBetaDecayStableZ(Afin)
                protons = intType
                neutrons = 0

        # Ajout des particules de sortie à la liste de sortie
        # (certaines parties sont commentées car des fonctions spécifiques ne sont pas définies ici)
        pass

        # Mise à jour du nombre de particules interagissantes
        input_particle.SetIntMult(intType)

        # Retourne la liste des particules générées lors de la propagation
        return output"""
    
    def GetProcess(self, input_particle, electrons=None):
        A = input_particle.GetMassNum()  # Masse atomique
        Z = input_particle.GetCharge()   # Charge
        E = input_particle.GetEprod() / 1.e9  # Énergie en GeV

        # Calcul du taux de processus
        rate = self.TotalRate(A, Z, 0, E, 0, None) * self.dtdz_com

        # Génération d'un nombre aléatoire selon une distribution logarithmique uniforme
        logu = math.log(random.random())

        # Initialisation des variables
        z_old, Gam_old, rate_old, logpSurv_old = 0, 0, 0, 0
        proc = eProcType.Earth   # Processus de propagation TEST
        rd, step = 0.01, 0.01

        if input_particle.GetType() in [eParticleType.Neutrino_e, eParticleType.Neutrino_mu, eParticleType.Neutrino_tau]:
            # Traitement des neutrinos
            while self.test_earth_geometry() != eWhere.OUTOFCUTS:
                self.build_process(input_particle)
                self.prob_process(step)
                rd = self.rand_()
                if self._p_cc > rd:
                    self.int_cc()
                    self._n_bang += 1
                    print("neutrino CC :","z =", self._z,"; E =", self._E,"; running length =", self._running_length)
                    break
                elif self._p_cc + self._p_nc > rd:
                    self.int_nc()
                else:
                    self.increment_running_length(step)
                    self._z += step
            if self.test_earth_geometry() == eWhere.OUTOFCUTS:
                input_particle.SetType(eParticleType.Null)
            else:
                input_particle.SetEint(self._E * 1.e9)
                input_particle.SetEprod(self._E * 1.e9)
                input_particle.SetCharge(1)
                input_particle.SetIntMult(0)
                input_particle.SetType(eParticleType.Tau)
                if (self.test_earth_geometry() == eWhere.INEARTH_AND_WITHINCUTS or self.test_earth_geometry() == eWhere.OUTEARTH_AND_WITHINCUTS):
                    input_particle.SetIntMult(-1)
        elif input_particle.GetType() == eParticleType.Tau:
            # Traitement des taus
            while self.test_earth_geometry() == eWhere.INEARTH:
                self.build_process(input_particle)
                self.prob_process(step)
                rd = self.rand_()
                if self._p_decay > rd:
                    self.decay()
                    print("tau decay","z =", self._z, "; running length =", self._running_length)
                    self._n_bang += 1
                    break
                elif self._p_decay + self._p_cc > rd:
                    self.int_cc()
                    self._n_bang += 1
                    print("tau CC", self._z)
                    break
                elif self._p_decay + self._p_cc + self._p_nc > rd:  
                    self.int_nc()
                else:   
                    self.increment_running_length(step)
                    self._z += step
                    
                    XSI = 4 / (self._E * 1E-9) ** 0.15 
                    self._E = self._E * math.exp(-step / XSI)
            if input_particle.GetType() == eParticleType.Tau and ((self.test_earth_geometry() == eWhere.INEARTH_AND_WITHINCUTS or self.test_earth_geometry() == eWhere.OUTEARTH_AND_WITHINCUTS)):
                input_particle.SetEint(self._E * 1.e9)
                input_particle.SetEprod(self._E * 1.e9)
                input_particle.SetIntMult(-1)
                input_particle.SetType(eParticleType.Tau)
                input_particle.SetCharge(1)
            else:
                input_particle.SetEint(self._E * 1.e9)
                input_particle.SetEprod(self._E * 1.e9)
                input_particle.SetIntMult(0)
                input_particle.SetType(eParticleType.Neutrino_tau)
                input_particle.SetCharge(0)
        return input_particle

    def build_process(self, input_particle):
        # Fonction pour calculer les processus de propagation en fonction du type de particule
        if input_particle.GetType() == eParticleType.Neutrino_e or input_particle.GetType() == eParticleType.Neutrino_mu or input_particle.GetType() == eParticleType.Neutrino_tau :
            # Calcul des sections efficaces pour la diffusion des neutrinos
            R = math.sqrt(self._x * self._x + self._y * self._y + self._z * self._z)
            dens = self.earth_density(R)  # Densité de la Terre en cm^-3
            eps = math.log10(self._E)
            y = math.log(eps + 1.826)
            sig_cc = pow(10., -17.31 - 6.406 * y + 1.431 * y * y - 17.91 / y)  # Section efficace de diffusion des CC en cm^2
            sig_nc = pow(10., -17.31 - 6.448 * y + 1.431 * y * y - 18.61 / y)  # Section efficace de diffusion des NC en cm^2
            # Calcul de la longueur libre moyenne pour les CC et les NC en km
            self._mfp_cc = pow(dens * sig_cc, -1.) * 1E-5 # in km
            self._mfp_nc = pow(dens * sig_nc, -1.) * 1E-5 # in km
            self._mfp_dec = 49 * self._E * 1E-9 # Inutilisé pour les neutrinos
        else:
            # Calcul de la longueur libre moyenne pour la décroissance des autres particules en km
            self._mfp_dec = 49 * self._E * 1E-9 # 49 km for 1.0 EeV
            R = math.sqrt(self._x * self._x + self._y * self._y + self._z * self._z)
            dens = self.earth_density(R)  # Densité de la Terre en cm^-3
            if self._E <= 0:
                print(self._E)
            eps = math.log10(self._E)
            y = math.log(eps + 1.826)
            sig_cc = pow(10., -17.31 - 6.406 * y + 1.431 * y * y - 17.91 / y)  # Section efficace de diffusion des CC en cm^2
            sig_nc = pow(10., -17.31 - 6.448 * y + 1.431 * y * y - 18.61 / y)  # Section efficace de diffusion des NC en cm^2
            # Calcul de la longueur libre moyenne pour les CC et les NC en km
            self._mfp_cc = pow(dens * sig_cc, -1.) * 1E-5
            self._mfp_nc = pow(dens * sig_nc, -1.) * 1E-5

    def prob_process(self, step):
        # Calcul des probabilités de processus pour un pas de propagation donné
        self._p_decay = step / self._mfp_dec  # Probabilité de décroissance
        self._p_cc = step / self._mfp_cc  # Probabilité de collision CC
        self._p_nc = step / self._mfp_nc  # Probabilité de collision NC

    def TotalRate(self, A, Z, j, E, z, type=None):
        # Calcul du taux total de processus pour une particule donnée
        self._E = E
        self._x, self._y, self._z = self.get_x(), self.get_y(), self.get_z()

        R = math.sqrt(self._x * self._x + self._y * self._y + self._z * self._z)
        dens = self.earth_density(R)  # Densité de la Terre en cm^-3
        eps = math.log10(self._E)
        y = math.log(eps + 1.826)
        sig_cc = pow(10., -17.31 - 6.406 * y + 1.431 * y * y - 17.91 / y)  # Section efficace de diffusion des CC en cm^2
        sig_nc = pow(10., -17.31 - 6.448 * y + 1.431 * y * y - 18.61 / y)  # Section efficace de diffusion des NC en cm^2
        # Calcul de la longueur libre moyenne pour les CC et les NC en km
        self._mfp_cc = pow(dens * sig_cc, -1.) * 1E-5
        self._mfp_nc = pow(dens * sig_nc, -1.) * 1E-5
        print("mfp_cc =", self._mfp_cc, "; mfp_nc =", self._mfp_nc, "; x =", self._x, "; y =", self._y,"; z =", self._z,"; E =",self._E)
        return self._mfp_cc + self._mfp_nc

    def int_cc(self):
        # Interaction par collision CC
        self._E = self._E * (1.002 - pow((.235 + .7654 * self.rand_()), 4.292))

    def int_nc(self):
        # Interaction par collision NC
        self._E = self._E * (1.002 - pow((.235 + .7654 * self.rand_()), 4.292))

    def PropagatePion(self, input_particle):
        # Fonction pour propager un pion
        output = []  # Initialisation de la liste des particules de sortie
        charge = input_particle.GetCharge()  # Récupération de la charge du pion
        E = input_particle.GetEprod()  # Récupération de l'énergie du pion
        nBr = input_particle.GetBranch() + 1  # Récupération du nombre de branches
        r = random.random()  # Génération d'un nombre aléatoire entre 0 et 1

        if charge == 0:
            # Si le pion est neutre
            output.append(Particle(r * E, nBr, eParticleType.Photon))
            output.append(Particle((1. - r) * E, nBr, eParticleType.Photon))
            print("photon ", r * E, "; photon ", (1. - r) * E)
        else:
            # Si le pion est chargé
            E_nu = (1. - (mmu * mmu) / (mpi * mpi)) * E * r
            output.append(Particle(E - E_nu, nBr, eParticleType.Muon))
            output.append(Particle(E_nu, nBr, eParticleType.Neutrino_mu))
            print("muon ", E - E_nu, "; neutrino ", E_nu)

        input_particle.SetIntMult(1002)  # Définition du multiplicateur d'interaction
        return output  # Retourne la liste des particules de sortie

    def PropagateMuon(self, input_particle):
        # Fonction pour propager un muon
        charge = input_particle.GetCharge()  # Récupération de la charge du muon
        E = input_particle.GetEprod()  # Récupération de l'énergie du muon
        nBr = input_particle.GetBranch() #+ 1  # Récupération du nombre de branches
        gamma = E / mmu  # Facteur de Lorentz

        # Initialisation des quantités dans le référentiel du centre de masse
        E1, E2, E3, p1, p2, p3 = 0, 0, 0, 0, 0, 0

        # Calcul des énergies et impulsions des particules dans le référentiel du centre de masse
        Enumax = (mmu - me * me / mmu) / 2.
        while True:
            E1 = random.uniform(0, Enumax-me)
            E2 = random.uniform(0, Enumax-me)
            if E1 + E2 < Enumax:
                E1 = Enumax - E1
                E2 = Enumax - E2 
            E3 = mmu - E1 - E2
            p1 = E1
            p2 = E2
            p3 = math.sqrt(E3 * E3 - me * me)
            if E3 >= me and p3 <= p1 + p2 and p1 <= p2 + p3 and p2 <= p3 + p1:
                break

        # Calcul des angles entre les neutrinos et avec la ligne de visée
        c12 = (p3 * p3 - p1 * p1 - p2 * p2) / (2. * p1 * p2)
        s12 = math.sqrt(1. - c12 * c12)
        c1 = random.uniform(-1., 1.)
        s1 = math.sqrt(1. - c1 * c1)
        ca = math.cos(random.uniform(0, 2 * math.pi))
        c2 = c12 * c1 - s12 * s1 * ca

        # Calcul des énergies des neutrinos et de l'électron dans le référentiel du laboratoire
        Enu1 = gamma * (E1 + p1 * c1)
        Enu2 = gamma * (E2 + p2 * c2)
        Ee = E - Enu1 - Enu2

        input_particle.SetIntMult(1003)  # Définition du multiplicateur d'interaction

        # Création et ajout des particules de sortie à la liste
        output = [
            Particle(Enu1, nBr, eParticleType.Neutrino_mu),
            Particle(Enu2, nBr, eParticleType.Neutrino_e),
            Particle(Ee, nBr, eParticleType.Electron)
        ] 
        print("electron ", Ee, "; neutrino ", Enu1, "; neutrino ", Enu2)
        return output  # Retourne la liste des particules de sortie

    def decay(self):
        # Fonction pour simuler la désintégration du tauon
        global read, _tauola_fraction

        if not read:
            # Lecture du fichier contenant les fractions de désintégration du tauon
            iname = "ener_nutau.taudecay"
            with open(iname, 'r') as file:
                _tauola_fraction = [float(line.strip()) for line in file.readlines()]
            read = True

        # Sélection aléatoire d'une fraction de désintégration
        ind = int(self.rand_() * len(_tauola_fraction))
        self._E *= _tauola_fraction[ind]  # Mise à jour de l'énergie du tauon
