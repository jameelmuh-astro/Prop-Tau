#!/usr/local/bin/python
# -*- coding: utf-8 -*-
"""
File : Particle.py 
For tests before implementing in Prop_Tau
"""
import math
from enum import Enum

# Définition des types de particules
class eParticleType(Enum):
    Null = 0
    Proton = 1
    Neutron = 2
    Nucleus = 3
    Pion = 4
    Muon = 5
    Photon = 6
    Neutrino_e = 7
    Antineutrino_e = -7
    Neutrino_mu = 8
    Antineutrino_mu = -8
    Electron = 9
    Neutrino_tau = 10
    Antineutrino_tau = -10
    Tau = 11  # Lepton tauon

class Particle:
    def __init__(self, Energy=None, branch=None, type=eParticleType.Null, massnum=-1, charge=-1, Eprod=0, Zprod=0):
        """ Initialisation
        Args:
            Energy (float): Energy of the particle
            branch (int): Particle branch
            type (eParticleType): Particle type
            massnum (int): Mass number 
            charge (int): Particle charge
            Eprod (float): Product energy
            Zprod (float): Product charge
        """
        self.fmassnum = massnum
        self.fcharge = charge
        self.fIntMult = 0      # Multiplicity if there are multiple particles of same type
        self.fEprod = Eprod
        self.fZprod = Zprod
        self.fEint = -1.       # By default
        self.fZint = -1.       # By default
        self.fType = type
        
        if Energy is not None and branch is not None and type is not None:
            self.fEnergy = Energy
            self.fbranch = branch

    def SetMassNum(self, val):
        self.fmassnum = val
    
    def SetCharge(self, val):
        self.fcharge = val
    
    def SetBranch(self, val):
        self.fbranch = val
    
    def SetIntMult(self, val):
        self.fIntMult = val
    
    def SetType(self, val):
        self.fType = val
    
    def SetEprod(self, val):
        self.fEprod = val
    
    def SetZprod(self, val):
        self.fZprod = val
    
    def SetEint(self, val):
        self.fEint = val
    
    def SetZint(self, val):
        self.fZint = val

    def GetNeutrons(self):
        return self.fmassnum - self.fcharge
    
    def GetMassNum(self):
        return self.fmassnum
    
    def GetCharge(self):
        return self.fcharge
    
    def GetBranch(self):
        return self.fbranch
    
    def GetEprod(self):  
        if self.fEprod == 0:
            return self.fEnergy
        else:
            return self.fEprod
    
    def GetZprod(self):
        return self.fZprod
    
    def GetEint(self):
        return self.fEint
    
    def GetZint(self):
        return self.fZint
    
    def GetIntMult(self):
        return self.fIntMult
    
    def GetType(self):
        return self.fType
    
    def GetFlavor(self):
        """ To get the flavor  -> int
        """
        switch = {
            eParticleType.Neutrino_e: 1,
            eParticleType.Antineutrino_e: -1,
            eParticleType.Neutrino_mu: 2,
            eParticleType.Antineutrino_mu: -2,
            eParticleType.Muon: 2 * -self.fcharge,
            eParticleType.Electron: 1 * -self.fcharge,
            eParticleType.Neutrino_tau: 3,
            eParticleType.Antineutrino_tau: -3,
            eParticleType.Tau: 3 * -self.fcharge
            
        }
        return switch.get(self.fType, 0)

    def GetName(self):
        """ To get the name of the particle based on the type -> str
        """
        switch = {
            eParticleType.Null: "null",
            eParticleType.Proton: "proton",
            eParticleType.Neutron: "neutron",
            eParticleType.Nucleus: "nucleus",
            eParticleType.Pion: "pion",
            eParticleType.Muon: "muon",
            eParticleType.Photon: "photon",
            eParticleType.Neutrino_e: "electron neutrino",
            eParticleType.Antineutrino_e: "electron antineutrino",
            eParticleType.Neutrino_mu: "muon neutrino",
            eParticleType.Antineutrino_mu: "muon antineutrino",
            eParticleType.Electron: "positron" if self.fcharge > 0 else "electron",
            eParticleType.Neutrino_tau: "tau neutrino",
            eParticleType.Antineutrino_tau: "tau antineutrino",
            eParticleType.Tau: "tauon" if self.fcharge < 0 else "antitauon"
        }
        return switch.get(self.fType, "null")

    def GetID(self):
        """
        Obtient l'ID de la particule en fonction de son type.
        
        Returns:
            int: ID de la particule.
        """
        switch = {
            eParticleType.Null: 0,
            eParticleType.Proton: 2212,
            eParticleType.Neutron: 2112,
            eParticleType.Nucleus: 1000000000 + 10000 * self.fcharge + 10 * self.fmassnum,
            eParticleType.Pion: 211 * self.fcharge if self.fcharge != 0 else 111,
            eParticleType.Muon: 13 * -self.fcharge,
            eParticleType.Photon: 22,
            eParticleType.Neutrino_e: 12,
            eParticleType.Antineutrino_e: -12,
            eParticleType.Neutrino_mu: 14,
            eParticleType.Antineutrino_mu: -14,
            eParticleType.Electron: 11 * -self.fcharge,
            eParticleType.Neutrino_tau: 16,
            eParticleType.Antineutrino_tau: -16,
            eParticleType.Tau: 15 * -self.fcharge
        }
        return switch.get(self.fType, 0)
