#!/usr/local/bin/python
# -*- coding: utf-8 -*-
"""
File : NucModel.py 
"""

import sys
import random
import math
import numpy as np
from Particle import eParticleType,Particle
import Propagate
from Constants import *
import Output
from Convert import *

# Constantes
epsIBL_min = 0
epsIBL_max = 0

# Variables globales
fZ = []  # Liste des charges des éléments
fA = []  # Liste des numéros atomiques des éléments
k_com = 0  # Variable pour la charge
Gam_com = 0  # Variable pour la vitesse
z_com = 0  # Variable pour la coordonnée
s_com = 0  # Variable pour la distance
dtdz_com = 0  # Variable pour le pas de temps

class NucModel:
    def __init__(self, model, LossType, Stoch, Bdecay, PairProd):
        global epsIBL_min, epsIBL_max
        self.fnelements = 0  # Nombre d'éléments
        self.fType = model  # Type de modèle
        self.fStoch = Stoch  # Stochastique
        self.fBdecay = Bdecay  # Décroissance β
        self.fPairProd = PairProd  # Production de paires
        global fLosses
        fLosses = LossType  # Type de pertes
        if fLosses == 1:
            epsIBL_min = 3.3e-3  # Énergie minimale IR/V/UV (eV)
            epsIBL_max = 1.e1  # Énergie maximale IR/V/UV (eV)
        elif fLosses == 2:
            epsIBL_min = 2.e-3  # Énergie minimale IR/V/UV (eV)
            epsIBL_max = 1.  # Énergie maximale IR/V/UV (eV)
        elif fLosses == 3:
            epsIBL_min = 1.26e-3
            epsIBL_max = 11.5
        elif fLosses in [4, 5, 6]:
            epsIBL_min = pow(10, logeps_EBL[0])  # Énergie minimale IR/V/UV (eV)
            epsIBL_max = pow(10, logeps_EBL[neps_EBL - 1])  # Énergie maximale IR/V/UV (eV)
        elif fLosses == 7:
            epsIBL_min = pow(10, logeps_EBLg[0])  # Énergie minimale IR/V/UV (eV)
            epsIBL_max = pow(10, logeps_EBLg[neps_EBLg - 1])  # Énergie maximale IR/V/UV (eV)

        if model == eSalomon:
            self.fnelements = 51  # Nombre d'éléments
            self.fe_1 = e_1  # Énergie minimale de photodésintégration (MeV)
            self.fe_max = e_max  # Énergie maximale de photodésintégration (MeV)
            for i in range(self.fnelements):
                fZ.append(ZZ[i])  # Ajout des charges
                fA.append(AA[i])  # Ajout des numéros atomiques
        elif model in [eGauss, eBreitWigner]:
            self.fnelements, self.fe_1, self.fe_max = map(float, input().split())  # Entrée des données
            for i in range(self.fnelements):
                ZZ, AA, eth_1, eth_2, e0_1, csi_1, delta_1, e0_2, csi_2, delta_2, zita = map(float, input().split())
                fZ.append(ZZ)  # Ajout des charges
                fA.append(AA)  # Ajout des numéros atomiques
        elif model in [eBreitWigner2, eGauss2]:
            self.fnelements, self.fe_1, self.fe_max = map(float, input().split())  # Entrée des données
            for i in range(self.fnelements):
                ZZ, AA, nuclt, nuclh1, nuclx1, nuclw1, nuclh2, nuclx2, nuclw2, nuclc, alpht, alphh1, alphx1, alphw1, alphh2, alphx2, alphw2, alphc = map(float, input().split())
                fZ.append(ZZ)  # Ajout des charges
                fA.append(AA)  # Ajout des numéros atomiques

    def __del__(self):
        pass  # Destructeur de la classe NucModel

    # Méthode pour valider un objet Particle
    # def Validate(self, input):
    #     pass

    # Méthode pour tester la charge primaire
    def TestPrimary(self, mass):
        if self.fType == eOnlyProton and mass == 1:  # Si le type est 'eOnlyProton' et la masse est 1
            return 1  # Retourne 1 (tous les particules sont transformés en protons)

        found = self.binary(mass)  # Recherche binaire pour trouver la charge de l'élément
        return found

    # Méthode de recherche binaire pour trouver la charge de l'élément
    def binary(self, mass):
        first = 0
        last = self.fnelements

        while first <= last:
            mid = (first + last) // 2
            if mass > self.fA[mid]:
                last = mid - 1
            elif mass < self.fA[mid]:
                first = mid + 1
            else:
                return self.fZ[mid]

        return -1  # Retourne -1 si la charge n'est pas trouvée

    # Méthode pour obtenir une particule aléatoire
    def GetRandomParticle(self, E, z):
        index = random.randint(0, self.fnelements - 1)  # Sélection aléatoire d'un index
        A = self.fA[index]  # Numéro atomique
        Z = self.fZ[index]  # Charge

        # Retourne une particule avec les paramètres donnés
        return Particle(A, Z, 0, E, z, eProton if A == 1 else eNucleus)  # Si A == 1, alors c'est un proton, sinon c'est un noyau
    
    ######################
    def PropagateParticle(self, input):
        # Affichage des informations sur la particule en cours de propagation
        print("... propagating ", input.GetName(), " (", input.GetMassNum(), ", ", input.GetCharge(), ") ", " Ecurr ", input.GetEprod(), " from z ", input.GetZprod())

        output = []  # Liste des particules de sortie
        # Sélection de la méthode de propagation en fonction du type de particule
        if input.GetType() == eProton:
            if self.fStoch <= 0:
                self.DetermProton(input)  # Si la méthode est déterministe, appel à la fonction DetermProton
            else:
                output = self.PropagateNucleus(input)  # Sinon, propagation de la particule comme un noyau
        elif input.GetType() in [eNeutron, eNucleus]:
            output = self.PropagateNucleus(input)  # Propagation d'un neutron ou d'un noyau
        elif input.GetType() == ePion:
            output = self.PropagatePion(input)  # Propagation d'un pion
        elif input.GetType() == ePhoton:
            output = self.PropagatePhoton(input)  # Propagation d'un photon
        elif input.GetType() == eMuon:
            output = self.PropagateMuon(input)  # Propagation d'un muon
        elif input.GetType() in [eNeutrino_e, eAntineutrino_e, eNeutrino_mu, eAntineutrino_mu]:
            output = self.PropagateNeutrino(input)  # Propagation d'un neutrino
        elif input.GetType() == eElectron:
            output = self.PropagateElectron(input)  # Propagation d'un électron
        elif input.GetType() == eNull:
            pass  # Type de particule nul

        # Affichage des informations sur la particule après la propagation
        if input.GetIntMult() > 1000:
            print("........ ", input.GetName(), " decays into ", input.GetIntMult() - 1000, " particles at z=", input.GetZint())
        elif input.GetIntMult() == 1000:
            print("........ ", input.GetName(), " produced with E=", input.GetEprod(), " at z=", input.GetZprod())
        elif input.GetIntMult() > 0:
            assert len(output) > 0
            print("...... int. type ", input.GetIntMult(), " (", input.GetMassNum(), ", ", input.GetCharge(), ") -> (", output[0].GetMassNum(), ", ", output[0].GetCharge(), ") zfin=", input.GetZint(), " Efin=", input.GetEint())
        elif input.GetZint() <= 0.:
            input.SetZint(0)
            print("........... ", input.GetName(), " (", input.GetMassNum(), ", ", input.GetCharge(), ") reaches Earth with E = ", input.GetEint())

        return output

    # Méthode pour la propagation déterministe d'un proton
    def DetermProton(self, input):
        Gamdet = self.evolveN(input.GetZprod(), input.GetGprod(), 0.)
        input.SetZint(0.)
        input.SetGint(Gamdet)
        input.SetIntMult(0)
        
    # Méthode pour obtenir le nombre de protons d'un noyau stable par désintégration bêta
    def GetBetaDecayStableZ(self, A):
        # Tableau des nombres de protons pour les noyaux stables par désintégration bêta
        Zstable = [0, 1,1,2,2,2, 3,3,4,4,5, 5,6,6,7,7, 8,8,8,9,10, 10,10,11,12,12, 12,13,14,14,14, 15,16,16,16,17, 16,17,18,19,18, 19,20,20,20,21, 22,22,22,22,22, 23,24,24,24,25, 26,26,26,27,28]
        Z = self.binary(A)  # Recherche binaire du nombre de protons pour le noyau A
        if Z > 0:
            return Z
        else:
            return Zstable[A]  # Retourne le nombre de protons du noyau stable par désintégration bêta

    # Méthode pour la propagation d'un noyau
    def PropagateNucleus(self, input):
        Acurr = input.GetMassNum()  # Nombre de masse actuel
        Zcurr = input.GetCharge()   # Nombre de protons actuel

        # Vérification si le noyau est stable par désintégration bêta
        if Zcurr != self.GetBetaDecayStableZ(Acurr):
            if self.fBdecay:
                return self.BetaDecay(input)  # Si la désintégration bêta est activée, effectuer la désintégration
            else:
                input.SetCharge(self.GetBetaDecayStableZ(Acurr))  # Sinon, définir le nombre de protons sur celui du noyau stable
                if input.GetType() == eNeutron:
                    input.SetType(eProton)  # Si la particule est un neutron, la changer en proton

        output = []  # Liste des particules de sortie

        # Vérification si le noyau est un proton et si la méthode est déterministe
        if input.GetType() == eProton and self.fStoch <= 0:
            self.DetermProton(input)
            return output

        nBr = input.GetBranch() + 1  # Numéro de la branche pour la particule
        proc = None  # Type de processus de propagation
        electrons = []  # Liste des électrons produits

        # Détermination du processus de propagation en fonction du nombre de masse et de la possibilité de production de paires
        if Acurr == 5 or Acurr == 8:
            input.SetZint(input.GetZprod())
            input.SetGint(input.GetGprod())
            input.SetIntMult(1002)
            if Acurr == 8:
                proc = eAlphaProd
            if Acurr == 5:
                proc = eNucleonProd
        else:
            proc = self.GetProcess(input, electrons if self.fPairProd else None)

        zfin = input.GetZint()  # Position finale
        Gfin = input.GetGint()  # Énergie cinétique finale

        # Initialisation des variables pour les particules de sortie
        Afin, Zfin, intType, protons, neutrons = None, None, None, None, None

        # Traitement du processus de propagation
        if proc == eEarth:
            return electrons
        elif proc == eAlphaProd:
            Afin = Acurr - 4
            Zfin = Zcurr - 2
            intType = 4
            protons = 0
            neutrons = 0
        elif proc in [eNucleonProd, eHadron, eHadron_IR]:
            Afin = Acurr - 1
            if Acurr == 3 or (Acurr != 5 and random.uniform(0., Acurr) < Zcurr):
                protons = 1
                neutrons = 0
                Zfin = Zcurr - 1
            else:
                protons = 0
                neutrons = 1
                Zfin = Zcurr
            intType = 1
        elif proc == eDisi:
            intType = input.GetIntMult()
            Afin = Acurr - intType
            # Vérification si le noyau est instable et doit se fragmenter
            if Afin > 1 and self.GetJNuc(Afin) < 0:
                jnuc = -self.GetJNuc(Afin)
                Afin = self.fA[jnuc]
                intType = Acurr - Afin
            # Calcul de la composition finale en protons et neutrons
            if Acurr == 9 and Afin == 4:
                proc = eAlphaProd
                intType = 1
                Zfin = 2
                protons = 0
                neutrons = 1
            else:
                Zfin = Zcurr
                for i in range(intType):
                    if random.uniform(0., Acurr - i) < Zfin:
                        Zfin -= 1  # Perte d'un proton, sinon perte d'un neutron
                if Afin > 1:
                    if Zfin == 0:
                        Zfin += 1
                    if Zfin == Afin:
                        Zfin -= 1  # Éviter les noyaux comme 2He
                protons = Zcurr - Zfin
                neutrons = intType - protons

        Epion = 0  # Énergie du pion
        pioncharge = 0  # Charge du pion
        if proc in [eHadron, eHadron_IR]:
            # Calcul de l'énergie du pion
            Epion = self.Photopion(zfin, input.GetGint(), proc)
            # Décision aléatoire de la charge du pion
            if random.uniform(0., 1.) < 1. / 3.:
                if protons:
                    protons -= 1
                    neutrons += 1
                    pioncharge = +1
                    print("....... stacking pion+, z=", zfin, ", E=", Epion)
                else:
                    protons += 1
                    neutrons -= 1
                    pioncharge = -1
                    print("....... stacking pion-, z=", zfin, ", E=", Epion)
            else:
                print("....... stacking pion0, z=", zfin, ", E=", Epion)
            if proc == eHadron:
                print("\n")
            if proc == eHadron_IR:
                print(" (from IBL)\n")

        if Afin == 1:
            if Zfin == 1:
                protons += 1
            else:
                neutrons += 1
            Afin = 0

        # Traitement des cas spéciaux pour les noyaux instables
        if not self.fBdecay:
            if Afin > 1:
                Zfin = self.GetBetaDecayStableZ(Afin)
                protons = intType
                neutrons = 0

        # Ajout des particules de sortie à la liste output
        if Afin != 0:
            output.append(Particle(Afin, Zfin, nBr, Afin * mN * Gfin, zfin, eNucleus))
        else:
            if Acurr > 1:
                print("Nucleus fully fragmented\n")

        if proc == eAlphaProd:
            print("....... stacking 4He, z=", zfin, ", E=", 4 * mN * Gfin, '\n')
            output.append(Particle(4, 2, nBr, 4 * mN * Gfin, zfin, eNucleus))

        if protons or neutrons:
            print("....... stacking ", protons, " p + ", neutrons, " n, z=", zfin, ", E=", mN * Gfin - Epion, '\n')

        for i in range(protons):
            output.append(Particle(1, 1, nBr, mN * Gfin - Epion, zfin, eProton))
        for i in range(neutrons):
            output.append(Particle(1, 0, nBr, mN * Gfin - Epion, zfin, eNeutron))
        if proc in [eHadron, eHadron_IR]:
            output.append(Particle(0, pioncharge, nBr, Epion, zfin, ePion))
        output.extend(electrons)
        input.SetIntMult(intType)  # Définition du nombre de particules après interaction
        return output

    def GetJNuc(self, ANuc):
        """
        Méthode pour obtenir l'index de l'élément nucléaire dans fA correspondant à un nombre de masse donné.

        Arguments :
            ANuc (int) : Nombre de masse de l'élément nucléaire recherché.

        Returns :
            jNuc (int) : Index de l'élément nucléaire dans fA.
                        Si l'élément nucléaire est trouvé, jNuc est positif et représente l'index dans fA.
                        Sinon, jNuc est négatif et représente la position où l'élément nucléaire aurait été inséré dans fA.
                        Si ANuc est inférieur à 0 ou supérieur au premier élément de fA, la fonction retourne 0.
        """
        jNuc = 0  # Initialisation de l'index de l'élément nucléaire
        i = 0     # Initialisation de l'itérateur

        # Vérifie si ANuc est dans la plage valide ou égal à 1 (qui correspond à l'hydrogène)
        if ANuc > self.fA[0] or ANuc < 0:
            return 0
        if ANuc == 1:
            return -1

        # Parcours de fA pour trouver l'index correspondant à ANuc
        for i in range(self.fnelements):
            if self.fA[i] <= ANuc:
                break

        # Vérifie si ANuc correspond exactement à un élément dans fA
        if self.fA[i] == ANuc:
            jNuc = i  # L'élément nucléaire a été trouvé, jNuc correspond à son index dans fA
        else:
            jNuc = -i  # L'élément nucléaire n'a pas été trouvé, jNuc correspond à l'index d'insertion dans fA

    def GetProcess(self, input, electrons=None):
        """
        Méthode pour obtenir le type de processus pour une particule donnée.

        Arguments :
            input (Particle) : Instance de la classe Particle représentant la particule d'entrée.
            electrons (list) : Liste optionnelle pour stocker les électrons générés pendant le processus.

        Returns :
            proc (eProcType) : Type de processus pour la particule d'entrée.
        """
        # Obtention du nombre de masse et de la charge de la particule d'entrée
        A = input.GetMassNum()
        Z = input.GetCharge()

        # Obtention de l'index de l'élément nucléaire dans fA correspondant à A
        J = self.GetJNuc(A)

        # Obtention de la largeur de la désintégration et de la coordonnée z de la production
        Gam = input.GetGprod()
        z = input.GetZprod()

        # À ce stade, vous pouvez utiliser les valeurs obtenues pour déterminer le type de processus.
        # Cependant, cela dépend de la logique métier spécifique de votre application, qui n'est pas fournie dans le code fourni.

        # La valeur de retour proc doit être un élément de l'énumération eProcType, représentant le type de processus.
        # Vous pouvez implémenter cette logique en fonction des valeurs de A, Z, J, Gam et z.

        # Par exemple :
        # Si A == 1 et Z == 1, cela pourrait indiquer un processus de production d'hydrogène.
        # Si J est négatif, cela pourrait indiquer un processus de désintégration.

        # Si des électrons sont générés pendant le processus, vous pouvez les stocker dans la liste des électrons.
        if electrons is not None:
            # Ajouter les électrons générés à la liste des électrons, si nécessaire
            pass

        # À la fin, retournez le type de processus déterminé
        return proc
        return jNuc

    def NEW_method(self, input, electrons=None):
        """
        Méthode pour effectuer une nouvelle opération.

        Arguments :
            input (Particle) : Instance de la classe Particle représentant la particule d'entrée.
            electrons (list) : Liste optionnelle pour stocker les électrons générés pendant le processus.

        Returns :
            proc (eProcType) : Type de processus pour la particule d'entrée.
        """
        # Obtention des valeurs B et Lcoh
        B = self.fB[0]
        Lcoh = self.fB[1]

        # Calcul du nombre de pas en fonction de z
        nstep = 200 * z + 20

        # Définition des valeurs initiales et finales de z
        zOri = z
        zEnd = 0

        # Calcul du pas de temps comobile
        dtdz_com = 1 / (H0 * (1 + z) * sqrt(om_m * pow(1 + z, 3) + om_l))

        # Calcul du taux total de réaction
        rate = TotalRate(A, Z, J, Gam, z, None) * dtdz_com

        # Initialisation des variables pour le calcul de la survie logarithmique
        logpSurv = 0
        thetasumq = 0
        Ttot = 0
        rad = 180 / 3.14159
        pi = 3.14159
        delta_fin = 0

        # Initialisation des variables pour la boucle
        logu = log(gRandom.Rndm())
        z_old, Gam_old, rate_old, logpSurv_old = 0, 0, 0, 0
        proc = eEarth

        # Boucle principale
        for j in range(1, nstep + 1):
            z_old = z
            Gam_old = Gam
            rate_old = rate
            logpSurv_old = logpSurv
            t = float(j) / nstep
            z = (1 - t) * zOri + t * zEnd
            if j == nstep:
                z = 0  # Correction des erreurs d'arrondi

            # Calcul du pas de temps comobile
            dtdz_com = 1 / (H0 * (1 + z) * sqrt(om_m * pow(1 + z, 3) + om_l))

            # Calcul de theta et th
            theta = GetTheta(z, z_old, A, Z, Gam_old, B, Lcoh)
            th = gRandom.Uniform(0., theta)
            thetasumq += theta * rad * theta * rad

            # Calcul de dtdz
            dtdz = Getdtdzeff(z, z_old, A, Z, Gam_old, B, Lcoh, th)
            dtdz_com = dtdz

            # Mise à jour de Ttot
            Ttot -= dtdz * (1 + z) * (z - z_old)

            # Évolution de Gam
            Gam = evolveA0(A, z_old, Gam_old, z)

            # Calcul du taux total de réaction
            rate = TotalRate(A, Z, J, Gam, z, None) * dtdz_com

            # Calcul de logpSurv
            logpSurv += (z - z_old) * (rate + rate_old) / 2.

            # Vérification de la survie
            if logpSurv <= logu:
                z -= (logu - logpSurv) / (logpSurv_old - logpSurv) * (z - z_old)
                Gam = evolveA0(A, z_old, Gam_old, z)
                TotalRate(A, Z, J, Gam, z, proc)

            # Génération des électrons
            if electrons:
                K = 2 * me / (A * mN + 2 * me)
                npairs_avg = log(Gam_old / (1 + z_old) / Gam * (1 + z)) / K
                npairs = gRandom.Poisson(npairs_avg)
                nBr = input.GetBranch() + 1
                for i in range(npairs):
                    Ee = K * A * exp(gRandom.Uniform(log(Gam), log(Gam_old))) * mN / 2
                    zprod = gRandom.Uniform(z, z_old)
                    electrons.append(Particle(0, -1, nBr, Ee, zprod, eElectron))
                    electrons.append(Particle(0, +1, nBr, Ee, zprod, eElectron))

            # Vérification du type de processus
            if proc != eEarth:
                break

        # Calcul de delta_fin
        delta_fin = sqrt(thetasumq)

        # Mise à jour des propriétés de la particule
        input.SetTtot(Ttot)
        input.Setdelta_fin(delta_fin)
        # Détermination du nombre de disintégrations
        ndis = 0
        if proc == eEarth:
            assert z == 0
            ndis = 0
        elif proc in [eHadron, eHadron_IR, eNucleonProd]:
            assert z > 0
            ndis = 1
        elif proc == eAlphaProd:
            assert z > 0
            ndis = 4
        elif proc == eDisi:
            assert z > 0
            bdis = SecInY * SpeedOfLight * m2cm / (Mpc2cm * losseAdisi(A, z, Gam, J, ndis))
            # case eDecay:
            #     ndis = 1003
            #     break

        # Mise à jour des propriétés de la particule
        input.SetZint(z)
        input.SetGint(Gam)
        input.SetIntMult(ndis)
        
        return proc

    def evolveA0(self, A_0, z_in, Gam_in, z_out):
        """
        Évolution avec le décalage vers le rouge de Gamma et A pour les noyaux.

        Arguments :
            A_0 (int) : Nombre initial de masse du noyau.
            z_in (float) : Décalage vers le rouge initial.
            Gam_in (float) : Valeur initiale de Gamma.
            z_out (float) : Décalage vers le rouge final.

        Returns :
            Gam_out (float) : Valeur finale de Gamma après l'évolution.
        """
        # Initialisation des variables
        k_com = self.GetJNuc(A_0)
        eps = 1e-9
        h1 = 1e-4
        hmin = 0.0

        # Définition des variables pour l'ODE
        Lny0 = log(Gam_in)
        nok, nbad = 0, 0

        # Résolution de l'ODE avec la méthode d'intégration odeint
        Lny = odeint(Lny0, z_in, z_out, eps, h1, hmin, nok, nbad, self.derivsA0)

        # Calcul de Gam_out à partir de Lny
        Gam_out = exp(Lny)

        return Gam_out
    
    def losseAdisi(self, Acurr, z, Gam, k0, multi):
        """
        Méthode qui calcule les pertes d'énergie de photodésintégration pour les noyaux exprimées en 1/y.
        Tous les milieux environnants (CMB+IR/V/UV) sont pris en compte.
        """

        # Si le noyau est un proton, les pertes d'énergie sont nulles
        if Acurr == 1:
            return 0.0

        beta_disi = 0.0

        # Définition des constantes globales
        self.Gam_com = Gam
        self.z_com = z
        self.k_com = k0
        self.s_com = 0
        self.gModel = self

        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        
        # Pertes d'énergie de photodésintégration (CMB + IRV/UV)

        lxmin1 = math.log(eth_1[k0])
        lxmax1 = math.log(eth_2[k0])
        sum1 = qgaus(fdisi, lxmin1, lxmax1)
        
        lxmin2 = math.log(eth_2[k0])
        lxmax2 = math.log(e_1)
        sum2 = qgaus(fdisi, lxmin2, lxmax2)
        
        lxmin3 = math.log(e_1)
        lxmax3 = math.log(e_max)
        sum3 = qgaus(fdisi, lxmin3, lxmax3)

        # Somme totale des pertes d'énergie
        sum_total = sum1 + sum2 + sum3
        
        # Si la somme est négative ou nulle, les pertes d'énergie sont nulles
        if sum_total <= 0:
            sum_total = 0.0
            if multi < 0:
                return sum_total
        
        # Réinitialisation de multi
        multi = 0
        
        # Sélection de la multiplicité en fonction de la somme des pertes d'énergie
        urange = random.random()
        umult = random.random()
        
        if urange <= sum1 / sum_total:
            multi = 1
        elif urange > (sum1 / sum_total) and urange <= (sum1 + sum2) / sum_total:
            self.s_com = 1
            sum21 = qgaus(fdisi, lxmin2, lxmax2)
            multi = 2
            if umult < sum21 / sum2:
                multi = 1
        else:
            dim = int(self.BRdisi(Acurr, 0))
            cum = 0
            for n in range(1, dim + 1):
                cum += self.BRdisi(Acurr, n)
                if umult < cum:
                    multi = n
                    break
            if Acurr == 9:
                multi = 5
        
        beta_disi = sum_total
        
        return beta_disi

    
    @staticmethod
    def derivsA0(z, Lny):
        """
        Fonction qui calcule le côté droit des équations différentielles pour les pertes d'énergie.

        Args:
            z (float): Valeur du décalage vers le rouge.
            Lny (float): Valeur du logarithme naturel de l'énergie.

        Returns:
            float: La valeur du côté droit de l'équation différentielle.
        """
        # Déclaration des variables globales utilisées dans la fonction
        global Stoch, dtdz_com, k_com, fA

        # Calcul de Gamma à partir du logarithme naturel de l'énergie
        Gam = exp(Lny)

        # Calcul du taux de variation du temps cosmique par rapport au décalage vers le rouge (z)
        dtdz = 1.0 / (H0 * (1 + z) * sqrt(om_m * pow((1 + z), 3) + om_l))

        # Utilisation du taux de variation du temps cosmique calculé précédemment
        # (cela peut être défini ailleurs dans le code)
        dtdz = dtdz_com

        # Calcul du coefficient beta en fonction du mode Stochastique ou Non-Stochastique
        if Stoch:
            beta = losseA0pair(z, Gam)
        else:
            beta = losseA0(z, Gam)

        # Modification du coefficient beta si certaines conditions sont remplies
        if k_com >= 0 and fA[k_com] > 1 and fA[k_com] < 56:
            beta *= (fA[k_com] / 4.0)

        # Calcul du côté droit de l'équation différentielle
        return dtdz * beta + 1 / (1 + z)
    
    @staticmethod
    def integIBL(eps, z):
        """
        Fonction qui calcule l'intégrale de n_IBL(e) / (2e^2) de e à +inf.

        Args:
            eps (float): Limite inférieure de l'intégrale.
            z (float): Valeur du décalage vers le rouge.

        Returns:
            float: Valeur de l'intégrale.
        """
        In_IR = 0.0

        # Réduction de l'epsilon si les pertes sont évaluées dans le référentiel comobile
        if fLosses == 3:
            eps /= (1 + z)

        # Conditions sur les valeurs d'epsilon
        if eps < epsIBL_max:
            if z > 5.0:
                z = 5.0
            if eps < epsIBL_min:
                eps = epsIBL_min

            # Calcul de l'intégrale en fonction du type de perte (fLosses)
            if fLosses == 1:
                # Utilisation du spectre de Malkan & Stecker (2006)
                logeps = log10(eps)
                logIn_IR = -splin2(z_IR, lE_IR, lint_IR, y2a, nz_IR, nE_IR, z, logeps)
                In_IR = pow(10.0, logIn_IR)
            elif fLosses == 2:
                # Approximation par une loi de puissance
                gam_IR = 2.5  # Indice de la loi de puissance
                N0_IR = 2.0 * 5.207e-2 / (gam_IR + 1)  # Normalisation
                sum_ = pow(eps, -(gam_IR + 1)) - pow(epsIBL_max, -(gam_IR + 1))
                if z <= 1.4:
                    In_IR = N0_IR * pow(1 + z, 3.1) * sum_
                else:
                    In_IR = N0_IR * pow(1 + 1.4, 3.1) * sum_
            elif fLosses == 3:
                # Kneise (2004), avec mise à l'échelle en z comme dans CRPropa (?)
                scale = splint(z_, scale_, scale2a, nz, z) / pow(1 + z, 3)
                I = pow(10, splint(logE_, logI_, logI2a, nE, log10(eps)))
                In_IR = (1 + z) * I * scale
            elif fLosses in [4, 5, 6]:
                # Dominguez (2011)
                logI_EBL = (
                    logImid_EBL if fLosses == 4 else
                    logIlow_EBL if fLosses == 5 else
                    logIupp_EBL
                )
                logeps = log10(eps)
                logIn_IR = -splin2(z_EBL, logeps_EBL, logI_EBL, y2a, nz_EBL, neps_EBL, z, logeps)
                In_IR = pow(10.0, logIn_IR) * pow(1 + z, 3)
            elif fLosses == 7:
                # Gilmore (2012)
                logeps = log10(eps)
                logIn_IR = -splin2(z_EBLg, logeps_EBLg, logI_EBLg, y2a, nz_EBLg, neps_EBLg, z, logeps)
                In_IR = pow(10.0, logIn_IR) * pow(1 + z, 3)

        return In_IR
    
    @staticmethod
    def inv_integIBL(In_IR, z):
        """
        Fonction qui calcule l'inverse de l'intégrale de n_IBL(e) / (2e^2) de e à +inf.

        Args:
            In_IR (float): Valeur de l'intégrale.
            z (float): Valeur du décalage vers le rouge.

        Returns:
            float: Valeur de l'énergie.
        """
        eps = 0.0

        # Limitation de z à 5
        if z > 5.0:
            z = 5.0

        # Calcul de l'énergie en fonction du type de perte (fLosses)
        if fLosses == 1:
            # Utilisation du spectre de Malkan & Stecker (2006)
            logIn_IR = log10(In_IR)
            logeps = splin2_inv(z_IR, lE_IR, lint_IR, y2a, nz_IR, nE_IR, z, -logIn_IR)
            eps = pow(10.0, logeps)
        elif fLosses == 2:
            # Approximation par une loi de puissance
            gam_IR = 2.5  # Indice de la loi de puissance
            N0_IR = 2.0 * 5.207e-2 / (gam_IR + 1)  # Normalisation
            if z <= 1.4:
                sum_ = In_IR / N0_IR * pow(1 + z, -3.1) + pow(epsIBL_max, -(gam_IR + 1))
            else:
                sum_ = In_IR / N0_IR * pow(1 + 1.4, -3.1) + pow(epsIBL_max, -(gam_IR + 1))
            eps = pow(sum_, -1.0 / (gam_IR + 1))
        elif fLosses == 3:
            # Kneise (2004), avec mise à l'échelle en z comme dans CRPropa (?)
            scale = splint(z_, scale_, scale2a, nz, z) / pow(1 + z, 3)
            I = In_IR / (1 + z) / scale
            if I <= 0:
                return (1 + z) * epsIBL_max
            logeps = splin2_inv(z_, scale_, scale2a, nz, z, -log10(I))
            eps = (1 + z) * pow(10, logeps)
        elif fLosses in [4, 5, 6]:
            # Dominguez (2011)
            logI_EBL = (
                logImid_EBL if fLosses == 4 else
                logIlow_EBL if fLosses == 5 else
                logIupp_EBL
            )
            logIn_IR = log10(In_IR / pow(1 + z, 3))
            logeps = splin2_inv(z_EBL, logeps_EBL, logI_EBL, y2a, nz_EBL, neps_EBL, z, -logIn_IR)
            eps = pow(10.0, logeps)
        elif fLosses == 7:
            # Gilmore (2012)
            logIn_IR = log10(In_IR / pow(1 + z, 3))
            logeps = splin2_inv(z_EBLg, logeps_EBLg, logI_EBLg, y2a, nz_EBLg, neps_EBLg, z, -logIn_IR)
            eps = pow(10.0, logeps)

        return eps
    
    @staticmethod
    def fdisi(lx):
        """
        Fonction qui calcule les pertes d'énergie par photodésintégration pour les noyaux.

        Args:
            lx (float): ln(eps'/MeV).

        Returns:
            float: Pertes d'énergie par photodésintégration en y^-1.
        """
        epsprime = exp(lx)  # en MeV
        sigma = 0.0  # en mb
        sigma_1, sigma_2, sigma_3 = 0.0, 0.0, 0.0
        gModel.sigma_disi(k_com, epsprime, sigma_1, sigma_2, sigma_3)

        # Calcul de la somme des sections efficaces en fonction de s_com
        if s_com == 0:
            sigma = sigma_1 + sigma_2 + sigma_3
        elif s_com == 1:
            sigma = sigma_1
        elif s_com == 2:
            sigma = sigma_2
        elif s_com == 3:
            sigma = sigma_3

        eps = epsprime / (2 * Gam_com) * 1e6  # en eV
        integ = integCMB(eps, z_com) + integIBL(eps, z_com)  # en 1/(y MeV^2 mb)

        return epsprime * epsprime / (Gam_com * Gam_com) * sigma * integ  # en y^-1
    
    @staticmethod
    def fpion_IR(lx):
        """
        Fonction qui calcule les pertes d'énergie par production de pions pour le rayonnement infrarouge.

        Args:
            lx (float): ln(eps'/MeV).

        Returns:
            float: Pertes d'énergie par production de pions en y^-1.
        """
        epsprime = exp(lx)  # en MeV
        sigma = 0.0  # en mb
        sigma_pion(epsprime, sigma)

        eps = epsprime / (2 * Gam_com) * 1e6  # en eV
        integ = integIBL(eps, z_com)

        return epsprime * epsprime / (Gam_com * Gam_com) * sigma * integ  # en y^-1
    
    #########################################

    def BRdisi(self, A, imult):
        # Définition des limites inférieures des intervalles A
        Alow = [1, 9, 10, 23]
        # Définition des limites supérieures des intervalles A
        Aup = [4, 9, 22, 56]
        # Nombre de valeurs dans chaque intervalle A
        dim = [2, 1, 6, 15]

        # Facteurs de branchement pour chaque A
        fBR = [0.8, 0.2,    # He
               1.0,         # Be
               0.1, 0.3, 0.1, 0.1, 0.2, 0.2,    # B,C,N,O,F,Ne
               0.1, 0.35, 0.1, 0.05, 0.15,       # Na,...,Fe
               0.045, 0.040, 0.035, 0.030, 0.025,
               0.020, 0.018, 0.015, 0.012, 0.010]

        indx = -1
        ptr = 0
        # Recherche de l'indice de l'intervalle A dans lequel A se trouve
        for i in range(4):
            if A >= Alow[i] and A <= Aup[i]:
                indx = i
                break
            ptr += dim[i]

        mult = 0
        # Si l'indice de l'intervalle est valide, détermine le nombre de valeurs dans cet intervalle
        if indx >= 0:
            mult = dim[indx]

        # Si imult est inférieur ou égal à 0, retourne le nombre de valeurs dans l'intervalle
        if imult <= 0:
            return float(mult)
        else:
            # Sinon, retourne le facteur de branchement correspondant à l'indice imult dans l'intervalle A
            return fBR[ptr + imult - 1]

    def sigma_disi(self, k, Ene, sigma_1, sigma_2, sigma_3):
        """
        Sous-routine qui calcule la section efficace de photodésintégration,
        nécessite (Z, A) et l'énergie des photons (MeV) en entrée.
        La sortie est en mbarn et est divisée en sigma_1: émission d'un seul nucléon (g, n) et (g, p);
        sigma_2: émission de deux nucléons (g, 2n), (g, 2p) et (g, np);
        sigma_3: section efficace normalisée constante.
        Nous avons suivi les méthodes de Stecker et Salamon ApJ 512 (1999) 521.
        """
        # Initialisation des variables de sortie
        sigma_1[0] = 0.0
        sigma_2[0] = 0.0
        sigma_3[0] = 0.0
        
        # Vérification des conditions limites pour l'énergie des photons
        if Ene >= e_max:
            return
        
        # Vérification de la validité de l'indice k
        assert k >= 0
        
        # Récupération de Z et A à partir de fA et fZ
        A = fA[k]
        Ze = fZ[k]
        
        # Calcul du facteur de normalisation
        Sigma_d = 60.0 * ((A - Ze) * Ze) / float(A)  # mb*MeV
        
        if A <= 4:
            assert A == 2 or A == 3 or A == 4
            k = nA + 1 - A
            
            # Calcul de sigma_1
            if Ene >= eth_1[k] and Ene <= e_1:
                W1 = W(eth_1[k], e0_1[k], delta_1[k])
                sigma_1[0] = (csi_1[k] * Sigma_d / W1) * \
                    math.exp(-2.0 * ((Ene - e0_1[k]) / delta_1[k]) ** 2)
            else:
                sigma_1[0] = 0.0
            
            # Calcul de sigma_2
            if Ene >= eth_2[k] and Ene <= e_1:
                W2 = W(eth_2[k], e0_2[k], delta_2[k])
                sigma_2[0] = (csi_2[k] * Sigma_d / W2) * \
                    math.exp(-2.0 * ((Ene - e0_2[k]) / delta_2[k]) ** 2)
            else:
                sigma_2[0] = 0.0
            
            # Calcul de sigma_3
            if Ene >= e_1 and Ene <= e_max:
                sigma_3[0] = zita[k] * Sigma_d / (e_max - e_1)
            else:
                sigma_3[0] = 0.0
            
            if fType == eBreitWigner2 or fType == eGauss2:
                sigma_1[0] += 2 * sigma_2[0] + 1.2 * sigma_3[0]
                sigma_2[0] = 0.0
                sigma_3[0] = 0.0
            return
        
        if fType == eBreitWigner2:
            sigma_2[0] = 0.0
            if Ene <= fe_1:
                sigma_1[0] = fNuclh1[k] / (1.0 + ((Ene - fNuclx1[k]) / fNuclw1[k]) ** 2) + \
                    fNuclh2[k] / (1.0 + ((Ene - fNuclx2[k]) / fNuclw2[k]) ** 2)
                sigma_3[0] = fAlphh1[k] / (1.0 + ((Ene - fAlphx1[k]) / fAlphw1[k]) ** 2) + \
                    fAlphh2[k] / (1.0 + ((Ene - fAlphx2[k]) / fAlphw2[k]) ** 2)
            elif Ene <= fe_max:
                sigma_1[0] = fNuclc[k]
                sigma_3[0] = fAlphc[k]
            else:
                sigma_1[0] = 0.0
                sigma_3[0] = 0.0
            return
        
        if fType == eGauss2:
            sigma_2[0] = 0.0
            if Ene <= fe_1:
                sigma_1[0] = fNuclh1[k] * math.exp(-((Ene - fNuclx1[k]) / fNuclw1[k]) ** 2)
                sigma_3[0] = fAlphh1[k] * math.exp(-((Ene - fAlphx1[k]) / fAlphw1[k]) ** 2)
            elif Ene <= fe_max:
                sigma_1[0] = fNuclc[k]
                sigma_3[0] = fAlphc[k]
            else:
                sigma_1[0] = 0.0
                sigma_3[0] = 0.0
            return
        
        if (Ene >= feth_1[k]) and (Ene <= fe_1):
            if fType == eBreitWigner and A > 4:
                sigma_1[0] = fcsi_1[k] / (1.0 + ((Ene - fe0_1[k]) / fdelta_1[k]) ** 2)
            else:
                W1 = W(feth_1[k], fe0_1[k], fdelta_1[k])
                sigma_1[0] = (fcsi_1[k] * Sigma_d / W1) * \
                    math.exp(-2.0 * ((Ene - fe0_1[k]) / fdelta_1[k]) ** 2)
        else:
            sigma_1[0] = 0.0
        
        if (Ene >= feth_2[k]) and (Ene <= fe_1):
            if fType == eBreitWigner and A > 4:
                sigma_2[0] = fcsi_2[k] / (1.0 + ((Ene - fe0_2[k]) / fdelta_2[k]) ** 2)
            else:
                W2 = W(feth_2[k], fe0_2[k], fdelta_2[k])
                sigma_2[0] = (fcsi_2[k] * Sigma_d / W2) * \
                    math.exp(-2.0 * ((Ene - fe0_2[k]) / fdelta_2[k]) ** 2)
        else:
            sigma_2[0] = 0.0
        
        if (Ene >= fe_1) and (Ene <= fe_max):
            sigma_3[0] = fzita[k] * Sigma_d / (fe_max - fe_1)
        else:
            sigma_3[0] = 0.0

    @staticmethod
    def sigma_pion(Ene):
        """
        Calcule la section efficace de production de pions.

        Args:
            Ene (float): Énergie du photon en NRF (MeV).

        Returns:
            float: Section efficace en mb.
        """
        if Ene >= epsthr:
            # Calcul de la variable s
            s = (mN / 1e9) * (mN / 1e9) + 2 * (mN / 1e9) * (Ene / 1e3)  # in GeV^2
        
            if s <= 5.:
                # Données tabulées de section efficace
                s_tab = [1.160, 1.165, 1.195, 1.250, 1.305, 1.375, 1.455, 1.470, 1.485, 1.500,
                         1.625, 1.680, 1.750, 1.795, 1.820, 1.875, 1.950, 1.980, 2.000, 2.050,
                         2.100, 2.150, 2.200, 2.250, 2.300, 2.350, 2.400, 2.450, 2.500, 2.625,
                        2.750, 2.875, 3.000, 3.250, 3.500, 4.000, 4.500, 5.000]
                sigma_tab = [.0000, .0003, .0502, .1279, .1952, .3173, .4970, .5229, .5414, .5317,
                             .2799, .2233, .1851, .1719, .1676, .1811, .2014, .2115, .2189, .2253,
                              .2343, .2491, .2690, .2874, .2915, .2752, .2490, .2257, .2112, .2051,
                              .2166, .2109, .1873, .1620, .1564, .1493, .1395, .1359]
                # Interpolation spline
                sigma = splint(s_tab, sigma_tab, y2a, len(s_tab), s)
            else:
                # Extrapolation à haute énergie
                logsqrts = [.3495, .5044, 1., 2., 3., 4., 4.548]
                logsigma = [-.867, -.915, -.941, -.851, -.684, -.514, -.423]
                # Interpolation spline
                sigma = 10 ** splint(logsqrts, logsigma, y2a, len(logsqrts), 0.5 * math.log10(s))
        
            if sigma < 0.:
                sigma = 0.
        else:
            sigma = 0.
    
        return sigma

    @staticmethod
    def losseA0pair(z, Gam):
        """
        Calcule les pertes d'énergie par production de paires pour les noyaux.

        Args:
            z (float): Facteur de redshift.
            Gam (float): Facteur de Lorentz.

        Returns:
            float: Pertes d'énergie par production de paires (1/y).
        """
        # Calcul de l'énergie correspondante
        E = mN * Gam * (1 + z)  # Énergie du proton correspondante (eV)
        beta_pair = 0.
        lbeta = 0.
    
        if E <= 1.e17:
            beta_pair = 0.
        else:
            if E <= 1.e24:
                lbeta = splint(lEne_ee, lbeta0_ee, y2a, nE_ee, math.log10(E))
            else:
                lbeta = -11.91234 - 0.61714286 * (math.log10(E) - 24.)
            beta_pair = (1 + z) ** 3 * 10 ** lbeta
    
        return beta_pair
    
    @staticmethod
    # Fonction pour les pertes d'énergie par production de paires et production de pions (approximation continue)
    def losseA0(z, Gam):
        """
        Calcule les pertes d'énergie par production de paires et de pions (approximées comme continues).

        Args:
            z (float): Facteur de redshift.
            Gam (float): Facteur de Lorentz.

        Returns:
            float: Pertes d'énergie par production de paires et de pions.
        """
        E = mN * Gam * (1 + z)  # Calcul de l'énergie correspondante
        beta = 0.
        if E > 1.e17:  # Seulement si l'énergie est suffisamment grande
            if not y2a:
                y2a = [splint(LogE, lbeta0, y2a, nE_p, math.log10(E))]
            lbeta = splint(LogE, lbeta0, y2a, nE_p, math.log10(E))
            beta = (1 + z) ** 3 * 10 ** lbeta
        return beta
    
    #@staticmethod
    # Fonction d'intégrale pour les pertes d'énergie dans l'infrarouge (non traduite)
    # def fdisi_IR(lx):
    #     """
    #     Fonction d'intégrale pour les pertes d'énergie dans l'infrarouge.
    #     """
    #     x = math.exp(lx)
    #     lxx = math.log10(x)
    #     if not y2a:
    #         y2a = [splie2(lE_IR, z_IR, ln_IR, nE_IR, nz_IR, y2a)]
    #     lnn = splin2(lE_IR, z_IR, ln_IR, y2a, nE_IR, nz_IR, lxx, z_com)
    #     ln = pow(10.0, lnn) / x
    #     f_IR = x * (ln) / x / x
    #     return f_IR

    def HadrRate(self, Gam, z):
        """
        Taux d'interaction pour la photoproduction de pions par des protons, en année^{-1}
        (celui pour les noyaux est A fois celui-ci)
        Seulement CMB.

        Args:
            Gam (float): Facteur de Lorentz.
            z (float): Facteur de redshift.

        Returns:
            float: Taux d'interaction.
        """
        rate = 0.
        if self.fStoch > 0 and Gam * mN * (1. + z) > 1.e19:
            if not y2a:
                y2a = [splint(LogE_pi, lbeta_pi, y2a, nE_pi, math.log10(Gam * mN * (1. + z)))]
            lrate = splint(LogE_pi, lbeta_pi, y2a, nE_pi, math.log10(Gam * mN * (1. + z)))
            rate = pow(1. + z, 3.) * pow(10., lrate)
        return rate

    def HadrRate_IR(self, Gam, z):
        """
        Taux d'interaction pour la photoproduction de pions par des protons, en année^{-1}
        (celui pour les noyaux est A fois celui-ci)
        Seulement IBL.

        Args:
            Gam (float): Facteur de Lorentz.
            z (float): Facteur de redshift.

        Returns:
            float: Taux d'interaction.
        """
        rate = 0.
        epsmin = epsthr
        epsmax = 2.e-6 * Gam * epsIBL_max
        if self.fLosses == 3:
            epsmax *= 1 + z
        if self.fStoch == 2 and epsmax > epsmin:
            if self.fLosses == 1:
                if z > 5.:
                    z = 5.
                if not y2a:
                    y2a = [splin2(lG_IR, z_IR, hrate_IR, y2a, nG_IR, nz_IR, math.log10(Gam), z)]
                rate = pow(10., splin2(lG_IR, z_IR, hrate_IR, y2a, nG_IR, nz_IR, math.log10(Gam), z))
            elif self.fLosses in [2, 3, 4, 5, 6, 7]:
                self.z_com = z
                self.Gam_com = Gam
                rate = qgaus(fpion_IR, math.log(epsmin), math.log(epsmax))
        return rate

    def NuclRate(self, j, Gam, z):
        """
        Taux d'interaction pour les noyaux.

        Args:
            j (int): Indice.
            Gam (float): Facteur de Lorentz.
            z (float): Facteur de redshift.

        Returns:
            float: Taux d'interaction.
        """
        self.gModel = self
        self.s_com = 1.
        self.k_com = j
        self.z_com = z
        self.Gam_com = Gam
        return qgaus(fdisi, math.log(self.fNuclt[j]), math.log(fe_1)) + qgaus(fdisi, math.log(fe_1), math.log(fe_max))

    def AlphRate(self, j, Gam, z):
        """
        Taux d'interaction pour les alpha.

        Args:
            j (int): Indice.
            Gam (float): Facteur de Lorentz.
            z (float): Facteur de redshift.

        Returns:
            float: Taux d'interaction.
        """
        self.gModel = self
        self.s_com = 3.
        self.k_com = j
        self.z_com = z
        self.Gam_com = Gam
        return qgaus(fdisi, math.log(self.fAlpht[j]), math.log(fe_1)) + qgaus(fdisi, math.log(fe_1), math.log(fe_max))

    def TotalRate(self, A, Z, j, Gam, z, type):
        """
        Taux d'interaction total.

        Args:
            A (int): Nombre de nucléons.
            Z (int): Nombre de protons.
            j (int): Indice.
            Gam (float): Facteur de Lorentz.
            z (float): Facteur de redshift.
            type (eProcType): Type de processus.

        Returns:
            float: Taux d'interaction total.
        """
        ndisDummy = -1
        rate_had = A * self.HadrRate(Gam, z)
        rate_had_IR = A * self.HadrRate_IR(Gam, z)
        rate_dis = 0
        rate_nuc = 0
        rate_alp = 0

        if self.fType == eBreitWigner2 or self.fType == eGauss2:
            if A > 1:
                rate_nuc = nFactor * self.NuclRate(j, Gam, z)
            if A > 4:
                rate_alp = aFactor * self.AlphRate(j, Gam, z)
        else:
            rate_dis = losseAdisi(A, z, Gam, j, ndisDummy)

        rate_tot = rate_had + rate_had_IR + rate_dis + rate_nuc + rate_alp

        if type is not None:
            u = random.uniform(0, rate_tot)
            if u < rate_had:
                type = eHadron
            elif u < rate_had + rate_had_IR:
                type = eHadron_IR
            elif u < rate_had + rate_had_IR + rate_dis:
                type = eDisi
            elif u < rate_had + rate_had_IR + rate_dis + rate_nuc:
                type = eNucleonProd
            else:
                type = eAlphaProd

        return rate_tot

    def PropagatePion(self, input):
        """
        Propage un pion.

        Args:
            input (Particle): Entrée.

        Returns:
            list: Liste des particules résultantes.
        """
        output = []
        charge = input.GetCharge()
        E = input.GetEprod()
        z = input.GetZprod()
        nBr = input.GetBranch() + 1
        assert E >= 0. and z >= 0.
        r = random.random()

        if charge == 0:
            output.append(Particle(0, 0, nBr, r * E, z, ePhoton))
            output.append(Particle(0, 0, nBr, (1. - r) * E, z, ePhoton))
            print("photon", r * E, "; photon", (1. - r) * E)
        else:
            E_nu = (1. - (mmu * mmu) / (mpi * mpi)) * E * r
            output.append(Particle(0, charge, nBr, E - E_nu, z, eMuon))
            output.append(Particle(0, 0, nBr, E_nu, z, eNeutrino_mu if charge > 0 else eAntineutrino_mu))
            print("muon", E - E_nu, "; neutrino", E_nu)

        input.SetEint(E)
        input.SetZint(z)
        input.SetIntMult(1002)
        return output

    def PropagateMuon(self, input):
        """
        Propage un muon.

        Args:
            input (Particle): Entrée.

        Returns:
            list: Liste des particules résultantes.
        """
        charge = input.GetCharge()
        E = input.GetEprod()
        z = input.GetZprod()
        nBr = input.GetBranch() + 1
        assert E >= 0. and z >= 0.

        gamma = E / mmu
        Enumax = (mmu - me * me / mmu) / 2.

        while True:
            E1 = random.uniform(0, Enumax)
            E2 = random.uniform(0, Enumax)
            if E1 + E2 < Enumax:
                E1 = Enumax - E1
                E2 = Enumax - E2
            E3 = mmu - E1 - E2
            p1 = E1
            p2 = E2
            p3 = math.sqrt(E3 * E3 - me * me)
            if E3 >= me and p3 <= p1 + p2 and p1 <= p2 + p3 and p2 <= p3 + p1:
                break

        c12 = (p3 * p3 - p1 * p1 - p2 * p2) / (2. * p1 * p2)
        s12 = math.sqrt(1. - c12 * c12)
        c1 = random.uniform(-1., 1.)
        s1 = math.sqrt(1. - c1 * c1)
        ca = math.cos(random.uniform(0, 2 * math.pi))
        c2 = c12 * c1 - s12 * s1 * ca

        Enu1 = gamma * (E1 + p1 * c1)
        Enu2 = gamma * (E2 + p2 * c2)
        Ee = E - Enu1 - Enu2

        input.SetEint(E)
        input.SetZint(z)
        input.SetIntMult(1003)

        output = [
            Particle(0, 0, nBr, Enu1, z, eAntineutrino_mu if charge > 0 else eNeutrino_mu),
            Particle(0, 0, nBr, Enu2, z, eNeutrino_e if charge > 0 else eAntineutrino_e),
            Particle(0, charge, nBr, Ee, z, eElectron)
        ]
        print("electron", Ee, "; neutrino", Enu1, "; neutrino", Enu2)
        return output

    def PropagateNeutrino(self, input):
        """
        Propage un neutrino.

        Args:
            input (Particle): Entrée.

        Returns:
            list: Liste des particules résultantes.
        """
        output = []
        E = input.GetEprod()
        z = input.GetZprod()
        input.SetEint(E / (1. + z))
        input.SetZint(0.)
        input.SetIntMult(0)
        return output

    def PropagatePhoton(self, input):
        """
        Propage un photon.

        Args:
            input (Particle): Entrée.

        Returns:
            list: Liste des particules résultantes.
        """
        output = []
        E = input.GetEprod()
        zOri = input.GetZprod()
        input.SetEint(E)
        input.SetZint(zOri)
        input.SetIntMult(1000)
        return output

    def PropagateElectron(self, input):
        """
        Propage un électron.

        Args:
            input (Particle): Entrée.

        Returns:
            list: Liste des particules résultantes.
        """
        output = []
        E = input.GetEprod()
        zOri = input.GetZprod()
        input.SetEint(E)
        input.SetZint(zOri)
        input.SetIntMult(1000)
        return output

    def Photopion(self, z, Gam, proc):
        """
        Calcul de l'énergie d'un pion photoproduit par un proton.

        Args:
            z (float): Facteur de redshift.
            Gam (float): Facteur de Lorentz.
            proc (eProcType): Type de processus.

        Returns:
            float: Énergie du pion photoproduit.
        """
        # Restituisce l'energia di un pione fotoprodotto da un protone a energia Gam*M e redshift z
        # (L'energia del protone sarà quella iniziale meno quella del pione.)
        angle = 0.0
        E = Gam * mN
        if proc == eHadron_IR:
            s = sample_sIBL(E, z, angle)
        else:
            eps = (1. + z) * sample_eps(E * (1. + z))
            s = sample_s(eps, E, angle)
        r_s = math.sqrt(s)
        p_star = math.sqrt((s - (mpi + mN) * (mpi + mN)) * (s - (mpi - mN) * (mpi - mN))) / (2. * r_s)
        E_star = (s - mN * mN + mpi * mpi) / (2. * r_s)
        cos_theta = random.uniform(-1., 1.)  # supp. distribuz. isotropa nel CM
        gamm = E / r_s  # da CM a Lab
        assert gamm * (E_star + p_star * cos_theta) <= mN * Gam
        return gamm * (E_star + p_star * cos_theta)

    def sample_sIBL(self, E, z, angle):
        """
        Echantillonne l'énergie CoM de l'interaction photohadronique sur l'IBL.

        Args:
            E (float): Energie.
            z (float): Redshift.
            angle (float): Angle.

        Returns:
            float: Energie de l'interaction.
        """
        # samples CoM energy of photohadronic interaction on the IBL
        # (complicated, maybe I'll write something about how it works later)
        # Armando di Matteo,  2 May 2013
        I, ph, eps, s = 0.0, 0.0, 0.0, 0.0
        Gam = E / mN
        eps_min = max(epsthr / (2.e-6 * Gam), epsIBL_min)
        I_max = integIBL(eps_min, z)
        I_min = max(integIBL(epsIBL_max, z), 1e-4 * I_max)
        N = 8
        I_ = np.zeros(N + 1)
        phi_ = np.zeros(N + 1)
        for i in range(N + 1):
            I_[i] = ((N - i) * I_min + i * I_max) / N
            phi_[i] = phi(mN * mN + 4 * E * inv_integIBL(I_[i], z))
        area_cum = np.zeros(N + 1)
        area_cum[0] = 0.
        for i in range(1, N + 1):
            area_cum[i] = area_cum[i - 1] + phi_[i - 1] * (I_[i] - I_[i - 1])
        while True:
            u = random.uniform(0., area_cum[N])
            i_min, i_max = 0, N
            while i_max - i_min > 1:
                i_mid = (i_min + i_max) // 2
                if u > area_cum[i_mid]:
                    i_min = i_mid
                else:
                    i_max = i_mid
            I = random.uniform(I_[i_min], I_[i_max])
            ph = random.uniform(0., phi_[i_min])
            eps = inv_integIBL(I, z)
            s = phi_inv(ph)
            if s <= mN * mN + 4 * E * eps:
                break
        angle[0] = math.acos(1 - (s - mN * mN) / (2 * E * eps))
        return s

    def GetDecay(self, A, Z):
        """
        Obtient la désintégration.

        Args:
            A (int): Nombre de masse.
            Z (int): Charge.

        Returns:
            float: Désintégration.
        """
        N = A - Z
        Q = np.nan
        if Z < 32 and N < 37:
            Q = qval[Z][N]
        if np.isnan(Q):
            Zstable = self.GetBetaDecayStableZ(A)
            mp, mn, me = 938.2720e6, 939.5653e6, 0.5110e6
            aC, aA, aP = 0.697e6, 23.285e6, 12.0e6
            if Z == Zstable:
                Q = 0.
            elif Z < Zstable:  # beta-minus decay
                Q = mn - mp - me - (2 * Z + 1) * pow(A, -1. / 3.) * aC + 4. * aA * (A - 2 * Z - 1) / A
                if A % 2:
                    if Z % 2:  # odd-odd to even-even
                        Q += 2 * aP * pow(A, -0.5)
                    else:  # even-even to odd-odd
                        Q -= 2 * aP * pow(A, -0.5)
            else:  # beta-plus decay
                Q = mp - me - mn + (2 * Z - 1) * pow(A, -1. / 3.) * aC - 4. * aA * (A - 2 * Z + 1) / A
                if A % 2:
                    if Z % 2:  # odd-odd to even-even
                        Q += 2 * aP * pow(A, -0.5)
                    else:  # even-even to odd-odd
                        Q -= 2 * aP * pow(A, -0.5)
                Q = -Q
        return Q

    def BetaDecay(self, input):
        """
        Désintégration bêta.

        Args:
            input (Particle): Entrée.

        Returns:
            list: Liste des particules résultantes.
        """
        A = input.GetMassNum()
        Z = input.GetCharge()
        Gam = input.GetGprod()
        z = input.GetZprod()
        nBr = input.GetBranch() + 1
        Q = self.GetDecay(A, Z)
        charge = -1  # beta-minus decay
        if Q < 0.:  # beta-plus decay
            Q = -Q
            charge = +1
        max_val = math.sqrt((me + Q) * (me + Q) - me * me) * (me + Q) * Q * Q
        while True:
            Ee = random.uniform(me, me + Q)
            u = random.uniform(0., max_val)
            if u <= math.sqrt(Ee * Ee - me * me) * Ee * (me + Q - Ee) * (me + Q - Ee):  # neglecting nucleus recoil and EM effects
                break
        Enumax = (mmu - me * me / mpi) / 2.
        while True:
            E1 = random.uniform(0, Enumax)
            E2 = random.uniform(0, Enumax)
            if E1 + E2 < Enumax:
                E1 = Enumax - E1
                E2 = Enumax - E2
            E3 = mmu - E1 - E2
            p1 = E1
            p2 = E2
            p3 = math.sqrt(E3 * E3 - me * me)
            if E3 >= me and p3 <= p1 + p2 and p1 <= p2 + p3 and p2 <= p3 + p1:
                break
        E1 = Gam * (E1 * random.uniform(0., 2.))  # (1+cos(theta))
        E2 = Gam * (E2 * random.uniform(0., 2.))  # in lab frame
        En = input.GetEprod() - E1 - E2
        input.SetGint(Gam)
        input.SetZint(z)
        input.SetIntMult(1003)
        output = []
        output.append(Particle(A, Z - charge, nBr, En, z, eNucleus if A > 1 else eProton))
        output.append(Particle(0, 0, nBr, E2, z, eNeutrino_e if charge > 0 else eAntineutrino_e))
        output.append(Particle(0, charge, nBr, E1, z, eElectron))
        print("nucleus", En, "; electron", E1, "; neutrino", E2)

        return output

# Tableau des valeurs qval[A][Z] contenant les énergies de désintégration [en MeV] pour différentes combinaisons de nombre de masse (A) et de numéro atomique (Z).
# Les valeurs NaN indiquent qu'elles ne sont pas définies.
qval = [
    [float('nan'), 782347, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [0, 0, 18591, 2.348e+07, 2.151e+07, 2.427e+07, 2.303e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), 0, 0, 0, 3.5083e+06, 1.1193e+07, 1.0651e+07, 1.5985e+07, 1.576e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [-1.374e+07, -2.29e+07, -290000, 0, 0, 1.60052e+07, 1.36066e+07, 2.0444e+07, 2.0623e+07, 2.502e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), -2.632e+07, -4.288e+06, -861890, 0, 0, 555900, 1.1506e+07, 1.1708e+07, 1.669e+07, 1.629e+07, 2.083e+07, 2.06e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), -2.523e+07, -1.21e+07, -1.79798e+07, -1.068e+06, 0, 0, 1.33689e+07, 1.34372e+07, 2.0644e+07, 1.9099e+07, 2.339e+07, 2.273e+07, 2.74e+07, 2.694e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), -1.2173e+07, -1.64948e+07, -3.64795e+06, -1.9824e+06, 0, 0, 156476, 9.7717e+06, 8.01e+06, 1.3167e+07, 1.181e+07, 1.656e+07, 1.579e+07, 2.071e+07, 2.124e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), -2.31e+07, -1.365e+07, -1.73381e+07, -2.22047e+06, 0, 0, 1.04207e+07, 8.68e+06, 1.3896e+07, 1.2527e+07, 1.797e+07, 1.719e+07, 2.275e+07, 2.378e+07, 2.847e+07, 2.906e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), -1.471e+07, -1.7767e+07, -5.14394e+06, -2.7542e+06, 0, 0, 0, 4.8223e+06, 3.8149e+06, 8.11e+06, 6.49e+06, 1.128e+07, 1.151e+07, 1.617e+07, 1.744e+07, 2.003e+07, 2.062e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -2.465e+07, -1.392e+07, -1.5417e+07, -2.76051e+06, -1.6552e+06, 0, 7.02453e+06, 5.6842e+06, 1.0818e+07, 8.48e+06, 1.351e+07, 1.338e+07, 1.784e+07, 1.786e+07, 2.198e+07, 2.224e+07, 2.58e+07, 2.545e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -1.3316e+07, -1.4509e+07, -4.4435e+06, -3.23883e+06, 0, 0, 0, 4.37581e+06, 2.4666e+06, 7.25e+06, 7.292e+06, 1.259e+07, 1.223e+07, 1.539e+07, 1.474e+07, 1.819e+07, 1.821e+07, 2.111e+07, 2.036e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -1.887e+07, -1.1175e+07, -1.389e+07, -3.5476e+06, -2.8423e+06, 0, 5.51545e+06, 3.835e+06, 9.352e+06, 9.069e+06, 1.4029e+07, 1.3284e+07, 1.7272e+07, 1.587e+07, 2.002e+07, 2e+07, 2.395e+07, 2.343e+07, 2.653e+07, 2.603e+07, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -2.011e+07, -1.0723e+07, -1.3095e+07, -4.7855e+06, -4.0561e+06, 0, 0, 0, 2.61001e+06, 1.8318e+06, 7.596e+06, 6.962e+06, 1.1736e+07, 1.011e+07, 1.342e+07, 1.174e+07, 1.628e+07, 1.564e+07, 1.93e+07, 1.895e+07, 2.217e+07, 2.094e+07, float('nan'), float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -1.521e+07, -1.858e+07, -1.2243e+07, -1.38766e+07, -4.2767e+06, -4.00427e+06, 0, 4.64236e+06, 3.6797e+06, 8.561e+06, 7.995e+06, 1.302e+07, 1.196e+07, 1.702e+07, 1.423e+07, 1.826e+07, 1.653e+07, 2.012e+07, 1.947e+07, 2.383e+07, 2.214e+07, 2.524e+07, float('nan'), float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -1.398e+07, -1.7e+07, -1.0812e+07, -1.274e+07, -5.066e+06, -4.81236e+06, 0, 0, 0, 1.49188e+06, 224310, 5.845e+06, 4.601e+06, 1.05e+07, 7.77e+06, 1.241e+07, 1.069e+07, 1.48e+07, 1.357e+07, 1.884e+07, 1.75e+07, 2.093e+07, 2.074e+07, float('nan'), float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -2.124e+07, -1.505e+07, -1.812e+07, -1.1667e+07, -1.4334e+07, -4.9424e+06, -4.2324e+06, 0, 1.71048e+06, 248500, 5.374e+06, 3.9886e+06, 1.0413e+07, 7.9e+06, 1.21e+07, 1.029e+07, 1.476e+07, 1.374e+07, 1.862e+07, 1.773e+07, 2.122e+07, 2.116e+07, 2.481e+07, float('nan')],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -1.5e+07, -1.826e+07, -1.123e+07, -1.379e+07, -6.138e+06, -5.3962e+06, 0, 0, 0, 167180, 1.14222e+06, 4.86517e+06, 2.937e+06, 6.64e+06, 4.69e+06, 8.29e+06, 7.24e+06, 1.22e+07, 1.111e+07, 1.511e+07, 1.541e+07, 1.852e+07, 1.79e+07],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -2.248e+07, -1.63e+07, -1.851e+07, -1.198e+07, -1.2686e+07, -5.5826e+06, -5.49201e+06, 0, 709680, 0, 4.9165e+06, 3.442e+06, 7.48e+06, 5.76e+06, 9.51e+06, 7.84e+06, 1.244e+07, 1.141e+07, 1.501e+07, 1.539e+07, 1.901e+07, 1.844e+07, 2.181e+07],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -1.564e+07, -1.836e+07, -1.113e+07, -1.16193e+07, -6.0626e+06, -5.9661e+06, 0, -813870, 0, 565000, 1.50469e+06, 2.4916e+06, 599000, 4.583e+06, 3.14e+06, 6.838e+06, 5.7e+06, 9.79e+06, 8.41e+06, 1.217e+07, 1.085e+07, 1.421e+07, 1.32e+07, 1.66e+07],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -2.262e+07, -1.615e+07, -1.69e+07, -1.1879e+07, -1.2805e+07, -6.14746e+06, -5.91386e+06, 0, 1.31107e+06, 0, 3.52552e+06, 1.815e+06, 5.66e+06, 4.204e+06, 7.717e+06, 6.644e+06, 1.209e+07, 1.097e+07, 1.422e+07, 1.386e+07, 1.631e+07, 1.59e+07, 1.849e+07, 1.785e+07],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -1.543e+07, -1.836e+07, -1.1077e+07, -1.1812e+07, -6.111e+06, -5.96025e+06, 0, 0, 0, 653500, 0, 2.8069e+06, 1.125e+06, 5.414e+06, 3.82e+06, 7.108e+06, 6.129e+06, 1.118e+07, 1.015e+07, 1.314e+07, 1.3e+07, 1.557e+07, 1.482e+07, 1.712e+07, 1.669e+07],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -2.173e+07, -1.608e+07, -1.684e+07, -1.2288e+07, -1.2495e+07, -6.262e+06, -5.94445e+06, 0, 0, 1.15521e+06, 2.25099e+06, 0, 4.176e+06, 2.615e+06, 6.554e+06, 5.72e+06, 1.06e+07, 9.392e+06, 1.189e+07, 1.163e+07, 1.407e+07, 1.359e+07, 1.58e+07, 1.521e+07],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -1.531e+07, -1.86e+07, -1.1075e+07, -1.272e+07, -6.127e+06, -5.93225e+06, 0, 0, 0, 0, 0, 0, 3.3184e+06, 1.562e+06, 5.797e+06, 4.167e+06, 7.714e+06, 6.448e+06, 1.125e+07, 1.025e+07, 1.259e+07, 1.242e+07, 1.479e+07, 1.429e+07],
    [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), -2.235e+07, -1.616e+07, -1.706e+07, -1.2387e+07, -1.2828e+07, -6.243e+06, -5.91995e+06, 0, 0, 0, 0, 0, 0, 3.0886e+06, 1.249e+06, 5.559e+06, 3.873e+06, 7.617e+06, 6.223e+06, 1.086e+07, 9.678e+06, 1.198e+07, 1.189e+07, 1.419e+07]
]





