"""
File : Main.py 
For tests before implementing in Prop_Tau
"""
VERSION = 'TEST_TAU'
REVISION = 'idk, maybe 1 or 2 ?'

import os
import sys
import time
import numpy as np

from Particle import eParticleType, Particle
import Propagate
import Output
from Constants import *

# Define global variables
Nevts = 100 # Number of events to be processed
Seed = 65539  # Random seed for reproducibility
OutType = 0  # Output type
Eini = 1.e18  # Initial energy
NatEarth = [0] * 56  # Placeholder for Nucleus data
Photons, Electrons, Positrons = 0, 0, 0  # Counters for photons, electrons, and positrons
NuatEarth = [0] * 6  # Counters for different types of neutrinos reaching Earth

def Follow(mystack, mymodel, myout):
    """
    Perform particle tracking.

    Parameters:
    - mystack: list of Particle objects representing the particle stack
    - mymodel: Propagate object responsible for particle propagation
    - myout: Output object for streaming particle data
    
    """
    it = mystack.pop(0)    # Pop the first particle from the stack
    if Particle.GetType(it) != eParticleType.Null:
        products = mymodel.PropagateParticle(it)    # Propagate the particle
        if products:        # If there are products, add them to the stack
            mystack.extend(products)  
        myout.StreamNuc(it)  # Stream nuclear data
        if it.GetIntMult() == 1000: # TEST
            UpdateParticleCounters(it)

def UpdateParticleCounters(particle):
    """
    Update particle counters based on the type of particle.

    Parameters:
    - particle: Particle object representing the particle to be counted
    
    """
    global Photons, Electrons, Positrons, NuatEarth
    if particle.GetType() == eParticleType.Photon:
        Photons += 1
    elif particle.GetType() == eParticleType.Electron:
        if particle.GetCharge() > 0:
            Positrons += 1
        else:
            Electrons += 1
    UpdateNuatEarth(particle)

def UpdateNuatEarth(particle):
    """
    Update neutrino counters reaching Earth based on the type of particle.

    Parameters:
    - particle: Particle object representing the neutrino particle
    
    """
    global NuatEarth
    if particle.GetZint() == 0:
        nu_type = particle.GetType()
        if nu_type == eParticleType.Neutrino_e:
            NuatEarth[0] += 1
        elif nu_type == eParticleType.Antineutrino_e:
            NuatEarth[1] += 1
        elif nu_type == eParticleType.Neutrino_mu:
            NuatEarth[2] += 1
        elif nu_type == eParticleType.Antineutrino_mu:
            NuatEarth[3] += 1
        elif nu_type == eParticleType.Neutrino_tau:
            NuatEarth[4] += 1
        elif nu_type == eParticleType.Antineutrino_tau:
            NuatEarth[5] += 1
            
def GetOptions(argc, argv):
    """
    Get command line options.

    Parameters:
    - argc: int, number of command line arguments
    - argv: list of str, command line arguments

    Returns:
    bool, indicating whether the options were successfully obtained
    """
    global Nevts, Seed, Eini, OutType
    if len(argv) < 2:
        print("USAGE: python Main.py -N [number of events] -s [random seed] -e [initial energy] -o [output type]")
        return False
    # Parse options...
    return True

def PrintSummary():
    """
    Print a summary of the results.

    Parameters:
    None

    Returns:
    None
    """
    global Photons, Electrons, Positrons, NuatEarth
    print(" --------------------------------------------------")
    print(" photons from pi0 decay:     ", Photons)
    print(" electrons produced:         ", Electrons)
    print(" positrons produced:         ", Positrons)
    print(" --------------------------------------------------")
    print("       Neutrino Type     Number reaching Earth")
    print(" --------------------------------------------------")
    print(" electron neutrinos:          ", NuatEarth[0])
    print(" electron antineutrinos:      ", NuatEarth[1])
    print(" muon neutrinos:              ", NuatEarth[2])
    print(" muon antineutrinos:          ", NuatEarth[3])
    print(" tau neutrinos:               ", NuatEarth[4])
    print(" tau antineutrinos:           ", NuatEarth[5])
    print(" --------------------------------------------------")
    total_nu = sum(NuatEarth)
    print("         Total          ", total_nu)

if __name__ == "__main__":
    # Check if command line options were successfully obtained
    if not GetOptions(len(sys.argv), sys.argv):
        sys.exit(1)

    # Set the random seed for reproducibility
    np.random.seed(Seed)

    # Initialize the Propagate object for particle propagation
    model = Propagate.Propagate() # TEST

    # Print initial energy of nu tau and prepare output file
    print("** Initial energy of nu tau is:", Eini)
    filename = f"Output-N{Nevts}_E{Eini:.1e}.txt"
    myout = Output.Output(filename, OutType)
    print("** Output file is:", filename)
    print("** Events to be generated:", Nevts)

    t_tot = 0  # Initialize total time elapsed
    for event in range(Nevts):
        branchxev = 0  # Initialize branch events counter
        print("** Event", event, "****************************************")
        start = time.time()  # Record start time of the event processing

        # Create primary particle (nu tau) and inject it into the system
        primary = Particle(Eini, 0, eParticleType.Neutrino_tau)
        myout.Injec(event, primary)
        stack = [primary]  # Initialize the particle stack with the primary particle
        while stack:
            Follow(stack, model, myout)  # Perform particle tracking
            branchxev += 1  # Increment branch events counter

        end = time.time()  # Record end time of the event processing
        cpu_time_used = end - start  # Calculate elapsed CPU time for the event
        print("** event ended - elapsed time (s)", cpu_time_used)
        t_tot += cpu_time_used  # Accumulate total time elapsed
        CurrentSeed = np.random.get_state()
        myout.StreamEv(cpu_time_used, branchxev, CurrentSeed)  # Stream event data
        myout.Write()  # Write event data to the output file


    print("  ++++  Succesfully processed", Nevts, "events ++++")
    print(" --------------------------------------------------")
    PrintSummary()  # Print a summary of the results
    print("\nTotal time elapsed:", t_tot)  # Print total time elapsed
    myout.file.close()  # Close the output file