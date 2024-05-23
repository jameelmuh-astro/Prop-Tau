import numpy as np
from Particle import eParticleType, Particle
import Propagate
import Output
from Constants import *

import numpy as np
from Particle import eParticleType, Particle
import Propagate
import Output
from Constants import *

class FluoDetector:
    """
    Model to simulate a fluorescence detector of the Pierre Auger Observatory.
    """

    def __init__(self, position, detection_threshold, noise_level):
        """
        Initializes a fluorescence detector.

        :param position: Position of the detector (tuple of coordinates x, y, z).
        :param detection_threshold: Detection threshold in terms of number of photons.
        :param noise_level: Background noise level (average number of photons of noise per unit of time).
        """
        self.position = position
        self.detection_threshold = detection_threshold
        self.noise_level = noise_level
        self.detected_events = []

    def detect(self, input_particle, model=None):
        """
        Method to detect photons coming from a particle if the number of received photons is above the threshold.

        - input_particle: Instance of the Particle class.
        :param propagation_model: Function or method to model the propagation of photons.
        :return: True if the particle is detected, otherwise False.
        """
        photons_received = Propagate.PropagatePhoton(input_particle, model)
        signal = photons_received - self.noise_level
        if signal >= self.detection_threshold:
            self.detected_events.append((particle, photons_received))
            return True
        return False

    def propagate_photons(self, particle, propagation_model):
        """
        Simulates the propagation of photons from the particle to the detector.

        :param particle: Instance of the Particle class.
        :param propagation_model: Function or method to model the propagation of photons.
        :return: Number of photons received by the detector.
        """
        distance = self.calculate_distance(particle)
        photons_emitted = propagation_model(particle)
        attenuation_factor = self.attenuation(distance)
        photons_received = photons_emitted * attenuation_factor
        return photons_received

    def calculate_distance(self, particle):
        """
        Calculates the distance between the detector and a given particle.

        :param particle: Instance of the Particle class.
        :return: Distance between the detector and the particle.
        """
        dx = self.position[0] - particle.position[0]
        dy = self.position[1] - particle.position[1]
        dz = self.position[2] - particle.position[2]
        return np.sqrt(dx**2 + dy**2 + dz**2)

    def attenuation(self, distance):
        """
        Calculates the attenuation factor of photons as a function of distance.

        :param distance: Distance between the source of photons and the detector.
        :return: Attenuation factor (between 0 and 1).
        """
        # Simple model of inverse square attenuation, can be refined with atmospheric factors.
        return 1 / (4 * np.pi * distance**2)

    def summary(self):
        """
        Returns a summary of the detected events.

        :return: List of detected events.
        """
        return self.detected_events


# Example of photon propagation model (to be refined according to the real model)
def example_propagation_model(particle):
    """
    Simple model of photon propagation from a particle.

    :param particle: Instance of the Particle class.
    :return: Number of photons emitted by the particle.
    """
    return particle.energy * 1e3  # Conversion of energy to number of photons


# Example usage:
# Assuming that the Particle class and its instances are already defined in your project.

# Creation of a detector at a given position with a detection threshold and a noise level.
detector = FluoDetector(position=(1000, 2000, 3000), detection_threshold=1000, noise_level=300)

# Example particle (to be defined according to your Particle class)
# particle = Particle(position=(0, 0, 0), energy=1e19, type='tau')

# Detection of the particle with the propagation model
# detected = detector.detect(particle, example_propagation_model)

# Get a summary of the detected events
# print(detector.summary())
