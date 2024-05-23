import numpy as np
from Particle import eParticleType, Particle
import Propagate
import io

class Output:
    def __init__(self, filename, type):
        """
        Initialize the Output object with a filename and type, and set up the initial data structure.

        Parameters:
        - filename (str): The name of the file where the output data will be saved.
        - type (str): The type of output, which could influence how data is processed or formatted.

        The data attribute is initialized with a dictionary containing keys for 'nuc', 'ev', and 'summary'.
        Each key maps to a list or a nested dictionary designed to hold specific types of data:
        - 'nuc': A list intended to store information about nuclear particles.
        - 'ev': A list intended to store event-related data.
        - 'summary': A dictionary with keys 'nuc', 'pho', 'neu', and 'ele', each mapping to a list intended to store summarized information about nuclear particles, photons, neutrinos, and electrons, respectively.
        """
        self.filename = filename
        self.type = type
        self.data = {
            'nuc': [],
            'ev': [],
            'summary': {
                'nuc': [],
                'pho': [],
                'neu': [],
                'ele': []
            }
        }
        self.buffer = ""  # Initialize the buffer as an empty string
        self.file = open(self.filename, 'w')  # Open the file in write mode

    """def __del__(self):
        # Convert any ndarray objects to lists before serializing
        serialized_data = self.data
        for key, value in serialized_data.items():
            if isinstance(value, np.ndarray):
                serialized_data[key] = value.tolist()""" # TEST
    
    def __del__(self):
        self.file.close()  # Close the file when the Output object is deleted
    
    def ndarray_to_list(self, item):
        """
        Recursively convert numpy ndarray objects to lists.
        """
        if isinstance(item, np.ndarray):
            return item.tolist()
        elif isinstance(item, dict):
            return {key: self.ndarray_to_list(value) for key, value in item.items()}
        elif isinstance(item, list):
            return [self.ndarray_to_list(element) for element in item]
        else:
            return item
        
    def Write(self):
        self.file.write(self.buffer)  # Write the buffer to the file
        self.buffer = ""  # Clear the buffer
        
        # Convert any ndarray objects to lists before serializing
        #serialized_data = self.ndarray_to_list(self.data)
                    
        #with io.open(self.filename, 'w', encoding='utf-8') as f:   # TEST (In __del__, this is maybe not necessary
        #    f.write(json.dumps(serialized_data, ensure_ascii=False, indent=4))

    def StreamNuc(self, mypart):
        self.data['nuc'].append({
            'branch': mypart.GetBranch(),
            'intmult': mypart.GetIntMult(),
            'Acurr': mypart.GetMassNum(),
            'Zecurr': mypart.GetCharge(),
            'Flav': mypart.GetFlavor(),
            'EOri': mypart.GetEprod(),
            'EEnd': mypart.GetEint(),
            'zEnd': mypart.GetZprod()
        })

    def StreamEv(self, timexev, branxev, cseed):
        self.data['ev'].append({
            'timexev': timexev,
            'branxev': branxev,
            'seed': cseed
        })

    def Injec(self, evt, primary):
        self.data['summary']['nuc'].append({
            'event': evt,
            'Energy_i': primary.GetEprod(),
            'A_i': primary.GetMassNum(),
            'Z_i': primary.GetCharge()
        })

    def Earth(self, product):
        if product.GetType() in [eParticleType.Proton, eParticleType.Neutron, eParticleType.Nucleus]:
            self.data['summary']['nuc'].append({
                'E_nuc': product.GetEint(),
                'A_nuc': product.GetMassNum(),
                'Z_nuc': product.GetCharge()
            })
        elif product.GetType() in [eParticleType.Neutrino_e, eParticleType.Antineutrino_e, eParticleType.Neutrino_mu, eParticleType.Antineutrino_mu]:
            self.data['summary']['neu'].append({
                'E_neu': product.GetEint(),
                'F_neu': product.GetFlavor()
            })
        elif product.GetType() in [eParticleType.Photon, eParticleType.Electron]:
            self.data['summary']['pho'].append({
                'E_pho': product.GetEprod(),
                'z_pho': product.GetZprod()
            })
            self.data['summary']['ele'].append({
                'Z_ele': product.GetCharge(),
                'E_ele': product.GetEprod(),
                'z_ele': product.GetZprod()
            })