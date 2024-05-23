import numpy as np

# Topping
class AndesGeometry:
    def __init__(self, filepath):
        """
        Initialise la classe AndesGeometry avec les données de topographie.

        Parameters:
            filepath (str): Chemin vers le fichier contenant les données de topographie.
        """
        self.filepath = filepath
        self.topography = None
        self.load_topography()

    def load_topography(self):
        """
        Charge les données de topographie à partir du fichier spécifié lors de l'initialisation.
        """
        try:
            # Charger les données [longitude, latitude, élévation] depuis le fichier
            self.topography = np.loadtxt(self.filepath, delimiter=' ')
        except Exception as e:
            print(f"Erreur lors du chargement des données de topographie : {e}")

    def get_elevation(self, longitude, latitude):
        """
        Retourne l'élévation à une position spécifique basée sur les coordonnées de longitude et latitude.

        Parameters:
            longitude (float): Longitude de la position.
            latitude (float): Latitude de la position.

        Returns:
            float: Élévation à la position spécifiée, ou 0 si les coordonnées sont en dehors de la plage des 
            coordonnées fournies dans le fichier de données de topographie
        """
        if self.topography is not None:
            min_longitude = np.min(self.topography[:, 0])
            max_longitude = np.max(self.topography[:, 0])
            min_latitude = np.min(self.topography[:, 1])
            max_latitude = np.max(self.topography[:, 1])

            if longitude < min_longitude or longitude > max_longitude or latitude < min_latitude or latitude > max_latitude:
                return 0.0
            else:
                distances = np.sqrt((self.topography[:,0] - longitude)**2 + (self.topography[:,1] - latitude)**2)
                nearest_index = np.argmin(distances)
                return self.topography[nearest_index, 2]
        else:
            print("Aucune topographie chargée.")
            return None

    def print_topography_info(self):
        """
        Affiche des informations de base sur les données de topographie chargées.
        """
        if self.topography is not None:
            print(f"Topographie chargée avec {self.topography.shape[0]} points.")
        else:
            print("Aucune topographie chargée.")

if __name__ == "__main__":
    # Demander à l'utilisateur de saisir la latitude et la longitude
    latitude = float(input("Entrez la latitude : "))
    longitude = float(input("Entrez la longitude : "))
    
    # Initialisez un objet AndesGeometry en lui passant le chemin du fichier de données de topographie
    geometry = AndesGeometry('DEM_Andes.txt')
    
    # Appelez la méthode get_elevation avec des coordonnées de longitude et latitude pour obtenir l'élévation à cette position
    elevation = geometry.get_elevation(longitude, latitude)
    print(f"Élévation : {elevation}")

    # Appelez la méthode print_topography_info pour afficher des informations de base sur les données de topographie chargées
    geometry.print_topography_info()