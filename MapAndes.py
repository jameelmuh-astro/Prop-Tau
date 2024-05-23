import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D

# Lecture des données à partir d'un fichier .txt
data = np.loadtxt('DEM_Andes.txt', delimiter=' ')
auger_data = np.loadtxt('Auger_coord.txt', delimiter=',')
longitudes = data[:, 0]
latitudes = data[:, 1]
elevations = data[:, 2]

# Conversion des élévations en km
elevations_km = elevations / 1000.0

# Créer une grille régulière de longitudes et latitudes
num_points = 100  # Nombre de points dans chaque dimension pour la grille
lon_grid = np.linspace(min(longitudes), max(longitudes), num_points)
lat_grid = np.linspace(min(latitudes), max(latitudes), num_points)
lon_grid, lat_grid = np.meshgrid(lon_grid, lat_grid)

# Interpoler les élévations à la grille régulière
elev_grid = griddata((longitudes, latitudes), elevations_km, (lon_grid, lat_grid), method='cubic')

# Projection des coordonnées de l'Observatoire Pierre Auger sur la surface
auger_lon = auger_data[:, 1]
auger_lat = auger_data[:, 0]
auger_elev = griddata((longitudes, latitudes), elevations_km, (auger_lon, auger_lat), method='cubic')

# Création de la figure et d'un axe 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Création de la surface 3D lisse
surface = ax.plot_surface(lon_grid, lat_grid, elev_grid, cmap='terrain', edgecolor='none')

# Ajout du contour de l'Observatoire Pierre Auger comme un hexagone rouge
ax.plot(auger_lon, auger_lat, auger_elev, 'r-')

# Ajout d'une grille
ax.grid(True)

# Étiquettes des axes avec longitudes et latitudes
ax.set_xlabel('Longitude (°)')
ax.set_ylabel('Latitude (°)')
ax.set_zlabel('Elevation (km)')

# Ajout d'une barre de couleur pour les élévations
cbar = fig.colorbar(surface, ax=ax, shrink=0.5, aspect=5)
cbar.set_label('Élévation (km)')

# Affichage de la figure
plt.show()
