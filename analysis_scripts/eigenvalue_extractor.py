# ----
# Calculate the eigenvalues of the CR Energy tensor
# from RAMSES output snapshots
# ----

# imports
import numpy as np
import matplotlib.pyplot as plt
import yt

# load in the data
filename = "/scratch/gpfs/ms0821/sampson_ramses/delta_test/output_00062/info_00062.txt"
ds = yt.load(filename)

# Grab the cr energy density and x,y,z values
my_sphere = ds.sphere("c", (20.0, "pc"))
cr_data   = my_sphere.to_dataframe(("ramses","hydro_CRegy_01"))['hydro_CRegy_01']
x_values  = my_sphere.to_dataframe(("ramses","x"))['x']
y_values  = my_sphere.to_dataframe(("ramses","y"))['y']
z_values  = my_sphere.to_dataframe(("ramses","z"))['z']

# calculate the moment of inertia tensor
I_Ecr = np.zeros((3, 3), dtype=float)

for idx, rho_cr in enumerate(cr_data):
    
    # column 1
    I_Ecr[0,0] += rho_cr * ((y_values[idx])**2 + (z_values[idx])**2) 
    I_Ecr[0,1] -= rho_cr * ((x_values[idx]) * (y_values[idx])) 
    I_Ecr[0,2] -= rho_cr * ((x_values[idx]) * (z_values[idx])) 
    
    # column 2
    I_Ecr[1,1] += rho_cr * ((x_values[idx])**2 + (z_values[idx])**2) 
    I_Ecr[1,2] -= rho_cr * ((y_values[idx]) * (z_values[idx])) 
    
    # column 3
    I_Ecr[2,2] += rho_cr * ((x_values[idx])**2 + (y_values[idx])**2) 
    
# symmetric properties
I_Ecr[1,0] = I_Ecr[0,1] 
I_Ecr[2,0] = I_Ecr[0,2] 
I_Ecr[2,1] = I_Ecr[1,2] 

# calculate the eigenvalues
eigenvalues, eigenvectors = np.linalg.eig(I_Ecr)
