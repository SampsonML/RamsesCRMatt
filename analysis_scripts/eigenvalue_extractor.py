# ----
# Calculate the eigenvalues of the CR Energy tensor
# from RAMSES output snapshots
# ----

# imports
import numpy as np
import matplotlib.pyplot as plt
import yt

def inertiaTensor(filename):
    """
    Calculates the inertia tensor for a given field
    from RAMSES data outputs
    
    Inputs
    ------
    filename: string -- path to the data
    
    Returns
    -------
    I_Ecr: 3 by 3 numpy array || the inertial tensor of the field
    eigs: numpy array         || the eigenvalues
    evecs: tuple              || the eigenvectors
    """

    # load the file in 
    ds = yt.load(filename)

    # grab the cr energy density and x,y,z values
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
    eigs, evecs = np.linalg.eig(I_Ecr)
    
    return I_Ecr, eigs, evecs
