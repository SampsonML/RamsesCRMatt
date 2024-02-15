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
    t = np.array( ds.current_time.in_units('Myr') )
    
    # grab the cr energy density and x,y,z values
    data = ds.all_data()
    cr_data = np.array( data[("ramses","hydro_CRegy_01")] ) 
    x_values = np.array( data[("ramses","x")] )
    y_values = np.array( data[("ramses","y")] )
    z_values = np.array( data[("ramses","z")] )

    # calculate the moment of inertia tensor
    I_Ecr = np.zeros((3, 3), dtype=float)

    for idx, rho_cr in enumerate(cr_data):

        # column 1
        I_Ecr[0,0] += rho_cr * ((x_values[idx]) * (x_values[idx])) 
        I_Ecr[0,1] += rho_cr * ((x_values[idx]) * (y_values[idx])) 
        I_Ecr[0,2] += rho_cr * ((x_values[idx]) * (z_values[idx])) 

        # column 2
        I_Ecr[1,1] += rho_cr * ((y_values[idx]) * (y_values[idx])) 
        I_Ecr[1,2] += rho_cr * ((y_values[idx]) * (z_values[idx])) 

        # column 3
        I_Ecr[2,2] += rho_cr * ((z_values[idx]) * (z_values[idx])) 

    # symmetric properties
    I_Ecr[1,0] = I_Ecr[0,1] 
    I_Ecr[2,0] = I_Ecr[0,2] 
    I_Ecr[2,1] = I_Ecr[1,2] 

    # calculate the eigenvalues
    eigs, evecs = np.linalg.eig(I_Ecr)
    
    return I_Ecr, eigs, evecs, t
