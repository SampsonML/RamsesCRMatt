# ----
# Calculate the eigenvalues of the CR Energy tensor
# from RAMSES output snapshots
# ----

# imports
import numpy as np
import matplotlib.pyplot as plt
import yt

def cr_spatial_covariance(filename, group=1, normalize=True):
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
    group_name = "hydro_CRegy_0" + str(group)
    ds = yt.load(filename)
    t = np.array( ds.current_time.in_units('s') )
    
    # grab the cr energy density and x,y,z values
    data = ds.all_data()
    cr_data = np.array( data[("ramses", group_name)] ) 
    x_values = np.array( data[("ramses","x")] )
    y_values = np.array( data[("ramses","y")] )
    z_values = np.array( data[("ramses","z")] )
    
    # Cell volume
    dx = np.array( data[("ramses","dx")] )
    dy = np.array( data[("ramses","dx")] )
    dz = np.array( data[("ramses","dx")] )
    E_tot = 0
    
    # max values
    code_len =  np.array(ds.domain_width)
    box_len = np.array(ds.length_unit) * code_len
    x_mid = box_len[0] / 2
    y_mid = box_len[1] / 2
    z_mid = box_len[2] / 2

    # calculate the moment of inertia tensor
    I_Ecr = np.zeros((3, 3), dtype=float)

    for idx, rho_cr in enumerate(cr_data):

        # column 1
        I_Ecr[0,0] += rho_cr * ((x_values[idx] - x_mid) * (x_values[idx]- x_mid)) # / dx[idx] * dy[idx] * dz[idx]
        I_Ecr[0,1] += rho_cr * ((x_values[idx] - x_mid) * (y_values[idx]- y_mid)) # / dx[idx] * dy[idx] * dz[idx]
        I_Ecr[0,2] += rho_cr * ((x_values[idx] - x_mid) * (z_values[idx]- z_mid)) # / dx[idx] * dy[idx] * dz[idx]

        # column 2
        I_Ecr[1,1] += rho_cr * ((y_values[idx]- y_mid) * (y_values[idx]- y_mid)) # / dx[idx] * dy[idx] * dz[idx]
        I_Ecr[1,2] += rho_cr * ((y_values[idx]- y_mid) * (z_values[idx]- z_mid)) # / dx[idx] * dy[idx] * dz[idx]

        # column 3
        I_Ecr[2,2] += rho_cr * ((z_values[idx]- z_mid) * (z_values[idx]- z_mid)) # / dx[idx] * dy[idx] * dz[idx]
        
        E_tot += rho_cr 

    # symmetric properties
    I_Ecr[1,0] = I_Ecr[0,1] 
    I_Ecr[2,0] = I_Ecr[0,2] 
    I_Ecr[2,1] = I_Ecr[1,2] 
    
    
    # Divide by total energy
    if normalize == True:
        E = cr_data * (dx[0] * 3.24078e-19) **3
        I_Ecr /= np.sum(E)

    # calculate the eigenvalues
    eigs, evecs = np.linalg.eig(I_Ecr)
    
    return I_Ecr, eigs, evecs, t
