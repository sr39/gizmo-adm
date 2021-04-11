################################################################################
###### This is an example script to generate HDF5-format ICs for GIZMO
######  The specific example below is obviously arbitrary, but could be generalized
######  to whatever IC you need.
################################################################################
################################################################################

## load libraries we will use
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py as h5py

# the main routine. this specific example builds an N-dimensional box of gas plus
#   a collisionless particle species, with a specified mass ratio. the initial
#   gas particles are distributed in a uniform lattice; the initial collisionless
#   particles laid down randomly according to a uniform probability distribution
#   with a specified random velocity dispersion


gamma_eos = 5./3. # polytropic index of ideal equation of state the run will assume
rd = 50 # kpc (size of gas cloud)
G = 4.3 * 10**4 # gravitational constant in code units
rho0 = 7.9*(10**(-4)) # NFW DM halo density parameter
rs = 15 # kpc, scale length of DM halo (NFW)
mg = 1/20.0 # mass fraction of gas to halo
M200 = 100 # 10^10 Msol (the virial mass of the dm halo)
N = 20000 # Initial number of particles randomly thrown into the cubic box
kt = 2.083 * 10**(-57) # In energy units, corresponding to temp of 30,000K
white_noise_factor = (1/10.)
velocity_factor = 0.9
ax = 1.2 # x rescaling
ay = 1.2 # y-rescaling
az = 0.8 # z-rescaling
random_seed = 1	# seed to initialise random assignment

def M(r):
    """
    Function that returns the NFW halo mass given the global parameters specified above.
    Input: radial distance in kpc
    Output: halo mass contained in radius r (in 10^10 Msol)
    """
    return 4 * np.pi * (rs**3) * rho0 * ( -r/(r + rs) + np.log(1 + r/rs) )


def make_IC():
    '''
    This is an example subroutine provided to demonstrate how to make HDF5-format
    ICs for GIZMO. The specific example here is arbitrary, but can be generalized
    to whatever IC you need
    '''
    fname='/tigress/sr39/gizmo_adm_testing/gas_cloud_elliptical_ADM_IC.hdf5'; # output filename

    #####################
    # initial positions #
    #####################
    # I randomly throw particle in rd^3 box. Then I select those that lie in the sphere
    qinit = np.random.rand(N,3)*2*rd - np.ones((N,3))*rd
    q = []
    for qi in qinit[:]:
        if np.dot(qi,qi) <= rd**2:
            q.append(qi)
    q = np.asarray(q)

    # The final number of gas particles
    Ngas = q.shape[0]

    # Rescaling the positions into an ellipsoid
    q[:,0] = ax*q[:,0]
    q[:,1] = ay*q[:,1]
    q[:,2] = az*q[:,2]

    """
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.scatter(q[:,0],q[:,1],q[:,2],marker='.',s=1)
    plt.show()
    """

    #######################
    # initial velocities
    #######################
    # Here, I'll give them circular velocities according to M(r) with some gaussian white noise
    # The particle orbits will be going through the origin with upward pointing angular momentum vector,
    # so the angular momenta will be random but will, on average, point upwards
    # The white noise is proportional to the velocities of the particles
    x = q[:,0]
    y = q[:,1]
    z = q[:,2]
    r = np.sqrt(x**2+y**2+z**2)
    rmatrix = np.zeros((Ngas,3))
    rmatrix[:,0] = r; rmatrix[:,1] = r; rmatrix[:,2] = r
    s = np.sqrt(x**2+y**2)

    # Angular momentum vector
    lx = (1/r)*(-np.sqrt(r**2 - s**2)*(x/s)*(z/np.abs(z)))
    ly = (1/r)*(-np.sqrt(r**2 - s**2)*(y/s)*(z/np.abs(z)))
    lz = (s/r)

    l = np.zeros((Ngas,3))
    l[:,0] = lx
    l[:,1] = ly
    l[:,2] = lz

    # v = omega cross r
    velocity_n = np.cross(l,q/rmatrix,axis=1)

    # Finally, giving the velocities the right magnitude and some white noise
    vvals = np.sqrt(G* M(r)/r)
    vmatrix = np.zeros((Ngas,3))
    vmatrix[:,0] = vvals; vmatrix[:,1] = vvals; vmatrix[:,2] = vvals;

    v = velocity_factor*vmatrix * (velocity_n + np.random.rand(Ngas,3)*white_noise_factor)


    ############################
    # Other initial conditions #
    ############################
    gamma_eos = 5./3. # polytropic index of ideal equation of state the run will assume
    # set the initial magnetic field in x/y/z directions (here zero).
    #  these can be overridden (if constant field values are desired) by BiniX/Y/Z in the parameterfile
    bx_g=0.*q[:,0]; by_g=0.*q[:,0]; bz_g=0.*q[:,0];
    # set the particle masses. Here we set it to be a list the same length, with all the same mass
    #   since their space-density is uniform this gives a uniform density, at the desired value
    mv_g=mg*M200/(Ngas) + 0.*q[:,0]
    # set the initial internal energy per unit mass. recall gizmo uses this as the initial 'temperature' variable
    #  this can be overridden with the InitGasTemp variable (which takes an actual temperature)
    uv_g=(1/(gamma_eos-1.))*(kt/mv_g[0]) + 0.*q[:,0]
    # set the gas IDs: here a simple integer list
    id_g=np.arange(1,Ngas+1)

    # Set the ADM gas types here. I'll set all the gas particles to ADM particles. Need to have ADM flag enabled on GIZMO
    id_adm = np.ones(Ngas,int)

    # now we get ready to actually write this out
    #  first - open the hdf5 ics file, with the desired filename
    file = h5py.File(fname,'w')

    # set particle number of each type into the 'npart' vector
    #  NOTE: this MUST MATCH the actual particle numbers assigned to each type, i.e.
    #   npart = np.array([number_of_PartType0_particles,number_of_PartType1_particles,number_of_PartType2_particles,
    #                     number_of_PartType3_particles,number_of_PartType4_particles,number_of_PartType5_particles])
    #   or else the code simply cannot read the IC file correctly!
    #
    npart = np.array([Ngas,0,0,0,0,0]) # we have gas and particles we will set for type 3 here, zero for all others

    # now we make the Header - the formatting here is peculiar, for historical (GADGET-compatibility) reasons
    h = file.create_group("Header");
    # here we set all the basic numbers that go into the header
    # (most of these will be written over anyways if it's an IC file; the only thing we actually *need* to be 'correct' is "npart")
    h.attrs['NumPart_ThisFile'] = npart; # npart set as above - this in general should be the same as NumPart_Total, it only differs
                                         #  if we make a multi-part IC file. with this simple script, we aren't equipped to do that.
    h.attrs['NumPart_Total'] = npart; # npart set as above
    h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
    h.attrs['MassTable'] = np.zeros(6); # these can be set if all particles will have constant masses for the entire run. however since
                                        # we set masses explicitly by-particle this should be zero. that is more flexible anyways, as it
                                        # allows for physics which can change particle masses
    ## all of the parameters below will be overwritten by whatever is set in the run-time parameterfile if
    ##   this file is read in as an IC file, so their values are irrelevant. they are only important if you treat this as a snapshot
    ##   for restarting. Which you shouldn't - it requires many more fields be set. But we still need to set some values for the code to read
    h.attrs['Time'] = 0.0;  # initial time
    h.attrs['Redshift'] = 0.0; # initial redshift
    h.attrs['BoxSize'] = 1.0; # box size
    h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
    h.attrs['Omega0'] = 1.0; # z=0 Omega_matter
    h.attrs['OmegaLambda'] = 0.0; # z=0 Omega_Lambda
    h.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
    h.attrs['Flag_Sfr'] = 0; # flag indicating whether star formation is on or off
    h.attrs['Flag_Cooling'] = 0; # flag indicating whether cooling is on or off
    h.attrs['Flag_ADM'] = 0; # flag indicating whether ADM physics is enabled or not
    h.attrs['Flag_StellarAge'] = 0; # flag indicating whether stellar ages are to be saved
    h.attrs['Flag_Metals'] = 0; # flag indicating whether metallicity are to be saved
    h.attrs['Flag_Feedback'] = 0; # flag indicating whether some parts of springel-hernquist model are active
    h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
    h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs
    ## ok, that ends the block of 'useless' parameters

    # Now, the actual data!
    #   These blocks should all be written in the order of their particle type (0,1,2,3,4,5)
    #   If there are no particles of a given type, nothing is needed (no block at all)
    #   PartType0 is 'special' as gas. All other PartTypes take the same, more limited set of information in their ICs

    # start with particle type zero. first (assuming we have any gas particles) create the group
    p = file.create_group("PartType0")
    # write it to the 'Coordinates' block
    p.create_dataset("Coordinates",data=q)
    # write it to the 'Velocities' block
    p.create_dataset("Velocities",data=v)
    # write particle ids to the ParticleIDs block
    p.create_dataset("ParticleIDs",data=id_g)
    # write particle masses to the Masses block
    p.create_dataset("Masses",data=mv_g)
    # write internal energies to the InternalEnergy block
    p.create_dataset("InternalEnergy",data=uv_g)
    # combine the xyz magnetic fields into a matrix with the correct format
    b=np.zeros((Ngas,3)); b[:,0]=bx_g; b[:,1]=by_g; b[:,2]=bz_g;
    # write magnetic fields to the MagneticField block. note that this is unnecessary if the code is compiled with
    #   MAGNETIC off. however, it is not a problem to have the field there, even if MAGNETIC is off, so you can
    #   always include it with some dummy values and then use the IC for either case
    p.create_dataset("MagneticField",data=b)

    # write adm particle type to dataset
    p.create_dataset("ADMType",data=id_adm)

    # no PartType1 for this IC
    # no PartType2 for this IC

    # now assign the collisionless particles to PartType3. note that this block looks exactly like
    #   what we had above for the gas. EXCEPT there are no "InternalEnergy" or "MagneticField" fields (for
    #   obvious reasons).
    """
    p = file.create_group("PartType3")
    q=np.zeros((Ngrains,3)); q[:,0]=xv_d; q[:,1]=yv_d; q[:,2]=zv_d;
    p.create_dataset("Coordinates",data=q)
    q=np.zeros((Ngrains,3)); q[:,0]=vx_d; q[:,1]=vy_d; q[:,2]=vz_d;
    p.create_dataset("Velocities",data=q)
    p.create_dataset("ParticleIDs",data=id_d)
    p.create_dataset("Masses",data=mv_d)
    p.create_dataset("PICParticleType",data=type_d)
    """
    # no PartType4 for this IC
    # no PartType5 for this IC

    # close the HDF5 file, which saves these outputs
    file.close()
    # all done!

make_IC()
