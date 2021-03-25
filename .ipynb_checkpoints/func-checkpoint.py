import gizmo_analysis as gizmo
import utilities as ut
import numpy as np
import matplotlib.pyplot as plt
import copy
import os
import glob


from IPython.utils import io

def x(ls):
    return ls['position'][:,0]

def y(ls):
    return ls['position'][:,1]

def z(ls):
    return ls['position'][:,2]

def vx(ls):
    return ls['velocity'][:,0]

def vy(ls):
    return ls['velocity'][:,1]

def vz(ls): 
    return ls['velocity'][:,2]

def Rgc(ls) :
    return np.sqrt(x(ls)*x(ls) + y(ls)*y(ls) + z(ls)*z(ls))

def R(ls) :
    return np.sqrt(x(ls)*x(ls) + y(ls)*y(ls))

def v(ls) :
    return np.sqrt(vx(ls)*vx(ls) + vy(ls)*vy(ls) + vz(ls)*vz(ls))

def vR(ls) :
    vR1 = vx(ls)*x(ls)/R(ls)
    vR1 += vy(ls)*y(ls)/R(ls)
    return vR1

def vphi(ls) :
    vphi1 = -vx(ls)*y(ls)/R(ls)
    vphi1 += vy(ls)*x(ls)/R(ls)
    return vphi1

def vr(ls) :
    vr1 = vR(ls)*R(ls)/Rgc(ls)
    vr1 += vz(ls)*z(ls)/Rgc(ls)
    return vr1

def vtheta(ls) :
    vtheta1= vR(ls)*z(ls)/Rgc(ls)
    vtheta1+= -vz(ls)*R(ls)/Rgc(ls)
    return vtheta1

def create_galaxy_arrays(simulation_directory, n_dm1, n_star1, n_dm2, n_star2, time_array):

    gal1=[]
    gal2=[]

    in_output_directory = glob.glob(simulation_directory + "/output/snapshot_*.hdf5")
    
    with io.capture_output() as captured:  #suppresses output

        for i in range(len(in_output_directory)-1):
            

            #Read in the particles from the snapshot
            
            #------some changes made to the package------#
            #assign_hosts_rotation: set distance_max=5 and mass_percent=90 and age_percent=None
            #in gizmo_io ReadClass, make sure to associate ID 2 with 'star' not 'dark2'
            #in utilities/basic/particle, make sure to replace part.snapshot['time.hubble'] with 1.0 when calculating velocity differences
            #-------------------------------------------#
            
            part = gizmo.io.Read.read_snapshots(['dark', 'star'], 'time', time_array[i], simulation_directory,
                                                properties=['position', 'mass', 'id', 'velocity', 'potential'],
                                                assign_hosts=True,  
                                                assign_hosts_rotation=True, 
                                                assign_orbits=True)

            #assigns orbital properties to the dark and star particles
            #positions and velocities relative to the COM and COV can now be found using host.distance, host.velocity
            gizmo.io.Read.assign_orbits(part, species=['dark', 'star'])

            #Rotate position and velocity coordinates, cartesian coordinates
            #When calculating the moment of inertia, it just looks at the particles within 5 kpc of the center
            pos_dark = ut.coordinate.get_coordinates_rotated(part['dark']['host.distance'], part.host['rotation'][0]) 
            vel_dark = ut.coordinate.get_coordinates_rotated(part['dark']['host.velocity'], part.host['rotation'][0]) 

            pos_star = ut.coordinate.get_coordinates_rotated(part['star']['host.distance'], part.host['rotation'][0])
            vel_star = ut.coordinate.get_coordinates_rotated(part['star']['host.velocity'], part.host['rotation'][0])

            #Transform from cartesian to cylindrical coordinates
            pos_cyl_dark = ut.coordinate.get_positions_in_coordinate_system(pos_dark, system_from='cartesian', system_to='cylindrical')
            vel_cyl_dark = ut.coordinate.get_velocities_in_coordinate_system(vel_dark, pos_dark, system_from='cartesian', system_to='cylindrical')

            pos_cyl_star = ut.coordinate.get_positions_in_coordinate_system(pos_star, system_from='cartesian', system_to='cylindrical')
            vel_cyl_star = ut.coordinate.get_velocities_in_coordinate_system(vel_star, pos_star, system_from='cartesian', system_to='cylindrical')

            #Transform from cartesian to spherical coordinates
            pos_sph_dark = ut.coordinate.get_positions_in_coordinate_system(pos_dark, system_from='cartesian', system_to='spherical')
            vel_sph_dark = ut.coordinate.get_velocities_in_coordinate_system(vel_dark, pos_dark, system_from='cartesian', system_to='spherical')

            pos_sph_star = ut.coordinate.get_positions_in_coordinate_system(pos_star, system_from='cartesian', system_to='spherical')
            vel_sph_star = ut.coordinate.get_velocities_in_coordinate_system(vel_star, pos_star, system_from='cartesian', system_to='spherical')
            
            #we will need to divide up the dark/star particles between galaxy 1 and galaxy 2, as per true origin
            gal1_dict = copy.deepcopy(part)
            gal2_dict = copy.deepcopy(part)
            
            for key in ['host.distance', 'host.distance.total', 'host.distance.norm', 'host.velocity', 'host.velocity.total',
                       'host.velocity.tan', 'host.velocity.rad', 'host.velocity.norm', 'host.velocity.ratio']:
                del gal1_dict['dark'][key]
                del gal1_dict['star'][key]
                del gal2_dict['dark'][key]
                del gal2_dict['star'][key]
                             
            #the following is the selection criteria for each
            haloID_DM = part['dark']['id'] <= n_dm1
            satID_DM = part['dark']['id'] > n_dm1 + n_star1

            haloID_star = part['star']['id'] <= n_dm1 + n_star1
            satID_star = part['star']['id'] > n_dm1 + n_star1 + n_dm2

            for GAL_DICT in [gal1_dict, gal2_dict]:

                if GAL_DICT == gal1_dict:
                    ID_DM = haloID_DM
                    ID_STAR = haloID_star
                    flag1 = True
                elif GAL_DICT == gal2_dict:
                    ID_DM = satID_DM
                    ID_STAR = satID_star
                    flag1 = False

                #Dark Matter Halo

                GAL_DICT['dark']['position'] = pos_dark[ID_DM]
                GAL_DICT['dark']['velocity'] = vel_dark[ID_DM]
                GAL_DICT['dark']['mass'] = part['dark']['mass'][ID_DM]
                GAL_DICT['dark']['id'] = part['dark']['id'][ID_DM]
                GAL_DICT['dark']['position.cyl'] = pos_cyl_dark[ID_DM]
                GAL_DICT['dark']['velocity.cyl'] = vel_cyl_dark[ID_DM]
                GAL_DICT['dark']['position.sph'] = pos_sph_dark[ID_DM]
                GAL_DICT['dark']['velocity.sph'] = vel_sph_dark[ID_DM]
                GAL_DICT['dark']['potential'] = part['dark']['potential'][ID_DM]

                #Stellar Disk

                GAL_DICT['star']['position'] = pos_star[ID_STAR]
                GAL_DICT['star']['velocity'] = vel_star[ID_STAR]
                GAL_DICT['star']['mass'] = part['star']['mass'][ID_STAR]
                GAL_DICT['star']['id'] = part['star']['id'][ID_STAR]
                GAL_DICT['star']['position.cyl'] = pos_cyl_star[ID_STAR]
                GAL_DICT['star']['velocity.cyl'] = vel_cyl_star[ID_STAR]
                GAL_DICT['star']['position.sph'] = pos_sph_star[ID_STAR]
                GAL_DICT['star']['velocity.sph'] = vel_sph_star[ID_STAR]
                GAL_DICT['star']['potential'] = part['star']['potential'][ID_STAR]

                if flag1 == True:  #if gal1_dict
                    gal1.append(GAL_DICT)
                elif flag1 == False: #if gal2_dict
                    gal2.append(GAL_DICT)   
                    
    return gal1, gal2
    
    
def apply_mask(ls, func, func_min_val, func_max_val):

    newls=[]

    for i in range(len(ls)):
    
        #mask = (func(ls[i]['dark'])<func_max_val) and (func(ls[i]['dark'])>func_max_val)
        mask_d = [(x > func_min_val and x < func_max_val) for x in func(ls[i]['dark'])]
        cut1_d = {k: v[mask_d] for k, v in ls[i]['dark'].items()}
        
        mask_s = [(x > func_min_val and x < func_max_val) for x in func(ls[i]['star'])]
        cut1_s = {k: v[mask_s] for k, v in ls[i]['star'].items()}
        
        
        cut2 = {'dark': cut1_d, 'star': cut1_s}
    
        newls.append(cut2)
    return newls