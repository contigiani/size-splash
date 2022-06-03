import h5py
import numpy as np
import sim_tools as st
from pathlib import Path

cluster_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 21, 22, 24, 25, 28, 29]

#Load galaxy positions at z=0
# data/fullgalaxies/*pos
for i in cluster_indices:
    print(i)
    simdir = 'Hydrangea/CE-'+str(i)+'/HYDRO/highlev/'
    destdir = 'data/galaxies/'+str(i)



    my_file = Path(simdir + "GalaxyCoordinates10Myr.hdf5")
    if my_file.is_file():
        ff = h5py.File(simdir + "GalaxyCoordinates10Myr.hdf5")
        pos = np.array(ff['InterpolatedPositions'])

        np.save(destdir+"pos.npy", pos)

exit()



# Write cluster accretion rate in 
# data/snap/hydro/*Gdynp5.npy - for hydro sims
# data/snap/dm/*Gdynp5.npy  - for Dark matter only sims
for z_label, z_index, z_index2, z1, z2 in zip(["", "0.5", "1"],
                                              [29, 17, 11],
                                              [21, 13, 6],
                                              [0, 0.474, 1.017],
                                              [0.297, 0.795, 1.993]):
    for i in cluster_indices:
        simdir = 'Hydrangea/CE-'+str(i)+'/HYDRO'
        subdir = st.form_files(simdir, isnap = z_index2, types = 'sub')
        M200m_1, _, _ =  st.eagleread(subdir, 'FOF/Group_M_Mean200', astro=True)
        subdir = st.form_files(simdir, isnap = z_index, types = 'sub')
        M200m_2, _, _ =  st.eagleread(subdir, 'FOF/Group_M_Mean200', astro=True)


        Gdyn = np.log(M200m_2[0]/M200m_1[0])/(np.log(1/(1+z1))-np.log(1/(1+z2)))

        write_dir = "data/snap/hydro"+z_label+"/"+str(i)+''


        np.save(write_dir+'Gdyn', Gdyn)
exit()



# Load subhalo information at high-redshift
for z, z_label, z_index in zip([0.474, 1.017], ["0.5", "1"], [17, 11]):
    for i in cluster_indices:
        print(i)
        simdir = 'Hydrangea/CE-'+str(i)+'/HYDRO/'
        destdir = 'data/sub/hydro'+z_label+'/'+str(i)

        subdir = st.form_files(simdir, isnap = z_index, types = 'sub')

        m200m_FoF, _, _ = st.eagleread(subdir, 'FOF/Group_M_Mean200', astro=True)
        #PotCen_FoF, _, _ = st.eagleread(subdir, 'Subhalo/CentreOfMass', astro=True)
        PotCen_FoF, _, _ = st.eagleread(subdir, 'FOF/GroupCentreOfPotential', astro=True)
        R200m_FoF, _, _ = st.eagleread(subdir, 'FOF/Group_R_Mean200', astro=True)

        numOfSubhaloes = st.eagleread(subdir, 'FOF/NumOfSubhalos', astro = False)


        cluster_mass = m200m_FoF[0]
        cluster_centre = PotCen_FoF[0]
        cluster_radius = R200m_FoF[0]
        cluster_richness = numOfSubhaloes[0]

        np.save(destdir+'info', [cluster_mass, cluster_centre, cluster_radius, cluster_richness])
exit()




# Load subhalo positions and masses at high-redshift
for z_label, z_index in zip(["0.5", "1"], [17, 11]):
    for i in cluster_indices:
        print(i)
        simdir = 'Hydrangea/CE-'+str(i)+'/HYDRO/'
        destdir = 'data/sub/hydro'+z_label+'/'+str(i)

        subdir = st.form_files(simdir, isnap = z_index, types = 'sub')

        mass, _, _ = st.eagleread(subdir, 'Subhalo/Mass', astro = True)
        pos, _, _ = st.eagleread(subdir, 'Subhalo/CentreOfMass', astro = True)
        #vel, _, _ = st.eagleread(subdir, 'Subhalo/Velocity', astro = True)

        np.save(destdir+'mass', mass)
        np.save(destdir+'pos', pos)
        #np.save(destdir+'vel', vel)


exit()






#Load galaxy positions at high-redshift
for z_label, z_index in zip(["0.5", "1"], [-13, -19]):
    for i in cluster_indices:
        print(i)
        simdir = 'Hydrangea/CE-'+str(i)+'/HYDRO/highlev/'
        destdir = 'data/galaxies'+z_label+'/'+str(i)

        ff = h5py.File(simdir + "FullGalaxyTables.hdf5")

        np.save(destdir+"mstar.npy", ff['Mstar'][:, z_index])

        ff = h5py.File(simdir + "GalaxyPositionsSnap.hdf5")

        np.save(destdir+"pos.npy", ff['Centre'][:, z_index, :])
        np.save(destdir+"vel.npy", ff['Velocity'][:, z_index, :])

