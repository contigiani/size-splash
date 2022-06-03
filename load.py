import sim_tools as st
import numpy as np
from pdb import set_trace


cluster_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 21, 22, 24, 25, 28, 29]



# Write mass and position of stellar particles in 
# data/stars/*coord.npy
# data/stars/*mass.npy 
for i in cluster_indices:

    labelr ='/HYDRO/'


    simdir = '/net/quasar/data3/Hydrangea/CE-'+str(i)+labelr+''
    subdir = st.form_files(simdir, isnap = 29, types = 'sub')
    PotCen_FoF, _, _ = st.eagleread(subdir, 'FOF/GroupCentreOfPotential', astro=True)
    cluster_center = PotCen_FoF[0]



    write_dir = "/net/zwin/data2/contigiani/splashletter/data/stars/"+str(i)+''
    simdir = '/net/quasar/data3/Hydrangea/CE-'+str(i)+labelr+''
    subdir = st.form_files(simdir, isnap = 29, types = 'snap')



    loadpaths = ['PartType4/Coordinates','PartType4/Mass']

    X, _, _ = st.eagleread(subdir, 'PartType4/Coordinates', astro=True)
    M, _, _ = st.eagleread(subdir, 'PartType4/Mass', astro=True)

    idx = np.sqrt(np.sum(np.abs(X - cluster_center)**2,axis=1)) < 0.3

    np.save(write_dir+"coord", X[idx])
    np.save(write_dir+"mass", M[idx])



# Write particle position, mass and velocity for the different particle types 
# data/snap/hydro/*pos, vel, mass [1-4]
# data/snap/dm/*pos, vel, mass [1-4]
for i in cluster_indices:
    for labelw, labelr in zip(['hydro', 'dm'], ['/HYDRO/', '/DM/']):

        simdir = '/net/quasar/data3/Hydrangea/CE-'+str(i)+labelr+''
        subdir = st.form_files(simdir, isnap = 29, types = 'sub')
        PotCen_FoF, _, _ = st.eagleread(subdir, 'FOF/GroupCentreOfPotential', astro=True)
        R200c_FoF, _, _ = st.eagleread(subdir, 'FOF/Group_R_Mean200', astro=True)
        cluster_center = PotCen_FoF[0]
        cluster_radius = R200m_FoF[0]

        write_dir = "data/snap/"+labelw+"/"+str(i)+''
        simdir = '/net/quasar/data3/Hydrangea/CE-'+str(i)+labelr+''
        subdir = st.form_files(simdir, isnap = 29, types = 'snap')


        try:
            m_dm = st.m_dm(subdir, issnapdir=True)
        except OSError:
            continue
        except TypeError:
            continue

        np.save(write_dir+"mdm", m_dm)


        loadpaths = ['PartType0/Coordinates', 'PartType0/Velocity', \
                      'PartType0/Mass', \
                      'PartType4/Coordinates', 'PartType4/Velocity', \
                      'PartType4/Mass', \
                      'PartType5/Coordinates', 'PartType5/Velocity', \
                      'PartType5/Mass', \
                      'PartType1/Coordinates', 'PartType1/Velocity']
        writepaths = ['pos0', 'vel0', 'mass0', 'pos4', 'vel4', 'mass4', 'pos5', 'vel5', 'mass5' 'pos1', 'vel1']

        for j, [loadpath, writepath] in enumerate(zip(loadpaths, writepaths)):
            buff, _, _ = st.eagleread(subdir, loadpath, astro=True)
            if(j%3==0):
                idx = abs(buff - cluster_center)/cluster_radius < 3.
            np.save(write_dir+writepath, buff[idx])

exit()



# Write subhalo information in
# data/sub/hydro/*mass, pos, vel, info
# data/sub/dm/*mass, pos, vel, info
def return_data(subdir):
    m200m_FoF, _, _ = st.eagleread(subdir, 'FOF/Group_M_Mean200', astro=True)
    PotCen_FoF, _, _ = st.eagleread(subdir, 'FOF/GroupCentreOfPotential', astro=True)
    R200m_FoF, _, _ = st.eagleread(subdir, 'FOF/Group_R_Mean200', astro=True)


    #firstSubhalo = st.eagleread(subdir, 'FOF/FirstSubhaloID', astro = False)
    numOfSubhaloes = st.eagleread(subdir, 'FOF/NumOfSubhalos', astro = False)
    #idxs = np.arange(firstSubhalo[0], firstSubhalo[0]+numOfSubhaloes[0])


    cluster_mass = m200m_FoF[0]
    cluster_centre = PotCen_FoF[0]
    cluster_radius = R200m_FoF[0]
    cluster_richness = numOfSubhaloes[0]

    mass, _, _ = st.eagleread(subdir, 'Subhalo/Mass', astro = True)
    pos, _, _ = st.eagleread(subdir, 'Subhalo/CentreOfMass', astro = True)
    vel, _, _ = st.eagleread(subdir, 'Subhalo/Velocity', astro = True)

    #return mass[idxs], pos[idxs], vel[idxs], cluster_mass, cluster_centre
    return mass, pos, vel, cluster_mass, cluster_centre, cluster_radius, cluster_richness


cluster_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 21, 22, 24, 25, 28, 29]
for i in cluster_indices:
    for labelw, labelr in zip(['hydro, dm'], ['/HYDRO/', '/DM/']):
        write_dir = "data/sub/"+labelw+"/"+str(i)+''
        simdir = '/net/quasar/data3/Hydrangea/CE-'+str(i)+labelr+''
        subdir = st.form_files(simdir, isnap = 29, types = 'sub')


        mass, pos, vel, cluster_mass, cluster_centre, cluster_radius, cluster_richness = return_data(subdir)

        np.save(write_dir+'mass', mass)
        np.save(write_dir+'pos', pos)
        np.save(write_dir+'vel', vel)
        np.save(write_dir+'info', [cluster_mass, cluster_centre, cluster_radius, cluster_richness])
