import h5py
import numpy as np
import sim_tools as st
import glob
import os.path
cluster_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 21, 22, 24, 25, 28, 29]
# Hubble constant
h = 0.6777


'''
Example script used to generate a 3D mass profile as a function of radius r

'''


binLims = np.geomspace(1e-3, 20, 400)
binLims = np.concatenate([[0], binLims])

for CEidx in cluster_indices:
    for z, z_index, zlabel, correction in zip([0.474, 1.017], [17, 11], ["snapshot_017_z000p474", "snapshot_011_z001p017"], [0., 0.]):
        #z = 0.
        #z_index = 29
        #zlabel = 'snapshot_029_z000p000'

        zdir = 'data/'+zlabel+'/'
        a = 1./(1.+z)
        volume =4./3.*np.pi*(binLims[1:]**3.-binLims[:-1]**3.)


        MassProfile3D = np.zeros([len(binLims)-1, 4])
        DensityProfile3D = np.copy(MassProfile3D)

        simdir = 'Hydrangea/CE-'+str(CEidx)+'/HYDRO/'
        subdir = st.form_files(simdir, isnap = z_index, types = 'sub')
        PotCen_FoF, _, _ = st.eagleread(subdir, 'FOF/GroupCentreOfPotential', astro=True)
        #COM_SUB, _, _ = st.eagleread(subdir, 'Subhalo/CentreOfMass', astro=True)
        cluster_centre = PotCen_FoF[0] # in Physical coordinates
        cluster_centre = cluster_centre + correction
        simdir = simdir+zdir

        for filename in glob.glob(simdir+"snap*.hdf5"):
            print(filename)

            f = h5py.File(filename)

            for j, pt in enumerate([0, 1, 4, 5]):

                X = np.array(f['PartType'+str(pt)+'/Coordinates'])/h*a - cluster_centre
                r = np.sqrt((X**2.).sum(axis=1))


                if(j == 1):
                    m = f['Header'].attrs['MassTable'][1]*1e10/h
                    temph, _ =  np.histogram(r, bins=binLims)
                    temph = temph*m
                else:
                    m = np.array(f['PartType'+str(pt)+'/Mass'])*1e10/h
                    temph, _ = np.histogram(r, bins=binLims, weights=m)
                MassProfile3D[:, j] += np.cumsum(temph)
                DensityProfile3D[:, j] += temph/volume

        np.save("cache/massprofile"+str(CEidx)+str(z), MassProfile3D)
        np.save("cache/densityprofile"+str(CEidx)+str(z), DensityProfile3D)
        np.save("cache/bins"+str(CEidx)+str(z), binLims)
