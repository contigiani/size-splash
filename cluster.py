import h5py
import numpy as np
import profiles
import emcee
from scipy.interpolate import interp1d


# Hard coded information about the Hydrangea clusters, see tables in Bahe+ 2017
cluster_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 21, 22, 24, 25, 28, 29]
TXs = [1.81, 2.18, 1.90, 2.01, 1.79, 1.75, 2.28, 2.48, 2.03, 2.60, 3.00, 2.66, 3.28, 3.20, 3.99, 3.56, 3.06, 4.22, 5.09, 6.04, 5.31, 7.95, 6.29, 7.66]
LXs = [43.234, 43.459, 43.227, 43.114, 42.888, 43.240, 43.518, 43.657, 43.518, 43.620, 43.890, 44.080,
   43.942, 44.117, 44.321, 44.037, 43.877,
   44.402, 44.540, 44.618, 44.412, 45.021,
   44.407, 44.942]

cs = [5.3, 3.7, 6.1, 6.5, 4.2, 6.7, 3.6, 4.6, 4.5, 5.2, 4.8, 6.5, 4.7, 6.3, 2.5, 2.2, 7.0, 4.8, 3.3, 4.4, 5.0, 2.5, 3.7, 1.8]
R200ms = [1.74, 1.63, 1.63, 1.84, 1.89, 1.90, 2.16, 2.17, 2.12, 2.36, 2.21, 2.34, 2.49, 2.52, 2.66,2.73, 2.84, 3.03, 3.34, 3.72, 3.61, 3.87, 4.06, 4.61]
M200ms = [14.24, 14.15, 14.15, 14.31, 14.34, 14.35, 14.52, 14.53, 14.49, 14.63, 14.55, 14.63, 14.71, 14.72, 14.79, 14.83, 14.88, 14.96, 15.09, 15.23, 15.19, 15.28, 15.34, 15.51]



# Class that stores information about individual clusters and produces various profiles
class cluster:
    def __init__(self, idx, zlabel=''):
        self.zlabel = zlabel
        label='hydro'+zlabel

        self.idx = int(idx)
        self.CEidx = cluster_indices[self.idx]
        self.TX = TXs[self.idx]
        self.LX = LXs[self.idx]
        self.c = cs[self.idx]
        self.R200m = R200ms[self.idx]
        self.M200m = M200ms[self.idx]



        self.r_sp, self.M_sp = np.load("data/fit"+zlabel+"/"+str(self.CEidx)+"dm.npy")
        self.r_sp_gal, self.M_sp_gal, self.M_sp_gal_stellar = np.load("data/fit"+zlabel+"/"+str(self.CEidx)+"gal.npy")
        self.r_sp_sub, self.M_sp_sub = np.load("data/fit"+zlabel+"/"+str(self.CEidx)+"sub.npy")


        info, Gdyn = np.load('data/sub/'+label+'/'+str(self.CEidx)+'info.npy',allow_pickle=True), \
                    np.load('data/snap/'+label+'/'+str(self.CEidx)+'Gdyn.npy')

        if(zlabel!=''):
            self.R200m = info[2]
            self.M200m = info[0]
        self.Gdyn = Gdyn
        self.R200c = info[2]
        self.M200c = info[0]
        self.posc = info[1]
        # Best fit parameters after profile fit
        self.dmbfparams = np.load('data/fit'+zlabel+'/chains/'+str(self.CEidx)+"dmbfparams.npy")
        self.galbfparams = np.load('data/fit'+zlabel+'/chains/'+str(self.CEidx)+"galbfparams.npy")
        self.subbfparams = np.load('data/fit'+zlabel+'/chains/'+str(self.CEidx)+"subbfparams.npy")
        return



    def get_profile(self, bins=np.linspace(0.5, 10, 30), sub=False, label='hydro', masscut=-1.):
        # 3d profile

        posc = self.posc
        rbins = ((bins[1:]**3. + bins[:-1]**3.) / 2.)**(1./3.)
        r_min = bins.min()
        r_max = bins.max()

        if(sub):
            mass, pos = np.load('data/sub/'+label+self.zlabel+'/'+str(self.CEidx)+'mass.npy'), \
                            np.load('data/sub/'+label+self.zlabel+'/'+str(self.CEidx)+'pos.npy')
        else:
            mstar, pos = np.load('data/galaxies'+self.zlabel+'/'+str(self.CEidx)+'mstar.npy'), \
                        np.load('data/galaxies'+self.zlabel+'/'+str(self.CEidx)+'pos.npy')



        r = np.sqrt( ( (pos-posc)**2.).sum(axis=1))

        if(sub):
            r = r[(r>r_min) & (r<r_max) & (np.log10(mass)>masscut)]
        else:
            r = r[(r>r_min) & (r<r_max)& (mstar>8.)]

        volume = 4./3. * np.pi *(bins[1:]**3. - bins[:-1]**3.)


        Sigma = np.histogram(r, bins=bins)[0]/volume

        return rbins, Sigma


    def get_ps(self, psrbins= np.linspace(np.log10(1.), np.log10(10), 100), psvbins= np.linspace(-1000., 1000., 100), sub=True):
        #Phase space distribution (r, v_r)

        label = 'hydro'

        velc = np.load('data/sub/hydro'+self.zlabel+'/'+str(self.CEidx)+'vel.npy')
        posc = self.posc
        velc = np.average(velc, axis=0)


        if(sub):
            mass, pos, vel = np.load('data/sub/hydro'+self.zlabel+'/'+str(self.CEidx)+'mass.npy'), \
                            np.load('data/sub/hydro'+self.zlabel+'/'+str(self.CEidx)+'pos.npy'), \
                            np.load('data/sub/hydro'+self.zlabel+'/'+str(self.CEidx)+'vel.npy')
        else:
            mstar, pos, vel = np.load('data/galaxies'+self.zlabel+'/'+str(self.CEidx)+'mstar.npy'), \
                        np.load('data/galaxies'+self.zlabel+'/'+str(self.CEidx)+'pos.npy'), \
                        np.load('data/galaxies'+self.zlabel+'/'+str(self.CEidx)+'vel.npy')


        r = np.sqrt( ( (pos-posc)**2.).sum(axis=1))
        vr = np.sum( (vel-velc) * (pos-posc), axis=1)/r

        if(sub):
            vr = vr[(np.log10(mass)>-1)]
            r = r[(np.log10(mass)>-1)]
        else:
            vr = vr[(mstar>8.)]
            r = r[(mstar>8.)]


        hist = np.histogram2d(np.log10(r), vr, bins=[psrbins, psvbins], density=True)[0]

        return hist, r, vr





    def get_2d_subs(self, Rmin, Rmax, dir_idx=0):
        # return subhalo positions in a 2d plane, dir_idx specifies which plane to use.

        label = 'hydro'+self.zlabel
        mass, pos = np.load('data/sub/'+label+'/'+str(self.CEidx)+'mass.npy'), \
                        np.load('data/sub/'+label+'/'+str(self.CEidx)+'pos.npy')

        if(dir_idx==0):
            dir_idxs = [1,2]
        if(dir_idx==1):
            dir_idxs = [0,2]
        if(dir_idx==2):
            dir_idxs = [0,1]

        posc = self.posc
        pos = (pos-posc)

        x1, x2 = pos[:, dir_idxs].T
        R = np.sqrt(x1**2.+x2**2.)

        idx = (R>Rmin) & (R<Rmax) & (np.log10(mass)>-1)
        return x1[idx], x2[idx], mass[idx]

    def get_2d_stars(self, Rmin=0, Rmax=0.03, dir_idx=0):
        # return star positions in a 2d plane, dir_idx specifies which plane to use.

        X = np.load('data/stars'+self.zlabel+'/'+str(self.CEidx)+'coord.npy')
        M = np.load('data/stars'+self.zlabel+'/'+str(self.CEidx)+'mass.npy')

        if(dir_idx==0):
            dir_idxs = [1,2]
        if(dir_idx==1):
            dir_idxs = [0,2]
        if(dir_idx==2):
            dir_idxs = [0,1]

        x1, x2 = (X[:, dir_idxs] - self.posc[dir_idxs]).T
        R = np.sqrt(x1**2.+x2**2.)

        idx = (R<Rmax) & (R>Rmin)

        return x1[idx], x2[idx], M[idx]

    def get_rotation(self, dir_idx=0, type=''):
        # return rotated positions according to some definition of rotation (defined by type)
        r200 = self.R200m

        theta = np.load('data/direction'+self.zlabel+'/'+str(self.CEidx)+str(dir_idx)+type+'.npy')

        c, s = np.cos(theta), np.sin(theta)
        x1, x2, _ = self.get_2d_subs(0, r200*10., dir_idx=dir_idx)

        return x1*c-x2*s, x1*s+x2*c

    def get_rotated_dist2d(self, dir_idx, type=0):
        # Use rotated subhalo positions to get a 2d histogram of the density distribution

        x, y = self.get_rotation(dir_idx, type) # Rotation is along y axis, swap for comfort

        hist, xe, ye = np.histogram2d(y/self.R200m, x/self.R200m, bins=np.linspace(-5, 5, 300))

        return hist, xe, ye
