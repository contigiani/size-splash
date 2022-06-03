from astropy import units as u
import numpy as np
import scipy.integrate as integrate
from scipy.integrate import simps
from astropy import constants as const

'''
    NFW profile, hardcoded cosmology
'''
def concentration(Mvir, z):
    # Mass-concentration from Duffy+ 2008
    return 7.85*( Mvir.to(u.Msun).value/2./1e12 )**(-0.081) * (1.+z)**(-0.71)

def rvir(Mvir, z):
    #Virial overdensity from Bryan&Norman 1998, hardcoded cosmology
    Om = 0.3
    OL = 1-Om
    H0 = 67*u.km/u.s/u.Mpc

    O_z = Om*(1.+z)**3./(Om*(1.+z)**3. + OL)
    Delta = (18*np.pi*np.pi + 82*(O_z-1.)-39*(O_z-1.)**2.)/O_z

    rho_bg = 3. * H0**2. * Om * (1+z)**3. / 8. / np.pi / const.G
    return  (3*Mvir/4./np.pi/Delta/rho_bg )**(1./3.)


def Mvir_to_M200m(Mvir, z, also_r=False, c=None):
    #Finds M200m starting from Mvir, hardcoded cosmology
    Om = 0.3
    OL = 1-Om
    H0 = 70*u.km/u.s/u.Mpc

    if(c is None):
        c = concentration(Mvir, z)
    rv = rvir(Mvir, z)
    rs = rv/c

    f = np.log(1.+c) - c/(1.+c)

    r = np.linspace(0, 20, 1001)*u.Mpc
    M_inside = Mvir/f * (np.log(1.+r/rs) - r/(rs+r) )

    rho_bg = 3. * H0**2. * Om * (1+z)**3. / 8. / np.pi / const.G
    # Finds r such that M_vir == 4pi/3 * 200rho_bg * r^3
    idx = (np.abs(Mvir - 4.*np.pi/3. * 200. * rho_bg * r**3.) ).argmin()
    if(also_r):
        return M_inside[idx], r[idx] # returns r_200m too!
    else:
        return M_inside[idx]


class NFW:
    '''
        Simple class for NFW profile, no astropy.units for faster computation.
        Everything returned is in combinations of Msun and Mpc.
        arXiv:astro-ph/9908213v1 for lensing expressions
    '''
    def __init__(self, Mvir, z, c=None):
        # computes concentration from mass-concentration relation if not specified
        if(c is None):
            c = concentration(Mvir, z)
        rv = rvir(Mvir, z)
        rs = rv/c
        f = np.log(1.+c) - c/(1.+c)

        self.rs = rs.to('Mpc')
        self.rv = rv.to('Mpc')
        self.c = c
        self.rs_dm = rs.to('Mpc').value
        self.prefactor = (Mvir/f/4./np.pi).to('Msun').value

    def rho(self, r):
        # Msun/Mpc/Mpc/Mpc
        rs = self.rs_dm
        return self.prefactor/r/(r+rs)**2.

    def D_Sigma(self, r):
        # Msun/Mpc/Mpc
        rs = self.rs_dm
        x =r/rs

        result = np.zeros(r.size)

        result[r<rs] =\
            8 * np.arctanh( np.sqrt( (1.-x[r<rs])/(1.+x[r<rs]) ) )/ x[r<rs]**2. / np.sqrt(1.-x[r<rs]**2.) +\
            4./x[r<rs]**2. * np.log(x[r<rs]/2.) - 2./(x[r<rs]**2. -1.) +\
            4 * np.arctanh( np.sqrt( (1.-x[r<rs])/(1.+x[r<rs]) ) )/(x[r<rs]**2. -1.)/np.sqrt(1.-x[r<rs]**2.)
        result[r==rs] = (10./3. + 4.*np.log(0.5))
        result[r>rs] =\
            8 * np.arctan( np.sqrt( (x[r>rs]-1.)/(1.+x[r>rs]) ) )/x[r>rs]**2./np.sqrt(x[r>rs]**2.-1.) +\
            4./x[r>rs]**2. * np.log(x[r>rs]/2.) - 2./(x[r>rs]**2. -1.) +\
            4 * np.arctan( np.sqrt( (x[r>rs]-1.)/(1.+x[r>rs]) ) )/(x[r>rs]**2.-1.)**(3./2.)

        return self.prefactor/self.rs_dm/self.rs_dm * result


    def Sigma(self, r):
        # Msun/Mpc/Mpc
        rs = self.rs_dm
        x =r/rs

        result = np.zeros(r.size)

        result[r<rs] =\
            2./(x[r<rs]**2.-1.)*(1.-2./np.sqrt(1.-x[r<rs]**2.)*np.arctanh( np.sqrt((1.-x[r<rs])/(1.+x[r<rs])) ) )
        result[r==rs] = 2./3.
        result[r>rs] =\
            2./(x[r>rs]**2. -1.) * (1.-2./np.sqrt(x[r>rs]**2.-1.)*np.arctan(np.sqrt((x[r>rs]-1.)/(1.+x[r>rs]))))
        return self.prefactor/self.rs_dm/self.rs_dm  * result

    def D_Sigma_reduced(self, r, Sigma_crit):
        result = np.zeros(r.size)
        return self.D_Sigma(r)/(1.-self.Sigma(r)/Sigma_crit)


'''
    DK14 profile, see Diemer&Kravtsov 2014, Baxter+ 2017
'''


def DK14(R, params, h_max=40., R_min=0.05, mode=0, epsr=0.001, hp=False):
    '''
        Return the excess surface density at R

        Parameters
        ----------
        mode : int
            One of 0, 1, 2, 3. Corresponding to:
            0: yes infall, yes f_trans (requires 8 parameters)
            1: no infall, no f_trans (requires 5 parameters)
        hp : bool
            forces high precision quad integral with relative error epsr

    '''
    sigmatemp = lambda x: x*DK14_S(x, params, h_max, mode, hp)
    result = np.zeros(len(R))
    mu = integrate.quad(sigmatemp, 0.0001, R_min, epsrel=epsr)[0]

    integral = 0
    for i in range(len(R)):
        if(hp):
            integral = integrate.quad(sigmatemp, R_min, R[i], epsrel=epsr)[0]
        else:
            if(i==0):
                int_x = np.linspace(R_min, R[i], 10)
            else:
                int_x = np.linspace(R[i-1], R[i], 10)
            integral += simps(int_x*DK14_S(int_x, params, h_max, mode), int_x)
        result[i] = 2.*(mu+integral)/R[i]/R[i] - DK14_S(R[i], params, h_max, mode, hp)

    return result

def DK14_S(r, params, h_max=40, mode=0, hp=False):
    '''
        Returns the line of sight integral of the 3D profile between -hmax and +hmax, evaluated at r.
    '''
    r = np.atleast_1d(r)
    result = np.zeros(len(r))

    for i in range(len(r)):
        if(hp):
            rho_temp = lambda h: rho(np.sqrt(r[i]**2.+h**2.), params, mode)
            integral = integrate.quad(rho_temp, 0., h_max)[0]
        else:
            int_h = np.linspace(0, h_max, 200)
            integral = simps(rho(np.sqrt(r[i]**2.+int_h**2.), params, mode), int_h)
        result[i] = 2*integral

    return result

def rho(r, params, mode=0):
    if(mode == 0):
        rho_s, r_s, logalpha, r_t, logbeta, loggamma, rho_0, s_e  = params
        return rho_Ein(r, [rho_s, r_s, logalpha])*f_trans(r, [r_t, logbeta, loggamma])+rho_infall(r, [rho_0, s_e])
    if(mode == 1):
        if(len(params)==8):
            rho_s, r_s, logalpha, r_t, logbeta, loggamma, rho_0, s_e  = params
        else:
            rho_s, r_s, logalpha, rho_0, s_e = params
        return rho_Ein(r, [rho_s, r_s, logalpha])#+rho_infall(r, [rho_0, s_e])

def rho_Ein(r, params):
    rho_s, r_s, logalpha = params
    alpha = (10.**logalpha)
    return rho_s * np.exp(-2./alpha*( (r/r_s)**alpha-1.) )

def f_trans(r, params):
    r_t, logbeta, loggamma = params
    beta = (10.**logbeta)
    gamma = (10.**loggamma)
    return (1.+(r/r_t)**beta)**(-gamma/beta)

def rho_infall(r, params):
    rho_0, s_e = params
    r_0 = 1.5
    return rho_0*(r/r_0)**(-s_e)
