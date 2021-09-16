
# import necessary modules

# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const


# Code created by G. Besla & E. Patel

class NFW_Z:

    def __init__(self, Mv, cosmo):
        """ Initialize the class with the current Virial mass of the halo 
        Functions now account for the redshift you want 
        input: virial mass in Msun (mass enclosed within Rvir, which is the radius at which the dark matter
        density is DeltaVir*avg dark matter density of the universe ). """
        
        # get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc*u.km**2/u.Msun/u.s**2).value  # kpc km^2/M/S^2
        #print(self.G)
          
        # initialize the virial mass global variable    
        self.Mvir = Mv
        
        if cosmo == "Jeon": # Cosmology Same as Jeon+2021
            self.h0 = 0.71
            self.omegaM0 = 0.265
            self.omegaL0 = 1- self.omegaM0
            
        elif cosmo == "TNG": # TNG: The Next Generation
            self.h0 = 0.6774
            self.omegaM0 = 0.3089
            self.omegaL0 = 0.6911
        
        self.omegaR0 = 8.24e-5 
        self.omegaK0 = 0
        self.Ho = self.h0*100
        
    def HubbleParameter(self, z):
        # Function that defines the Hubble Parameter as a function of redshift
        # Input:   Redshift z 
        # Returns: The Hubble parameter at the given redshift in units of km/s/Mpc        
        
        # FILL THIS IN 
        M = self.omegaM0*(1+z)**3
        R = self.omegaR0*(1+z)**4
        L = self.omegaL0
        K = self.omegaK0*(1+z)**2
        
        return  self.Ho*np.sqrt(M+R+L+K)
    
    def hz(self,z):
        return self.HubbleParameter(z)/100
    
    def OmegaM_Z(self,z):
        # Function that defines the matter density parameter as a function of redshift
        # Input:  Redshift z . Can be an array
        # Output:  Matter Density Parameter at the given redshift.
        
        # FILL THIS IN
        return self.omegaM0*(1+z)**3*self.Ho**2/self.HubbleParameter(z)**2
    
    def OmegaR_Z(self,z):
        # Function that defines the radiation density parameter as a function of redshift
        # Input:  Redshift z . Can be an array
        # Output:  Radiation Density Parameter at the given redshift.
        
        # FILL THIS IN
        return self.omegaR0*(1+z)**4*self.Ho**2/self.HubbleParameter(z)**2
    
    
    def OmegaL_Z(self,z):
        # Function that defines the dark energy density parameter as a function of redshift
        # Input:  Redshift z . Can be an array
        # Output:  Dark Energy Density Parameter at the given redshift.
        
        # FILL THIS IN
        return self.omegaL0*self.Ho**2/self.HubbleParameter(z)**2
    
    def rho_crit(self,z):
        # Msun/kpc^3 = 3H^2/8piG    
        # R200 defined where density is 200*rho_crit
        # units:  (km/s/Mpc)**2 /  ( kpc (km/s)**2 / Msun ) ==> need to divide by Mpc^2 * kpc --> kpc^3 
        return self.HubbleParameter(z)**2*3/8/np.pi/self.G/1000**2
    
    def rho_vir(self,z):  
        # Msun/ kpc^3
        # virial density is defined by Deltavir*OmegaM*rhocrit , where OmegaM*rhocrit = avg density
        # NOT SURE THIS MAKES SENSE YET
        return self.delta_vir(self.OmegaM_Z(z))*self.OmegaM_Z(z)*self.rho_crit(z)
    
    def rho_mean(self,z):
        # Mean density of the halo is this: Mvir*3/4/pi/Rvir^3. But this is not rho_vir
        # Msun/kpc^3
        return self.Mvir*3/4/np.pi/self.r_vir(z)
    
    
    
    """ To compare against subfind  
    Note that we define the virial radius Rvir of a FOF-halo as the radius
of a sphere which is centered on the most-bound particle of the
group and has an overdensity 200 with respect to the critical
density. We take the enclosed mass 
mvir = 100H^2Rvir^3/G as the
virial mass, and we define the virial velocity as Vvir^2 = G mvir/Rvir.""" 
    
    
    def c_vir(self,z):
        ### NEED TO ADDJUST THIS FOR REDSHIFT 
        # Concentration parameter for halo of a given virial mass
        # taken from Klypin 2011 on Bolshoi simulations, equation 10
        a = self.Mvir*self.hz(z)/ 1e12
        return 9.60*(a)**(-0.075)

    
    def delta_vir(self,z):
        # Overdensity to define virial radius (OmegaM is a function of z)
        # delta_c taken from Bryan and Norman (1998)

        x = self.OmegaM_Z(z) - 1.
        deltac = 18*np.pi**2 + 82*x -39*x**2
        return deltac/self.OmegaM_Z(z)

    
    def r_vir(self,z):
        # virial radius. Where the density of the halo is equal to DeltaVir*AvgDensity of the universe
        # taken from van der Marel 2012, equation A1 
        #("THE M31 VELOCITY VECTOR. II. RADIAL ORBIT TOWARD THE MILKY WAY AND IMPLIED LOCAL GROUP MASS")
        # MIGHT NEED TO ADJUST THIS FOR REDSHIFT 
        a = 206./self.hz(z)
        b = self.delta_vir(z) * self.OmegaM_Z(z) / 97.2
        c = self.Mvir * self.hz(z)/1e12
        return a * b**(-1./3.) * c**(1./3.)

    def M_vir(self,r_vir,z):
        # virial mass. equation derived by Binh. 
        a = 206./self.hz(z)
        b = self.delta_vir(z) * self.OmegaM_Z(z) / 97.2
        c = self.hz(z)/1e12
        return (r_vir/a)**3 * b * (1/c)    
    
    def v_vir(self,z):
        # Circular speed at the virial radius (in km/s)
        rvir = self.r_vir(z)
        return np.sqrt(self.G*self.Mvir/rvir)
    
    
    def r_s(self, z, c=False):
        # Scale length for the NFW profile (in kpc)
        if c:
            return self.r_vir(z)/c
        else: 
            c = self.c_vir(z)
            return self.r_vir(z)/c
    
    
    def f(self,x):
        a = np.log(1+x) 
        b = x/(1+x)
        return a - b
    
    
    
    
    def mass(self, z, r, c=False):
        """NFW mass enclosed as a function of r
        Input: r = Galactocentric distance (kpc)
        c = concentration - Can take concentration as given (cvir) or give it a value
        """
        if c:
            cvir = c
        else:
            cvir = self.c_vir(z)
        
        x = r/self.r_s(z,c=cvir)
        
        return self.Mvir*self.f(x)/self.f(cvir)

    
    
    def rho(self,z,r,c=False):
        """NFW density profile as a function of r
        Input: r = Galactocentric distance (kpc)
        c = concentration - Can take concentration as given (cvir) or give it a value
        """
        if c:
            cvir = c
        else:
            cvir = self.c_vir(z)
        x = r/self.r_s(z,c=cvir)
        rho_s = self.Mvir/(4.*np.pi*self.r_s(z,c=cvir)**3.*self.f(cvir))
        rho_N = rho_s/(x*(1+x)**2.)
        return rho_N
    
    def rho_meanR(self,z,r):
        # Mean density of the entire halo is this: Mvir*3/4/pi/Rvir^3. 
        # Msun/kpc^3
        return self.mass(z,r)*3/4/np.pi/self.r_vir(z)
    
    def v_max(self,z,c=False):
        """ Maximal circular speed (km/s);  occurs at rmax = 2.163*(r_s) 
        Input: r = Galactocentric distance (kpc)
        c = concentration - Can take concentration as given (cvir) or give it a value
        """
        if c:
            cvir = c
        else:
            cvir = self.c_vir(z)
            
        return 0.465*self.v_vir*np.sqrt(cvir/self.f(cvir))
    
  
      
  
