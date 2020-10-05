import numpy as np
import numpy.ma as ma

import matplotlib.pyplot as plt

from astropy.io import fits
import scipy.interpolate as interp
import scipy.signal as signal
from scipy.optimize import curve_fit

import time


def load_sim_data(sim_name=None, data_dir=None, plot_color="k",
                  pointers=None,
                  little_h=1.,boxsize=100.,
                  kpc_to_mpc=False, to_pos_hcorr=False,
                  delog_sfr=False, ssfr_to_sfr=False,
                  to_logm=False, to_mass_hcorr=False, mass_units=1,
                  to_haew_abs_neg=False, 
                  columns=None, **kwargs):
    """
    loads the sim data into a dictionary for use with sightline functions
    assumes full path to each file is data_dir+pointer
    looks for seperate pointers for positions, mass, spectra, and spec outputs
    """
    ## spectra
    if pointers.get("pointer_to_spectra") == None:
        specs = None
        spec_wave = None
    else:
        spec_wave = np.genfromtxt(data_dir+pointers["pointer_to_specwav"])
        if pointers["pointer_to_spectra"][-3:] == 'npy':
            specs = np.load(data_dir+pointers["pointer_to_spectra"])
        else:
            specs = np.genfromtxt(data_dir+pointers["pointer_to_spectra"])
    
    ## positions
    if pointers.get("pointer_to_posvel") == None:
        posvel = None
        print(sim_name+": missing position data!")
    else:
        if pointers["pointer_to_posvel"][-3:] == 'csv':
            posvel_dat = np.genfromtxt(data_dir+pointers["pointer_to_posvel"], delimiter=',',dtype=None)
        else:
            posvel_dat = np.genfromtxt(data_dir+pointers["pointer_to_posvel"])

        posvel = posvel_dat[:,columns["posvel_col_start"]:(columns["posvel_col_start"]+6)]
        if kpc_to_mpc==True:
            posvel /= 1000.
            
        if to_pos_hcorr == True:
            posvel *= little_h

            
    ## mass
    if pointers.get("pointer_to_mass") is not None:
        if pointers["pointer_to_mass"] == pointers["pointer_to_posvel"]:
            mass = posvel_dat[:,columns["mass_col"]]
        else:
            mass = np.genfromtxt(data_dir+pointers["pointer_to_mass"])[:,columns["mass_col"]]
            
        if to_logm==True and to_mass_hcorr==False:
            mass = np.log10(mass*mass_units)
        elif to_logm==False and to_mass_hcorr==False:
            mass = mass*mass_units
        elif to_logm==True and to_mass_hcorr==True:
            mass = np.log10(mass/little_h*mass_units)
        elif to_logm==False and to_mass_hcorr==True:
            mass = np.log10(10**mass/little_h*mass_units)
        else:
            print(sim_name+": uh oh, something's up with the mass (check your mass flags)")
            mass = None

            
    ## sfr
    if ssfr_to_sfr == False:
        if pointers.get("pointer_to_sfr") is not None:
            if pointers["pointer_to_sfr"] == pointers["pointer_to_posvel"]:
                sfr = posvel_dat[:,columns["sfr_col"]]
            else:
                if pointers.get("pointer_to_sfr")[-3:] == 'csv':
                    sfr = np.genfromtxt(data_dir+pointers["pointer_to_sfr"], delimiter=',',skip_header=1)[:,columns["sfr_col"]]
                else:
                    sfr = np.genfromtxt(data_dir+pointers["pointer_to_sfr"])[:,columns["sfr_col"]]

            if delog_sfr==True:
                sfr = 10**sfr
            
        else:
            sfr = None
            
    elif ssfr_to_sfr == True:
        if pointers.get("pointer_to_ssfr") is not None:
            if pointers["pointer_to_ssfr"] == pointers["pointer_to_posvel"]:
                ssfr = posvel_dat[:,columns["ssfr_col"]]
            else:
                if pointers.get("pointer_to_ssfr")[-3:] == 'csv':
                    ssfr = np.genfromtxt(data_dir+pointers["pointer_to_ssfr"], delimiter=',',skip_header=1)[:,columns["ssfr_col"]]
                else:
                    ssfr = np.genfromtxt(data_dir+pointers["pointer_to_ssfr"])[:,columns["ssfr_col"]]
                    
            sfr = ssfr * 10**mass

    ## censat
    if pointers.get("pointer_to_censat") is not None:
        if pointers["pointer_to_censat"] == pointers["pointer_to_posvel"]:
            censat = posvel_dat[:,columns["censat_col"]]
        else:
            if pointers.get("pointer_to_censat")[-3:] == 'csv':
                censat = np.genfromtxt(data_dir+pointers["pointer_to_censat"], delimiter=',',skip_header=1)[:,columns["censat_col"]]
            else:
                censat = np.genfromtxt(data_dir+pointers["pointer_to_censat"])[:,columns["censat_col"]]
    else:
        censat = None
    

    ## analyzed spectra quantities
    if pointers["pointer_to_specdat"]== None:
        simdat = None
        d400 = None
        haew = None
        rmag = None
        gmag = None
        print(sim_name+": Missing analyzed spectra data!")
    else:        
        simdat = np.genfromtxt(data_dir+pointers["pointer_to_specdat"], dtype=None)
        
        if columns["d400_col"] is not None:
            d400=simdat[:,columns["d400_col"]]
        else:
            d400=None
            
        if columns["rmag_col"] is not None:
            rmag=simdat[:,columns["rmag_col"]]
        else:
            rmag=None

        if columns["gmag_col"] is not None:
            gmag=simdat[:,columns["gmag_col"]]
        else:
            gmag=None
            
        ## haew
        if to_haew_abs_neg == True and columns["haew_col"] is not None:
            haew = -1*simdat[:,columns["haew_col"]]
        elif to_haew_abs_neg == False and columns["haew_col"] is not None:
            haew = simdat[:,columns["haew_col"]]
        else:
            haew = None

        if posvel is not None:
            ngals = len(posvel)
            if len(d400) != len(posvel):
                print(sim_name+": phot data and position data are not the same length!")
        else:
            ngals = None
    
    data_dict = {"sim_name":sim_name,
                 "plot_color":plot_color,
                 "little_h":little_h,
                 "boxsize":boxsize,
                 "posvel":posvel, 
                 "mass":mass,
                 "sfr":sfr,
                 "censat":censat,
                 "haew":haew,
                 "d400":d400,
                 "rmag":rmag,
                 "gmag":gmag,
                 "spectra":specs,
                 "wavelength":spec_wave,
                 "ngals":ngals}
    
    return data_dict
    
def get_dn4000(wave,spec):
    interp_spec = interp.interp1d(wave,spec)
    blue_wav = np.linspace(3850,3950,100)
    red_wav = np.linspace(4000,4100,100)
    d4000 = np.sum(interp_spec(red_wav)) / np.sum(interp_spec(blue_wav))
    return d4000

def get_HAEW(wave, spec):

    spec1 = signal.medfilt(spec, 151)
    spec_em = spec - spec1

    bandw_Ha = np.logical_and(wave > 6554.6, wave < 6574.6)
    bandw_Ha_blueside = np.logical_and(wave > 6483.0, wave < 6513.0)
    bandw_Ha_redside = np.logical_and(wave > 6623.0, wave < 6653.0)

    spec_cont_Ha = spec1[bandw_Ha] + np.median(np.concatenate((spec_em[bandw_Ha_blueside],spec_em[bandw_Ha_redside]), axis = 0))

    HaEW = tsum(wave[bandw_Ha],np.divide((spec_cont_Ha - spec[bandw_Ha]), spec_cont_Ha))

    return -1*HaEW

def tsum(xin, yin):
    tsum = np.sum(np.abs((xin[1:]-xin[:-1]))*(yin[1:]+yin[:-1])/2. )
    return tsum
    
    
    
class Simulation(object):
    """a class containing all necessary data and functions 
    to calculate the 
    for multiple sightlines for a single simulation"""
    
    
    def __init__(self, run_params, iterable=(), **kwargs):
        """initializes the simulation class object
        expects a dictionary of run_params and 
        a dictionary of load_params. 
        creates the massive_box parameter is mass and posvel are included"""
        self.__dict__.update(iterable, **kwargs)
        self.run_params = run_params
        self.data = load_sim_data(**kwargs)
        self.obs = {}
        if self.data["mass"] is not None and self.data["posvel"] is not None:
            self.generate_massive_box(self.run_params["host_mass_thresh"])
            
            
    def convert_posvel_to_spherical(self, pos_type="galaxies", cartesian_origin=(0,0,0)):
        """sets r, theta, phi for all positions given a cartesian x,y,z origin
        assumes (0,0,0) if nothing is returned
        positions is data["posvel"] by default, but can be changed to massive_box
        by specifying pos_type="host" """
        
        if pos_type == "galaxies":
            xx = self.data["posvel"][:,0]
            yy = self.data["posvel"][:,1]
            zz = self.data["posvel"][:,2]
            self.obs.update(rr = np.sqrt((xx-cartesian_origin[0])**2+(yy-cartesian_origin[1])**2+(zz-cartesian_origin[2])**2))
            self.obs.update(theta = np.arccos((zz-cartesian_origin[2])/self.obs["rr"]))
            self.obs.update(phi = np.arctan2(yy-cartesian_origin[1],xx-cartesian_origin[0]))
        elif pos_type == "host":
            xx = self.massive_box[:,0]
            yy = self.massive_box[:,1]
            zz = self.massive_box[:,2]
            host_rr = np.sqrt((xx-cartesian_origin[0])**2+(yy-cartesian_origin[1])**2+(zz-cartesian_origin[2])**2)
            host_theta = np.arccos((zz-cartesian_origin[2])/host_rr)
            host_phi = np.arctan2(yy-cartesian_origin[1],xx-cartesian_origin[0])
            self.obs.update(proj_host_pos = (host_rr, host_theta, host_phi))
        else:
            print("unknown position key word. currently accepts 'galaxies' or 'host'")
        
        
    def generate_massive_box(self, mass_threshold=10.39):
        """uses the periodic bounding conditions to create a larger volume
        to search for host galaxies"""
        
        offsets = np.array(((0,100,0),(0,0,100),(0,-100,0),(0,0,-100),
                    (0,100,-100),(0,-100,100),(0,100,100),(0,-100,-100),
                    (100,100,0),(100,0,100),(100,-100,0),(100,0,-100),
                    (100,100,-100),(100,-100,100),(100,100,100),(100,-100,-100),
                    (-100,100,0),(-100,0,100),(-100,-100,0),(-100,0,-100),(-100,100,-100),
                    (-100,-100,100),(-100,100,100),(-100,-100,-100),(-100,0,0),(100,0,0)))/100.
        
        bounding_box = self.data["posvel"][self.data["mass"] > mass_threshold][:,:3]
        addon = self.data["posvel"][self.data["mass"] > mass_threshold][:,:3]
        for i in range(26):
            bounding_box = np.append(bounding_box,
                                     addon+offsets[i]*self.data["boxsize"]/self.data["little_h"],
                                     axis=0)
        self.massive_box = bounding_box
        self.massive_mass = np.ravel([self.data["mass"][self.data["mass"] > mass_threshold]]*27)
        

    def calculate_dhost_true(self,):
        """calculates the 3d distance to the nearest massive neighbor for all galaxies
        not meaningful for galaxies above mass_thresh set in generate massive box"""
        
        self.dhost_true = np.zeros((self.data["ngals"]))
        for i in range(self.data["ngals"]):
            gal_pos = self.data["posvel"][i,:3]
            dist_to_lmgal = np.sqrt(np.sum((gal_pos - self.massive_box)**2,axis=1))
            self.dhost_true[i] = dist_to_lmgal[(dist_to_lmgal > 0)].min()
            

    def calculate_qfrac_true(self, mass_bins, dhost_bins):
        """calculate the number of quenched and all galaxies per dhost and mass bins
        also saves input mass and dhost bins (will override with subsequent runs of calc_qfrac)"""
        
        self.qfrac_true = np.zeros((len(mass_bins)-1,len(dhost_bins)-1,2))
        dhost_distribution = self.dhost_true
        for i in range(len(mass_bins)-1):
            for j in range(len(dhost_bins)-1):
                n_all = len(self.data["mass"][(self.data["mass"] > mass_bins[i]) & (self.data["mass"] < mass_bins[i+1]) & 
                                 (self.dhost_true > dhost_bins[j]) & (self.dhost_true < dhost_bins[j+1])])
                n_q = len(self.data["mass"][(self.data["mass"] > mass_bins[i]) & (self.data["mass"] < mass_bins[i+1]) & 
                               (self.dhost_true > dhost_bins[j]) & (self.dhost_true < dhost_bins[j+1]) & 
                               (self.data["d400"] > 0.6 + 0.1*self.data["mass"]) & (self.data["haew"] < 2)])
                self.qfrac_true[i][j][0] = n_q
                self.qfrac_true[i][j][1] = n_all
        self.mass_bins_true = mass_bins
        self.mass_bins_true = dhost_bins
        
    
    def generate_observer_location(self,dist_away=5):
        """randomly generates a location DIST_AWAY from 
        a randomly selected face of the box"""
        boxsize = self.data["boxsize"]
        n1 = np.random.rand()
        n2 = np.random.rand()
        side = np.random.randint(6)
        if side == 0:
            origin = (n1*boxsize, n2*boxsize, -1*dist_away)
        elif side == 1:
            origin = (n1*boxsize, n2*boxsize, boxsize+dist_away)
        elif side == 2:
            origin = (n1*boxsize, -1*dist_away, n2*boxsize)
        elif side == 3:
            origin = (n1*boxsize, boxsize+dist_away, n2*boxsize)
        elif side == 4:
            origin = (-1*dist_away, n1*boxsize, n2*boxsize)
        elif side == 5:
            origin = (boxsize+dist_away, n1*boxsize, n2*boxsize)
        self.obs.update(origin=origin)
    
    
    def create_single_sightline(self, dist_away=5):
        """goes through the process of creating a single sightline
        observer is on a spherical shell dist_away Mpc beyond box radius
        saves all data to the obs dict"""
        
        ### first place the observer outside the box
        self.generate_observer_location()
        
        ### galaxies and massive box from observer's perspective
        self.convert_posvel_to_spherical(cartesian_origin=self.obs["origin"], pos_type="galaxies")
        self.convert_posvel_to_spherical(cartesian_origin=self.obs["origin"], pos_type="host")
        
        self.obs.update(rmag_app = self.data["rmag"]+5*np.log10(self.obs["rr"]*10**5))
        
        ### dhost calc
        self.calculate_dhost_proj()
        print("Dhost calculated")
        
        if self.run_params["add_noise_to_spectra"] == True:
            self.add_noise_to_spectra()
            if self.run_params["add_gr_noise"] == True:
                print("gr color + rmag noise added")
            else:
                print("Rmag noise added")
            
    
    def calculate_vmax(self, n_points=50):
        """calculates the vmax corrections for each galaxy 
        from the observer's perspective. depends on the apparent mag limit"""
        
        xx = np.linspace(0,self.data["boxsize"],n_points)
        yy = np.linspace(0,self.data["boxsize"],n_points)
        zz = np.linspace(0,self.data["boxsize"],n_points)
        xx_g, yy_g, zz_g = np.meshgrid(xx,yy,zz)
        xx_g, yy_g, zz_g = np.meshgrid(xx,yy,zz)
        box_rr = np.ravel(np.sqrt((xx_g-self.obs["origin"][0])**2+(yy_g-self.obs["origin"][1])**2+(zz_g-self.obs["origin"][2])**2))

        mag_lim = self.run_params["obs_mag_lim"]
        lim_dist = (10**((mag_lim-self.data["rmag"])/5.))/10**5

        vmax = np.zeros((self.data["ngals"]))
        for i in range(self.data["ngals"]):
            vmax[i] = np.float(box_rr[box_rr < lim_dist[i]].shape[0])/(n_points)**3
        self.obs.update(vmax=vmax)

        
    def calculate_dhost_proj(self):
        """calculate dhost in the project frame. 
        assumes that both the galaxies and the massive_box 
        have been converted into the spherical frame"""
        
        dhost = np.zeros((self.data["ngals"]))
        for i in range(self.data["ngals"]):
            ang_dist = ((self.obs["theta"][i]-self.obs["proj_host_pos"][1])**2+
                        (self.obs["phi"][i]-self.obs["proj_host_pos"][2])**2)**0.5
            possible_dhosts = (ang_dist*self.obs["rr"][i])[(ang_dist < 7/self.obs["rr"][i]) & 
                                                           (np.abs(self.obs["proj_host_pos"][0] - self.obs["rr"][i]) < 10) & 
                                                           (self.obs["proj_host_pos"][0] != self.obs["rr"][i])]
            if len(possible_dhosts) > 0:
                dhost[i] = possible_dhosts.min()
            else:
                dhost[i] = 7.
        self.obs.update(dhost_proj = dhost)

        
    def add_noise_to_spectra(self):
        """add noise to the spectra, remeasures dn4000 and haew 
        and then remeasure haew and dn4000 for noisy_spectra in obs
        only looks at galaxies above apparent mag limit
        only does rmag noise at the moment
        eventually this will do different things based on the gr_noise_model flag
        """
        
        if self.run_params["add_noise_to_spectra"] == True and self.run_params["add_gr_noise"] == False:
            rmag_app = self.obs["rmag_app"]
            rmag_bins = self.run_params["rmag_bins"]
            rbins_mask = self.run_params["rbins_mask"]
            binned_snr = self.run_params["sdss_binned_spectra"]
            sdss_wavelength = self.run_params["sdss_wavelength"]

            noisy_haew = np.zeros((self.data["ngals"]))
            noisy_d400 = np.zeros((self.data["ngals"]))

            digi_rmag = np.digitize(rmag_app,rmag_bins)

            noisy_spec = np.zeros((self.data["ngals"], len(sdss_wavelength)))
            for i in range(self.data["ngals"]):
                if rmag_app[i] < self.run_params["obs_mag_lim"]:
                    if digi_rmag[i] == 0:
                        new_digi = np.arange(1,len(rmag_bins))[rbins_mask][0]
                    elif digi_rmag[i] == len(rmag_bins):
                        new_digi = np.arange(1,len(rmag_bins))[rbins_mask][-1]
                    elif rbins_mask[digi_rmag[i]-1] == False:
                        new_digi = np.arange(1,len(rmag_bins))[rbins_mask][np.abs(digi_rmag[0]-1-np.arange(1,len(rmag_bins))[rbins_mask]).argmin()]
                    else:
                        new_digi = digi_rmag[i]

                    interp_spec = np.interp(sdss_wavelength, self.data["wavelength"],self.data["spectra"][i])
                    noisy_spec[i] = np.random.normal(loc=interp_spec, scale=interp_spec/binned_snr[new_digi-1])
                    noisy_haew[i] = get_HAEW(sdss_wavelength, noisy_spec[i])
                    noisy_d400[i] = get_dn4000(sdss_wavelength, noisy_spec[i])
                    
            self.obs.update(noisy_wavelength = sdss_wavelength)
            self.obs.update(noisy_spectra = noisy_spec)
            self.obs.update(noisy_haew = noisy_haew)
            self.obs.update(noisy_d400 = noisy_d400)
                    
        elif self.run_params["add_noise_to_spectra"] == True and self.run_params["add_gr_noise"] == True:
            rmag_app = self.obs["rmag_app"]
            gr_color = self.data["gmag"] - self.data["rmag"]
            
            bin_center_points = self.run_params["sdss_gr_bin_list"]
            binned_snr = self.run_params["sdss_gr_binned_snr"]
            sdss_wavelength = self.run_params["sdss_wavelength"]
            
            noisy_spec = np.zeros((self.data["ngals"], len(sdss_wavelength)))
            noisy_haew = np.zeros((self.data["ngals"]))
            noisy_d400 = np.zeros((self.data["ngals"]))
            
            for i in range(self.data["ngals"]):
                if rmag_app[i] < self.run_params["obs_mag_lim"]:
                    bin_id = (np.sum(((bin_center_points - (rmag_app[i],gr_color[i]))/(0.5,0.1))**2,axis=1)**0.5).argmin()
                    interp_spec = np.interp(sdss_wavelength, self.data["wavelength"], self.data["spectra"][i])
                    scale = np.abs(interp_spec/binned_snr[bin_id])
                    scale[np.isfinite(scale) == False] = 0
                    noisy_spec[i] = np.random.normal(loc=interp_spec, scale=scale)    
                    noisy_haew[i] = get_HAEW(sdss_wavelength, noisy_spec[i])
                    noisy_d400[i] = get_dn4000(sdss_wavelength, noisy_spec[i])
                
            self.obs.update(noisy_wavelength = sdss_wavelength)
            self.obs.update(noisy_spectra = noisy_spec)
            self.obs.update(noisy_haew = noisy_haew)
            self.obs.update(noisy_d400 = noisy_d400)
            
        else:
            self.obs.update(noisy_wavelength = self.data["wavelength"])
            self.obs.update(noisy_spectra = self.data["spectra"])
            self.obs.update(noisy_haew = self.data["haew"])
            self.obs.update(noisy_d400 = self.data["d400"])
        
        
    def calculate_qfrac_proj(self, mass_bins, dhost_bins):
        """calculate the number of quenched and total galaxies
        in the project frame given the mass and dhost bins
        assumes that you have already run create_sightline"""
        
        self.mass_bins_proj = mass_bins
        self.dhost_bins_proj = dhost_bins
        
        mag_lim = self.run_params["obs_mag_lim"]
        sb_lim = self.run_params["obs_sb_lim"]
        
        self.calculate_surface_brightness()

        mass = self.data["mass"]
        rmag_app = self.obs["rmag_app"]
        mu_app = self.obs["mu_rmag"]
        dhost_distribution = self.obs["dhost_proj"]
        
        digi_mass = np.digitize(mass, mass_bins)            
        digi_dhost= np.digitize(dhost_distribution, dhost_bins)
        
        if self.run_params["add_noise_to_spectra"]==True:
            haew = self.obs["noisy_haew"]
            d400 = self.obs["noisy_d400"]
        else:
            haew = self.data["haew"]
            d400 = self.data["d400"]
            
            
        self.obs["qfrac_proj_vmax"] = np.zeros((len(mass_bins)-1,len(dhost_bins)-1, 2))
        self.obs["qfrac_proj_novm"] = np.zeros((len(mass_bins)-1,len(dhost_bins)-1, 2))

        self.calculate_vmax()
        print("Vmax calculated")
        
        vmax_weights = self.obs["vmax"]
        
        for j in range(len(mass_bins)-1):
            for i in range(len(dhost_bins)-1):
                self.obs["qfrac_proj_vmax"][j][i][0] = np.sum(1/vmax_weights[(digi_mass == j+1) & (digi_dhost == i+1) & 
                                                           (haew < 2) & (d400 > 0.6+0.1*mass) & 
                                                           (rmag_app < mag_lim) & (mu_app < sb_lim)])
                self.obs["qfrac_proj_vmax"][j][i][1] = np.sum(1/vmax_weights[(digi_mass == j+1) & (digi_dhost == i+1) & 
                                                           (rmag_app < mag_lim) & (mu_app < sb_lim)])

                self.obs["qfrac_proj_novm"][j][i][0] = len(mass[(digi_mass == j+1) & (digi_dhost == i+1) & 
                                                           (haew < 2) & (d400 > 0.6+0.1*mass) & 
                                                           (rmag_app < mag_lim) & (mu_app < sb_lim)])
                self.obs["qfrac_proj_novm"][j][i][1] = len(mass[(digi_mass == j+1) & (digi_dhost == i+1) & 
                                                           (rmag_app < mag_lim) & (mu_app < sb_lim)])
            
            
    
    def create_multiple_sightlines(self, mass_bins, dhost_bins):
        """creates multiple sightlines and measures the qfrac
        currently overwrites the obs dictionary for each sightline
        """
        tstart = time.time() 
        n_sightlines = self.run_params["n_sightlines"]
        
        self.quenched_fractions = np.zeros((n_sightlines, len(mass_bins)-1, len(dhost_bins)-1, 2))
        for z in range(n_sightlines):
            ttstart = time.time() 
            print("Sightline "+np.str(z))
            self.create_single_sightline()
            self.calculate_qfrac_proj(mass_bins, dhost_bins)
            if self.run_params["vmax_correction"] == True:
                self.quenched_fractions[z] = self.obs["qfrac_proj_vmax"]
            else:
                self.quenched_fractions[z] = self.obs["qfrac_proj_novm"]
            
            ndur = time.time() - ttstart
            print("Sightline "+np.str(z)+" complete in {0}s".format(ndur))
            print("\n")
            
        if self.run_params["save_qfrac"] == True:
            np.save(self.__dict__["data_dir"]+self.data["sim_name"]+"_qfrac",self.quenched_fractions)
            
        ndur = time.time() - tstart
        print(np.str(n_sightlines)+" sightlines complete in {0}s".format(ndur))
        
        
        
    def calculate_galaxy_sizes(self):
        """calculates galaxy sizes with scatter"""
        Mstar = self.data["mass"]
        massSizes = np.loadtxt(self.__dict__["data_dir"]+'robs_gama.dat', skiprows = 1, usecols=[0,2])

        sigma = np.random.normal(loc = 0.0, scale= 0.5, size = Mstar.size)
        inter_size = interp.interp1d(massSizes[:,0], massSizes[:,1])

        Msize = np.ones((Mstar.shape))
        for i in range(len(Mstar)):
            if Mstar[i] > 11.0:
                Msize[i] = np.abs(np.random.normal(loc = inter_size(11.0), scale = 1.0) + sigma[i])

            else:
                if Mstar[i] < 8.35:
                    Msize[i] = np.abs(np.random.normal(loc = inter_size(8.35), scale = 0.5))
                else:
                    Msize[i] = np.abs(inter_size(Mstar[i]) + sigma[i])
        
        self.data["size"] = Msize
        
        
    def calculate_surface_brightness(self):
        """assumes to be used within single_sightline"""
        
        self.calculate_galaxy_sizes()
        size_arcsec = self.data["size"]/(self.obs["rr"]*1000)*206265
        
        rmag_flux_Jy = 10**(self.obs["rmag_app"]/(-2.5))*3631
        half_light_mag = -2.5*np.log10(rmag_flux_Jy / 2 / 3631) 
        
        half_light_surface_brightness =  half_light_mag + 2.5*np.log10(size_arcsec**2)
        self.obs["mu_rmag"] = half_light_surface_brightness
    