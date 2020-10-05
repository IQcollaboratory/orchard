# orchard

![Orchard logo](https://github.com/IQcollaboratory/orchard/blob/master/assets/orchard_logo.png)

A major goal of the IQ Collaboratory is the apples to apples comparison of isolated & quiescent galaxies across observations and simulations. To do so, we have developed `Orchard`, a framework for creating mock surveys from simulated volumes.

`Orchard` is currently set up to create SDSS-like surveys from large volume, hydrodynamical simulations. However, the framework is flexible and can be adjusted based on the desired input (semi-analytic models, zoom-in simulations) and output (other surveys such LSST or DESI).

This package is currently intended for IQ internal use. If you are interested in adapting `Orchard` for your own work, please get in touch with Claire Dickey (http://clairedickey.com).

## input

| Required inputs           | Units | Definition      |
|---------------------------|-------|-----------------|
| ``posvel``          | Mpc,  km/s   | Positions (X, Y, Z) and velocities (Vx, Vy, Vz) for each galaxy in the simulation box |
| ``mass``            | M⊙    | Total stellar mass for each galaxy |
| ``spectra``         |       | Synthetic spectra generated for each galaxy |
| ``wavelength``      | Å     | Corresponding wavelength array for synthetic spectra |
| ``boxsize``         | Mpc   | Simulation box size |
| ``little_h``        |       | Value of little h used in the simulation |
| ``gmag``, ``rmag``  |       | SDSS g,r absolute magnitudes derived from the synthetic spectra | 
| ``d4000``           |       | Dn4000 index measured from the synthetic spectra |
| ``haew``            | Å     | Ha equivalent width measured from the synthetic spectra |
| ``ngals``           | N     | Number of resolved galaxies in the simulation box |

| Orchard parameters        | Units | Definition      |
|---------------------------|-------|-----------------|
| ``dhost_iso_thresh`` | Mpc | The threshold in projected distance to the nearest potential host galaxy beyond which galaxies are considered isolated |
| ``host_mass_thresh``      | log10(M⊙)    | The stellar mass above which galaxies are potential host galaxies for determining isolation |
| ``obs_mag_lim``     | apparent mag   | SDSS r band magnitude threshold below which galaxies are not observed in the mock survey |
| ``obs_sb_lim`` | mag / arcsec^2 | Surface brightness threshold below which galaxies are not observed in the mock survey |
| ``noise_model``     |   |  Orchard currently has two models for SDSS-like noise; one based on the apparent r magnitude and based on r magnitude and g-r color |
| ``n_sightlines`` | N | The number of mock survey sightlines to be generated for the simulation box |

## output

| Outputs          | Units | Definition      |
|------------------|-------|-----------------|
| ``noisy_spec`` | | Noise-added spectra for all galaxies which are observable in the mock survey |
| ``noisy_d4000``   |     | Dn4000 measured from the noise-added spectra, used to select quiescent galaxies |
| ``noisy_haew``   | Å     | Ha EW measured from the noise-added spectra, used to select quiescent galaxies |
| ``dhost_proj`` | Mpc | Projected distance to the nearest potential host galaxy, used to select isolated galaxies |
| ``quiescent_fraction`` | | The quiescent fraction of isolated galaxies measured from the mock survey |

## example simulation initialization

```python
tng_load_params = {"sim_name":"TNG","plot_color":"#49DF67",
                   ### where the data lives
                   "data_dir":data_dir,
                   "pointers":{
                       "pointer_to_mass":"IQSFSdataTNG-MstarPosVelSat.txt",
                       "pointer_to_posvel":"IQSFSdataTNG-MstarPosVelSat.txt",
                       "pointer_to_specdat":"analyzespectra_TNG_ranincl_medcontflux_sidebands.dat",
                       "pointer_to_spectra":"tng_specs.npy",
                       "pointer_to_specwav":"FSPS_wave_Full.txt"},
                   ### sim params
                   "little_h":0.6774, "boxsize":110.7,
                   ### data column indices
                   "columns":{
                       "mass_col":1, "posvel_col_start":2, 
                       "haew_col":21, "d400_col":19, 
                       "rmag_col":5, "gmag_col":6,},}
                       
run_params = {"dhost_iso_thresh":1.5, ## Mpc thresh for determining isolation
              "host_mass_thresh":10.39, ## log(mass) thresh for determining what's a host gal
              "obs_mag_lim":17.7, ## only galaxies brighter than this get counted
              "obs_sb_lim":23.0, ## only galaxies brighter than this get counted
              "mass_bins":np.arange(8,11.1,0.25), ## stellar mass bins for quiescent fraction measurement
              "noise_model":"sdss_gr_model", ## select noise model
              "n_sightlines":25 ## how many sightlines to run
             }
             
tng = Simulation(run_params,**tng_load_params)
tng.create_multiple_sightlines()
```

## flexibility

Currently, ``Orchard`` defines galaxies as quiescent based on Dn4000 and HaEW, and isolated based on projected host distance. Galaxies are selected as observable based on SDSS survey limits, and the noise models are derived from the noise characteristics of SDSS spectroscopy. All of these definitions and limits are easily adjustable, to either use different definitions of quiescence or isolation (e.g., UVJ or distance from the star forming sequence, or central/satellite) or to create synthetic versions of other surveys.
