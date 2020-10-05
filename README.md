# orchard

![Orchard logo](https://github.com/IQcollaboratory/orchard/blob/master/assets/orchard_logo.png)

A major goal of the IQ Collaboratory is the apples to apples comparison of isolated & quiescent galaxies across observations and simulations. To do so, we have developed `Orchard`, a framework for creating mock surveys from simulated volumes.

`Orchard` is currently set up to create SDSS-like surveys from large volume, hydrodynamical simulations. However, the framework is flexible and can be adjusted based on the desired input (semi-analytic models, zoom-in simulations) and output (other surveys such LSST or DESI).

This package is currently intended for IQ internal use. If you are interested in adapting `Orchard` for your own work, please get in touch with Claire Dickey (http://clairedickey.com).

## input

| Required Inputs           | Units | Definition      |
|---------------------------|-------|-----------------|
| ``posvel``         | Mpc,  km/s   | Positions (X, Y, Z) and velocities (Vx, Vy, Vz) for each galaxy in the simulation box |
| Galaxy velocities         | km/s  |
| Galaxy stellar mass       | M⊙    |
| Synthetic spectra         |       |
| Synthetic wavelength grid | Å     |
| Box size                  | Mpc   |
| h                         |       |
| Absolute photometry       | SDSS g,r bands|
| Dn4000                    |       |
| Ha EW                     | Å     |
| Number of galaxies        |       |

| Orchard parameters        | Units |
|---------------------------|-------|
| Host mass threshold       | M⊙    |
| Limiting r magnitude      | apparent mag   |
| Limiting surface brightness | mag / arcsec^2 |
| Noise model      | rmag or g-r + rmag  |

## output

| Outputs          | Units |
|------------------|-------|
| Noise-added spectra | |
| Noise-added Dn4000   |       |
| Noise-added Ha EW    | Å     |
| D_host (projected) | Mpc |
| Quiescent fraction | |
