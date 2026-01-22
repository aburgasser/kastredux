# kastredux
 Specialized Lick/KAST optical spectral reduction package

## Installation

`kastredux` can be installed from pip:

	pip install kastredux

or from git:

	git clone https://github.com/aburgasser/kastredux.git
	cd kastredux
	python -m setup.py install

It is recommended that you install in a conda environment to ensure the dependencies do not conflict with installations of other packages

	conda create -n kastredux python=3.13
	conda activate kastredux

`kastredux` uses the following external packages:

* `astropy`: https://www.astropy.org/
* `matplotlib`: https://matplotlib.org/
* `numpy<2.0`: https://numpy.org/
* `pandas`: https://pandas.pydata.org/
* `scipy`: https://scipy.org/

## Resources

`kastredux` comes equipped with several spectral files of use in calibration and analysis of low mass stellar data:

### Flux standards

The code includes spectrophotometric flux standards from 
[Oke (1990)](https://ui.adsabs.harvard.edu/abs/1990AJ.....99.1621O),
[Hamuy et al. (1992)](https://ui.adsabs.harvard.edu/abs/1992PASP..104..533H) , and
[Hamuy et al. (1994)](https://ui.adsabs.harvard.edu/abs/1994PASP..106..566H);
see ESO's pages on [Oke standards](https://www.eso.org/sci/observing/tools/standards/spectra/okestandards.html) and [Hamuy standards](https://www.eso.org/sci/observing/tools/standards/spectra/hamuystandards.html).

### Spectral standards

* SDSS K and M dwarf. subdwarf, and giant spectral templates from 
[Kesseli et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJS..230...16K) and 
[Kesseli et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019AJ....157...63K). 
* SDSS M and L dwarf spectral templates from 
[Bochanski et al. (2007)](https://ui.adsabs.harvard.edu/abs/2007AJ....133..531B) and
[Schmidt et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014PASP..126..642S).
* Low metallicity M subdwarf optical standards from 
[Lepine et al. (2007)](https://ui.adsabs.harvard.edu/abs/2007ApJ...669.1235L). 
* L dwarf optical standards from
[Kirkpatrick et al. (1999)](https://ui.adsabs.harvard.edu/abs/1999ApJ...519..802K). 
* T dwarf optical standards from
[Burgasser et al. (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...594..510B).

### Telluric standards

OBAFGKM standards from the [Pickles (1998) library](https://ui.adsabs.harvard.edu/abs/1998PASP..110..863P).

### Sample spectra

The code also includes reduced spectra from Kast for the following sources in `*.fits` format:

* [G 233-42](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=G+233-42&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id): M5 dwarf, G=14.0
* [PM J0102+5254](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=%409627359&Name=PM%20J01020%2b5254&submit=submit): M3 dwarf, G=12.5
* [LP 555-25](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=%404008479&Name=LP%20%20555-25&submit=submit): M6.5 dwarf, G=15.8
* [TOI-2084](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=%4017843263&Name=TOI-2084&submit=submit): M2 dwarf, G=14.4, data from [Barkaoui et al. (2003)](https://ui.adsabs.harvard.edu/abs/2023A%26A...677A..38B/abstract)
* [LP 389-13](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=LP+389-13&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id): M2 dwarf, G=11.6, data from [Rothermich et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024AJ....167..253R/abstract)

## Reduction

`kastredux` includes a semi-automated pipeline for reducing optical spectra with the Kast RED and BLUE cameras for these dispersers (those in bold are fully tested)

* RED: 300/4230, 300/7500, 600/3000, 600/5000, **600/7500**, 830/8460, 1200/5000
* BLUE: 452/3306, **600/4310**, 830/3460
* LDSS-3: VPH-RED (experimental)

The main steps conducted by the reduction programs include:

* Creation of reduction script from fits files
* Creation of bias, flat field, and masking arrays
* Calibration of science images
* Determination of wavelength solution from arc lamps
* Spectral dispersion tracing using bright stars
* Boxcar and optimal extraction of 1D spectrum with background subtraction
* Flux calibration with spectrophotometric standard
* Telluric correction and secondary flux calibration with G2 or A0 flux standard

These steps are conducted using the following calls:

### Generate instructions

The first step for reduction is to use the acquired data to generate an instruction file, and to edit the file to ensure sources are assigned correctly.  Assuming all of the data are contained in a folder `data`:
 
	import kastredux as kr
	kr.makeInstructions("data",savelog=True)

This will create a `reductions` folder parallel to data generate two instruction files: `input.txt` and `input_blue.txt`, which include lines indicating the relevant files for flats, biases, and arcs; and tab-delimited lines providing extraction information for the flux calibrator (`FLUXCAL`), telluric stars (`TELLURIC`) and science targets (`SOURCE`). Baseline parameters are preset, and these can be modified or additional source extractions added as needed.

The keyword `savelog=True` also generates the CSV files  `log_[date]_RED.csv` and  `log_[date]_BLUE.csv` generated from the fits files which includes guesses on what sources are arcs, flats, biases, science, telluric standard, and flux calibrators.

At this stage, it is helpful to review the instruction files and make edits, including:

* File numbers (`FILES=###-###`), useful if there were initial "test" exposures
* Source names (`NAME=`)
* Flux calibrator (`FLUXCAL`) - note that at least one catalogued flux calibrator must be included
* Initial guess for the center of the spectral trace (`CENTER=###`) 
* Extraction window around center trace (`WINDOW=###`) 
* Regions around center trace to sample the background (`BACK=NNN,MMM`)
* Telluric star assignment (`TELLURIC=`)

You can check the initial centering of spectral traces by inspecting spatial profiles:

	kr.profileCheck("data/input.txt",verbose=True)

In addition to listing the source extraction centers (`verbose=True`), this function will produce a series of PDF files `diagnostic_profile_[NAME]_[RED/BLUE].pdf` which can be used to inspect spatial source centering.

<img width="300" alt="profile" align="middle" src="https://github.com/user-attachments/assets/7fe7f0f0-5ce2-4081-88f9-e2406d63e40e" />

Finally, there are additional keywords that can be added to the SOURCE lines to help with extraction:

* `RECENTER` = True (default): set to False to use the initial center guess as a fixed number, useful for faint sources or source in close proximity to each other
* `APPLY_TRACE` = False (default): set to True to use the profile trace from the Telluric standard,  useful for faint sources or source in close proximity to each other

### Run bulk extraction code

When you are happy with the instruction file, the full reduction can be called with the command:

	kr.reduce(instructions="data/input.txt",reset=True,verbose=True)

The keywords `reset=True` means it will overwrite prior extractions, while `verbose=True` gives detailed feedback in the output. 

Assuming everything goes well, the following diagnostic plots will be generated to check the reductions:

* `diagnostic_wavecal_[NAME]_[RED/BLUE].pdf`: two-panel plot showing the difference between arc line location and wavelength calibration fit (including RMS in Angstroms and km/s), and the arc spectrum as a function of wavelength with fit lines labeled.

<img width="300" alt="wavecal" align="middle" src="https://github.com/user-attachments/assets/c44c6462-2c22-440f-9ccf-5c30f218908e" />

* `diagnostic_trace_[NAME]_[RED/BLUE].pdf`: two-panel plot showing the trace (peak count pixel) as a function of X and Y pixel coordinate and the trace fit, and difference between peak pixels and trace fit

<img width="300" alt="trace" align="middle" src="https://github.com/user-attachments/assets/afa2e524-10bc-4670-a1ab-d907a4b36a24" />

* `diagnostic_extraction_[NAME]_[RED/BLUE].pdf`: five-panel plot showing the spatial profile, extracted count rate, signal-to-noise, background count rate, and 2D image around source trace

<img width="300" alt="extraction" align="middle" src="https://github.com/user-attachments/assets/db541611-cf45-438e-bb7c-6651a8f8d4a5" />

* `diagnostic_fluxcal_[RED/BLUE].pdf`: two-panel plot showing the ratio of calibrated to observed count rate for the flux calibrator and correction fit, and comparing the calibrated, observed, and corrected flux calibrator spectra

<img width="300" alt="fluxcal" align="middle" src="https://github.com/user-attachments/assets/d3f0d0d1-dbe4-45b2-8baf-5d75bee01f15" />

* `diagnostic_telluric_[NAME]_[RED/BLUE].pdf`: two-panel plot illustrating the observed telluric spectrum with telluric regions masked, and the telluric correction spectrum.

<img width="300" alt="telluric" align="middle" src="https://github.com/user-attachments/assets/aad2f1ce-cabc-45f4-86f5-123cfe721097" />

* `diagnostic_reflux_[NAME]_[RED/BLUE].pdf`: two-panel plot showing the ratio of observed and model telluric spectrum and corresponding correction function, and comparing the model, observed, and corrected telluric spectra

<img width="300" alt="reflux" align="middle" src="https://github.com/user-attachments/assets/3807eb2a-cd2d-4a9a-9430-81f2f5a72cd7" />

* `kast[RED/BLUE]_[NAME]_[DATE].pdf`: Final spectrum

<img width="300" alt="spectrum" align="middle" src="https://github.com/user-attachments/assets/0d732265-64f8-4686-9d77-16a846214c2d" />

It is recommended to evaluate the diagnostic plots and adjust the instruction files as needed.

The output fits file include:

* `bias_[RED/BLUE].fits`: bias file
* `flat_[RED/BLUE].fits`: normalized flatfield file
* `mask_[RED/BLUE].fits`: mask file
* `kast[RED/BLUE]_[NAME]_[DATE].fits`: Final spectrum

There are also a series of pickle (`*.pkl`) files containing intermediate data products.

### Extraction step-by-step

It is also possible to conduct reductions step-by-step if more control over the process is desired with the following steps:

1. Set up necessary information

Start with import statements and the variables needed for your reduction; the *f variables correspond to image frame numbers

	import kastredux as kr
	import numpy as np
	import os

	data_folder = "data"
	reduction_folder = "reduction"
	camera = "RED"
	grating = "600/7500"
	prefix = "r"
	darkf1,darkf2 = 1000,1010
	baisf1,biasf2 = 1000,1010
	arcf = 1021
	scif1,scif2 = 1030,1031
	sciname = "RedStar"
	flxf = 1030
	flxname = "Hiltner600"
	tellf = 1040
	
2. Generate calilbration frames

First make the bias frame

	files = ["{}{}.fits".format(prefix,int(n)) for n in np.arange(darkf1,darkf2)]
	bias_out = os.path.join(reduction_folder,"bias{}.fits",format(camera))
	bias, _ = kr.makeBias(files,folder=data_folder,mode=camera,output=bias_out)

Then make the normalized flat field frame

	files = ["{}{}.fits".format(prefix,int(n)) for n in np.arange(baisf1,biasf2)]
	flat_out = os.path.join(reduction_folder,"flat{}.fits",format(camera))
	flat, _ = kr.makeBias(files,bias,folder=data_folder,mode=camera,output=flat_out)

Then make the mask frame from the bias and flatfield

	mask_out = os.path.join(reduction_folder,"mask{}.fits",format(camera))
	mask = kr.makeMask(bias,flat,mode=camera,output=mask_file)
	flatc = kr.maskClean(flat,mask,replace=1.)

Finally generate the wavelength calibration

	arc,_ = kr.readFiles("{}{}.fits".format(prefix,int(arcf)),folder=data_folder,mode=camera)
	diagplot = "diagnostic_wavecal.pdf"
	wavecal = kr.waveCalibrateArcs(arc,,dispersion=grating,mode=camera,middle=True,plot_file=diagplot)

3. Generate flux calibration

Use your flux calibrator observation to make the flux correction function

	im,hd = kr.readFiles("{}{}.fits".format(prefix,int(flxf)),folder=data_folder,mode=camera)
	imr,var = kr.reduceScienceImage(im,bias,flat,mask,hd=hd)
	cntr = kr.findPeak(imr)
	trace = kr.traceDispersion(imr,cntr=cntr,window=10,method='maximum')
	imrect = kr.rectify(imr,trace)
	varrect = kr.rectify(var,trace)
	maskrect = kr.rectify(mask,trace)
	flatrect = kr.rectify(flat,trace)
	arcrect = kr.rectify(arc,trace)
	cntr = kr.findPeak(imrect,cntr=cntr,window=10)
	diagplot = "diagnostic_extraction_fluxstd.pdf"
	flxsp = kr. extractSpectrum(imrect,cntr=cntr,var=varrect,mask=maskrect,src_wnd=10,bck_wnd="20,40",method="boxcar",plot_file=diagplot)
	wavecal_new = kr.waveCalibrateArcs(arcrect,cntr=cntr,prior=wavecal,mode=camera)
	flxsp.applyWaveCal(wavecal_new)
	diagplot = "diagnostic_fluxcal.pdf"
	fluxcal = kr.fluxCalibrate(flxsp,flxname,fit_order=fit_order,fit_scale=flux_fit_scale,fit_range=[6000,9000],plot_file=diagplot

4. Generate telluric correction

Use your G2V telluric star (if obtained) to make the second-order flux correction and telluric correction functions

	im,hd = kr.readFiles("{}{}.fits".format(prefix,int(tellf)),folder=data_folder,mode=camera)
	imr,var = kr.reduceScienceImage(im,bias,flat,mask,hd=hd)
	cntr = kr.findPeak(imr)
	trace = kr.traceDispersion(imr,cntr=cntr,window=10,method='maximum')
	imrect = kr.rectify(imr,trace)
	varrect = kr.rectify(var,trace)
	maskrect = kr.rectify(mask,trace)
	flatrect = kr.rectify(flat,trace)
	arcrect = kr.rectify(arc,trace)
	cntr = kr.findPeak(imrect,cntr=cntr,window=10)
	diagplot = "diagnostic_extraction_tellstd.pdf"
	tellsp = kr. extractSpectrum(imrect,cntr=cntr,var=varrect,mask=maskrect,src_wnd=10,bck_wnd="20,40",method="boxcar",plot_file=diagplot)
	wavecal_new = kr.waveCalibrateArcs(arcrect,cntr=cntr,prior=wavecal,mode=camera)
	tellsp.applyWaveCal(wavecal_new)
	tellsp.applyFluxCal(fluxcal)
	diagplot = "diagnostic_telluric.pdf"
	tellcorr = kr.telluricCalibrate(tellsp,spt="G2V",plot_file=diagplot)
	diagplot = "diagnostic_reflux_tellstd.pdf"
	tellfluxcorr = kr.fluxReCalibrate(tellsp,spt="G2V",plot_file=diagplot)

5. Extract science spectrum

This case uses the trace from the telluric standard, and applies flux calibration and telluric correction

	files = ["{}{}.fits".format(prefix,int(n)) for n in np.arange(scif1,scif2)]
	ims,_ = kr.readFiles(files,folder=data_folder,mode=camera)
	im = crRejectCombine(ims,verbose=verbose)
	imr,var = kr.reduceScienceImage(im,bias,flat,mask,hd=flxhd)
	imrect = kr.rectify(imr,trace)
	varrect = kr.rectify(var,trace)
	maskrect = kr.rectify(mask,trace)
	flatrect = kr.rectify(flat,trace)
	arcrect = kr.rectify(arc,trace)
	cntr = kr.findPeak(imrect,cntr=cntr,window=10)
	diagplot = "diagnostic_extraction_science.pdf"
	scisp = kr. extractSpectrum(imrect,cntr=cntr,var=varrect,mask=maskrect,src_wnd=10,bck_wnd="20,40",method="boxcar",plot_file=diagplot)
	wavecal_new = kr.waveCalibrateArcs(arcrect,cntr=cntr,prior=wavecal,mode=camera)
	scisp.applyWaveCal(wavecal_new)
	scisp.applyFluxCal(fluxcal)
	scisp.applyFluxCal(tellfluxcorr)
	scisp.applyTelluricCal(tellcorr)

## Analysis

`kastredux` comes with several routines for analyzing optical spectra of low-mass stars and brown dwarfs. These routines operate on an `Spectrum` class object that contains the spectral data and allows for various spectral operations.

### Spectrum class

The kastredux Spectrum class is the primary data object for spectral data, and is similar to the [astropy specutils](https://specutils.readthedocs.io/en/stable/) class. In addition to arrays for wavelength, flux, uncertainty, variance and masking, the internal functions for the Spectrum class include:

* spectral math: built in functions are provided to add, subtract, multiple, and divide spectra, accounting for the appropriate wavelength solution and uncertainty propagation
* `scale(val)`: scale the flux and variance by a constant value
* `sample([w1,w2],method="median")`: sample the spectrum in a specified wavelength range using the specified statistic
* `normalize([w1,w2])`: scale the spectrum based on the maximum flux in a specified wavelength range
* `trim([w1,w2])`: trim the spectrum to the specified wavelength range
* `shift(val)`: shift the spectrum by a constant wavelength or velocity
* `smooth(width)`: apply a smoothing profile of a given pixel width
* `cleanCR()`: cleans discrepant pixels in the spectrum
* `applyMask(mask)`: applies a pixel mask, where the mask array specifies pixels to exclude as either True or 1
* `maskWave([w1,w2])`: mask pixels in a specified wavelength range 
* `applyWaveCal(wavecal)`: applies the wavelength calibration computed in `kr.waveCalibrateArcs()`
* `applyFluxCal(fluxcal)`: applies the flux calibration computed in `kr.fluxCalibrate()`
* `applyTelluricCal(tellcal)`: applies the telluric correction computed in `kr.telluricCalibrate()`
* `redden(val)`: applies a reddening to a spectrum using the [Cardelli, Clayton, and Mathis (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C/abstract) model
* `reset()`: reset the Spectrum object to its original read in state
* `convertWave(unit)`: convert wavelength array to the given wavelength unit 
* `convertFlux(unit)`: convert flux and uncertainty arrays to the given flux unit 
* `plot()`: plots the spectrum for visualization 
* `toFile(file)`: saves the spectrum to a file, including fits and tab-delimited ascii files

Spectrum class objects can be initiated through a file name or specifying wavelength, flux, and uncertainty arrays:

	import kastredux as kr
	sp = kr.Spectrum('kast_spectrum.fits")

	import astropy.unit as u
	wave_unit = u.Angstrom # default
	flux_unit = u.erg/u.s/u.cm/u.cm/u.Angstrom # default
	sp = kr.Spectrum(wave=[array]*wave_unit,flux=[array]*flux_unit,unc=[array]*flux_unit),
	
### Analysis routines
The `kastredux' analysis routines are defined primarily for late-type stars and brown dwarfs (M, L, and T dwarfs). These are the functions currently defined that operate on the Spectrum class objects (sp):

* `compareSpectra(sp1,sp2)`: compares two spectra using a defined statistic
* `classifyTemplate(sp)`: compares a spectrum to defined templates initiated using `initializeStandards()`
* `measureIndex(sp,ranges,sample="median",method="ratio")`: measures a spectral index defined by wavelengths specified in the ranges array, where sample describes how the spectral flux is measured and method describes how the fluxes are combined
* `measureIndexSet(sp,ref="lepine2003)`: measures a predefined set of spectral indices; call as `measureIndexSet(info=True)` to obtain a current list of spectral indices
* `classifyIndices(sp,ref="lepine2003)`: computes a spectral type from a predefined set of spectral indices; call as `classifyIndices(info=True)` to obtain a current list of spectral indices
* `measureEW(sp,w0)`: computes the equivalent width of a feature centered at wavelength w0
* `measureEWElement(sp,element)`: computes the equivalent widths for transitions of a given element
* `measureEWSet(sp,ref='mann2013')`: computes the equivalent widths for a predefined set of lines 
* `metallicity(sp,ref='mann2013')`: determines the metallicity from an empirical calibration of the zeta index
* `chiFactor(sp,ref='schmidt2014')`: computes the chi correction factor, and if desired the relative Halpha to bolometric luminosity using a predefined spectral-type based calibration

The function `kr.theWorks(sp)` runs all of these analysis routines together


## Citing the code

If you use this code in your research, publications, or presentations, please include the following citation:

	Burgasser (2026). aburgasser/kastredux (vXXX). Zenodo. https://doi.org/10.5281/zenodo.18333308

or in bibtex:

	@software{adam_burgasser_2026_18333308,
	  author       = {Adam Burgasser},
	  title        = {aburgasser/kastredux: vXXX},
	  month        = jan,
	  year         = 2026,
	  publisher    = {Zenodo},
	  version      = {v1.1},
	  doi          = {10.5281/zenodo.18333308},
	  url          = {https://doi.org/10.5281/zenodo.18333308},
	}

 where (vXXX) corresponds to the version used.






