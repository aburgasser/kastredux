.. KASTREDUX reduction instructions


KASTREDUX Reduction Instructions
================================

.. _`reduce()`: api.html#kastredux.reduce
.. _`readInstructions()`: api.html#kastredux.readInstructions


Overview of reduction procedures
--------------------------------

KASTREDUX performs a standard suite of reduction procedures for optical spectroscopy, optimized for point source extraction and calibration. The following are the steps that KASTREDUX will undertake:

  * Generation of an observing log based on provided data
  * Generation of bias, mask and normalized flat field frames
  * Basic image reduction: bias-subtraction, flat-field calibration, masking, re-orientation, variance calculation
  * Cosmic ray rejection for multiple exposures
  * Tracing of curved spectral dispersions and image rectification
  * Source extraction using various spatial profiles with background subtraction
  * Wavelength calibration using arc emission lamps and line identifications from `NIST <https://physics.nist.gov/PhysRefData/ASD/lines_form.html>`_
  * Flux calibration using NOAO spectrophotometric standards from `Massey et al. (1998) <https://ui.adsabs.harvard.edu/abs/1988ApJ...328..315M/abstract>`_, `Massey & Gornwall (1990) <https://ui.adsabs.harvard.edu/abs/1990ApJ...358..344M/abstract>`_, `Hamuy et al. (1992) <https://ui.adsabs.harvard.edu/abs/1992PASP..104..533H/abstract>`_ and `Hamuy et al. (1994) <https://ui.adsabs.harvard.edu/abs/1994PASP..106..566H/abstract>`_
  * Telluric absorption calibration using A0V or G2V stars

Capabilities still to be implemented:

  * Stitching together of red and blue orders
  * Combining of spectral data
  * Refined wavelength calibration using OH emission and telluric absorption

Downloading and organizing your data
------------------------------------

It is recommended that you set up a folder for the date of the observations, then two subfolders for the data, reduction, and final spectra:

:: 

  [DATE]/
    data/
    reduction/
    final/

KAST data can be downloaded from the Mt Hamilton archive `<https://mthamilton.ucolick.org/data/>`_. Expand the gzip tar file into the ``DATE/data`` folder


Preparing the input script
--------------------------

KASTREDUX takes as an input an ascii file that directs how the data is to be processed. The default name of this file is ``input.txt``, but this can be changed through the command line call. Each line contains a command word and either a value or set of keyword=value pairs, separated by tabs. An example input script is the following:

:: 

  MODE  RED
  DATE  20190918
  DATA_FOLDER /Users/adam/data/kast/190920/data/
  REDUCTION_FOLDER  /Users/adam/data/kast/190920/reduction/
  OUTPUT_FOLDER  /Users/adam/data/kast/190920/final/
  BIAS  FILES=8002-8032
  FLAT  FILES=8069-8098
  ARC_SHALLOW FILES=8000  LAMPS=AR,HG,NE
  ARC_DEEP  FILES=8001  LAMPS=AR,HG,NE
  FLUXCAL NAME=FEIGE110 FILES=2080  CENTER=335
  SOURCE  NAME=J0112+1502 FILES=8052-8053 TELLURIC=HD12846  CENTER=336
  SOURCE  NAME=J0338+0634 FILES=8063  TELLURIC=HD157170 CENTER=338
  SOURCE  NAME=J0347+0417 FILES=8061-8062 TELLURIC=HD157170 CENTER=334
  SOURCE  NAME=J1856-0822 FILES=8037-8038 TELLURIC=HD85333  CENTER=337
  SOURCE  NAME=J1856-0822B  FILES=8037-8038 TELLURIC=HD85333  CENTER=397  APPLY_TRACE=True
  SOURCE  NAME=J1856-0822C  FILES=8037-8038 TELLURIC=HD85333  CENTER=285
  SOURCE  NAME=J1856-0822D  FILES=8037-8038 TELLURIC=HD85333  CENTER=300
  SOURCE  NAME=J1856-0822E  FILES=8037-8038 TELLURIC=HD85333  CENTER=227
  SOURCE  NAME=J1928+1301 FILES=8041-8042 TELLURIC=HD17082  CENTER=335
  SOURCE  NAME=J1928+1301B  FILES=8041-8042 TELLURIC=HD17082  CENTER=409  WINDOW=5  BACK=8,15
  SOURCE  NAME=J1928+1301C  FILES=8041-8042 TELLURIC=HD17082  CENTER=431  WINDOW=5  BACK=8,15
  SOURCE  NAME=J1928+1301D  FILES=8041-8042 TELLURIC=HD17082  CENTER=262
  TELLURIC  NAME=HD12846  FILES=8055  SPT=G2V
  TELLURIC  NAME=HD17082  FILES=8043  SPT=G0V
  TELLURIC  NAME=HD85333  FILES=8040  SPT=A0V
  TELLURIC  NAME=HD157170 FILES=8034  SPT=A0V
  OBSERVERS Roman, Christian, Adam


Commands for Data Organization and Calibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DATA_FOLDER
  is the full path to where the data is; i.e., ``[DATE]/data`` in the architecture given above.

REDUCTION_FOLDER
  is the full path to where the reduction products should go; i.e., ``[DATE]/reduction`` in the architecture given above.

DDATE
  is the 8-digit observing date in the format YYYYMMDD (e.g., 2018 November 12 = 20181112), used primarily for output filename.

MODE
  is either ``RED`` or ``BLUE``; necessary (for the moment) for image rotation.

.. NOTE::
  currently only ``RED`` is supported!

BIAS
  takes the keyword ``FILES=###-###``, where the two numbers are the first and last frame numbers of the bias files.

FLAT
  takes the keyword ``FILES=###-###``, where the two numbers are the first and last frame numbers of the dome flat field files.

ARC, ARC_SHALLOW, and ARC_DEEP
  each of these takes two keywords: ``FILES=###-###``, where the two numbers are the first and last frame numbers of the arc files; and ``LAMPS=[list]`` where [list] is a comma-delimited list of lamps used. ARC is the default required command, but ARC_SHALLOW and ARC_DEEP can be used to specific short-exposure and long-exposure lamps, useful for bringing up faint arc lines.


Commands for Spectral Reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SOURCE
  contains the extraction parameters for science target; the primary keywords are ``NAME=[name]`` (name of source) and ``FILES=###-###`` (start and end file numbers); can also take ``TELLURIC=[name]`` (name of telluric source; see note below) and additional extraction keywords listed below.

FLUXCAL
  contains the extraction parameters for flux calibrator; takes the keywords ``NAME=[name of source]`` and ``FILES=###-###`` (start and end file numbers); can also take ``TELLURIC=[name]`` (name of telluric source; see note below) and additional extraction keywords listed below.

.. NOTE::
  ``NAME`` must correspond to one of the stored flux calibrators, currently FIEGE110 and HILTNER600. 

.. NOTE::
  If no ``FLUXCAL`` is provided, the program will look for the file ``cal_flux.pkl`` in the REDUCTION_FOLDER, which may have been produced on another night in the run [AJB: IS THIS TRUE?]

TELLURIC
  extraction parameters for telluric absorption calibrator, typically a G2V or A0V star; takes the keywords ``NAME=[name]`` (name of source), ``FILES=###-###`` (start and end file numbers) and ``SPT=[G2V or A0V]`` and additional extraction keywords listed below.

.. NOTE::
  the ``TELLURIC`` keyword value for the FLUXCAL and SOURCE commands must be a ``NAME`` specified in one of the TELLURIC commands. If the ``TELLURIC`` keyword is not included in FLUXCAL and SOURCE commands, or the telluric source is not found, no telluric correction will be done


Additional Extraction Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SOURCE, FLUXCAL and TELLURIC commands also take the following optional keywords to refine the extraction:
  
    * ``CENTER=###`` - center row for source
    * ``WINDOW=###`` - width of extraction in pixels
    * ``BACK=###,###`` - start and end radii for background in pixels
    * ``APPLY_TRACE=True/False`` - in SOURCE command, apply the telluric star trace (default = True)
    * ``RECENTER=True/False`` - recenter the dispersion directly from the image (default = True)
    * ``PROFILE=[source,telluric,boxcar]`` - specifies the spatial profile to use for extraction (default = source)

Additional Commands
~~~~~~~~~~~~~~~~~~~

Additional commands may be included in the input file, and will be added to the headers of the output fits files. These each take the form [COMMAND]<tab>[VALUE]. Some useful commands include:

OUTPUT_FOLDER
  full path to where the final spectral data products should go; i.e., ``DATE/final`` in the architecture given above

OBSERVERS
  specifies the names of the observers



Providing inputs through the program call
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inputs can also be provided through the call to `reduce()`_; see below


Calling the program reduce()
----------------------------

Once you are satisfied with your input file, the reduction is initiated through python starting from your reduction folder:

>>> import kastredux
>>> redux = kastredux.reduce()

If you have named your input file something else, or you are starting from a different folder, then the reduction command can be entered as:

>>> redux = kastredux.reduce(instructions='<full path to>/<command file>')

The ``redux`` output of the reduction script is a dictionary containing all of the reduction outputs. You can continue a prior reduction by feeding this dictionary back into a `reduce()`_ call:

>>> redux = kastredux.reduce(redux=redux)

.. NOTE::
  this functionality is still in development

You can feed in a previously generated parameter file generated from `readInstructions()`_ as:

>>> parameters = kastredux.readInstructions('<full path to>/<command file>')
>>> redux = kastredux.reduce(parameters=parameters)

You can also feed in a previously generated bias, flat, or mask file:

>>> redux = kastredux.reduce(bias_file='bias.fits')
>>> redux = kastredux.reduce(flat_file='flatRED.fits')
>>> redux = kastredux.reduce(mask_file='maskRED.fits')

The wavelength and flux calibration steps produce output dictionaries saved as pickle (*.pkl) files; these can be fed back into a `reduce()`_ call as follows:

>>> redux = kastredux.reduce(cal_wave_file='cal_wave.pkl')
>>> redux = kastredux.reduce(cal_flux_file='cal_flux.pkl')

It is assumed these files reside in the ``REDUCTION_FOLDER``. If they are not found, the associated reduction steps are recomputed.

You can also set some of the analysis parameters in the command line:

>>> redux = kastredux.reduce(src_wnd=5) # set the source extraction window to 5 pixels
>>> redux = kastredux.reduce(bck_wnd=[15,20]) # set the background extraction annulus to be 15-20 pixels


Outputs
-------

`reduce()`_ runs through the steps listed above, outputs intermediate and final analysis files into the ``REDUCTION_FOLDER``, and returns a dictionary containing all outputs. Specific files generated include:

bias.fits
  Median-combined bias frame

flat[MODE].fits
  Flat field fits file generated for MODE = RED or BLUE

mask[MODE].fits
  Mask fits file generated for MODE = RED or BLUE

cal_wave.pkl
  Pickle file containing a dictionary with all wavlength calibration parameters

cal_flux.pkl
  Pickle file containing a dictionary with all flux calibration parameters

cal_tell_[NAME].pkl
  Pickle file containing a dictionary with all telluric absorption correction parameters for telluric source [NAME]

diagnostic_wavecal.pdf
  PDF figure showing diagnostics for the arc lamp wavelength calibration, including the line identifications and fit residuals

diagnostic_fluxcal.pdf
  PDF figure showing diagnostics for the flux calibration, including the fit to the flux correction and the corrected flux calibrator spectrum

diagnostic_telluic_[NAME].pdf
  PDF figure showing diagnostics for the telluric absorption correction for telluric calibrator [NAME], including the continuum fit and the correction function

diagnostic_trace_[NAME].pdf
  PDF figure showing diagnostics for tracing the dispersion of source [NAME], including the trace fit and residuals

diagnostic_profile_[NAME].pdf
  PDF figure showing diagnostics for the spatial profile for source [NAME], useful for revising the centering and extraction/background windows of the extraction

diagnostic_extraction_[NAME].pdf
  PDF figure showing diagnostics for the extraction of source [NAME], including the spatial extraction profile, the raw extraction (count rate vs pixel), signal-to-noise, and background count rate

kast[MODE]_[NAME]_[DATE].fits
  Fits file containing the final spectrum of source [NAME] in mode [RED], consisting of four arrays specifying the wavelength, flux density, flux density uncertainty, and background
  
kast[MODE]_[NAME]_[DATE].txt
  Ascii file containing the final spectrum of source [NAME] in mode [RED]. The header parameters are given as commented lines (preceded by ``#``), followed by wavelength, flux density, flux density uncertainty, and background as tab-delimited columns
  
kast[MODE]_[NAME]_[DATE].pdf
  PDF figure showing the final spectrum of source [NAME] in mode [RED]
  
reduction_structure.pkl
  Pickle file containing a dictionary with all of the reduction outputs (same as ``redux`` dictionary delivered as output to `reduce()`_).


*Contents*

.. toctree::
   :maxdepth: 3

   index
   api
    
*Search*

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

