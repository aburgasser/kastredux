.. KASTREDUX documentation master file

KASTREDUX: KAST optical spectral reduction package
==================================================

KASTREDUX is a python-based spectral reduction package for the Lick/KAST spectrograph.


Installation and Dependencies
-----------------------------

KASTREDUX can be cloned from the github site `https://github.com/aburgasser/kastredux <https://github.com/aburgasser/kastredux>`_. which is updated on a regular basis. More detailed instructions are on the `installation <installation.html>`_ page. 

KASTREDUX has the following external dependencies:
  * `astropy <http://www.astropy.org/>`_
  * `matplotlib <http://matplotlib.org/>`_
  * `numpy <http://www.numpy.org/>`_
  * `pandas <http://pandas.pydata.org/>`_
  * `scipy <https://www.scipy.org/>`_

To make sure the code can access associated data files, you need to set the code path, which can be done in the following ways:

  * include the path to the kastredux top-level folder in your system ``PATH`` environmental variable 
  * include the path to the kastredux top-level folder in your ``PYTHONPATH`` environmental variable 
  * set a new environmental variable called ``KASTREDUX_PATH`` to the kastredux top-level folder


KASTREDUX has been tested on python 3.6.x. 

KASTREDUX has not yet reached v1.0, so bugs are likely. Please help us squish them by 
sending bug reports to aburgasser@ucsd.edu or start an issue on the github site <https://github.com/aburgasser/kastredux/issues>`_.

Using KASTREDUX
---------------

.. _`Spectrum`: splat.html?highlight=spectrum#the-splat-spectrum-object
.. _`getSpectrum()`: api.html#splat.getSpectrum
.. _`fluxCalibrate()`: api.html#splat.Spectrum.fluxCalibrate
.. _`plot()`: api.html#splat.Spectrum.plot
.. _`plotSpectrum()`: api.html#splat.plot.plotSpectrum
.. _`measureIndex()`: api.html#splat.measureIndex
.. _`measureIndexSet()`: api.html#splat.measureIndexSet
.. _`classifyGravity()`: api.html#splat.classifyGravity
.. _`classifyByXXX`: api.html#spectral-classification
.. _`compareSpectra()`: api.html#splat.compareSpectra
.. _`modelFitMCMC()`: api.html#splat.model.modelFitMCMC


Download and organize your data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is recommended that you set up a folder for the date of the observations, then two subfolders for the data, reduction, and final spectra:

DATE/
  data/
  reduction/
  final/

KAST data can be downloaded from the Mt Hamilton archive <https://mthamilton.ucolick.org/data/>`_.

Save your data into the DATE/data folder

Set up the input script
^^^^^^^^^^^^^^^^^^^^^^^

KASTREDUX takes as an input an ascii file that directs how the data is to be processed. Each line contains a command word and either a value or set of keyword=value pairs, separated by tabs. An example input script is the following:

:: 

  # Data for 20 Sep 2019
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
  SOURCE  NAME=J1856-0822B  FILES=8037-8038 TELLURIC=HD85333  CENTER=397
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

Required Keywords
~~~~~~~~~~~~~~~~~

MODE
  either ``RED`` or ``BLUE``

.. NOTE::
  currently only ``RED`` is supported

DATE
  8-digit observing date YYYYMMDD, used for filename outputs

DATA_FOLDER
  full path to where the data is; i.e., "DATE/data" in the architecture given above

REDUCTION_FOLDER
  full path to where the reduction products should go; i.e., "DATE/reduction" in the architecture given above

BIAS
  takes the keyword FILES=###-### where the two numbers are the first and last frame numbers of the bias files

FLAT
  takes the keyword FILES=###-### where the two numbers are the first and last frame numbers of the dome flat field files

ARC, ARC_SHALLOW, ARC_DEEP
  takes two keywords: FILES=###-###, where the two numbers are the first and last frame numbers of the arc files; and LAMPS=[list] where [list] is a comma-delimited list of lamps used. ARC is the default parameter, but ARC_SHALLOW and ARC_DEEP can be used to specific short-exposure and long-exposure lamps, useful for bringing up faint arc lines

FLUXCAL
  extraction parameters for flux calibrator; takes the keywords NAME=[name of source] and FILES=###-### (start and end file numbers). 
  Can optionally take keywords TELLURIC=[name of telluric source], CENTER=### (center row for source), WINDOW=### (width of extraction in pixels) and BACK=###,### (start and end radii for background in pixels)

.. NOTE::
  NAME must correspond to one of the stored flux calibrators, currently FIEGE110 and HILTNER600. 

.. NOTE::
  TELLURIC keyword value must be a name specific in one of the TELLURIC commands. If the TELLURIC keyword is not included, or the telluric source is not listed, no telluric correction will be done

.. NOTE::
  If no FLUXCAL is provided, the program will look for the file ``cal_flux.pkl`` in the REDUCTION_FOLDER, which may have been produced on another night in the run [IS THIS TRUE?]

SOURCE
  extraction parameters for science target; takes the keywords NAME=[name of source] and FILES=###-### (start and end file numbers). 
  Can optionally take keywords TELLURIC=[name of telluric source], CENTER=### (center row for source), WINDOW=### (width of extraction in pixels) and BACK=###,### (start and end radii for background in pixels)

.. NOTE::
  TELLURIC keyword value must be a name specific in one of the TELLURIC commands. If the TELLURIC keyword is not included, or the telluric source is not listed, no telluric correction will be done


Optional Keywords
~~~~~~~~~~~~~~~~~

OUTPUT_FOLDER
  full path to where the final spectral data products should go; i.e., "DATE/final" in the architecture given above

TELLURIC
  extraction parameters for telluric absorption calibrator, typically a G2V or A0V star; takes the keywords NAME=[name of source], FILES=###-### (start and end file numbers) and SPT=[G2V or A0V]
  Can optionally take CENTER=### (center row for source), WINDOW=### (width of extraction in pixels) and BACK=###,### (start and end radii for background in pixels)

OBSERVERS
  optional keyword specifying the names of the observers


Command line inputs
^^^^^^^^^^^^^^^^^^^



Acknowledgements
----------------


*Contents*

.. toctree::
   :maxdepth: 3

   quickstart
   reduction
   api
   
*Search*

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

