API
===============================================

Spectrum Class
--------------

.. autoclass:: kastredux.core.Spectrum
	:members:


Routines
--------

Imaging data reduction
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: kastredux.core.readInstructions
.. autofunction:: kastredux.core.makeBias
.. autofunction:: kastredux.core.makeFlat
.. autofunction:: kastredux.core.makeMask
.. autofunction:: kastredux.core.reduceScienceImage

Spectral data reduction
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: kastredux.core.findPeak
.. autofunction:: kastredux.core.traceDispersion
.. autofunction:: kastredux.core.extract2D
.. autofunction:: kastredux.core.spatialProfile
.. autofunction:: kastredux.core.extractSpectrum
.. autofunction:: kastredux.core.waveCalibrateArcs
.. autofunction:: kastredux.core.fluxCalibrate
.. autofunction:: kastredux.core.telluricCalibrate
.. autofunction:: kastredux.core.profileCheck
.. autofunction:: kastredux.core.reduce

Imaging data manipulation
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: kastredux.core.makeLog
.. autofunction:: kastredux.core.readKastFiles
.. autofunction:: kastredux.core.combineImages
.. autofunction:: kastredux.core.crRejectCombine
.. autofunction:: kastredux.core.maskClean
.. autofunction:: kastredux.core.crRejectCombine

Spectrum manipulation
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: kastredux.core.readSpectrum
.. autofunction:: kastredux.core.compareSpectra

Other helper functions
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: kastredux.core.checkDict
.. autofunction:: kastredux.core.isUnit
.. autofunction:: kastredux.core.isNumber
.. autofunction:: kastredux.core.numberList
.. autofunction:: kastredux.core.typeToNum



*Contents*

.. toctree::
   :maxdepth: 3

   :index:
   :reduction:
   :api:

*Search*


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

