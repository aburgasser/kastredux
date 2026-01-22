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

To make sure the code can access associated data files, you need to set the code path, which can be done in the following wayss:

  * include the path to the kastredux top-level folder in your system ``PATH`` environmental variable 
  * include the path to the kastredux top-level folder in your ``PYTHONPATH`` environmental variable 
  * set a new environmental variable called ``KASTREDUX_PATH`` to the kastredux top-level folder


KASTREDUX has been tested on python 3.6.x. 

KASTREDUX has not yet reached v1.0, so bugs are likely. Please help us squish them by 
sending bug reports to aburgasser@ucsd.edu or start an issue on the github site `<https://github.com/aburgasser/kastredux/issues>`_.


Acknowledgements
----------------
TBD


*Contents*

.. toctree::
   :maxdepth: 3

   :reduction:
   :api:
   
*Search*

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

