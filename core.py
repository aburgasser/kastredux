# -*- coding: utf-8 -*-
from __future__ import print_function

# WORKING COPY OF KAST REDUCTION CODE
# THINGS THAT NEED FIXING
# - flux calibration isn't quite right - not returning a G star to a G star (esp in blue)
# - fix flux calibration "wiggle" at red end of blue dispersion
# - documentation
# - function docstrings
# - if output dictionary exists, won't do anything! -> can't used saved files
# - enlarge xticks yticks in output plots

# imports - internal
import copy
import os
import re
import shutil
import sys
import warnings

# imports - external
import matplotlib
import matplotlib.pyplot as plt
import glob
import numpy
import pandas
import pickle
from astropy.io import ascii, fits			# for reading in spreadsheet
from astropy.coordinates import SkyCoord	  # coordinate conversion
from astropy import units as u			# standard units
from astropy import constants as const		# physical constants in SI units
from scipy import stats, signal
from scipy.integrate import trapz		# for numerical integration
from scipy.interpolate import interp1d

if sys.version_info.major != 2 and sys.version_info.major != 3:
	raise NameError('\nKASTREDUX only works on Python 2.7 and 3.X\n')
if sys.version_info.major == 2:	 # switch for those using python 3
	import string


############################################################
# PROGRAM PACKAGE KEYWORDS
############################################################

NAME = 'kastredux'
VERSION = '2021.11.06'
__version__ = VERSION
GITHUB_URL = 'https://github.com/aburgasser/kastredux/'

#set the CODE_PATH, either from set environment variable or from PYTHONPATH or from sys.path
CODE_PATH = ''
if os.environ.get('{}_PATH'.format(NAME.upper())) != None:
	CODE_PATH = os.environ['{}_PATH'.format(NAME.upper())]
if os.environ.get('{}PATH'.format(NAME.upper())) != None and CODE_PATH == '':
	CODE_PATH = os.environ['{}PATH'.format(NAME.upper())]
# get from PYTHONPATH
if os.environ.get('PYTHONPATH') != None and CODE_PATH == '':
	path = os.environ['PYTHONPATH']
	for i in path.split(':'):
		if NAME in i or NAME.upper() in i:
			CODE_PATH = i
# get from system path
if CODE_PATH == '':
	checkpath = [NAME in r for r in sys.path]
	if max(checkpath):
		CODE_PATH = sys.path[checkpath.index(max(checkpath))]
	checkpath = [NAME.upper() in r for r in sys.path]
	if max(checkpath):
		CODE_PATH = sys.path[checkpath.index(max(checkpath))]
if CODE_PATH == '':
	print('Warning: could not set CODE_PATH variable from PYTHONPATH or system PATH environmental variables; some functionality may not work')
	CODE_PATH = './'

############################################################
# WELCOME!
############################################################

print('\n\nWelcome to the KASTredux reduction package!')
print('This package was developed by Adam Burgasser (aburgasser@ucsd.edu)')
print('You are currently using version {}'.format(VERSION))
#print('If you make use of any features of this toolkit for your research, please remember to cite the SPLAT paper:')
#print('\n{}; Bibcode: {}\n'.format(CITATION,BIBCODE))
print('If you make use of any spectra or models in this package, please remember to cite the original source.')
print('Please report any errors are feature requests to our github page, {}\n\n'.format(GITHUB_URL))

############################################################
# RESOURCE DATA
############################################################

FLUXCALFOLDER = CODE_PATH+'/resources/flux_standards/'
FLUXCALS = {
	'FEIGE66': {'FILE' : 'ffeige66.dat'},
	'FEIGE67': {'FILE' : 'ffeige67.dat'},
	'FEIGE110': {'FILE' : 'ffeige110.dat'},
	'HILTNER600': {'FILE' : 'fhilt600.dat'},
	'LTT7987': {'FILE' : 'fltt7987.dat'},
}

SPTSTDFOLDER = CODE_PATH+'/resources/spectral_standards/'
SPTSTDS = {}
PLOT_DEFAULTS = {'figsize': [6,4], 'fontsize': 16,
'color': 'k', 'ls': '-', 'alpha': 1, 
'background_color': 'grey', 'background_ls': '--', 'background_alpha': 1,
'comparison_color': 'm', 'comparison_ls': '--', 'comparison_alpha': 1,
'unc_color': 'grey', 'unc_ls': '--', 'unc_alpha': 1,
'zero_color': 'k', 'zero_ls': '--', 'zero_alpha': 1,
}

TELLSTDFOLDER = CODE_PATH+'/resources/telluric_standards/'

ERROR_CHECKING = False


############################################################
# KAST INSTRUMENT CONSTANTS AND FUNCTIONS
############################################################

DEFAULT_WAVE_UNIT = u.Angstrom
DEFAULT_FLUX_UNIT = u.erg/u.cm/u.cm/u.Angstrom/u.s

INSTRUMENT_MODES = {
	'RED': {'PREFIX': 'r', 'SUFFIX': '', 'ALTNAME': ['kast-red','r','rd','long','ir','nir'], 'NAME': 'KAST red', 'VERSION': 'kastr', 'ROTATE': True, 'TRIM': [],},
	'BLUE': {'PREFIX': 'b', 'SUFFIX': '', 'ALTNAME': ['kast-blue','b','bl','short','uv','vis'], 'NAME': 'KAST blue', 'VERSION': 'kastb', 'ROTATE': False, 'TRIM': [],},
	'LDSS3': {'PREFIX': 'ccd', 'SUFFIX': 'c1', 'ALTNAME': ['ldss-3','ldss-3c'], 'NAME': 'LDSS-3', 'VERSION': '', 'ROTATE': True, 'TRIM': [],},
}

DISPERSIONS = {
	'452/3306': {'MODE': 'BLUE', 'RESOLUTION': 1.41, 'LAM0': 5460.74},
	'600/3000': {'MODE': 'RED', 'RESOLUTION': 1.29, 'LAM0': 5460.74},
	'600/4310': {'MODE': 'BLUE', 'RESOLUTION': 1.02, 'LAM0': 4358.33},
	'830/3460': {'MODE': 'BLUE', 'RESOLUTION': 0.63, 'LAM0': 5460.74},
	'1200/5000': {'MODE': 'RED', 'RESOLUTION': 0.65, 'LAM0': 5460.74},
	'600/5000': {'MODE': 'RED', 'RESOLUTION': 1.3,  'LAM0': 7032.41},
	'600/7500': {'MODE': 'RED', 'RESOLUTION': 1.3, 'LAM0': 7032.41},
	'830/8460': {'MODE': 'RED', 'RESOLUTION': 0.94, 'LAM0': 7032.41},
	'300/4230': {'MODE': 'RED', 'RESOLUTION': 2.53, 'LAM0': 7032.41},
	'300/7500': {'MODE': 'RED', 'RESOLUTION': 2.53, 'LAM0': 7032.41},
	'VPH-RED': {'MODE': 'LDSS3', 'RESOLUTION': 1.18, 'LAM0': 9122.95},
}

CCD_HEADER_KEYWORDS = {
	'RED': {
		'MODE': 'VERSION', # kastr = red, kastb = blue
		'SLIT': 'SLIT_N',
		'DISPERSER': 'GRATNG_N',
		'DATE-OBS': 'DATE-OBS',
		'OBSERVER': 'OBSERVER',
		'OBJECT': 'OBJECT',
		'RA': 'RA',
		'DEC': 'DEC',
		'HA': 'HA',
		'AIRMASS': 'AIRMASS',
		'EXPTIME': 'EXPTIME',
	},
	'BLUE': {
		'MODE': 'VERSION', # kastr = red, kastb = blue
		'SLIT': 'SLIT_N',
		'DISPERSER': 'GRISM_N',
		'DATE-OBS': 'DATE-OBS',
		'OBSERVER': 'OBSERVER',
		'OBJECT': 'OBJECT',
		'RA': 'RA',
		'DEC': 'DEC',
		'HA': 'HA',
		'AIRMASS': 'AIRMASS',
		'EXPTIME': 'EXPTIME',
	},
	'LDSS3': {
		'MODE': 'INSTRUME',
		'SLIT': 'APERTURE',
		'DISPERSER': 'GRISM',
		'SPEED': 'SPEED',
		'GAIN': 'EGAIN',
		'RN': 'ENOISE',
		'DATE-OBS': 'DATE-OBS',
		'TIME-OBS': 'TIME-OBS',
		'OBSERVER': 'OBSERVER',
		'OBJECT': 'OBJECT',
		'RA': 'RA',
		'DEC': 'DEC',
		'HA': 'HA',
		'AIRMASS': 'AIRMASS',
		'EXPTIME': 'EXPTIME',
	},
}

# these will be replaced with the above merged dictionary
KAST_CCD_HEADER_KEYWORDS = {
	'MODE': 'VERSION', # kastr = red, kastb = blue
	'GRISM': 'GRISM_N',
	'BLUE_DISPERSION': 'GRISM_N',
	'GRATING': 'GRATNG_N',
	'RED_DISPERSION': 'GRATNG_N',
	'SLIT': 'SLIT_N',
	'AIRMASS': 'AIRMASS',
	'RA': 'RA',
	'DEC': 'DEC',
	'OBJECT': 'OBJECT',
	'DATE-OBS': 'DATE-OBS',
	'EXPTIME': 'EXPTIME',
}

LDSS3_CCD_HEADER_KEYWORDS = {
	'MODE': 'INSTRUME',
	'SLIT': 'APERTURE',
	'DISPERSER': 'GRISM',
	'SPEED': 'SPEED',
	'GAIN': 'EGAIN',
	'RN': 'ENOISE',
	'DATE-OBS': 'DATE-OBS',
	'TIME-OBS': 'TIME-OBS',
	'OBSERVER': 'OBSERVER',
	'OBJECT': 'OBJECT',
	'RA': 'RA',
	'DEC': 'DEC',
	'HA': 'HA',
	'AIRMASS': 'AIRMASS',
	'EXPTIME': 'EXPTIME',
}

CCD_PARAMETERS = {
	'RED-FAST': {'GAIN': 0.55, 'RN': 4.3},
	'RED-SLOW': {'GAIN': 1.9, 'RN': 3.7},
	'BLUE-FAST': {'GAIN': 1.3, 'RN': 6.5},
	'BLUE-SLOW': {'GAIN': 1.2, 'RN': 3.8},
	'LDSS3-SLOW': {'GAIN': 0.13, 'RN': 5.3},
	'LDSS3-FAST': {'GAIN': 1.45, 'RN': 7.2},
}


# extract image mode from image header
def kastRBMode(hdr,keyword='VERSION'):
	if keyword not in list(hdr.keys()):
		raise ValueError('Header does not contain red/blue mode keyword {}'.format(keyword))
	md = hdr[keyword].strip()
	mode = ''
	for r in list(INSTRUMENT_MODES.keys()):
		if INSTRUMENT_MODES[r]['VERSION'] == md: mode=r
	if mode=='': raise ValueError('Could not identify spectral image mode from keyword {}'.format(keyword))
	return mode

# extract read mode from image header
def kastReadMode(hdr,keyword_mode='VERSION',keyword_read='READ-SPD'):
# get mode
	mode = kastRBMode(hdr,keyword=keyword_mode)

# get read mode
	if keyword_read not in list(hdr.keys()):
		raise ValueError('Header does not contain read speed keyword {}'.format(keyword_read))
	rd = float(hdr[keyword_read])
	rdspd = ''
	if mode=='RED':
		if rd==20.: rdspd='FAST'
		elif rd==40.: rdspd='SLOW'
		else: raise ValueError('Read speed {} does not conform to a standard read speed for mode {}'.format(rd,mode))
	elif mode=='BLUE':
		if rd==40.: rdspd='FAST'
		elif rd==80.: rdspd='SLOW'
		else: raise ValueError('Read speed {} does not conform to a standard read speed for mode {}'.format(rd,mode))
	else:
		raise ValueError('Do not recognize mode {}'.format(mode))
	return rdspd

# extract gain from image header
def kastGain(hdr,keyword_mode='VERSION',keyword_read='READ-SPD'):
# get R/B mode
	mode = kastRBMode(hdr,keyword=keyword_mode)
# get read mode
	rdmode = kastReadMode(hdr,keyword_mode=keyword_mode,keyword_read=keyword_read)
	index = '{}-{}'.format(mode,rdmode)
	if index not in list(CCD_PARAMETERS.keys()):
		raise ValueError('Do not know how to interpret mode {} and read mode {}'.format(mode,rdmode))
	else:
		return CCD_PARAMETERS[index]['GAIN']

# extract read noise from image header
def kastRN(hdr,keyword_mode='VERSION',keyword_read='READ-SPD'):
# get R/B mode
	mode = kastRBMode(hdr,keyword=keyword_mode)
# get read mode
	rdmode = kastReadMode(hdr,keyword_mode=keyword_mode,keyword_read=keyword_read)
	index = '{}-{}'.format(mode,rdmode)
	if index not in list(CCD_PARAMETERS.keys()):
		raise ValueError('Do not know how to interpret mode {} and read mode {}'.format(mode,rdmode))
	else:
		return CCD_PARAMETERS[index]['RN']

# extract dispersion from image header
def kastDispersion(hdr,keyword_mode='VERSION',keyword_blue_dispersion='GRISM_N',keyword_red_dispersion='GRATNG_N'):
# get R/B mode
	mode = kastRBMode(hdr,keyword=keyword_mode)
# get key
	if mode=='RED': key = keyword_red_dispersion
	elif mode=='BLUE': key = keyword_blue_dispersion
	else: 
		raise ValueError('Do not know how to interpret dispersion mode {}'.format(mode))
# get header
	return kastHeaderValue(hdr,key)

# extract header parameter using look-up table as backup
def kastHeaderValue(hdr,keyword):
	if keyword in list(hdr.keys()): return hdr[keyword]
	elif keyword.upper() in list(hdr.keys()): return hdr[keyword.upper()]
	elif keyword in list(KAST_CCD_HEADER_KEYWORDS.keys()): return hdr[KAST_CCD_HEADER_KEYWORDS[keyword]]
	elif keyword.upper() in list(KAST_CCD_HEADER_KEYWORDS.keys()): return hdr[KAST_CCD_HEADER_KEYWORDS[keyword.upper()]]
	else:
		raise ValueError('Cannot find keyword {} in header or header lookup table'.format(keyword))


############################################################
# SPECTRUM CLASS
############################################################

class Spectrum(object):
	'''
	Container for spectrum object
	Includes wavelength, flux, variance, background, mask vectors; trace; and header
	Includes methods for combining spectrum objects together, reading/writing, conversion
	'''
	def __init__(self,**kwargs):
		core_attributes = {'instr': 'KAST red','name': 'Unknown source','wave': [],'flux': [],'unc': [],'variance':[],'background': [],'mask': [],'header': {},}
		default_units = {'wave': DEFAULT_WAVE_UNIT,'flux':DEFAULT_FLUX_UNIT,'unc': DEFAULT_FLUX_UNIT,'background': DEFAULT_FLUX_UNIT,'variance': DEFAULT_FLUX_UNIT**2}

# set inputs
		for k in list(core_attributes.keys()): setattr(self,k,core_attributes[k])
		for k in list(kwargs.keys()): setattr(self,k.lower(),kwargs[k])
		for k in ['wave','flux','unc','variance','background','mask']: 
			if not isinstance(getattr(self,k),numpy.ndarray): setattr(self,k,numpy.array(getattr(self,k)))
		if len(self.flux) == 0: raise ValueError('Spectrum object must be initiated with a flux array')
		if len(self.wave) == 0: self.wave = numpy.arange(len(self.flux))
		if len(self.unc) == 0: self.unc = numpy.array([numpy.nan]*len(self.flux))
#		if len(self.unc) == 0: self.unc = numpy.zeros(len(self.flux))
		if len(self.variance) == 0: self.variance = self.unc**2
		if len(self.background) == 0: self.background = numpy.zeros(len(self.flux))
		if len(self.mask) == 0: self.mask = numpy.zeros(len(self.flux))
# set units
		for k in list(default_units.keys()):
			if isUnit(getattr(self,k))==False: setattr(self,k,getattr(self,k)*default_units[k])
		for k in list(default_units.keys()):
			try:
				setattr(self,k,getattr(self,k)).to(default_units[k])
			except:
				pass
# clean up
		self.original = copy.deepcopy(self)
		self.history = ['{} spectrum of {} successfully loaded'.format(self.instr,self.name)]
		return

	def setbase(self):
		'''
		:Purpose: Sets the current state of spectrum as default, eliminates prior original
		'''
		self.original = copy.deepcopy(self)
		return

	def reset(self):
		'''
		:Purpose: Reset a spectrum to its original read-in state
		'''
		for k in list(self.original.__dict__.keys()):
			if k != 'history':
				try:
					setattr(self,k,getattr(self.original,k))
				except:
					pass

		self.history.append('Returned to original state')
		self.original = copy.deepcopy(self)
		return

	def clean(self):
		'''
		:Purpose: Cleans up spectrum elements to make sure they are properly configured
		'''
# set up units
		try: funit = self.flux.unit
		except: funit = DEFAULT_FLUX_UNIT
		try: wunit = self.wave.unit
		except: wunit = DEFAULT_WAVE_UNIT

# clean wavelength vector
		if isUnit(self.wave):
			self.wave = numpy.array(self.wave.value)*wunit
		else:
			self.wave = numpy.array(self.wave)*wunit

# clean flux vector
		for k in ['flux','unc','background']:
			if isUnit(getattr(self,k)):
				setattr(self,k,numpy.array(getattr(self,k).value)*funit)
			else:
				setattr(self,k,numpy.array(getattr(self,k))*funit)
# set variance
		self.variance = self.unc**2
# need to: 
		self.history.append('Spectrum cleaned')

		return

	def __copy__(self):
		'''
		:Purpose: Make a copy of a Spectrum object
		'''
		s = type(self)()
		s.__dict__.update(self.__dict__)
		return s

	def __repr__(self):
		'''
		:Purpose: A simple representation of the Spectrum object
		'''
		return '{} spectrum of {}'.format(self.instr,self.name)

	def __add__(self,other):
		'''
		:Purpose: A representation of addition for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

		:Output: a new Spectrum object equal to the spectral sum of the inputs
		'''
		try:
			other.wave = other.wave.to(self.wave.unit)
		except:
			raise ValueError('Cannot add spectra with wave units {} and {}'.format(self.wave.unit,other.wave.unit))
		try:
			other.flux = other.flux.to(self.flux.unit)
		except:
			raise ValueError('Cannot add spectra with flux units {} and {}'.format(self.flux.unit,other.flux.unit))

# establish the baseline wavelength grid
		out = copy.deepcopy(self)
		wave = numpy.array(copy.deepcopy(self.wave.value))
		wself = numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))
		out.wave = wave[wself]
		out.wave = out.wave*self.wave.unit
		out.mask = self.mask[wself]
#		wother = numpy.where(numpy.logical_and(other.wave.value <= numpy.nanmax(out.wave.value),other.wave.value >= numpy.nanmin(out.wave.value)))

# do the math
		for k in ['flux','background']:
			fself = interp1d(self.wave.value,getattr(self,k).value,bounds_error=False,fill_value=0.)
			fother = interp1d(other.wave.value,getattr(other,k).value,bounds_error=False,fill_value=0.)
			setattr(out,k,(fself(out.wave.value)+fother(out.wave.value))*(getattr(self,k)).unit)
# special for variance
		fself = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
		fother = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)
		if numpy.random.choice(numpy.isfinite(self.variance.value))==True and numpy.random.choice(numpy.isfinite(other.variance.value))==True: 
			out.variance = (fself(out.wave.value)+fother(out.wave.value))*self.variance.unit
		elif numpy.random.choice(numpy.isfinite(self.variance.value))==True:
			out.variance = fself(out.wave.value)*self.variance.unit
		elif numpy.random.choice(numpy.isfinite(other.variance.value))==True:
			out.variance = fother(out.wave.value)*self.variance.unit
		else:
			out.variance = numpy.array([numpy.nan]*len(out.wave))
		out.unc = out.variance**0.5

# update other information
		out.name = self.name+' + '+other.name
		out.history.append('Sum of {} and {}'.format(self.name,other.name))
		out.original = copy.deepcopy(out)
		return out
	
	def __sub__(self,other):
		'''
		:Purpose: A representation of addition for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

		:Output: a new Spectrum object equal to the spectral sum of the inputs
		'''
		try:
			other.wave = other.wave.to(self.wave.unit)
		except:
			raise ValueError('Cannot subtract spectra with wave units {} and {}'.format(self.wave.unit,other.wave.unit))
		try:
			other.flux = other.flux.to(self.flux.unit)
		except:
			raise ValueError('Cannot subtract spectra with flux units {} and {}'.format(self.flux.unit,other.flux.unit))

# establish the baseline wavelength grid
		out = copy.deepcopy(self)
		wave = numpy.array(copy.deepcopy(self.wave.value))
		wself = numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))
		out.wave = wave[wself]
		out.wave = out.wave*self.wave.unit
		out.mask = self.mask[wself]
#		wother = numpy.where(numpy.logical_and(other.wave.value <= numpy.nanmax(out.wave.value),other.wave.value >= numpy.nanmin(out.wave.value)))

# do the math
		for k in ['flux','background']:
			fself = interp1d(self.wave.value,getattr(self,k).value,bounds_error=False,fill_value=0.)
			fother = interp1d(other.wave.value,getattr(other,k).value,bounds_error=False,fill_value=0.)
			if k=='variance':
				setattr(out,k,(fself(out.wave.value)+fother(out.wave.value))*(getattr(self,k).unit))
			else:
				setattr(out,k,(fself(out.wave.value)-fother(out.wave.value))*(getattr(self,k).unit))
# special for variance
		fself = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
		fother = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)
		if numpy.random.choice(numpy.isfinite(self.variance.value))==True and numpy.random.choice(numpy.isfinite(other.variance.value))==True: 
			out.variance = (fself(out.wave.value)+fother(out.wave.value))*self.variance.unit
		elif numpy.random.choice(numpy.isfinite(self.variance.value))==True:
			out.variance = fself(out.wave.value)*self.variance.unit
		elif numpy.random.choice(numpy.isfinite(other.variance.value))==True:
			out.variance = fother(out.wave.value)*self.variance.unit
		else:
			out.variance = numpy.array([numpy.nan]*len(out.wave))
		out.unc = out.variance**0.5

# update other information
		out.name = self.name+' - '+other.name
		out.history.append('Difference of {} and {}'.format(self.name,other.name))
		out.original = copy.deepcopy(out)
		return out

	def __mul__(self,other):
		'''
		:Purpose: A representation of addition for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

		:Output: a new Spectrum object equal to the spectral sum of the inputs
		'''
		try:
			other.wave = other.wave.to(self.wave.unit)
		except:
			raise ValueError('Cannot multiply spectra with wave units {} and {}'.format(self.wave.unit,other.wave.unit))

# establish the baseline wavelength grid
		out = copy.deepcopy(self)
		wave = numpy.array(copy.deepcopy(self.wave.value))
		wself = numpy.where(numpy.logical_and(self.wave.value < numpy.nanmax(other.wave.value),self.wave.value > numpy.nanmin(other.wave.value)))
		out.wave = wave[wself]
		out.wave = out.wave*self.wave.unit
		out.mask = self.mask[wself]
#		wother = numpy.where(numpy.logical_and(other.wave.value <= numpy.nanmax(out.wave.value),other.wave.value >= numpy.nanmin(out.wave.value)))

# do the math
		for k in ['flux','background']:
			fself = interp1d(self.wave.value,getattr(self,k).value,bounds_error=False,fill_value=0.)
			fother = interp1d(other.wave.value,getattr(other,k).value,bounds_error=False,fill_value=0.)
			setattr(out,k,numpy.multiply(fself(out.wave.value),fother(out.wave.value))*(getattr(self,k).unit)*(getattr(other,k).unit))
			if k=='flux': flxs,flxo = fself,fother
# special for variance
		fself = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
		fother = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)
		if numpy.random.choice(numpy.isfinite(self.variance.value))==True and numpy.random.choice(numpy.isfinite(other.variance.value))==True: 
			out.variance = numpy.multiply(out.flux.value**2,((numpy.divide(flxs(out.wave.value),fself(out.wave.value),out=numpy.zeros_like(flxs(out.wave.value)), where=fself(out.wave.value)!=0)**2)+(numpy.divide(flxo(out.wave.value),fother(out.wave.value),out=numpy.zeros_like(flxo(out.wave.value)), where=fother(out.wave.value)!=0)**2)))*self.variance.unit*other.variance.unit
		elif numpy.random.choice(numpy.isfinite(self.variance.value))==True:
			out.variance = (numpy.multiply(flxo(out.wave.value),fself(out.wave.value))**2)*self.variance.unit*(other.flux.unit**2)
		elif numpy.random.choice(numpy.isfinite(other.variance.value))==True:
			out.variance = (numpy.multiply(flxs(out.wave.value),fother(out.wave.value))**2)*other.variance.unit*(self.flux.unit**2)
		else:
			out.variance = numpy.array([numpy.nan]*len(out.wave))
		out.unc = out.variance**0.5

# update other information
		out.name = self.name+' - '+other.name
		out.history.append('Product of {} and {}'.format(self.name,other.name))
		out.original = copy.deepcopy(out)
		return out
	

	def __div__(self,other):
		'''
		:Purpose: A representation of addition for Spectrum objects which correctly interpolates as a function of wavelength and combines variances

		:Output: a new Spectrum object equal to the spectral sum of the inputs
		'''
		try:
			other.wave = other.wave.to(self.wave.unit)
		except:
			raise ValueError('Cannot divide spectra with wave units {} and {}'.format(self.wave.unit,other.wave.unit))

# establish the baseline wavelength grid
		out = copy.deepcopy(self)
		wave = numpy.array(copy.deepcopy(self.wave.value))
		wave = wave[wave<=numpy.nanmax(other.wave.value)]
		wave = wave[wave>=numpy.nanmin(other.wave.value)]
		out.wave = wave*self.wave.unit
#		out.mask = self.mask[wself]
#		wother = numpy.where(numpy.logical_and(other.wave.value <= numpy.nanmax(out.wave.value),other.wave.value >= numpy.nanmin(out.wave.value)))

# do the math
		for k in ['flux','background']:
			fself = interp1d(self.wave.value,getattr(self,k).value,bounds_error=False,fill_value=0.)
			fother = interp1d(other.wave.value,getattr(other,k).value,bounds_error=False,fill_value=0.)
			setattr(out,k,numpy.divide(fself(out.wave.value),fother(out.wave.value),out=numpy.zeros_like(fself(out.wave.value)), where=fother(out.wave.value)!=0)*(getattr(self,k).unit)/(getattr(other,k).unit))
			if k=='flux': flxs,flxo = fself,fother
# special for variance
		fself = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
		fother = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)
		if numpy.random.choice(numpy.isfinite(self.variance.value))==True and numpy.random.choice(numpy.isfinite(other.variance.value))==True: 
			out.variance = numpy.multiply(out.flux.value**2,((numpy.divide(flxs(out.wave.value),fself(out.wave.value),out=numpy.zeros_like(flxs(out.wave.value)), where=fself(out.wave.value)!=0)**2)+(numpy.divide(flxo(out.wave.value),fother(out.wave.value),out=numpy.zeros_like(flxo(out.wave.value)), where=fother(out.wave.value)!=0)**2)))*self.variance.unit/other.variance.unit
		elif numpy.random.choice(numpy.isfinite(self.variance.value))==True:
			out.variance = (numpy.divide(fself(out.wave.value),flxo(out.wave.value),out=numpy.zeros_like(fself(out.wave.value)), where=flxo(out.wave.value)!=0)**2)*self.variance.unit/(other.flux.unit**2)
		elif numpy.random.choice(numpy.isfinite(other.variance.value))==True:
			out.variance = (numpy.divide(fother(out.wave.value),flxs(out.wave.value),out=numpy.zeros_like(fother(out.wave.value)), where=flxs(out.wave.value)!=0)**2)*self.variance.unit/(other.flux.unit**2)
		else:
			out.variance = numpy.array([numpy.nan]*len(out.wave))
		out.unc = out.variance**0.5

# update other information
		out.name = self.name+' - '+other.name
		out.history.append('Ratio of {} and {}'.format(self.name,other.name))
		out.original = copy.deepcopy(out)
		return out

	def __truediv__(self,other):
		return self.__div__(other)
	
	def scale(self,fact):
		'''
		Scale spectrum by a float value
		'''
		for k in ['flux','background','unc']:
			if k in list(self.__dict__.keys()): setattr(self,k,getattr(self,k)*fact)
		self.variance = self.unc**2
		self.history.append('Spectrum scaled by factor {}'.format(fact))
		return

	def normalize(self,rng=[]):
		'''
		Scale spectrum by a float value
		'''
		if len(rng)==0: 
			rng=[numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)]
		if isUnit(rng[0]): rng = [r.to(self.wave.unit).value for r in rng]
		if isUnit(rng): rng = rng.to(self.wave.unit).value
		if numpy.nanmin(rng) > numpy.nanmax(self.wave.value) or numpy.nanmax(rng) < numpy.nanmin(self.wave.value):
			print('Warning: normalization range {} is outside range of spectrum wave array: {}'.format(rng,[numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)]))
			return
		if numpy.nanmax(rng) > numpy.nanmax(self.wave.value): rng[1] = numpy.nanmax(self.wave.value)
		if numpy.nanmin(rng) < numpy.nanmin(self.wave.value): rng[0] = numpy.nanmin(self.wave.value)

		w = numpy.where(numpy.logical_and(self.wave.value>=numpy.nanmin(rng),self.wave.value<=numpy.nanmax(rng)))
		factor = numpy.nanmax(self.flux.value[w])
		if factor == 0.: 
			print('\nWarning: normalize is attempting to divide by zero; ignoring')
		elif numpy.isnan(factor) == True: 
			print('\nWarning: normalize is attempting to divide by nan; ignoring')
		else: 
			self.scale(1./factor)
			self.history.append('Spectrum normalized')
		return

	def smooth(self,scale,method='median'):
		'''
		Apply a smoothing profile
		currently this is hanning or median
		'''
		# smoothing function
		if method=='hanning': 
			smwin = signal.windows.hann(scale)
			for k in ['flux','unc','background']:
				setattr(self,k,(signal.convolve(getattr(self,k).value,smwin,mode='same')/numpy.sum(smwin))*getattr(self,k).unit)

		elif method=='median':
			xsamp = numpy.arange(0,len(self.wave)-scale+1,scale)
#			self.wave = numpy.array([self.wave.value[x+int(0.5*scale)] for x in xsamp])*self.wave.unit
			for k in ['flux','unc','background']:
				repl = []
				for x in xsamp: 
					if len(repl)<len(self.flux.value):
						repl.extend([numpy.nanmedian(getattr(self,k).value[x:x+scale])]*scale)
				setattr(self,k,numpy.array(repl)*getattr(self,k).unit)
#			setattr(self,k,numpy.array([numpy.nanmedian(getattr(self,k).value[x:x+scale]) for x in xsamp])*getattr(self,k).unit)
		else:
			print('Warning: cannot smooth using method {}, no change to spectrum'.format(method))
			return

		self.variance = self.unc**2
		self.history.append('Smoothed by {} using a scale of {} pixels'.format(method,scale))
		return


	def applyWaveCal(self,cal_wave):
		'''
		Apply wavelength calibration
		'''
# check for relevant parameters
		required = ['COEFF','UNIT']
		for r in required:
			if r not in list(cal_wave.keys()): raise ValueError('Required parameter {} not in wavelength calibration structure'.format(r))

		self.wave = numpy.polyval(cal_wave['COEFF'],numpy.arange(len(self.flux)))*cal_wave['UNIT']
		self.history.append('Wavelength calibration applied')
		if 'RMS' in list(cal_wave.keys()): self.history.append('RMS error on wavelength solution = {} {}'.format(cal_wave['RMS'],self.wave.unit))
		return

	def applyFluxCal(self,cal_flux,log=True,debug=False):
		'''
		Apply flux calibration
		'''
# check for relevant parameters
		required = ['COEFF','UNIT']
		for r in required:
			if r not in list(cal_flux.keys()): raise ValueError('Required parameter {} not in flux calibration structure'.format(r))

		if debug==True: 
			flx_old = copy.deepcopy(self.flux)
			print('applying flux calibration with coefficients {}'.format(cal_flux['COEFF']))
		
		for k in ['flux','unc','background']:
			if k in list(self.__dict__.keys()):
				if log==True:
					setattr(self,k,(getattr(self,k).value)*(10.**(numpy.polyval(cal_flux['COEFF'],self.wave.value)))*cal_flux['UNIT'])
				else:
					setattr(self,k,(getattr(self,k).value)*numpy.polyval(cal_flux['COEFF'],self.wave.value)*cal_flux['UNIT'])
			else: 
				if debug==True: print('cannot find parameter {} in spectrum'.format(k))
		self.variance = self.unc**2
		self.history.append('Flux calibration applied')
		if 'NAME' in list(cal_flux.keys()): self.history.append('Flux calibrator {} used'.format(cal_flux['NAME']))

		if debug==True:
			plt.plot(self.wave.value,self.flux.value,'m-')
			plt.plot(self.wave.value,flx_old.value*numpy.nanmedian(self.flux.value)/numpy.nanmedian(flx_old.value),'b-')
			plt.legend(['Corrected Spectrum','Prior Spectrum'])
			plt.xticks(fontsize=16)
			plt.yticks(fontsize=16)
			plt.xlim([numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)])
			plt.ylim([0,numpy.nanquantile(self.flux.value,0.9)*1.2])
			plt.xlabel('Wavelength (Ang)',fontsize=16)
			plt.ylabel('Apparent Flux density ({})'.format(self.flux.unit),fontsize=16)
			plt.tight_layout()
			plt.savefig('debug.pdf')
			plt.close()
			raise ValueError('check debug.pdf')

		return

	def applyTelluricCal(self,cal_tell):
		'''
		Apply telluric correction
		'''
# check for relevant parameters
		required = ['WAVE','CORRECTION']
		for r in required:
			if r not in list(cal_tell.keys()): raise ValueError('Required parameter {} not in telluric calibration structure'.format(r))

		f = interp1d(cal_tell['WAVE'],cal_tell['CORRECTION'],bounds_error=False,fill_value=0.)
		for k in ['flux','unc','background']:
			if k in list(self.__dict__.keys()):
				setattr(self,k,(getattr(self,k).value)*f(self.wave.value)*self.flux.unit)
		self.variance = self.unc**2
		self.history.append('Telluric correction applied')
		if 'NAME' in list(cal_tell.keys()): self.history.append('Telluric calibrator {} used'.format(cal_tell['NAME']))
		return


	def convertWave(self,wunit):
		'''
		Convert the wavelength to a new unit
		'''
		try:
			self.wave = self.wave.to(wunit)
			self.history.append('Spectrum wavelength converted to units of {}'.format(wunit))
		except:
			print('Cannot convert spectrum with wave units {} to {}'.format(self.wave.unit,wunit))
		return

	def convertFlux(self,funit):
		'''
		Convert the flux to a new unit
		'''
		for k in ['flux','unc','background']:
			if k in list(self.__dict__.keys()):
				try:
					setattr(self,k,(getattr(self,k).to(funit)))
				except:
					print('Cannot convert spectrum element {} into flux units {}'.format(k,funit))
		self.variance = self.unc**2
		self.history.append('Spectrum flux converted to units of {}'.format(funit))
		return

	# def clean(self,positive=True):
	# 	'''
	# 	Clean spectrum of nan, zero flux, and mask
	# 	'''
	# 	if 'mask' in list(self.__dict__.keys()): mask = copy.deepcopy(self.mask)
	# 	else: mask = numpy.zeros(len(self.wave))
	# 	if not isinstance(mask,numpy.ndarray): mask = numpy.array(mask)
	# 	mask[numpy.isnan(self.flux.value)==True] = 1
	# 	mask[numpy.isnan(self.unc.value)==True] = 1
	# 	if positive==True: mask[self.flux.value<=0] = 1
	# 	for k in ['wave','flux','unc','background','mask']:
	# 		if k in list(self.__dict__.keys()):
	# 			setattr(self,k,(getattr(self,k).value[mask==0])*getattr(self,k).unit)
	# 	self.variance = self.unc**2
	# 	self.history.append('Spectrum cleaned of {:.0f} pixels'.format(numpy.total(mask)))
	# 	return

	def plot(self,**kwargs):
		'''
		Plot the spectrum
		'''
		f,fig = plt.subplots(figsize=kwargs.get('figsize',[6,4]))
		fig.plot(self.wave.value,self.flux.value,c=kwargs.get('color','k'),ls=kwargs.get('ls','-'),alpha=kwargs.get('alpha',1))
		leg = [kwargs.get('label',self.name)]
		if kwargs.get('plot_background',False)==True:
			fig.plot(self.wave.value,self.background.value,c=kwargs.get('background_color','grey'),ls=kwargs.get('background_ls','-'),alpha=kwargs.get('background_alpha',1))
			leg.append('Background')
		fig.legend(leg)
		if kwargs.get('plot_uncertainty',True)==True:
			fig.plot(self.wave.value,self.unc.value,c=kwargs.get('unc_color','grey'),ls=kwargs.get('unc_ls','-'),alpha=kwargs.get('unc_alpha',1))
		if kwargs.get('plot_zero',True)==True:
			fig.plot(self.wave.value,numpy.zeros(len(self.wave)),c=kwargs.get('zero_color','k'),ls=kwargs.get('zero_ls','--'),alpha=kwargs.get('zero_alpha',1))
		fig.set_xlim(kwargs.get('xlim',[numpy.nanmin(self.wave.value),numpy.nanmax(self.wave.value)]))
		fig.set_ylim(kwargs.get('ylim',[numpy.nanmin(self.flux.value),numpy.nanmax(self.flux.value)]))
		fig.set_xlabel(kwargs.get('xlabel',r'Wavelength ({})'.format(self.wave.unit)),fontsize=kwargs.get('fontsize',16))
		fig.set_ylabel(kwargs.get('xlabel',r'Flux ({})'.format(self.flux.unit)),fontsize=kwargs.get('fontsize',16))
#		fig.set_xticks(fontsize=kwargs.get('fontsize',16))
#		fig.set_yticks(fontsize=kwargs.get('fontsize',16))
#		fig.show()
		return fig

	def toFile(self,file,overwrite=True,csv=False,delimiter='\t',save_header=True,save_noise=False,save_background=False,save_mask=False,comment='#',**kwargs):
		'''
		Exports a spectrum to a file
		'''
# what are we saving?
		output = [self.wave.value,self.flux.value,self.unc.value]
		labels = ['Wave ({})'.format(self.wave.unit),'Flux ({})'.format(self.flux.unit),'Uncertainty ({})'.format(self.unc.unit),]
		if save_background==True: 
			output.append(self.background.value)
			labels.append('Background ({})'.format(self.flux.unit))
		if save_mask==True: 
			output.append(self.mask)
			labels.append('Mask')
# determine which type of file
		ftype = file.split('.')[-1]
# fits file
		if (ftype == 'fit' or ftype == 'fits'):
			output = tuple(output)
			data = numpy.vstack(output)
			hdu = fits.PrimaryHDU(data)
			for k in list(self.header.keys()):
				if k.upper() not in ['SIMPLE','HISTORY','COMMENT','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND'] and k.replace('#','') != '': # and k not in list(hdu.header.keys()):
					hdu.header[k] = str(self.header[k])
			for k in list(self.__dict__.keys()):
				if isinstance(getattr(self,k),str) == True or isinstance(getattr(self,k),int) == True or isinstance(getattr(self,k),bool) == True or (isinstance(getattr(self,k),float) == True and numpy.isnan(getattr(self,k)) == False):
					hdu.header[k.upper()] = str(getattr(self,k))
			hdu.writeto(file,overwrite=overwrite)

# ascii file - by default tab delimited
		else:
			f = open(file,'w')
			if save_header == True:
				for k in list(self.header.keys()):
					if k.upper() not in ['HISTORY','COMMENT'] and k.replace('#','') != '':
						f.write('{}{} = {}\n'.format(comment,k.upper(),self.header[k]))
				for k in list(self.__dict__.keys()):
					if isinstance(getattr(self,k),str) == True or isinstance(getattr(self,k),int) == True or isinstance(getattr(self,k),bool) == True or (isinstance(getattr(self,k),float) == True and numpy.isnan(getattr(self,k)) == False):
						f.write('{}{} = {}\n'.format(comment,k.upper(),getattr(self,k)))
				lhead = '{}{}'.format(comment,labels[0])
				for l in labels[1:]: lhead=lhead+'{}{}'.format(delimiter,l)
				f.write('{}\n'.format(lhead))
			for i in range(len(self.wave.value)): 
				ln = '{}'.format(output[0][i])
				for j in range(1,len(labels)): ln=ln+'{}{}'.format(delimiter,output[j][i])
				f.write('{}\n'.format(ln))
			f.close()

		self.history.append('Spectrum saved to {}'.format(file))
		return

	def write(self,file,**kwargs): 
		self.toFile(file,**kwargs)
		return
		
	def save(self,file,**kwargs): 
		self.toFile(file,**kwargs)
		return
		


############################################################
# HELPER FUNCTIONS
# These functions are general purpose and used in set up
############################################################

def checkDict(ref,refdict,altref='altname',replace=[],verbose=ERROR_CHECKING):
	'''
	Purpose: 
		General usage program to check if a key is present in a dictionary, with the option to look through alternate names

	Required Inputs:
		:param ref: A string containing the reference for lumiosity/SpT relation, should be among the keys and alternate names in refdict
		:param refdict: dictionary containing empirical relation information

	Optional Inputs:
		None

	Output:
		A string containing SPLAT's default name for a given reference set, or False if that reference is not present

	Example:

	>>> import splat
	>>> print(splat.checkEmpiricalRelation('filippazzo',splat.SPT_LBOL_RELATIONS))
		filippazzo2015
	>>> print(splat.checkEmpiricalRelation('burgasser',splat.SPT_BC_RELATIONS))
		False
	'''
	output = False
	refc = copy.deepcopy(ref)

# check reference	
	if not isinstance(refc,str):
		return output
	if len(replace) > 0:
		for rep in replace:
			if isinstance(rep,list) == True and len(rep) > 0: refc = refc.replace(rep[0],rep[1])
	for k in list(refdict.keys()):
		if refc.lower()==k.lower(): output = k
		elif altref.lower() in [x.lower() for x in list(refdict[k].keys())]:
			try: 
				if refc.lower() in [x.lower() for x in list(refdict[k][altref.lower()])]: output = k
			except:
				if refc.lower() in [x.lower() for x in list(refdict[k][altref.upper()])]: output = k
		else: pass
	if output == False:
		if verbose: print('\nCould not find item {} in input dictionary; try: {}'.format(ref,list(refdict.keys())))

	return output


def isUnit(s):
	'''
	:Purpose: 

		Check if a quantity is an astropy unitted quantity
	
	:Required Input:

		:param s: quantity to be evaluated

	:Optional Input:

		None

	:Output:

		True if quantity is unitted, False otherwise

	:Usage:

		TBD

	'''	
	return isinstance(s,u.quantity.Quantity) or \
		isinstance(s,u.core.Unit) or \
		isinstance(s,u.core.CompositeUnit) or \
		isinstance(s,u.core.IrreducibleUnit) or \
		isinstance(s,u.core.NamedUnit) or \
		isinstance(s,u.core.PrefixUnit)

def isNumber(s):
	'''
	:Purpose: 

		Check if a quantity is a number or can be converted into a number
	
	:Required Input:

		:param s: quantity to be evaluated

	:Optional Input:

		None

	:Output:

		True if quantity is a number, False otherwise

	:Usage:

		TBD

	'''		
	s1 = copy.deepcopy(s)
	if isinstance(s1,bool): return False
	if isinstance(s1,u.quantity.Quantity): s1 = s1.value
	if isinstance(s1,float): return (True and not numpy.isnan(s1))
	if isinstance(s1,int): return (True and not numpy.isnan(s1))
	try:
		s1 = float(s1)
		return (True and not numpy.isnan(s1))
	except:
		return False
	


def makeLog(folder,output='log.csv',prefix='',suffix='',mode='RED',sort='DATETIME'):
	"""
	Generates a log dynamically from the original fits files and outputs to a csv, tsv, or excel file

	Parameters
	----------
	folder : str
		The full path to the folder containing the data

	output : str, default = 'log.csv'
		The filename for the output log file, which should include the full path. 
		The format of the file is determined by the file suffix, with three recognized options:
			* `.csv`: comma-delimited file
			* `.txt` or `.tsv`: tab-delimited file
			* `.htm` or `.html`: html file
			* `.tex`: latex table file
			* `.xls` or `.xslx`: excel-style file

	prefix : str, default = ''
		Prefix for data files; if not provided, pre-populated with instrument defaults

	suffix : str, default = ''
		Suffix for data files; if not provided, pre-populated with instrument defaults

	mode : str, default = 'kast'
		Instrument/mode that acquired data; currently red, blue, and ldss3 are recognized

	sort : str, default = 'DATETIME'
		Parameter to sort the log file by. The default variable `DATETIME` sorts by a combination of DATE-OBS and TIME-OBS

	Returns
	-------
	pandas Dataframe
		Dataframe for log; also written to file specified by `output` parameter

	Examples
	--------
		Here's an example that generates separate log files for red and blue data and writes these to html tables

		>>> import kastredux as kr
		>>> log_red = kr.makeLog('/data/kast/200108/data/',output='/data/kast/200108/log_red.html',mode='RED')
		>>> log_blue = kr.makeLog('/data/kast/200108/data/',output='/data/kast/200108/log_blue.html',mode='BLUE')
		>>> log_read.head()
					File	DATE-OBS	OBJECT	RA	DEC	HA	AIRMASS	EXPTIME	SLIT	DISPERSER	Note	TIME-OBS
		... 0	r1001.fits	2020-01-08	Arcs	23:34:15.67	37:15:31.6	00:00:04.55	1	3	0.5 arcsec	600/7500		00:33:07.37
		... 1	r1002.fits	2020-01-08	Arcs	23:35:10.74	37:15:31.4	00:00:04.55	1	30	0.5 arcsec	600/7500		00:33:35.64
		... 2	r1003.fits	2020-01-08	Bias	23:36:56.11	37:15:31.1	00:00:04.55	1	0	0.5 arcsec	600/7500		00:35:47.68
		... 3	r1004.fits	2020-01-08	Bias	23:37:31.15	37:15:31.0	00:00:04.55	1	0	0.5 arcsec	600/7500		00:36:26.59
		... 4	r1005.fits	2020-01-08	Bias	23:38:01.19	37:15:30.9	00:00:04.55	1	0	0.5 arcsec	600/7500		00:36:54.36

"""
# check for data
	if os.path.isdir(folder) == False: 
		print('Cannot find folder {}; no log created'.format(folder))
		return

# check instrument
	mode = checkDict(mode,INSTRUMENT_MODES)
	if mode == False: raise ValueError('Do not recognize instrument {}; use RED, BLUE or LDSS3'.format(instrument))
	hdkey = CCD_HEADER_KEYWORDS[mode]

# check files
	if prefix=='': prefix = INSTRUMENT_MODES[mode]['PREFIX']
	if suffix=='': suffix = INSTRUMENT_MODES[mode]['SUFFIX']
	files = glob.glob('{}/{}*{}.fits'.format(folder,prefix,suffix))
	if len(files) == 0: 
		print('Cannot find any {}*{}.fits data files in {}; no log created'.format(prefix,suffix,folder))
		return

# set up data frame
	dp = pandas.DataFrame()
	dp['File'] = [f.split('/')[-1] for f in files]

# fill in key log parameters
	logkeys = ['DATE-OBS','TIME-OBS','OBJECT','RA','DEC','HA','AIRMASS','EXPTIME','SLIT','DISPERSER']
	for c in logkeys: 
		if c in list(hdkey.keys()): dp[c] = ['']*len(files)
	dp['Note'] = ['']*len(files)
	for i,f in enumerate(files):
		hdu = fits.open(f)
		hdu.verify('silentfix')
		h = hdu[0].header
		hdu.close()
		for c in logkeys: 
			if c in list(hdkey.keys()): dp.loc[i,c] = h[hdkey[c]]

# for kast data, TIME-OBS is created from DATE-OBS
	if mode=='RED' or mode=='BLUE':
		dp['TIME-OBS'] = ['']*len(files)
		for i,date in enumerate(dp['DATE-OBS']):
			d,t = date.split('T')
			dp.loc[i,'DATE-OBS'] = d
			dp.loc[i,'TIME-OBS'] = t

# sort entries
	if sort in list(dp.columns): 
		dp.sort_values(sort,inplace=True)
		dp.reset_index(inplace=True,drop=True)
	if sort=='DATETIME':
		dp['DT'] = dp['DATE-OBS']+dp['TIME-OBS']		
		dp.sort_values('DT',inplace=True)
		dp.reset_index(inplace=True,drop=True)
		del dp['DT']

# save file
	if 'xls' in output.split('.')[-1]: dp.to_excel(output,index=False)	
	elif 'csv' in output.split('.')[-1]: dp.to_csv(output,index=False)	
	elif 'txt' in output.split('.')[-1]: dp.to_csv(output,index=False,sep='\t')
	elif 'htm' in output.split('.')[-1]: dp.to_html(output,index=False)
	elif 'tex' in output.split('.')[-1]: dp.to_tex(output,index=False)
	else: raise ValueError('Could not save to output file {}; try .xslx, .csv, or .txt file formats instead'.format(output))
	return dp

def numberList(numstr,sort=False):
	'''
	:Purpose: 

		Convert a string listing of numbers into an array of numbers
	
	:Required Input:

		:param numstr: string indicating number list, e.g., '45,50-67,69,72-90'

	:Optional Input:

		:param sort: set to True to sort output list (default = False)

	:Output:

		list of integers specified by string

	:Usage:
		>>> import kastredux
		>>> a = kastredux.numberList('45,50-67,69,72-90')
		>>> print(a)
			[45, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 69, 
			72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90]
	'''
# check
	if not isinstance(numstr,str): 
		raise ValueError('\nInput to numberList {} must be a string'.format(numstr))

	numlist = []
	tmp1 = numstr.replace(' ','')
	tmp2 = tmp1.split(',')
	for a in tmp2:
		tmp3 = a.split(';')
		for b in tmp3:
			tmp4 = b.split('-')
			if len(tmp4) > 1:
				numlist.extend(list(range(int(tmp4[0]),int(tmp4[1])+1)))
			else:
				numlist.append(int(tmp4[0]))
	
	if sort==True: numlist = sorted(numlist)
	return numlist


def typeToNum(input, prefix='',suffix='',verbose=ERROR_CHECKING):
	'''
	:Purpose: 

		Converts between string and numeric spectral types, with the option of specifying the class prefix/suffix and uncertainty tags

	:Required inputs: 

		:param inp: Spectral type to convert. Can convert a number or a string from 0.0 (K0) and 49.0 (Y9).

	:Optional inputs: 

		:param: error = '': flag to indicate magnitude of classification uncertainty; by default ':' for uncertainty > 1 subtypes and '::' for uncertainty > 2 subtype added as suffix to string output. Can also use `err`.
		:param: uncertainty = 0: numerical uncertainty of classification; can also use `unc`
		:param: subclass = 'dwarf': spectral class; options include:

			- *field* or *fld* or *alpha*: object is a field dwarf - no prefix/suffix to string output
			- *sd* or *subdwarf*: object is a subdwarf - 'sd' prefix to string output
			- *dsd* or *d/sd*: object is an intermediate subdwarf - 'd/sd' prefix to string output
			- *esd*: object is an extreme subdwarf - 'esd' prefix to string output
			- *usd*: object is an ultra subdwarf - 'usd' prefix to string output
			- *delta*: object is a extremely low surface gravity dwarf (~1 Myr) - 'delta' suffix to string output
			- *vlg* or *gamma* or *lowg*: object is a low surface gravity dwarf (~10 Myr) - 'gamma' suffix to string output
			- *intg* or *beta*: object is an intermediate surface gravity dwarf (~100 Myr) - 'beta' suffix to string output
			- *giant*: object is a giant with luminosity class III suffix added to string output
			- *subgiant*: object is a subgiant with luminosity class IV suffix added to string output
			- *supergiant*: object is a supergiant with luminosity class I suffix added to string output

		:param: metallicity_class = '': metallicity class of object, traditionally represented by 'sd','d/sd','esd','usd', and added on as prefix to string output. Can also use `lumclass`
		:param: luminosity_class = '': luminosity class of object traditionally represented by roman numerals (e.g., 'III') and added on as suffix to string output. Can also use `lumclass`
		:param: age_class = '': age class of object, traditionally one of 'alpha', 'beta', 'gamma', 'delta' and added on as suffix to string output (see subclass). Can also use 'ageclass'
		:param: color_class: color class of object, traditionally 'b' (for blue) or 'r' (for red), added as prefix to string output. Can also use 'colorclass'
		:param: peculiar = False: Set to True if object is peculiar, which adds a 'pec' suffix to string output
		:param: verbose = False: Set to True to provide more feedback

	:Outputs: 

		The number or string of a spectral type

	:Example:
		>>> import splat
		>>> print splat.typeToNum(30)
			T0.0
		>>> print splat.typeToNum('T0.0')
			30.0
		>>> print splat.typeToNum(27, peculiar = True, uncertainty = 1.2, lumclass = 'II')
			L7.0IIp:
		>>> print splat.typeToNum(50)
			Spectral type number must be between 0 (K0) and 49.0 (Y9)
			nan
	'''
# keywords
	spletter = 'OBAFGKMLTY'
	strip = ['esd','usd','sd','alpha','beta','gamma','I','V','+/-',':','pec','p']
	inp = copy.deepcopy(input)

# only works on individual inputs
	if isinstance(inp,list):
		raise ValueError('\nInput to typeToNum() must be a single element (string or number)')

# number -> spectral type
	if isNumber(inp):
		spind = int(abs(inp/10.))
		if spind < 0 or spind > len(spletter): 
			if verbose: print('Spectral type number must be between 0 ({}0) and {} ({}9)'.format(spletter[0],len(spletter)*10.-1.,spletter[-1]))
			return 'N/A'
		spdec = numpy.around(inp,1)-spind*10.
		
		return '{}{}{:3.1f}{}'.format(prefix,spletter[spind],spdec,suffix)

# spectral type -> number
	elif isinstance(inp,str):
		for s in strip: inp.replace(s,'')
		sptype = re.findall('[{}]'.format(spletter),inp.upper())
		outval = 0.

# specialty classes				
		if len(sptype) >= 1:
			outval = spletter.find(sptype[0])*10.
			spind = inp.find(sptype[0])+1
			if spind < len(inp):
				if inp.find('.') < 0:
					if isNumber(inp[spind]):
						outval = outval+float(inp[spind])
				else:
					try:
						outval = outval+float(inp[spind:spind+3])
						spind = spind+3
					except:
						if verbose: print('\nProblem converting input type {} to a numeric type'.format(inp))
						outval = numpy.nan
			return outval

		else:
			if verbose: print('\nOnly spectral classes {} are handled by typeToNum'.format(spletter))
			return numpy.nan

# none of the above - return the input
	else:
		if verbose: print('\nWarning: could not recognize format of spectral type {}\n'.format(inp))
		return inp



############################################################
# SPECTRUM MANIPULATION FUNCTIONS
# These functions are used in conjunctino with Spectrum structure
############################################################


def readSpectrum(filename,file_type='',delimiter='\s+',comment='#',columns=['wave','flux','unc','background','mask'],wave_unit=DEFAULT_WAVE_UNIT,flux_unit=DEFAULT_FLUX_UNIT,verbose=ERROR_CHECKING,**kwargs):
	'''
	:Purpose:

		Reads in a file and puts into a Spectrum structure

	:Required Inputs: 

		:param file: full path to data file

	:Optional Inputs: 

		:param file_type='': preset for delimiter for ascii files 
			* 'tab' = tab-delimited
			* 'csv' = comma-delimited
			* 'semicsv' = semi-colon delimited
			* 'ascii' = text with specified delimiter
			* 'excel' = excel spreadsheet
			* 'fit' or 'fits' = fits file
		:param delimiter='\s+': delimiter for ascii file (by default this is "white space")
		:param comment='#': character(s) prefixing a comment line
		:param columns=['wave','flux','unc','background','mask']: if not specified by a header, this is the assumed column order
		:param wave_unit=DEFAULT_WAVE_UNIT: default unit for wavelength
		:param flux_unit=DEFAULT_FLUX_UNIT: default unit for flux, uncertainty and background

		Additional parameters for setting properties in Spectrum object can be passed through **kwargs

	:Output: 

		Spectrum object

	:Usage: 

		TBD

	'''
# check file
	if os.path.exists(filename)==False:
		raise ValueError('Could not find file {}'.format(filename))
	ftype = (os.path.basename(filename).split('.'))[-1]


	if file_type=='space': delimiter = '\s+'
	if file_type=='tab': delimiter = '\t'
	if file_type=='csv': delimiter = ','
	if file_type=='semicsv': delimiter = ';'

# fits - only in special case at the moment
	if ftype in ['fit','fits'] or file_type in ['fits','fit']:
		with fits.open(os.path.normpath(filename),ignore_missing_end=True,do_not_scale_image_data=True) as data:
			data.verify('silentfix+ignore')
			if 'NAXIS3' in list(data[0].header.keys()):
				d = numpy.copy(data[0].data[0,:,:])
			else:
				d =  numpy.copy(data[0].data)
			header = data[0].header
		if header['NAXIS']==1:
			w0,w1 = None,None
			for x in ['CRVAL1']:
				if x in list(header.keys()): w0 = header[x]
			for x in ['CDELT1','CD1_1']:
				if x in list(header.keys()): w1 = header[x]
			if w0==None or w1==None: raise ValueError('Could not find the necessary keywords in fits header to generate wavelength scale')
			if 'linear' in header['CTYPE1'].lower(): wave = numpy.arange(len(d))*w1+w0
			else: raise ValueError('Cannot handle wavelength scale type {}'.format(header['CTYPE1']))
			sp = Spectrum(wave=wave*wave_unit,flux=d*flux_unit,header=header)
		else:
			sp = Spectrum(wave=d[0,:]*wave_unit,flux=d[1,:]*flux_unit,header=header)
			if len(d[:,0])>2:
				for i in range(2,numpy.min([len(columns),len(d[:,0])])):
					setattr(sp,columns[i],d[i,:])
		sp.clean()

	else:
# csv/tab/ascii
		if ftype in ['txt','dat','asc','csv'] or file_type in ['ascii','space','csv']:
			pd = pandas.read_csv(filename,delimiter=delimiter,comment=comment,names=columns,index_col=False)
# excel
		elif ftype=='xls' or ftype=='xlsx' or file_type=='excel':
			pd = pandas.read_excel(filename,comment=comment,names=columns,index_col=False)
		else:
			raise ValueError('Do not recognize file type {}'.format(filename))
#		for i,c in enumerate(list(pd.columns)): pd.rename(columns={c:column_order[i]},inplace=True)
		sp = Spectrum(wave=list(pd['wave'])*wave_unit,flux=list(pd['flux'])*flux_unit)
		for c in list(pd.columns):
			if c=='wave' or c=='flux': pass
			elif c=='unc' or c=='background': setattr(sp,c,list(pd[c])*flux_unit)
			else: setattr(sp,c,list(pd[c]))
		sp.clean()
		# 	sp.unc = list(pd['unc'])*flux_unit
		# if 'background' in list(pd.columns): sp.unc = list(pd['unc'])*flux_unit
		# if pd.ndim>2:
		# 	for i in range(2,pd.ndim):
		# 		setattr(sp,column_order[i],list(pd[column_order[i]]))

		# cols = list(pd.columns)
		# if isNumber(cols[0]):
		# 	pd = pandas.read_csv(filename,delimiter=delimiter,comment=comment,names=column_order[:len(cols)])
		# sp = Spectrum(wave=list(pd[cols[0]])*wave_unit,flux=list(pd[cols[1]])*flux_unit)
		# if pd.ndim>2:
		# 	for i in range(2,pd.ndim):
		# 		setattr(sp,column_order[i],list(pd[cols[i]]))

# other parameters
	for k in list(kwargs.keys()): setattr(sp,k.lower(),kwargs[k])
	return sp



############################################################
# KAST IMAGE DATA FUNCTIONS
# These functions allow for I/O with KAST data files, and image manipulation
############################################################



def readFiles(num,folder='./',mode='RED',prefix='',suffix='',rotate=False,trim=[],verbose=ERROR_CHECKING):
	'''
	:Purpose:

		Reads in a KAST fits file(s) and returns 2D image & header

	:Required Inputs: 

		:param num: file number

	:Optional Inputs: 

		:param folder='./': data folder containing data
		:param mode='RED': data mode (RED or BLUE)

	:Output: 

		If one file number provided, a 2D numpy array containing image and dictionary containing header
		If more than one file number provided, a 3D numpy array of [n files x n rows x n columns] and array of header dictionaries

	:Usage: 

		TBD

	'''
# variables
	images,hds = [],[]


# generate list of files depending on input format	
	nlist = copy.deepcopy(num)
	if not isinstance(nlist,list): nlist = [num]
	files = []
# files provided	
	if isinstance(nlist[0],str): 
		if os.path.exists(nlist[0]): files = nlist
		elif os.path.exists('{}/{}'.format(folder,nlist[0])): files = ['{}/{}'.format(folder,n) for n in nlist]
# file numbers syntax provided - need mode or prefix as input
		else:
			try:
				nlist = numberList(nlist[0])
			except: 
				raise ValueError('Could not located data files from input {}'.format(num))

# generate file names from numbers using instrument defaults
	if len(files) == 0:

# check instrument mode and set defaults
		if checkDict(mode,INSTRUMENT_MODES) == False: 
			print('Warning: unknown instrument mode {}'.format(mode))
			mode = ''
		else: 
			mode = checkDict(mode,INSTRUMENT_MODES)
			if prefix=='': prefix = INSTRUMENT_MODES[mode]['PREFIX']
			if suffix=='': suffix = INSTRUMENT_MODES[mode]['SUFFIX']
			rotate = INSTRUMENT_MODES[mode]['ROTATE']
			trim = INSTRUMENT_MODES[mode]['TRIM']
#		files = ['{}/{}{}.fits'.format(folder,prefix,str(n).zfill(4)) for n in nlist]
		files = ['{}/{}{}{}.fits'.format(folder,prefix,str(n),suffix) for n in nlist]

# read in files
	for f in files:
		if not os.path.exists(f): print('Warning: cannot find data file {}; skipping'.format(f))
		else:
			hdulist = fits.open(f)
			im = hdulist[0].data
			hdr = hdulist[0].header
			hdulist.close()
# try to infer header from file
			if mode=='':
				try: 
					mode = kastRBMode(hdr)
					rotate = INSTRUMENT_MODES[mode]['ROTATE']
					trim = INSTRUMENT_MODES[mode]['TRIM']
				except: 
					pass
# rotate and trim
			if rotate==True: im = numpy.rot90(im,k=1) # only one direction, may need to generalize more
			if len(trim) > 0: im = im[trim[0]:trim[1],:]
# subtract overscan: TBD
#	ims = subtractOverscan(im)
			images.append(im.astype(float))
			hds.append(hdr)
	if len(files) == 1: return images[0],hds[0]
	else: return images,hds


def combineImages(imarr,method='median',axis=0,sclip=5.,verbose=ERROR_CHECKING,**kwargs):
	'''
	:Purpose:

		Takes an array of 2D images and combines them with a give statistic

	:Required Inputs: 

		:param imarr: by default, a 3D array of N images, where N is length of ``axis''; 
		alternately this can be a string specifying the file numbers or 1D array of file numbers, 
		with images read in using readFiles()

	:Optional Inputs: 

		:param method='median': method for combining images (more than one of these can be included)
			* 'median': median combine
			* 'add': add all pixels
			* 'average: average value of all pixels
			* 'sigclip': sigma clipping, rejecting outliers N x std deviation, where N is specified by input ``sclip''
		:param axis=0.: axis upon which images are combined
		:param sclip=5.: number of standard deviations to apply signma clipping

		If reading in files, additional parameters for readFiles() can be included though **kwargs

	:Output: 

		Combined 2D image

	:Usage: 

		TBD

	'''
	images = copy.deepcopy(imarr)
	if isinstance(imarr,str) or len(numpy.shape(imarr))==1:
		try: 
			images,hds = readFiles(imarr,**kwargs)
		except: 
			raise ValueError('Cannot process image input {}'.format(imarr))
	if len(numpy.shape(images))==2:
		if verbose==True: print('Warning: image array is just a single image of dimensions {}; returning'.format(numpy.shape(images)))
		return images

	if 'sigclip' in method.lower(): std = numpy.nanstd(images,axis=axis)
	if 'median' in method.lower(): f = numpy.nanmedian
	elif 'add' in method.lower(): f = numpy.nansum
	elif 'average' in method.lower() or 'mean' in method.lower(): f = numpy.nanmean
	else:
		raise ValueError('Combine method {} has not be implemented; returning median'.format(method))
		return numpy.nanmedian(images,axis=axis)

	output = f(images,axis=axis)
	if 'sigclip' in method.lower(): 
		for i in len(output[:,0]):
			for j in len(output[0,:]):
# note: this only works if axis=0				
				sl = numpy.resize(images[:,i,j],len(images[:,i,j]))
				output[i,j] = f(sl[numpy.absolute(sl-output[i,j])<sclip*std[i,j]])
	return output


# combine images with cosmic ray rejection
def crRejectCombine(ims,sclip=1.,nfitcycle=5,rejlimit=0.01,scliplimit=10.,smooth=False,smooth_scale=3,verbose=ERROR_CHECKING):
	msks = []
	for i in range(len(ims)):
		diff = ims[i]-ims[numpy.mod(i+1,len(ims))]
		if smooth==True:
			diffs = scipy.ndimage.filters.uniform_filter(diff,int(smooth_scale))
			diffc = diff-diffs
		else: diffc = diff
# now do a rejection loop, limiting the fraction of rejected pixels by upping the sigma clip scale to a maximum value		
		rejflag = False
		sc = copy.deepcopy(sclip)
		while rejflag==False:
			msk = diff*0
			for j in range(nfitcycle):
				std,md = numpy.nanstd(diffc[msk==0]),numpy.nanmedian(diffc[msk==0])
				msk[(diffc-md)>(sc*std)] = 1
			if numpy.sum(msk)/numpy.size(msk) > rejlimit: sc=sc+1.
			else: rejflag=True
			if sc>scliplimit: 
				msk = diff*0
				rejflag=True
		msks.append(msk)
		if verbose==True: print('CR Reject on image {} rejected {} pixels, {}\% of total'.format(i,numpy.sum(msk),100.*numpy.sum(msk)/numpy.size(msk)))
	msks = numpy.array(msks)
	imc = numpy.nanmedian(ims,axis=0)
	imcc = copy.deepcopy(imc)
	mskc = numpy.nansum(msks,axis=0)
	for i in range(len(ims)):
		imcc[msks[i]==1] = ims[numpy.mod(i+1,len(ims))][msks[i]==1]
	imcc[mskc==len(ims)] = numpy.nanmedian(imcc)
	return imcc


# clean images using a mask
def maskClean(image,mask,replace=0.,verbose=ERROR_CHECKING):
	'''
	:Purpose:

		Applies a mask to clean an image

	:Required Inputs: 

		:param image: 2D array
		:param mask: 2D array where 0 = good and 1 = mask

	:Optional Inputs: 

		:param replace=0.: replace bad pixels with this value

	:Output: 

		Cleaned 2D image

	:Usage: 

		TBD

	'''
# process image input for flexibility
# string or list of file numbers  
	imout = copy.deepcopy(image)
	if not isinstance(imout,numpy.ndarray): imout = numpy.array(imout)
	imout[mask==1] = replace
	return imout

def subtractOverscan():
	'''
	Subtracts overscan bias voltage
	Input: image
	Output: image with overscan removed
	'''
	pass

def computeVariance():
	'''
	Produces a variance image based on basic noise model
	Input: image in units of e-
	Output: variance
	'''
	pass

############################################################
# KAST IMAGE REDUCTION FUNCTIONS
# These functions contain primary reduction steps for imaging data
############################################################


def readInstructions(file,comment='#',verbose=ERROR_CHECKING):
	'''
	:Purpose:

		Reads in an reduction instruction file and place parameters into a dictionary

	:Required Inputs: 

		:param file: full path filename to instruction file, which is a set of tab-separated key and value pairs

	:Optional Inputs: 

		:param comment='#': character(s) prefixing a line that is a comment

	:Output: 

		Dictionary containing input parameters

	:Usage: 

		The instruction file should have the format as follows:

			# A comment
			MODE	RED
			DATE	20190918
			DATA_FOLDER	/Users/adam/data/kast/190918/data/
			REDUCTION_FOLDER	/Users/adam/data/kast/190918/reduction/
			BIAS	FILES=2003-2032
			FLAT	FILES=2037-2066
			ARC_SHALLOW	FILES=2000	LAMPS=AR,HG,NE
			ARC_DEEP	FILES=2001	LAMPS=AR,HG,NE
			SOURCE	NAME=J1833+2225	FILES=2068-2069	TELLURIC=HD177082
			TELLURIC	NAME=HD177082	FILES=2070-2071	SPT=G2V
			FLUXCAL	NAME=FIEGE110	FILES=2080	TELLURIC=HD219833
			TELLURIC	NAME=HD219833	FILES=2082	SPT=G2V
			OBSERVERS	Roman, Christian, Adam

		The required keywords are 'DATA_FOLDER','MODE','BIAS','FLAT','ARC_SHALLOW','ARC_DEEP', and 'EXPORT'
		Note that sources and tellurics are matched through their name fields

		[example of usage]

	'''
	parameters = {}
	required = ['DATA_FOLDER','MODE','BIAS','FLAT','DATE']

	if not os.path.exists(file): raise ValueError('Cannot find instruction file {}'.format(file))
	f = open(file,'r')
	lines = f.readlines()
	f.close()
	for ln in lines:
		if ln[0]!=comment:
			parts = ln.split('\t')
			ref = parts[0].upper().strip()
			if ref=='MODE': 
				parameters[ref] = parts[1].upper().strip()
				if checkDict(parameters[ref],INSTRUMENT_MODES) == False: raise ValueError('Mode {} not defined'.format(parameters[ref]))
				else: parameters[ref]=checkDict(parameters[ref],INSTRUMENT_MODES)
			elif ref in ['BIAS','FLAT','ARC','ARC_SHALLOW','ARC_DEEP']:
				dt = {}
				for p in parts[1:]:
					kv = p.split('=')
					dt[kv[0].upper().strip()] = kv[1].upper().strip()
				if 'FILES' not in list(dt.keys()): print('Warning: FILES keyword must be included in instruction file line for {}'.format(ref))
				else: dt['FILES'] = numberList(dt['FILES'])
				dt['REUSE'] = False
				if 'FILE' in list(dt.keys()): dt['REUSE'] = True
				if 'BACK' in list(dt.keys()): 
					dt['BACK'] = dt['BACK'].replace('[','').replace(']','').replace('(','').replace(')','').split(',')
					dt['BACK'] = [int(x) for x in dt['BACK']]
				if 'WINDOW' in list(dt.keys()): dt['WINDOW'] = int(dt['WINDOW'])
				if 'METHOD' not in list(dt.keys()): dt['METHOD'] = 'optimal'
				parameters[ref] = dt
			elif ref in ['SOURCE','TELLURIC','FLUXCAL']:
				if ref not in (parameters.keys()): parameters[ref] = {}
				dt = {}				
				for p in parts[1:]:
					kv = p.split('=')
					if len(kv)>1: dt[kv[0].upper().strip()] = kv[1].upper().strip()
				if 'NAME' not in list(dt.keys()): raise ValueError('information list for {} requires NAME keyword; current lists gives {}'.format(ref,dt))
				if 'FILES' not in list(dt.keys()): raise ValueError('information list for {} requires FILES keyword; current lists gives {}'.format(ref,dt))
				dt['FILES'] = numberList(dt['FILES'])
				if 'BACK' in list(dt.keys()): 
					dt['BACK'] = dt['BACK'].replace('[','').replace(']','').replace('(','').replace(')','').split(',')
					dt['BACK'] = [int(x) for x in dt['BACK']]
				if 'WINDOW' in list(dt.keys()): dt['WINDOW'] = int(dt['WINDOW'])
				if 'METHOD' not in list(dt.keys()): dt['METHOD'] = 'optimal'
				if ref=='TELLURIC' and 'SPT' not in list(dt.keys()): dt['SPT'] = 'G2V'
				parameters[ref][dt['NAME']] = dt
			elif ref in ['WAVE_INITIAL']:
				dt = {}				
				for p in parts[1:]:
					kv = p.split('=')
					if len(kv)>1: dt[kv[0].upper().strip()] = [float(x) for x in kv[1].split(',')]
				if 'WAVE' not in list(dt.keys()): raise ValueError('information list for {} requires WAVE keyword; current lists gives {}'.format(ref,dt))
				if 'PIXEL' not in list(dt.keys()): raise ValueError('information list for {} requires PIXEL keyword; current lists gives {}'.format(ref,dt))
				parameters[ref] = dt
			else: parameters[ref] = parts[1].strip()

# place holder for reusing existing calibration files
	if 'BIAS_REUSE' in list(parameters.keys()): parameters['BIAS'] = {'FILES': '0-1'}
	if 'FLAT_REUSE' in list(parameters.keys()): parameters['FLAT'] = {'FILES': '0-1'}

	for r in required:
		if r not in list(parameters.keys()): print('WARNING: required parameter {} not in instruction file'.format(r))
	if 'REDUCTION_FOLDER' not in list(parameters.keys()): parameters['REDUCTION_FOLDER'] = parameters['DATA_FOLDER']+'/../reduction/'
	if 'OUTPUT_FOLDER' not in list(parameters.keys()): parameters['OUTPUT_FOLDER'] = parameters['REDUCTION_FOLDER']
	if 'ARC_SHALLOW' not in list(parameters.keys()) and 'ARC' in list(parameters.keys()): parameters['ARC_SHALLOW'] = parameters['ARC']
	if 'ARC_DEEP' not in list(parameters.keys()) and 'ARC' in list(parameters.keys()): parameters['ARC_DEEP'] = parameters['ARC']
	if 'WAVE_INITIAL' not in list(parameters.keys()): parameters['WAVE_INITIAL'] = {'WAVE': [], 'PIXEL': []}
	for r in ['ARC_SHALLOW','ARC_DEEP']:
		if r not in list(parameters.keys()): print('WARNING: required parameter {} not in instruction file'.format(r))

	return parameters



def makeBias(files,method='median',folder='./',mode='',prefix='',overwrite=True,verbose=ERROR_CHECKING,output='',input_headers=None):
	'''
	:Purpose:

		Creates a bias frame by median combining bias frame files

	:Required Inputs: 

		:param files: file numbers to combine, as either string sequence or 1D array

	:Optional Inputs: 

		:param method='median': method for combining files; see ``combineImages()''
		:param folder='./': data folder containing data
		:param mode='RED': data mode (RED or BLUE)

	:Output: 

		2D image array file for bias frame and header file from first image

	:Usage: 
		
		Note: currently this is only a pass through to ``combineImages()''; could include additional checks

	'''
	if isinstance(files[0],str) or isinstance(files[0],int) or isinstance(files[0],float): 
		images,hds = readFiles(files,folder=folder,mode=mode,prefix=prefix)
	elif isinstance(files[0],numpy.ndarray):
		images = copy.deepcopy(files)
		hds = copy.deepcopy(input_headers)
		if not isinstance(hds,list): raise ValueError('If passing an image array, be sure to pass the list of image headers in input_headers keyword')
		if len(hds) != len(images): raise ValueError('list of input headers is not the same length as list of images')
	else:
		raise ValueError('Input should be a list of files or numpy 2D arrays')

	if len(numpy.shape(images))==2:
		if verbose==True: print('Warning: image array is just a single image of dimensions {}; returning'.format(numpy.shape(images)))
		return images,hds

	im = combineImages(images,method=method)

	if output == '': output='bias_{}.fits'.format(mode)
	if verbose==True: print('Writing bias frame to {}'.format(output))
	try:
		hdu = fits.PrimaryHDU(data=im,header=hds[0])
		hdu.writeto(output,overwrite=overwrite)
	except:
		print('Warning: could not write to file {}'.format(output))

	return im,hds[0]


def makeFlat(files,bias,method='median',quantile=0.9,folder='./',mode='',prefix='',verbose=ERROR_CHECKING,overwrite=True,output='',input_headers=None):
	'''
	:Purpose:

		Creates a bias frame by median combining bias frame files

	:Required Inputs: 

		:param files: file numbers to combine, as either string sequence or 1D array
		:param bias: 2D bias image

	:Optional Inputs: 

		:param method='median': method for combining files; see ``combineImages()''
		:param quantile=0.9: quantile to use for normalizing flat
		:param folder='./': data folder containing data
		:param mode='RED': data mode (RED or BLUE)

	:Output: 

		2D image array file for flat field frame and header from first file

	:Usage: 
		
		TBD

	'''
	if isinstance(files[0],str) or isinstance(files[0],int) or isinstance(files[0],float): 
		images,hds = readFiles(files,folder=folder,mode=mode,prefix=prefix)
	elif isinstance(files[0],numpy.ndarray):
		images = copy.deepcopy(files)
		hds = copy.deepcopy(input_headers)
		if not isinstance(hds,list): raise ValueError('If passing an image array, be sure to pass the list of image headers in input_headers keyword')
		if len(hds) != len(images): raise ValueError('list of input headers is not the same length as list of images')
	else:
		raise ValueError('Input should be a list of files or numpy 2D arrays')

	if len(numpy.shape(images))==2:
		if verbose==True: print('Warning: image array is just a single image of dimensions {}; returning'.format(numpy.shape(images)))
		return images,hds

# read in bias if necessary
	if isinstance(bias,str): 
		if os.path.exists(bias): bias,bhd = readFiles(bias)
		elif os.path.exists(folder+'/'+bias): bias,bhd = readFiles(folder+'/'+bias)
		else: raise ValueError('Cannot find bias file {}'.format(bias))
		if kastRBMode(bhd) != kastRBMode(hds[0]): raise ValueError('Red/blue mode of bias frame is {} while that of flat field frames is {}'.format(kastRBMode(bhd),kastRBMode(hds[0])))

	images = [f-bias for f in images]
	flat = combineImages(images,method=method)
	im = flat/numpy.quantile(flat,0.9)

	if output == '': output='flat_{}.fits'.format(mode)
	if verbose==True: print('Writing flat field frame to {}'.format(output))
	try:
		hdu = fits.PrimaryHDU(data=im,header=hds[0])
		hdu.writeto(output,overwrite=overwrite)
	except:
		print('Warning: could not write to file {}'.format(output))

	return im,hds[0]


def makeMask(bias,flat,sclip=3.,flatmin=0.01,mode='',verbose=ERROR_CHECKING,overwrite=True,output=''):
	'''
	Creates a mask file: 0 = OK, 1 = BAD
	Input: flat, bias
	Output: mask image
	'''

# read in bias and flat if necessary
	if isinstance(bias,str): 
		if os.path.exists(bias): bias,bhd = readFiles(bias)
		elif os.path.exists(folder+'/'+bias): bias,bhd = readFiles(folder+'/'+bias)
		else: raise ValueError('Cannot find bias file {}'.format(bias))
		if mode=='': mode = kastRBMode(bhd)
	if isinstance(flat,str): 
		if os.path.exists(flat): flat,fhd = readFiles(flat)
		elif os.path.exists(folder+'/'+flat): flat,fhd = readFiles(folder+'/'+flat)
		else: raise ValueError('Cannot find flat file {}'.format(flat))
		if mode=='': mode = kastRBMode(fhd)
	if mode=='': 
		print('Warning: no red/blue mode provided for mask, so unclear which setting this is for')
		mode = 'UNKNOWN'

	mask = copy.deepcopy(bias)*0
# remove negative values	
	mask[flat<0.] = 1
	mask[bias<0.] = 1
# mask out pixels outside of illuminated region
	mask[flat<(flatmin*numpy.quantile(flat[mask==0],0.9))] = 1
#	profile = numpy.nanmedian(flat,axis=axis)
#	notorder = numpy.arange(len(profile))
#	notorder = notorder[profile<(0.01*numpy.nanmax(profile))]
#	print(notorder)
#	for n in notorder: mask[:,n] = 1
#	for n in notorder: print(n,mask[:10,n])

# mask out discrepant pixels
	mask[numpy.absolute(flat-numpy.nanmedian(flat[mask==0]))>sclip*numpy.nanstd(flat[mask==0])] = 1
	mask[numpy.absolute(bias-numpy.nanmedian(bias[mask==0]))>sclip*numpy.nanstd(bias[mask==0])] = 1
#	mask[redux['BIAS']>(5.*numpy.nanstd(redux['BIAS'])+numpy.nanmedian(redux['BIAS']))] = 1

	if output == '': output='mask_{}.fits'.format(mode)
	if verbose==True: print('Writing mask frame to {}'.format(output))
	try:
		hdu = fits.PrimaryHDU(data=mask)
		hdu.writeto(output,overwrite=overwrite)
	except:
		print('Warning: could not write to file {}'.format(output))

	return mask


def reduceScienceImage(image, bias, flat, mask=[], hd={}, mode='', rdmode='', gain=0., rn=0., exposure='EXPTIME', mask_image=True, folder='./', verbose=ERROR_CHECKING):
	'''
	:Purpose:

		Performs the following image reduction steps on a science image: 
		* subtracts bias
		* converts to e/s
		* computes variance
		* divides by normalized flat field
		* masks bad pixels

	:Required Inputs: 

		:param images: either full filename or file number of the science frame
		:param bias: 2D bias image
		:param flat: 2D flat image

	:Optional Inputs: 

		:param mask=None: set to a 2D mask image if masking is desired
		:param mask_image=True: set to True to mask clean a science image
		:param exposure='EXPTIME': either the header keyword for exposure time or exposure in seconds
		:param folder='./': data folder containing data
		:param mode='RED': data mode (RED or BLUE)

	:Output: 

		2D image array file for flat field frame

	:Usage: 
		
		TBD

	'''
# read in image if necessary
	im = copy.deepcopy(image)	
	if isinstance(image,str) or isinstance(image,int): 
		im,hd = readFiles(image,folder=folder,mode=mode)
# CCD keywords
	if mode=='':
		try: mode = kastRBMode(hd)
		except: raise ValueError('Cannot determine red/blue mode of this image; be sure to pass this value or the original header')
	if rdmode=='':
		try: rdmode = kastReadMode(hd)
		except: raise ValueError('Cannot determine read mode of this image; be sure to pass this value or the original header')
	if gain==0.:
		try: gain = kastGain(hd)
		except: raise ValueError('Cannot determine gain of this image; be sure to pass this value or the original header')
#		if verbose==True: print('Using gain of {:.1f}'.format(gain))
	if rn==0.:
		try: rn = kastRN(hd)
		except: raise ValueError('Cannot determine read noise of this image; be sure to pass this value or the original header')
#		if verbose==True: print('Using read noise of {:.1f}'.format(rn))
	if isinstance(exposure,str):
		expkey = copy.deepcopy(exposure) 
		try: exposure = kastHeaderValue(hd,expkey)
		except: raise ValueError('Could not identify exposure keyword {} in header; you may need to pass a numerical value or the original header'.format(exposure))
	if exposure==0.:
		raise ValueError('Module inferred exposure time of {}; try passing a different exposure keyword'.format(exposure))

# read in bias and flat if necessary
	if isinstance(bias,str): 
		if os.path.exists(bias): bias,bhd = readFiles(bias)
		elif os.path.exists(folder+'/'+bias): bias,bhd = readFiles(folder+'/'+bias)
		else: raise ValueError('Cannot find bias file {}'.format(bias))
		if kastRBMode(bhd) != kastRBMode(hd): raise ValueError('Red/blue mode of bias frame is {} while that of science frames is {}'.format(kastRBMode(bhd),kastRBMode(hd)))
	if isinstance(flat,str): 
		if os.path.exists(flat): flat,fhd = readFiles(flat)
		elif os.path.exists(folder+'/'+flat): flat,fhd = readFiles(folder+'/'+flat)
		else: raise ValueError('Cannot find flat file {}'.format(flat))
		if kastRBMode(fhd) != kastRBMode(hd): raise ValueError('Red/blue mode of flat frame is {} while that of science frames is {}'.format(kastRBMode(bhd),kastRBMode(hd)))

# process image
	if numpy.shape(bias)==numpy.shape(im): image_b = im-bias
	else: 
		print('WARNING! bias image is not the same shape ({}) as science image ({}); ignoring bias subtraction'.format(numpy.shape(bias),numpy.shape(im)))
		image_b = copy.deepcopy(im)
	image_e = image_b*gain
# mask flat 
	if len(mask)==0: mask = image_e*0
	if numpy.shape(flat)==numpy.shape(mask): fm = maskClean(flat,mask,replace=1.)
	else: 
		print('WARNING! flat image is not the same shape ({}) as mask image ({}); ignoring mask clean'.format(numpy.shape(flat),numpy.shape(mask)))
		fm = copy.deepcopy(flat)
# divide image by flat 
	if numpy.shape(fm)==numpy.shape(image_e): 
		image_f = image_e/(fm*exposure)
		variance = (image_e/fm+rn**2)/(exposure**2)
	else: 
		print('WARNING! flat image is not the same shape ({}) as science image ({}); ignoring flat fielding'.format(numpy.shape(fm),numpy.shape(image_e)))
		image_f = image_e/exposure
		variance = (image_e+rn**2)/(exposure**2)
# mask image
	if mask_image==True:
		if numpy.shape(image_f)==numpy.shape(mask):
			image_f = maskClean(image_f,mask,replace=numpy.nanmedian(image_f))	
			variance = maskClean(variance,mask,replace=numpy.nanmax(variance))	
		else: 
			print('WARNING! mask image is not the same shape ({}) as science image ({}); ignoring mask clean'.format(numpy.shape(mask),numpy.shape(image_f)))

	return image_f,variance


############################################################
# KAST SPECTRAL REDUCTION FUNCTIONS
# These functions contain primary reduction steps for spectral data
############################################################

def findPeak(im,rng=[],cntr=-1,window=50,trace_slice=[],method='maximum',verbose=ERROR_CHECKING,plot_file=''):
	'''
	Finds a source spatial peak in a median combined spectral stack of an image
	'''
# find center of profile trace
	if len(rng)==0 and cntr>=0 and window!=None: rng=[int(cntr-window),int(cntr+window)]
	if len(rng)==0: rng=[0,len(im[:,0])]
	if len(trace_slice)==0: trace_slice = [150,len(im[0,:])-150]
	sprofile = numpy.nanmedian(im[:,trace_slice[0]:trace_slice[1]],axis=1) 
	profile = numpy.nanmedian(im[rng[0]:rng[1],trace_slice[0]:trace_slice[1]],axis=1)

# only doing maximum now
#	if method=='maximum': 
	cntr = numpy.argmax(profile)+rng[0]

	if plot_file != '':
		yrng = [numpy.nanmedian(sprofile),numpy.nanmax(profile)]
		plt.clf()
		plt.figure(figsize=[8,4])
		x = numpy.arange(len(profile))+rng[0]
		plt.plot(x,profile,'k-')
		plt.plot([cntr,cntr],yrng,'m--')
		plt.legend(['Profile Cut','Center at {:.0f}'.format(cntr)])
		plt.plot(numpy.arange(len(sprofile)),sprofile,'k-')
		plt.plot([rng[0],rng[0]],yrng,'m--')
		plt.plot([rng[1],rng[1]],yrng,'m--')
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('X Pixel',fontsize=16)
		plt.ylabel('Counts',fontsize=16)
		plt.ylim(yrng)
		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			if verbose==True: print('Warning: could not save flux calibration diagnostics to file {}'.format(plot_file))
#		plt.clf()

# just using maximum; could do something smarter
#	if method=='maximum': 
	return cntr


def traceDispersion(im,cntr=-1,window=5,step_size=5,trace_slice=[],fit_order=3,fitcycle=5,sigclip=3,method='maximum',verbose=ERROR_CHECKING,plot_file=''):
	'''
	Traces dispersion

	NEED TO FIX REGIONS THAT ARE MASKED
	'''
# fix trace_slice if not provided
	if len(trace_slice)==0: trace_slice=[150,len(im[0,:])-150]
# if no center provided, find the peak
	if cntr<0: 
		cntr = findPeak(im,window=window,trace_slice=trace_slice,method=method)
	if method=='nofit' or method=='single': return numpy.zeros(len(im[0,:]))+cntr

# tracing peaks; find peak closest to stacked peak
	xs = numpy.arange(trace_slice[0]+step_size,trace_slice[1]-step_size,step_size)
	ys = []
	for x in xs:
		profile = numpy.nanmedian(im[int(cntr-window):int(cntr+window),int(x-step_size):int(x+step_size)],axis=1)
#		if method=='maximum': 
		ys.append(numpy.argmax(profile)+int(cntr-window))
	ys = numpy.array(ys)

# fit this to a polynomial through a fit cycle
	w = numpy.where(xs>0)
	for i in range(fitcycle):
		p = numpy.polyfit(xs[w],ys[w],fit_order)
		diff = ys-numpy.polyval(p,xs)
		w = numpy.where(numpy.absolute(diff)<sigclip*numpy.nanstd(diff))
		rms = numpy.nanstd(diff[w])
		if verbose==True: print('Dispersion fit pass {}: RMS={:.2f} pixels for {} peaks and fit order {}'.format(i,rms,len(xs[w]),fit_order))

# plot outcome
	if plot_file != '':
		plt.clf()
		plt.figure(figsize=[8,6])
		plt.subplot(211)
		plt.plot(xs,ys,'k.',alpha=0.5)
		plt.plot(xs[w],ys[w],'bo')
		xim = numpy.arange(len(im[0,:]))
		plt.plot(xim,numpy.polyval(p,xim))
		plt.legend(['Centers','Fit Centers','Trace Fit (o={})'.format(fit_order)])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('X Pixel',fontsize=16)
		plt.ylabel('Y Pixel',fontsize=16)
		plt.ylim([numpy.nanmin(ys[w]),numpy.nanmax(ys[w])])
		plt.subplot(212)
		plt.plot(xs,ys-numpy.polyval(p,xs),'k.')
		plt.plot(xs[w],ys[w]-numpy.polyval(p,xs[w]),'bo')
		plt.plot([0,numpy.nanmax(xim)],[0,0],'k--')
		plt.plot([0,numpy.nanmax(xim)],[rms,rms],'k--')
		plt.plot([0,numpy.nanmax(xim)],[-1.*rms,-1.*rms],'k--')
		plt.legend(['O-C','RMS={:.3f}'.format(rms)])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('X Pixel',fontsize=16)
		plt.ylabel('O-C',fontsize=16)
		plt.ylim([-3.*rms,3.*rms])
		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			print('Warning: could not save flux calibration diagnostics to file {}'.format(plot_file))
#		plt.clf()

	return numpy.polyval(p,numpy.arange(len(im[0,:])))


def extract2D(im,trace,start=0,window=5,plot_file=''):
	'''
	extracts a 2D spectrum centered on trace or central column
	THIS FUNCTION HAS BEEN SUPERSEDED BY RECTIFY
	'''
	stack = []
	tr = numpy.array([int(t) for t in trace])

	for i,t in enumerate(tr):
		stack.append(im[(t-int(window)):(t+int(window)+1),int(start)+i])

	return stack


def spatialProfile(im,cntr=-1,window=10,verbose=ERROR_CHECKING,plot_file=''):
	'''
	generates spatial profile for optimal extraction
	NOTE: assumes it is provided a rectified image
	'''
	profile = numpy.nanmedian(im,axis=1)
	if cntr<0: cntr = numpy.argmax(profile)
	cntr = int(cntr)
	window = int(window)
	profile = profile[(cntr-window):(cntr+window+1)]
	profile = profile-numpy.nanmin(profile)
	profile = profile/numpy.nanmax(profile)

	if plot_file != '':
		yrng = [0,1.2]
		plt.clf()
		plt.figure(figsize=[8,4])
		x = numpy.arange(len(profile))+cntr-0.5*(len(profile)-1)
		plt.plot(x,profile,'k-')
		plt.plot([cntr,cntr],yrng,'m--')
		plt.legend(['Profile Cut','Center at {:.0f}'.format(cntr)])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('X Pixel',fontsize=16)
		plt.ylabel('Normalized Profile',fontsize=16)
		plt.ylim(yrng)
		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			print('Warning: could not save profile plot to file {}'.format(plot_file))
#		plt.clf()
	return profile



def extractSpectrum(im,var=[],mask=[],method='optimal',cntr=-1,profile=[],src_wnd=5,bck_wnd=[10,20],subtract_background=True,verbose=ERROR_CHECKING,plot_file='',center_kwargs={},**kwargs):
	'''
	Purpose:
		Extracts spectrum from science image; can be done with pre-determined trace, 
		and with various options for extraction depending on use

	Required Inputs: 
		image: 

	Optional Inputs: 
		variance = []: variance image
		method = 'optimal':
		mask = []: 
		spatial_rng = []: 
		start = 0:
		cntr = None:
		src_wnd = 5:
		bck_wnd = [10,20]:
		trace = []:
		trace_slice = []:
		shift_trace = False:
		profile = []: 
		plot_file = '':

	Output: spectrum object of each source, which includes flux, variance, background
	'''

# assert integerness of all indices
	src_wnd = int(src_wnd)
	bck_wnd = [int(i) for i in bck_wnd]
	if bck_wnd[0] <= src_wnd: bck_wnd = [i+src_wnd-bck_wnd[0]+1 for i in bck_wnd]
#	if len(trace_slice)==0: trace_slice = [150,len(im[0,:])-150]
#	trace_slice = [int(i) for i in trace_slice]
	if len(mask) == 0: mask = im*0
	if len(var) == 0: var = im*0

# identify spatial location of source if not provided
	if cntr<=0.: cntr=findPeak(im,verbose=verbose,**center_kwargs)   
	cntr = int(cntr)

#	if len(trace)==0 and shift_trace==False:
#		trace = traceDispersion(im,cntr=cntr,window=src_wnd,trace_slice=trace_slice,method='maximum')
# or shift trace from prior extraction
#	elif len(trace)>0 and shift_trace==True:
#		trace = numpy.array(trace)+cntr-numpy.nanmedian(trace[trace_slice[0]:trace_slice[1]])
#	elif len(trace)==0 and shift_trace==True:
#		raise ValueError('You must provide a trace to shift it to the current source')
#	else: pass
#	trc = [int(t) for t in trace]

# some method presets
#	if verbose==True: print('input profile length =  {}'.format(len(profile)))
	if method.lower()=='arc':
		subtract_background = False
	if method.lower()=='optimal' and len(profile)==0:
		profile = spatialProfile(im,cntr=cntr,window=src_wnd)
	elif method.lower()=='boxcar' or len(profile)==0:
		profile = numpy.ones(2*src_wnd+1)
	else:
		if verbose==True: print('Extraction method {} unrecognized; using boxcar'.format(method))
		profile = numpy.ones(2*src_wnd+1)
# force agreement between source window and extraction profile so profile can be used as input
	if len(profile)>(2*src_wnd+1): 
		profile = profile[int(0.5*(len(profile)-1)-src_wnd):int(0.5*(len(profile)-1)+src_wnd+1)]
	if len(profile)<(2*src_wnd+1): 
		src_wnd = int(0.5*(len(profile)-1))
	if verbose==True: print('Using extraction method {}'.format(method))

# compute background and subtract from extraction region
	bck = []
	for i in range(len(im[0,:])):
		bck1 = im[(cntr-bck_wnd[1]):(cntr-bck_wnd[0]),i]
		bck2 = im[(cntr+bck_wnd[0]):(cntr+bck_wnd[1]),i]
		bck12 = numpy.append(bck1,bck2)
# TBD: what about variance of background?
		# vbck1 = var[(cntr-bck_wnd[1]):(cntr-bck_wnd[0]),i]
		# vbck2 = var[(cntr+bck_wnd[0]):(cntr+bck_wnd[1]),i]
		# vbck12 = numpy.append(vbck1,vbck2)
		mbck1 = mask[(cntr-bck_wnd[1]):(cntr-bck_wnd[0]),i]
		mbck2 = mask[(cntr+bck_wnd[0]):(cntr+bck_wnd[1]),i]
		mbck12 = numpy.append(mbck1,mbck2)
		bck12m = bck12[mbck12!=1]
		if len(bck12m) == 0: bck.append(numpy.nanmedian(bck12))
		else: bck.append(numpy.nanmedian(bck12m))
#		if numpy.isnan(bck[-1]): bck[-1] = numpy.nanmedian(bck12)
	bck = numpy.array(bck)

# subtract background from extraction region
	imsub = copy.deepcopy(im)
	if subtract_background==True: 
		for i in range(cntr-bck_wnd[1],cntr+bck_wnd[1]): imsub[i,:] = imsub[i,:]-bck


# extract spectrum and uncertainty
	flx,unc = [],[]
#	print(numpy.shape(imsub),cntr,src_wnd)
	for i in range(len(im[0,:])):
		src = imsub[(cntr-src_wnd):(cntr+src_wnd+1),i]
		vsrc = var[(cntr-src_wnd):(cntr+src_wnd+1),i]
		msrc = mask[(cntr-src_wnd):(cntr+src_wnd+1),i]
# mask profile
		prf = copy.deepcopy(profile)
		prf[msrc==1] = 0.
# catch if we mask all pixels		
		if numpy.nansum(prf)==0.: prf = copy.deepcopy(profile)
# extract fluxes & uncertainties - not entirely sure about the variance calculation...
		flx.append(numpy.nansum(src*prf)/numpy.nansum(prf))
#
# NOTE: THERE IS A PERSISTENT ERROR HERE
#
		unc.append((numpy.nansum(vsrc*prf)**0.5)/numpy.nansum(prf))
	flx = numpy.array(flx)
	unc = numpy.array(unc)

	sp = Spectrum(flux=flx,unc=unc,background=bck,**kwargs)

	if plot_file != '':
		x = numpy.arange(len(flx))
		spsm = copy.deepcopy(sp)
		spsm.smooth(15)
		sn = flx/unc

		plt.clf()
		plt.figure(figsize=[8,15])
		plt.subplot(511)
		slc = numpy.nanmedian(im,axis=1)
		slc = slc-numpy.nanmedian(slc)
		slc = slc/numpy.nanmax(slc)
		# xslc = numpy.arange(len(slc))+cntr-2*bck_wnd[1]
		# xref = int(numpy.nanmedian(x))
		# slc = numpy.nanmedian(im[int(trace[xref]-2*bck_wnd[1]):int(trace[xref]+2*bck_wnd[1]),int(xref-10):int(xref+10)],axis=1)

		xrng = [cntr-2*bck_wnd[1],cntr+2*bck_wnd[1]]
		yrng = [-0.05,1.2]
		plt.plot(numpy.arange(len(slc)),slc,'k-',alpha=0.5)
		plt.plot(numpy.arange(len(profile))+cntr-0.5*(len(profile)-1),profile,'b-')
		plt.legend(['Spatial Profile','Extraction Profile'])
		slc[:int(cntr-0.5*(len(profile)-1))] = numpy.nan
		slc[int(cntr+0.5*(len(profile)-1))] = numpy.nan
		slc = slc/numpy.nanmax(slc)
		plt.plot(numpy.arange(len(slc)),slc,'m-',alpha=0.5)
		# slcpc = slc[int(cntr-0.5*(len(profile)-1)):int(cntr+0.5*(len(profile)))]
		# slcpcx = slc[int(cntr-0.5*(len(profile)-1)):int(cntr+0.5*(len(profile)))]
		# slcpc = slcpc/numpy.nanmax(slcpc)
		# plt.plot(numpy.arange(len(slcpc))+cntr-0.5*(len(slcpc)-1),slcpc,'m-')
#		plt.plot(xslc[int(2*bck_wnd[1]-src_wnd):int(2*bck_wnd[1]+src_wnd)],slc[int(2*bck_wnd[1]-src_wnd):int(2*bck_wnd[1]+src_wnd)],'k-')
#		plt.plot(xslc[bck_wnd[1]:int(2*bck_wnd[1]-bck_wnd[0])],slc[bck_wnd[1]:int(2*bck_wnd[1]-bck_wnd[0])],'m-')
#		plt.plot(xslc[int(2*bck_wnd[1]+bck_wnd[0]):int(3*bck_wnd[1])],slc[int(2*bck_wnd[1]+bck_wnd[0]):int(3*bck_wnd[1])],'m-')
		for b in bck_wnd: 
			plt.plot([cntr-b,cntr-b],yrng,'m--')
			plt.plot([cntr+b,cntr+b],yrng,'m--')
		plt.plot([cntr-src_wnd,cntr-src_wnd],yrng,'k--')
		plt.plot([cntr+src_wnd,cntr+src_wnd],yrng,'k--')
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Pixel',fontsize=16)
		plt.ylabel('Count Rate per Pixel',fontsize=16)
		plt.xlim(xrng)
		plt.ylim(yrng)

		plt.subplot(512)
		plt.plot(x,flx,'k-')
		plt.plot(x,unc,'k--',alpha=0.5)
		plt.legend(['Flux','Uncertainty'])
		plt.plot(x,x*0.,'k--')
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Pixel',fontsize=16)
		plt.ylabel('Count Rate per Pixel',fontsize=16)
		plt.ylim([-2.*numpy.nanmedian(unc),1.5*numpy.nanquantile(spsm.flux.value,0.9)])

		plt.subplot(513)
		plt.plot(x,flx/unc,'k-')
		plt.plot(x,x*0.,'k--')
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Pixel',fontsize=16)
		plt.ylabel('Signal/Noise',fontsize=16)
		plt.ylim([0,1.5*numpy.nanquantile(sn,0.9)])

		plt.subplot(514)
		plt.plot(x,bck,'k-')
		plt.legend(['Background'])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Pixel',fontsize=16)
		plt.ylabel('Count Rate per Pixel',fontsize=16)
		plt.ylim([0.,1.5*numpy.nanquantile(spsm.background.value,0.9)])

		plt.subplot(515)
		plt.imshow(imsub,vmin=numpy.nanmedian(imsub)-3.*numpy.nanstd(imsub),vmax=numpy.nanmedian(imsub)+3.*numpy.nanstd(imsub))
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Pixel',fontsize=16)
		plt.ylabel('Pixel',fontsize=16)

		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			print('Warning: could not save flux calibration diagnostics to file {}'.format(plot_file))
#		plt.clf()
# generate a spectrum object
	return sp


def waveCalibrateArcs(arcim,deep=[],dispersion='',mode='',trace=[],prior={},fit_order=-1,sfit_order=2,verbose=ERROR_CHECKING,middle=True,resolution=0.,lam0=[],pixel0=[],cntr=0.,fitcycle=7,sclip=2.,plot_file=''):
	'''
	Wavelength calibration from arc lamp
	Input: arc, list of lines
	Output: dictionary containing pixel->wave conversion, fit diagnostics
	NEED TO INCLUDE HELIOCENTRIC CORRECTION EXPLICITLY
	'''
# check inputs
	if dispersion in list(DISPERSIONS.keys()):
		if resolution == 0.: resolution = DISPERSIONS[dispersion]['RESOLUTION']
		if len(lam0)==0.: lam0 = [DISPERSIONS[dispersion]['LAM0']]
	if len(prior)==0 and (resolution==0. or lam0 ==0.):
		raise ValueError('You must provide prior fit coefficients (prior={}), specify the dispersion (dispersion={}), or specify the resolution (resolution={}) and strongest line wavelength (lam0={})'.format(prior,dispersion,resolution,lam0))

# strongest line in red - THIS WILL NEED TO BE CHANGED TO ALLOW BLUE EXTRACTIONS
	if mode=='RED':
		extwidth = 30
		swindow = 30
		dwindow = 15
		deep_threshold = 0.003
		hehgcd = numpy.array([4358.33,4471.50,4678.16,4799.92,5015.68,5085.82,5460.74,5769.59,5790.65,5875.62,\
				6678.15,7065.19,8667.943,9122.97,9224.50,9354.22,9657.78,9784.50,10139.5])
		near = numpy.array([5852.49,5881.19,5944.83,5975.28,6030.00,6074.34,6096.16,6143.06,6163.59,6217.28,\
				6266.50,6304.79,6334.40,6382.99,6402.25,6506.53,6532.88,6598.95,6678.20,6717.04,6929.47,6965.43,\
				7032.41,7173.94,7245.17,7383.98,7438.90,7635.105,7723.76,7948.175,8115.31,8264.52,8300.33,8377.61,8418.43,8424.65,8495.36,\
	#		  8581.26,8654.45,8667.943,
				8780.6223,8865.7562,8919.5007,8988.58,9148.68,9300.85,9813.98,9373.28,9425.38,9459.21,9486.68,\
				9534.17,9665.43])
		alines = numpy.append(near,hehgcd)
		strong = numpy.array([5944.83,6143.06,6402.25,6506.52,6678.2,6929.47,7032.41,7245.17,7438.90,7635.105,8115.31,8377.61,8919.5007,8988.58,9122.9660])
		if fit_order < 0: fit_order = 6
	elif mode=='BLUE':
		extwidth = 30
		swindow = 30
		dwindow = 25
		deep_threshold = 0.0015
		hehgcd = numpy.array([3261.05,3341.08,3403.65,3466.55,3610.51,3650.15,3662.88,3888.65,4046.56,4077.83,\
				4358.33,4471.50,4678.16,4799.92,4921.93,5015.68,5085.82,5460.74,5769.59,5790.65])
		alines = hehgcd
		strong = numpy.array([3466.55,4046.56,4358.33,4678.16,4799.92,5085.82,5460.74,5875.62,5944.83])
		if fit_order < 0: fit_order = 4
	elif mode=='LDSS3':
		extwidth = 30
		swindow = 30
		dwindow = 30
		cntr = 250
		sfit_order = 3
		deep_threshold = 0.003
		alines = numpy.array([6096.1630,6163.5939,6217.2813,6266.495,6304.7892,6334.4279,6382.9914,6402.246,\
			6506.5279,6532.8824,6598.9529,6678.20,6717.0428,6929.468,7032.4127,7173.939,7245.167,7281.349,7383.98,\
			7503.867,7635.105,7723.8,8014.786,8115.311,8264.521,8300.325,8424.647,8521.441,8654.383,8782.1872,\
			9122.966,9224.498,9657.784])
		strong = numpy.array([6402.246,6678.20,7383.98,7503.867,7723.8,8014.786,8424.647,8521.441,9122.966,9657.784])
		if fit_order < 0: fit_order = 4
	else: raise ValueError('You must specify the red/blue mode (mode=RED or BLUE); mode={} was passed'.format(mode))

	cal_wave = {}

# extract arc traces from arc files
	if len(trace)==0:
		if cntr==0.: cntr = int(0.5*len(arcim[:,0]))
		trace = numpy.zeros(len(arcim[0,:]))+cntr
	tr = [int(t) for t in trace]
	ext = extractSpectrum(arcim,method='arc',subtract_background=False,cntr=cntr,trace=tr,src_wnd=extwidth,shift_trace=False,verbose=verbose)
	arctrace = ext.flux.value
	arctrace = arctrace-numpy.nanmin(arctrace)
	arctrace = arctrace/numpy.nanmax(arctrace)
	if len(deep)>0:
		ext = extractSpectrum(deep,method='arc',subtract_background=False,cntr=cntr,trace=tr,src_wnd=extwidth,shift_trace=False,verbose=verbose)
		arctrace_deep = ext.flux.value
		arctrace_deep = arctrace_deep-numpy.nanmin(arctrace_deep)
		arctrace_deep = arctrace_deep/numpy.nanmax(arctrace_deep)

	# 	arctrace = numpy.nanmedian(arcsh[(middle-extwidth):(middle+extwidth),:],axis=0)
	# 	arctrace_deep = numpy.nanmedian(arcdp[(middle-extwidth):(middle+extwidth),:],axis=0)
	# else:
	# 	arctrace = numpy.nanmedian(arcsh,axis=0)
	# 	arctrace_deep = numpy.nanmedian(arcdp,axis=0)

# if no prior provided, make initial guess from resolution, center to strongest line, and do 1st order fit
	if len(pixel0) == 0: pixel0=[numpy.argmax(arctrace)]
	pixel0 = [int(x) for x in pixel0]
	if 'COEFF' in list(prior.keys()):
		p1 = prior['COEFF']
		wave1 = numpy.polyval(p1,numpy.arange(len(arctrace)))
	else: 
		if len(lam0) == 0: raise ValueError('Problem: there is no initial estimate of central wavelength to perform wavelenth calibration')
		elif len(lam0) == 1: p0 = [resolution,lam0[0]-resolution*pixel0[0]]
		else: p0 = numpy.polyfit(pixel0,lam0,int(numpy.nanmin([sfit_order,len(lam0)-1])))
		wave0 = numpy.polyval(p0,numpy.arange(len(arctrace)))

# now center the strongest feature
# this could be done better using a cross-correlation method
		lines_y = strong[strong<(numpy.nanmax(wave0)-0.1*(numpy.nanmax(wave0)-numpy.nanmin(wave0)))]
		lines_y = lines_y[lines_y>(numpy.nanmin(wave0)+0.1*(numpy.nanmax(wave0)-numpy.nanmin(wave0)))]
		lines_x = []
		for n in lines_y:
			x0 = numpy.argmin(numpy.absolute(wave0-n))
			aslice = arctrace[(x0-swindow):(x0+swindow+1)]
			lines_x.append(x0+numpy.argmax(aslice)-swindow+1)
#		print(dispersion,resolution,lam0,numpy.argmax(arctrace),numpy.nanmin(wave0),numpy.nanmax(wave0),lines_y,lines_x)
#		sfit_order = numpy.nanmin([len(lines_x)-1,fit_order])
		p1 = numpy.polyfit(lines_x,lines_y,sfit_order)
		wave1 = numpy.polyval(p1,numpy.arange(len(arctrace)))
		rms = numpy.nanstd(lines_y-numpy.polyval(p1,lines_x))
		if verbose==True: print('Strong line pass: RMS={:.3f} Ang at {:.2f} Ang for {} lines and fit order {}; dRV = {:.1f} km/s'.format(rms,numpy.nanmedian(wave1),len(lines_y),sfit_order,3.e5*rms/numpy.nanmedian(wave1)))

	# plt.plot(wave1,arctrace)
	# raise ValueError

# now repeat for all of the lines 
	all_select = alines[alines<(numpy.nanmax(wave1)-0.02*(numpy.nanmax(wave1)-numpy.nanmin(wave1)))]
	all_select = all_select[all_select>(numpy.nanmin(wave1)+0.02*(numpy.nanmax(wave1)-numpy.nanmin(wave1)))]
#	all_select = alines[alines<(numpy.nanmax(wave1))]
#	all_select = all_select[all_select>(numpy.nanmin(wave1))]
# use deep exposure if you have it
	atrace = arctrace
	if len(deep)>0: atrace = arctrace_deep

	lines_y,lines_x = [],[]
	for n in all_select:
		x0 = numpy.argmin(numpy.absolute(wave1-n))
		aslice = atrace[(x0-dwindow):(x0+dwindow+1)]
		if numpy.nanmax(aslice)>deep_threshold:
			lines_y.append(n)
			lines_x.append(x0+numpy.argmax(aslice)-dwindow+1)
	lines_x = numpy.array(lines_x)
	lines_y = numpy.array(lines_y)
	w = numpy.where(lines_x>0)
# cycle through the fit	
	for i in range(fitcycle):
		sfit_order = numpy.nanmin([len(lines_x[w])-1,fit_order])
		p2 = numpy.polyfit(lines_x[w],lines_y[w],sfit_order)
		diff = lines_y-numpy.polyval(p2,lines_x)
		w = numpy.where(numpy.absolute(diff)<sclip*numpy.nanstd(diff[w]))
		rms = numpy.nanstd(diff[w])
		wave = numpy.polyval(p2,numpy.arange(len(atrace)))
		if verbose==True: print('Deep line pass {}: RMS={:.3f} Ang for {} lines and fit order {}, dRV = {:.1f} km/s'.format(i,rms,len(lines_y[w]),sfit_order,3.e5*rms/numpy.nanmedian(wave)))
		lines_yp,lines_xp = [],[]

# reselect lines with improved wavelength calibration
		all_select = alines[alines<(numpy.nanmax(wave)-0.02*(numpy.nanmax(wave)-numpy.nanmin(wave)))]
		all_select = all_select[all_select>(numpy.nanmin(wave)+0.02*(numpy.nanmax(wave)-numpy.nanmin(wave)))]
		for n in all_select:
			x0 = numpy.argmin(numpy.absolute(wave-n))
			aslice = atrace[(x0-dwindow):(x0+dwindow+1)]
#			print(n,numpy.nanmax(aslice),deep_threshold)
			if numpy.nanmax(aslice)>deep_threshold:
				lines_yp.append(n)
				lines_xp.append(x0+numpy.argmax(aslice)-dwindow+1)
		if len(lines_xp) > len(lines_x):
			lines_x = numpy.array(lines_xp)
			lines_y = numpy.array(lines_yp)
			w = numpy.where(lines_x>0)

	lines_x = lines_x[w]
	lines_y = lines_y[w]
	diff = lines_y-numpy.polyval(p2,lines_x)
	wave = numpy.polyval(p2,numpy.arange(len(atrace)))
	rms = numpy.nanstd(lines_y-numpy.polyval(p2,lines_x))
	if verbose==True: 
		print('Final deep line pass: RMS={:.3f} Ang at {:.2f} Ang for {} lines and fit order {}, dRV = {:.1f} km/s'.format(rms,numpy.nanmedian(wave),len(lines_y),sfit_order,3.e5*rms/numpy.nanmedian(wave)))
		print('Coefficients:')
		for coeff in p2: print('\t{}'.format(coeff))

# generate reporting structure	
	cal_wave['COEFF'] = p2
	cal_wave['WAVE'] = wave
	cal_wave['UNIT'] = u.Angstrom
	cal_wave['RMS'] = rms
	cal_wave['FITLINES'] = lines_y
	cal_wave['FITPIXELS'] = lines_x
	cal_wave['ARCTRACE'] = atrace

# plot results
	if plot_file!='':
		plt.clf()
		plt.figure(figsize=[12,8])
		plt.subplot(211)
		plt.plot(lines_x,diff,'bo')
		plt.legend(['RMS={:.2f} Ang\nRMS={:.0f} km/s'.format(rms,3.e5*rms/numpy.nanmedian(lines_y))])
		plt.plot([0,len(wave)],[0,0],'k--')
		plt.plot([0,len(wave)],[rms,rms],'k--')
		plt.plot([0,len(wave)],[-1.*rms,-1.*rms],'k--')
		plt.xlim([0,len(wave)])
		plt.ylim([-3.*rms,3.*rms])
		for i,x in enumerate(lines_x): 
			plt.text(x,diff[i],' {:.2f} \n {:.0f} '.format(lines_y[i],x),rotation=45.,horizontalalignment='left',fontsize=5)
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Pixel',fontsize=16)
		plt.ylabel('Line O-C (Ang)',fontsize=16)

		plt.subplot(212)
		yrng=[0.01,1.5]
		plt.semilogy(wave,arctrace+0.01)
		if len(deep)>0: 
			plt.semilogy(wave,arctrace_deep+0.01)
			plt.legend(['Shallow arc','Deep arc'])
		else: plt.legend(['Arc'])
		for x in lines_y: 
			plt.plot([x,x],yrng,'k--',alpha=0.5)
			plt.text(x,yrng[1],'{:.3f} '.format(x),rotation=90.,horizontalalignment='right',fontsize=5)
		for x in all_select:
			if x not in lines_y: 
				plt.plot([x,x],yrng,'k--',alpha=0.1)
				plt.text(x,yrng[1],'{:.3f} '.format(x),rotation=90.,horizontalalignment='right',fontsize=4,color='grey')
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Wavelength (Ang)',fontsize=16)
		plt.ylabel('Normalized flux density',fontsize=16)
		plt.ylim(yrng)
		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			print('Warning: could not save arc wavelength calibration diagnostics to file {}'.format(plot_file))
#		plt.clf()
	return cal_wave



def fluxCalibrate(fluxsp,name,fit_order=5,fit_cycle=10,sclip=3.,fit_range=[],fluxcalfolder=FLUXCALFOLDER,calwaveunit=DEFAULT_WAVE_UNIT,calfluxunit=DEFAULT_FLUX_UNIT,plot_file='',verbose=ERROR_CHECKING):
	'''
	Uses spectrum of flux calibrator to determine the flux calibration correction
	Input: flux cal spectrum, flux cal name
	Output: dictionary containing  observed to physical flux correction, fit diagnostics
	'''
	calscalefactor = 1.e-16
	if name.upper() not in list(FLUXCALS.keys()):
		raise ValueError('Flux standard {} is not present in flux calibration library; those sources are {}'.format(name,list(FLUXCALS.keys())))
	calfile = '{}/{}'.format(fluxcalfolder,FLUXCALS[name.upper()]['FILE'])
	if os.path.exists(calfile)==False:
		raise ValueError('Cannot find flux standard file {} please check your paths'.format(calfile))
# read in flux cal and compute ratio
	fcal = readSpectrum(calfile,flux_unit=calfluxunit,wave_unit=calwaveunit)
	ratio = fcal/fluxsp
	if len(fit_range) == 0: fit_range = [numpy.nanmin(fcal.wave.value),numpy.nanmax(fcal.wave.value)]

# iteratively fit LOG of ratio with rejection
	mask = numpy.zeros(len(ratio.wave))
	mask[ratio.flux<=0] = 1
	mask[ratio.wave.value<fit_range[0]] = 1
	mask[ratio.wave.value>fit_range[1]] = 1
	mask[numpy.isfinite(ratio.flux)==False] = 1
	for i in range(fit_cycle):
		p = numpy.polyfit(ratio.wave.value[mask==0],numpy.log10(ratio.flux[mask==0]),fit_order)
		diff = numpy.log10(ratio.flux)-numpy.polyval(p,ratio.wave.value)
		mask[numpy.absolute(diff)>(sclip*numpy.nanstd(diff[mask==0]))] = 1
	p[-1]=p[-1]+numpy.log10(calscalefactor)

# generate reporting structure	
	cal_flux = {}
	cal_flux['COEFF'] = p
	cal_flux['CORRECTION'] = 10.**numpy.polyval(p,fluxsp.wave.value)
	cal_flux['UNIT'] = calfluxunit
	cal_flux['FITBANDS'] = len(mask[mask==0])
	cal_flux['NAME'] = name

# plot results
	if plot_file!='':
		plt.clf()
		plt.figure(figsize=[8,8])
		plt.subplot(211)
		plt.semilogy(ratio.wave.value,ratio.flux.value,'k.',alpha=0.5)
		plt.semilogy(ratio.wave.value[mask==0],ratio.flux.value[mask==0],'bo')
		plt.semilogy(fluxsp.wave.value,(10.**numpy.polyval(p,fluxsp.wave.value))/calscalefactor,'m--')
		plt.legend(['Ratio','RMS={:.2f} dex'.format(numpy.nanstd(diff)),'Fit'])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlim([numpy.nanmin(fluxsp.wave.value),numpy.nanmax(fluxsp.wave.value)])
		plt.xlabel('Wavelength (Ang)',fontsize=16)
		plt.ylabel('Calibrator/Observed',fontsize=16)
		plt.subplot(212)
#		yrng=[0,1.2]
		plt.plot(fcal.wave.value,fcal.flux.value*calscalefactor,'k-')
		f = interp1d(fcal.wave.value,fcal.flux.value*calscalefactor,bounds_error=False,fill_value=0.)
		plt.plot(fluxsp.wave.value,fluxsp.flux.value*numpy.nanmedian(f(fluxsp.wave.value))/numpy.nanmedian(fluxsp.flux.value),'b-')
		plt.plot(fluxsp.wave.value,fluxsp.flux.value*(10.**numpy.polyval(p,fluxsp.wave.value)),'m-')
		plt.legend(['Calibration Spectrum of {}'.format(name),'Observed Spectrum (scaled)','Corrected Spectrum'])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlim([numpy.nanmin(fluxsp.wave.value),numpy.nanmax(fluxsp.wave.value)])
		plt.ylim([0,numpy.nanmax(f(fluxsp.wave.value))])
		plt.xlabel('Wavelength (Ang)',fontsize=16)
		plt.ylabel('Apparent Flux density ({})'.format(calfluxunit),fontsize=16)
		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			print('Warning: could not save flux calibration diagnostics to file {}'.format(plot_file))
#		plt.clf()
	return cal_flux


def fluxReCalibrate(tellsp,spt,fit_order=4,fit_cycle=10,sclip=3.,smooth=5,snlimit=30,fit_range=[],calfolder=TELLSTDFOLDER,stdfile='',plot_file='',verbose=ERROR_CHECKING):
	'''
	Uses spectrum of telluric calibrator to determine a second-order flux calibration correction
	Input: telluric cal spectrum, spt
	Output: dictionary containing  observed to physical flux correction, fit diagnostics

	TODO
	need to add wavelength shift around H lines, and smooth data and model to same resolution
	'''

	tspt = copy.deepcopy(spt).upper()
	if tspt[-1] != 'V': tspt=tspt+'V'

	cal_reflux = {}
	cal_reflux['SPT'] = tspt
	cal_reflux['UNIT'] = u.dimensionless_unscaled

	wv = tellsp.wave.value	
	flx = tellsp.flux.value
	mask = numpy.zeros(len(wv))

# establish fit range
# quick exit if spectrum is outside fit range
	if len(fit_range) == 0: fit_range = [numpy.nanmin(tellsp.wave.value),numpy.nanmax(tellsp.wave.value)]
	if numpy.nanmax(wv) < numpy.nanmin(fit_range) or numpy.nanmin(wv) > numpy.nanmax(fit_range): 
		if verbose==True: print('Fit range {} is outside spectral range {}'.format(fit_range,[numpy.nanmin(wv),numpy.nanmax(wv)]))
		cal_reflux['COEFF'] = [1.]
		cal_reflux['MASK'] = mask
		cal_reflux['WAVE'] = wv
		cal_reflux['CORRECTION'] = numpy.ones(len(wv))
		return cal_reflux

# spectral standards - currently only A0V and G2V
	stdfile = calfolder+'{}_pickles1998.txt'.format(tspt)
	# stds = {
	# 	'A0': calfolder+'A00_SDSS-kesseli2017.txt',
	# 	'G2': calfolder+'G20_SDSS-kesseli2017.txt',
	# }

	# stdfile = ''
	# for k in list(stds.keys()):
		# if (spt.upper())[:2] in k: stdfile = stds[k]
	# if stdfile=='':
	if os.path.exists(stdfile)==False:
		# if verbose==True: print('Telluric spectral type {} is not among of current telluric standard spectra: {}'.format(spt,list(stds.keys())))
		if verbose==True: print('Telluric spectral type {} or calibration file {} is not among of current telluric standards in {}; no change made'.format(spt,'{}_pickles1998.txt'.format(tspt),calfolder))
		cal_reflux['COEFF'] = [1.]
		cal_reflux['MASK'] = mask
		cal_reflux['WAVE'] = wv
		cal_reflux['CORRECTION'] = numpy.ones(len(wv))
		return cal_reflux
	fcal = readSpectrum(stdfile,flux_unit=tellsp.flux.unit,wave_unit=tellsp.wave.unit)

# compute, normalize, and smooth ratio of observed telluric to calibrator telluric
	ratio = tellsp/fcal
	calscalefactor = numpy.nanmedian(ratio.flux.value)
	ratio.scale(1./calscalefactor)
	ratio.smooth(smooth)

# force ends
#	ratio[:50] = numpy.nanmedian(ratio[50:150])
#	ratio[-50:] = numpy.nanmedian(ratio[-150:-50])

# regions of strong telluric absorption	
	tranges = [[5800,6000],[6230,6340],[6850,7050],[7140,7350],[7580,7700],[8080,8360],[8920,9800]]
# regions of G/A star line absorption	
	glines = [[6548,6578]]
	if 'A' in spt.upper():
		tranges = [[5800,6000],[6230,6340],[6850,7050],[7140,7350],[7580,7700],[8080,8360],[9050,9800]]
		glines = [[6535,6593],[8194,8214],[8480,8550],[8571,8621],[8638,8688],[8718,8778],[8830,8890],[8992,9032],[9196,9256],[9526,9576]]

# define mask
	mask[numpy.isnan(flx)==True] = 1
	mask[ratio.flux.value<=0.4] = 1
	mask[ratio.flux.value>=1.6] = 1
	mask[numpy.isfinite(ratio.flux.value)==False] = 1
	mask[numpy.isfinite(tellsp.flux.value/tellsp.unc.value)==False] = 1
	mask[tellsp.flux.value/tellsp.unc.value < snlimit] = 1
	for t in tranges: mask[numpy.logical_and(wv>=t[0],wv<=t[1])]=1
#	for t in glines: mask[numpy.logical_and(wv>=t[0],wv<=t[1])]=1
	mask[numpy.logical_or(wv<fit_range[0],wv>fit_range[1])]=1



# determine a flux correction in unmasked regions of ratio
	for i in range(fit_cycle):
		p = numpy.polyfit(ratio.wave.value[mask==0],1./ratio.flux.value[mask==0],fit_order)
		diff = ratio.flux.value-(1./numpy.polyval(p,ratio.wave.value))
		mask[numpy.absolute(diff)>(sclip*numpy.nanstd(diff[mask==0]))] = 1

# generate interpolated correction function defaulting to 1
	f = interp1d(ratio.wave.value[mask==0],ratio.flux.value[mask==0],bounds_error=False,fill_value=1.)
	correction = 1./f(wv)

# generate reporting structure	
	cal_reflux['COEFF'] = p
	cal_reflux['MASK'] = mask
	cal_reflux['WAVE'] = wv
	cal_reflux['CORRECTION'] = correction

# plot results
	if plot_file!='':
		xlim = [numpy.nanmin(tellsp.wave.value),numpy.nanmax(tellsp.wave.value)]
		plt.clf()
		plt.figure(figsize=[8,8])
		plt.subplot(211)
		plt.plot(ratio.wave.value,ratio.flux.value,'k-',alpha=0.5)
		plt.plot(ratio.wave.value[mask==0],ratio.flux.value[mask==0],'k-')
		plt.plot(wv,1./correction,'m-')
		plt.plot(wv,1./numpy.polyval(p,wv),'m--')
		plt.legend(['Ratio','RMS={:.2f} dex'.format(numpy.nanstd(diff)),'Correction','Fit'])
		plt.plot(xlim,[1,1],'k--')
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlim(xlim)
		plt.xlabel('Wavelength (Ang)',fontsize=16)
		plt.ylabel('Standard/Observed',fontsize=16)
		plt.subplot(212)
#		yrng=[0,1.2]
		plt.plot(fcal.wave.value,fcal.flux.value*calscalefactor,'k-')
		plt.plot(tellsp.wave.value,tellsp.flux.value,'b-')
		plt.plot(tellsp.wave.value,tellsp.flux.value*numpy.polyval(p,tellsp.wave.value),'m-')
		plt.legend(['{} Calibration Spectrum'.format(spt),'Observed Telluric Spectrum (scaled)','Corrected Telluric Spectrum'])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlim(xlim)
		plt.ylim([0,numpy.nanquantile(tellsp.flux.value,0.9)*1.2])
		plt.xlabel('Wavelength (Ang)',fontsize=16)
		plt.ylabel('Apparent Flux density ({})'.format(tellsp.flux.unit),fontsize=16)
		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			print('Warning: could not save flux calibration diagnostics to file {}'.format(plot_file))

	return cal_reflux


def telluricCalibrate(tellsp,fit_range=[6000,9200],fit_order=4,fitcycle=10,sclip=3.,spt='G2',plot_file='',verbose=ERROR_CHECKING):
	'''
	Uses spectrum of telluric calibrator to determine corrections to telluric absorption
	Note: this is ONLY for the red
	Input: telluric cal spectrum, spectral type
	Output: dictionary containing telluric absorption correction, fit diagnostics
	'''

# regions of strong telluric & absorption	
	tranges = [[6230,6340],[6850,7050],[7140,7350],[7580,7700],[8080,8360],[8920,9200],[9250,9800]]
# regions of G/A star line absorption	
	glines = [[6548,6578]]
	if 'A' in spt.upper():
		tranges = [[6230,6340],[6850,7050],[7140,7350],[7580,7700],[8080,8360],[9050,9200],[9250,9800]]
		glines = [[6535,6593],[8194,8214],[8480,8550],[8571,8621],[8638,8688],[8718,8778],[8830,8890],[8992,9032],[9196,9256],[9526,9576]]
# mask out telluric and star lines to get fit of continuum
	wv = tellsp.wave.value	
	flx = tellsp.flux.value

# quick exit if spectrum is outside telluric fit range (e.g., for Blue)
	if numpy.nanmax(wv) < numpy.nanmin(fit_range): 
		cal_tell = {}
		cal_tell['WAVE'] = wv
		cal_tell['CORRECTION'] = numpy.ones(len(wv))
		return cal_tell

# mask and fit continuum
	mask = numpy.zeros(len(wv))
	mask[numpy.isnan(flx)==True] = 1
	mask[flx<=0] = 1
#	print(len(mask),numpy.nansum(mask))
	for t in tranges: mask[numpy.logical_and(wv>=t[0],wv<=t[1])]=1
#	print(len(mask),numpy.nansum(mask))
	for t in glines: mask[numpy.logical_and(wv>=t[0],wv<=t[1])]=1
#	print(len(mask),numpy.nansum(mask))
	mask[numpy.logical_or(wv<fit_range[0],wv>fit_range[1])]=1
	fmask = copy.deepcopy(mask)
	#plt.plot(wv[msk==0],flx[msk==0])
	for i in range(fitcycle):
		if numpy.nansum(fmask) >= len(fmask): fmask = numpy.zeros(len(fmask))
		p = numpy.polyfit(wv[fmask==0],flx[fmask==0],fit_order)
		continuum = numpy.polyval(p,wv)
		diff = continuum-flx
		fmask[numpy.absolute(diff)>(sclip*numpy.nanstd(diff[fmask==0]))] = 1
	#plt.plot(wv[msk==0],continuum(msk==0))

# now compute correction as ratio of continuum/observed in telluric absorption regions
	correction = numpy.ones(len(wv))
	tmask = numpy.zeros(len(wv))
	for t in tranges: tmask[numpy.logical_and(wv>=t[0],wv<=t[1])]=1
	for t in glines: tmask[numpy.logical_and(wv>=t[0],wv<=t[1])]=0
	if numpy.nansum(tmask)>0:
		correction[numpy.logical_and(tmask==1,flx>0)] = continuum[numpy.logical_and(tmask==1,flx>0)]/flx[numpy.logical_and(tmask==1,flx>0)]

# outside fit_range no correction
	correction[wv>numpy.nanmax(fit_range)]=1
	correction[wv<numpy.nanmin(fit_range)]=1

# generate reporting structure	
	cal_tell = {}
	cal_tell['WAVE'] = wv
	cal_tell['CORRECTION'] = correction

# plot results
#	print(len(wv),len(flx),len(mask),numpy.nansum(mask))
	if plot_file!='':
		plt.clf()
		plt.figure(figsize=[8,8])
		plt.subplot(211)
		plt.plot(wv,flx,'k-',alpha=0.5)
		plt.plot(wv[tmask==0],flx[tmask==0],'b-',)
		plt.plot(wv,continuum,'m-')
		plt.legend(['Observed','Stellar ({})'.format(spt),'Continuum'])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Wavelength (Ang)',fontsize=16)
		plt.ylabel('Observed Flux Density',fontsize=16)
		plt.subplot(212)
#		yrng=[0,1.2]
		plt.plot(wv,correction,'k-')
		plt.legend(['Correction'])
		plt.xticks(fontsize=16)
		plt.yticks(fontsize=16)
		plt.xlabel('Wavelength (Ang)',fontsize=16)
		plt.ylabel('Correction factor',fontsize=16)
		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			print('Warning: could not save flux calibration diagnostics to file {}'.format(plot_file))
#		plt.clf()

	return cal_tell


def rectify(im,trace=[],cntr=-1,trim=[],verbose=ERROR_CHECKING,center_kwargs={},trace_kwargs={},save_image='',*args):
	'''
	Purpose:
		Rectify image using trace, so that dispersion axis lies along a single row
	'''

# if necessary, generate a new trace
	if len(trace)==0:
		if cntr<0: cntr=findPeak(im,**kwargs)   
		cntr = int(cntr)
		trace = traceDispersion(im,cntr=cntr,**kwargs)

# starting a middle of trace, slide columns so that rows align
	shift_trace = trace[int(len(trace)/2.)]-numpy.array(trace)
	shift_trace = [int(t) for t in shift_trace]
	imshift = copy.deepcopy(im)
	for i,sh in enumerate(shift_trace): imshift[:,i] = numpy.roll(imshift[:,i],sh)

# rectify other images - WORK IN PROGRESS
	# if len(args)>0:
	# 	argout = []
	# 	for a in args:
	# 		aout = copy.deepcopy(a)
	# 		for i,sh in enumerate(shift_trace): aout[:,i] = numpy.roll(a[:,i],sh)

	if save_image!='':
		hdu = fits.PrimaryHDU(imshift)
		try: hdu.writeto(save_image,overwrite=overwrite)
		except: print('Warning: could not write rectified image to file {}'.format(save_image))

#  WORK IN PROGRESS
	# if len(args)>0: 
	# 	out = copy.deepcopy(imshift)
	# 	out = []
	# 	for a in args: out.append(a)
	# 	return *out
	return imshift

def applyCalibrations():
	'''
	Applies wavelength, flux and/or telluric calibrations to spectrum
	Input: spectrum object
	Output: spectrum object with corrections applied
	'''
	pass

def profileCheck(instructions='',cntr=335,verbose=ERROR_CHECKING,trace_slice=[250,1250],**kwargs):
	'''
	Examines profiles from instruction file to prep the extraction axes
	'''
# first read instructions
	if instructions!='': 
		if os.path.exists(instructions)==False: raise ValueError('Cannot find instruction file {}; be sure full path is correct'.format(instructions))
		parameters = readInstructions(instructions)
	if cntr==0: cntr=int(0.5*len(im[:,0]))
	if len(trace_slice)==0: trace_slice = [150,len(im[0,:])-150]

# now do abbreviated "peak check" and plot out results
	for k in ['SOURCE','TELLURIC','FLUXCAL']:
	# for src in list(parameters['SOURCE'].keys()):
	# 	if verbose==True: print('Checking spatial profile of {}'.format(src))
	# 	im,hd = readFiles(parameters['SOURCE'][src]['FILES'][0],folder=parameters['DATA_FOLDER'],mode=parameters['MODE'])
	# 	if 'CENTER' in list(parameters['SOURCE'][src].keys()): cntr = int(parameters['SOURCE'][src]['CENTER'])
	# 	cntr = findPeak(im,cntr=cntr,window=30,trace_slice=trace_slice,plot_file='{}/diagnostic_profile_{}_{}.pdf'.format(parameters['REDUCTION_FOLDER'],src,parameters['MODE']))
	# 	if verbose==True: print('Best center = {}'.format(cntr))

		if k in list(parameters.keys()):
			if len(list(parameters[k].keys()))>0:
				for src in list(parameters[k].keys()):
					if verbose==True: print('Checking spatial profile of {}'.format(src))
					im,hd = readFiles(parameters[k][src]['FILES'][0],folder=parameters['DATA_FOLDER'],mode=parameters['MODE'])
					if 'CENTER' in list(parameters[k][src].keys()): cntr = int(parameters[k][src]['CENTER'])
					cntr = findPeak(im,cntr=cntr,window=30,trace_slice=trace_slice,plot_file='{}/diagnostic_profile_{}_{}.pdf'.format(parameters['REDUCTION_FOLDER'],src,parameters['MODE']))
					if verbose==True: print('Best center = {}'.format(cntr))

	return


# full reduction pipeline
def reduce(redux={},parameters={},instructions='input.txt',bias_file='',flat_file='',mask_file='',cal_wave_file='',cal_flux_file='',reset=False,verbose=ERROR_CHECKING,src_wnd=5,bck_wnd=[15,20],fit_order=5,overwrite=True,**kwargs):
	'''
	Full reduction package
	'''
	src_wnd0 = copy.deepcopy(src_wnd)
	bck_wnd0 = copy.deepcopy(bck_wnd)
	fit_order0 = copy.deepcopy(fit_order)
	if instructions!='': 
		if os.path.exists(instructions)==False: raise ValueError('Cannot find instruction file {}; be sure full path is correct'.format(instructions))
		parameters = readInstructions(instructions)

# read in parameters and check for mandatory ones
	if 'PARAMETERS' not in list(redux.keys()) or reset==True:
		redux['PARAMETERS'] = parameters
	required_parameters = ['BIAS','FLAT','DATA_FOLDER','REDUCTION_FOLDER','MODE']
	for r in required_parameters:
		if r not in list(redux['PARAMETERS'].keys()): raise ValueError('Input instructions do not include required parameter {}'.format(r))
# set some parameters based on instructions
	if 'FLUXCAL_REUSE' in list(redux['PARAMETERS'].keys()): 
		if os.path.exists(redux['PARAMETERS']['FLUXCAL_REUSE'])==True: cal_flux_file = redux['PARAMETERS']['FLUXCAL_REUSE']
		elif verbose==True: print('Could not find flux calibration data {}'.format(redux['PARAMETERS']['FLUXCAL_REUSE']))
	if 'BIAS_REUSE' in list(redux['PARAMETERS'].keys()): bias_file = redux['PARAMETERS']['BIAS_REUSE']
		# bias_file = redux['PARAMETERS']['BIAS_REUSE']
		# if os.path.exists(bias_file) == False: bias_file=redux['PARAMETERS']['REDUCTION_FOLDER']+'/'+bias_file
		# if os.path.exists(bias_file) == False: 
		# 	if verbose==True: print('\nCannot find bias file {}'.format(bias_file))
		# 	bias_file=''
		# else:
		# 	if verbose==True: print('\nReading in bias frame from file {}'.format(bias_file))
		# 	redux['BIAS'],redux['BIAS_HEADER'] = readFiles(bias_file,mode=redux['PARAMETERS']['MODE'],rotate=False)
	if 'FLAT_REUSE' in list(redux['PARAMETERS'].keys()): flat_file = redux['PARAMETERS']['FLAT_REUSE']
	if 'MASK_REUSE' in list(redux['PARAMETERS'].keys()): mask_file = redux['PARAMETERS']['MASK_REUSE']

# make sure necessary folders are in place
	if os.path.exists(redux['PARAMETERS']['DATA_FOLDER']) == False: raise ValueError('Could not located data folder {}; please check the path'.format(redux['PARAMETERS']['DATA_FOLDER']))
	if os.path.exists(redux['PARAMETERS']['REDUCTION_FOLDER']) == False: 
		if verbose==True: print('Creating reduction folder {}'.format(redux['PARAMETERS']['REDUCTION_FOLDER']))
		try:
			os.mkdir(redux['PARAMETERS']['REDUCTION_FOLDER'])
		except:
			raise ValueError('Could not create reduction folder {}; please check the path'.format(redux['PARAMETERS']['REDUCTION_FOLDER']))

# IMAGE ANALYSIS
# read or create bias frame
	if 'BIAS' not in list(redux.keys()) and bias_file!='':
		tmp = copy.deepcopy(bias_file)
		if os.path.exists(bias_file) == False: bias_file=redux['PARAMETERS']['REDUCTION_FOLDER']+'/'+bias_file
		if os.path.exists(bias_file) == False: 
			if verbose==True: print('\nCannot find bias file {} or {}'.format(tmp,bias_file))
			bias_file=''
		else: 
			if verbose==True: print('\nReading in bias frame from file {}'.format(bias_file))
			redux['BIAS'],redux['BIAS_HEADER'] = readFiles(bias_file,mode=redux['PARAMETERS']['MODE'],rotate=False)
	if 'BIAS' not in list(redux.keys()) or (reset==True and 'BIAS_REUSE' not in list(redux['PARAMETERS'].keys())) or bias_file=='':
		bias_file='{}/bias_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'])
		if verbose==True: print('\nReducing bias frames')
		redux['BIAS'],redux['BIAS_HEADER'] = makeBias(redux['PARAMETERS']['BIAS']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'],output=bias_file)

	# if reset==True and 'BIAS_REUSE' not in list(redux['PARAMETERS'].keys())
	# if 'BIAS' not in list(redux.keys()) or reset==True:
	# 	print(reset,bias_file)
	# 	if bias_file!='' and reset==False:
	# 		if verbose==True: print('\nReading in bias frame from file {}'.format(bias_file))
	# 		redux['BIAS'],redux['BIAS_HEADER'] = readFiles(bias_file,mode=redux['PARAMETERS']['MODE'],rotate=False)
	# 	else:

# read or create flat field frame		
	if 'FLAT' not in list(redux.keys()) and flat_file!='':
		tmp = copy.deepcopy(flat_file)
		if os.path.exists(flat_file) == False: flat_file=redux['PARAMETERS']['REDUCTION_FOLDER']+flat_file
		if os.path.exists(flat_file) == False: 
			if verbose==True: print('\nCannot find flat file {} or {}'.format(tmp,flat_file))
			flat_file=''
		else: 
			if verbose==True: print('\nReading in flat field frame from file {}'.format(flat_file))
			redux['FLAT'],redux['FLAT_HEADER'] = readFiles(flat_file,mode=redux['PARAMETERS']['MODE'],rotate=False)
	if 'FLAT' not in list(redux.keys()) or (reset==True and 'FLAT_REUSE' not in list(redux['PARAMETERS'].keys())) or flat_file=='':
		flat_file='{}/flat_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'])
		if verbose==True: print('\nReducing flat field frames')
		redux['FLAT'],redux['FLAT_HEADER'] = makeFlat(redux['PARAMETERS']['FLAT']['FILES'],redux['BIAS'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'],output=flat_file)

# read or create mask frame		
	if 'MASK' not in list(redux.keys()) and mask_file!='':
		tmp = copy.deepcopy(mask_file)
		if os.path.exists(mask_file) == False: mask_file=redux['PARAMETERS']['REDUCTION_FOLDER']+mask_file
		if os.path.exists(mask_file) == False: 
			if verbose==True: print('\nCannot find mask file {} or {}'.format(tmp,mask_file))
			mask_file=''
		else: 
			if verbose==True: print('\nReading in mask frame from file {}'.format(mask_file))
			redux['MASK'],redux['MASK_HEADER'] = readFiles(mask_file,mode=redux['PARAMETERS']['MODE'],rotate=False)
	if 'MASK' not in list(redux.keys()) or (reset==True and 'MASK_REUSE' not in list(redux['PARAMETERS'].keys())) or mask_file=='':
		mask_file='{}/mask_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'])
		if verbose==True: print('\nGenerating mask file')
		redux['MASK'] = makeMask(redux['BIAS'],redux['FLAT'],mode=redux['PARAMETERS']['MODE'],output=mask_file)
		if verbose==True: print('Masking {} bad pixels ({:.0f}%)'.format(int(numpy.nansum(redux['MASK'])),100.*numpy.nansum(redux['MASK'])/redux['MASK'].size))

# clean the flat
	if verbose==True: print('\nCleaning flat frames')
	redux['FLAT'] = maskClean(redux['FLAT'],redux['MASK'],replace=1.)
	if reset==True:
		if flat_file!='':
			hdu = fits.PrimaryHDU(data=redux['FLAT'],header=redux['FLAT_HEADER'])
			hdu.writeto(flat_file,overwrite=overwrite)

# SPECTRAL CALIBRATIONS
# initial arc solution
	if cal_wave_file != '':
		if os.path.exists(cal_wave_file) == False: cal_wave_file=redux['PARAMETERS']['REDUCTION_FOLDER']+cal_wave_file
		if os.path.exists(cal_wave_file) == False: cal_wave_file=''
	if cal_wave_file != '':
		try: redux['CAL_WAVE'] = pickle.load(open(cal_wave_file,'rb'))
		except: print('WARNING: could not read in prior wavelength calibration structure from {}'.format(cal_wave_file))
	if 'CAL_WAVE' not in list(redux.keys()) or reset==True or kwargs.get('reset_wavecal',False) == True:
		if verbose==True: print('\nDetermining baseline wavelength solution')
		arcsh,arcshhd = readFiles(redux['PARAMETERS']['ARC_SHALLOW']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
		arcdp,arcdphd = readFiles(redux['PARAMETERS']['ARC_DEEP']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
# check mode
		arcmode = kastRBMode(arcshhd)
		if arcmode!=redux['PARAMETERS']['MODE']: raise ValueError('Warning: arc image mode {} is not the same as reduction mode {}'.format(arcmode,redux['PARAMETERS']['MODE']))
		grating = kastDispersion(arcshhd)
		if grating not in list(DISPERSIONS.keys()): raise ValueError('Do not have parameters for grating {}'.format(grating))
		redux['CAL_WAVE'] = waveCalibrateArcs(arcsh,deep=arcdp,dispersion=grating,middle=True,lam0=redux['PARAMETERS']['WAVE_INITIAL']['WAVE'],pixel0=redux['PARAMETERS']['WAVE_INITIAL']['PIXEL'],mode=redux['PARAMETERS']['MODE'],plot_file='{}/diagnostic_wavecal_initial_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']))
# save this to pickle file
		f = open('{}/cal_wave_{}.pkl'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']),'wb')
		pickle.dump(redux['CAL_WAVE'],f)
		f.close()
		
# determine flux calibration
	if cal_flux_file != '':
		if os.path.exists(cal_flux_file) == False: cal_flux_file=redux['PARAMETERS']['REDUCTION_FOLDER']+cal_flux_file
		if os.path.exists(cal_flux_file) == False: cal_flux_file=''
	if cal_flux_file != '':
		try: redux['CAL_FLUX'] = pickle.load(open(cal_flux_file,'rb'))
		except: print('WARNING: could not read in prior flux calibration structure from {}'.format(cal_flux_file))
	if ('CAL_FLUX' not in list(redux.keys()) or reset==True or kwargs.get('reset_fluxcal',False) == True) and 'FLUXCAL' in list(redux['PARAMETERS'].keys()):
		ref = list(redux['PARAMETERS']['FLUXCAL'].keys())[0]
		if verbose==True: print('\nAnalyzing flux calibrator {}'.format(ref))
		ims,hds = readFiles(redux['PARAMETERS']['FLUXCAL'][ref]['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
		if len(redux['PARAMETERS']['FLUXCAL'][ref]['FILES']) > 1:
			im,hd = crRejectCombine(ims,verbose=verbose),hds[0]
		else: im,hd = ims,hds
# some parameters
		cntr=-1
		if 'CENTER' in list(redux['PARAMETERS']['FLUXCAL'][ref].keys()): cntr = int(redux['PARAMETERS']['FLUXCAL'][ref]['CENTER'])
		src_wnd = copy.deepcopy(src_wnd0)
		if 'WINDOW' in list(redux['PARAMETERS']['FLUXCAL'][ref].keys()): src_wnd = redux['PARAMETERS']['FLUXCAL'][ref]['WINDOW']
		bck_wnd = copy.deepcopy(bck_wnd0)
		if 'BACK' in list(redux['PARAMETERS']['FLUXCAL'][ref].keys()): bck_wnd = redux['PARAMETERS']['FLUXCAL'][ref]['BACK']
		fit_order = copy.deepcopy(fit_order0)
		if 'FIT_ORDER' in list(redux['PARAMETERS']['FLUXCAL'][ref].keys()): fit_order = redux['PARAMETERS']['FLUXCAL'][ref]['FIT_ORDER']
# reduce data, determine trace and extract spectrum
		imr,var = reduceScienceImage(im,redux['BIAS'],redux['FLAT'],redux['MASK'],hd=hd)
		cntr = findPeak(imr,cntr=cntr)
		trace = traceDispersion(imr,cntr=cntr,window=src_wnd,method='maximum',plot_file='{}/diagnostic_trace_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],redux['PARAMETERS']['MODE']))
		imrect = rectify(imr,trace)
		varrect = rectify(var,trace)
		maskrect = rectify(redux['MASK'],trace)
		cntr = findPeak(imrect,cntr=cntr,window=src_wnd,plot_file='{}/diagnostic_profile_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],redux['PARAMETERS']['MODE']))
#		profile = spatialProfile(imrect,cntr=cntr,window=src_wnd)
#		profile = numpy.ones(int(2*src_wnd+1))
		spflx = extractSpectrum(imrect,cntr=cntr,var=varrect,mask=maskrect,src_wnd=src_wnd,bck_wnd=bck_wnd,method=redux['PARAMETERS']['FLUXCAL'][ref]['METHOD'],plot_file='{}/diagnostic_extraction_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],redux['PARAMETERS']['MODE']))
		spflx.name = redux['PARAMETERS']['FLUXCAL'][ref]['NAME']
		spflx.header = hd
#		spflx = extractSpectrum(imr,var,mask=redux['MASK'],trace=trace,src_wnd=10)
		arcdp,arcdphd = readFiles(redux['PARAMETERS']['ARC_DEEP']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
		arcrect = rectify(arcdp,trace)
		arcrecal = waveCalibrateArcs(arcrect,cntr=cntr,prior=redux['CAL_WAVE'],mode=redux['PARAMETERS']['MODE'],plot_file='{}/diagnostic_wavecal_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],redux['PARAMETERS']['MODE']))
		spflx.applyWaveCal(arcrecal)
		if redux['PARAMETERS']['MODE'] == 'BLUE':
			redux['CAL_FLUX'] = fluxCalibrate(spflx,redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],fit_order=fit_order,fit_range=[3500.,5550.],plot_file='{}/diagnostic_fluxcal_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']))
		else:
			redux['CAL_FLUX'] = fluxCalibrate(spflx,redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],fit_order=fit_order,plot_file='{}/diagnostic_fluxcal_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']))
		redux['CAL_FLUX']['TRACE'] = trace
#		redux['CAL_FLUX']['PROFILE'] = profile
		spflx.applyFluxCal(redux['CAL_FLUX'],log=True)
		redux[redux['PARAMETERS']['FLUXCAL'][ref]['NAME']] = spflx
# save this to pickle file
		f = open('{}/cal_flux_{}.pkl'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']),'wb')
		pickle.dump(redux['CAL_FLUX'],f)
		f.close()
# plot  and save
		fig = spflx.plot(ylim=[0,1.2*numpy.quantile(spflx.flux.value,0.95)])
		fig.figure.savefig('{}/kast{}_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],redux['PARAMETERS']['DATE']))
		plt.close()
#		plt.clf()
		spflx.toFile('{}/kast{}_{}_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],redux['PARAMETERS']['DATE']))
		spflx.toFile('{}/kast{}_{}_{}.txt'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],redux['PARAMETERS']['FLUXCAL'][ref]['NAME'],redux['PARAMETERS']['DATE']))

# determine telluric calibrations
	if ('CAL_TELL' not in list(redux.keys()) or reset==True or kwargs.get('reset_tellcal',False) == True) and 'TELLURIC' in list(redux['PARAMETERS'].keys()):
		if verbose==True: print('\nComputing telluric corrections')
		redux['CAL_TELL'] = {}
		for tstar in list(redux['PARAMETERS']['TELLURIC'].keys()):
			ims,hds = readFiles(redux['PARAMETERS']['TELLURIC'][tstar]['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
			if len(redux['PARAMETERS']['TELLURIC'][tstar]['FILES']) > 1: im,hd = crRejectCombine(ims,verbose=verbose),hds[0]
			else: im,hd = ims,hds
# some parameters
			cntr=-1
			if 'CENTER' in list(redux['PARAMETERS']['TELLURIC'][tstar].keys()): cntr = int(redux['PARAMETERS']['TELLURIC'][tstar]['CENTER'])
			src_wnd = copy.deepcopy(src_wnd0)
			if 'WINDOW' in list(redux['PARAMETERS']['TELLURIC'][tstar].keys()): src_wnd = redux['PARAMETERS']['TELLURIC'][tstar]['WINDOW']
			bck_wnd = copy.deepcopy(bck_wnd0)
			if 'BACK' in list(redux['PARAMETERS']['TELLURIC'][tstar].keys()): bck_wnd = redux['PARAMETERS']['TELLURIC'][tstar]['BACK']
# reduce imaging data
			imr,var = reduceScienceImage(im,redux['BIAS'],redux['FLAT'],redux['MASK'],hd=hd)
#			imr,var = reduceScienceImage(im,redux['BIAS'],redux['FLAT'],redux['MASK'],hd['EXPTIME'])
# trace source and rectify
			cntr = findPeak(imr,cntr=cntr)
			trace = traceDispersion(imr,cntr=cntr,window=src_wnd,method='maximum',plot_file='{}/diagnostic_trace_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
			imrect = rectify(imr,trace)
			varrect = rectify(var,trace)
			maskrect = rectify(redux['MASK'],trace)
# find center
			cntr = findPeak(imrect,cntr=cntr,window=src_wnd,plot_file='{}/diagnostic_profile_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
# make profile
#			profile = spatialProfile(imrect,cntr=cntr,window=src_wnd)
#			profile = numpy.ones(int(2*src_wnd+1))
# extract spectrum
			spflx = extractSpectrum(imrect,cntr=cntr,var=varrect,mask=maskrect,src_wnd=src_wnd,bck_wnd=bck_wnd,method=redux['PARAMETERS']['TELLURIC'][tstar]['METHOD'],plot_file='{}/diagnostic_extraction_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
			spflx.name = tstar
			spflx.header = hd
# reidentify arc lines
			arcdp,arcdphd = readFiles(redux['PARAMETERS']['ARC_DEEP']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
			arcrect = rectify(arcdp,trace)
			arcrecal = waveCalibrateArcs(arcrect,trace=trace,prior=redux['CAL_WAVE'],mode=redux['PARAMETERS']['MODE'],plot_file='{}/diagnostic_wavecal_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
# apply wavelength calibration
			spflx.applyWaveCal(arcrecal)
#			spflx.applyWaveCal(redux['CAL_WAVE'])
# apply flux calibration
			spflx.applyFluxCal(redux['CAL_FLUX'],debug=False,log=True)
# compute telluric corection
#			print(tstar,redux['PARAMETERS']['TELLURIC'][tstar]['SPT'])
			redux['CAL_TELL'][tstar] = telluricCalibrate(spflx,spt=redux['PARAMETERS']['TELLURIC'][tstar]['SPT'],plot_file='{}/diagnostic_telluric_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
			redux['CAL_TELL'][tstar]['NAME'] = tstar
			redux['CAL_TELL'][tstar]['TRACE'] = trace
#			redux['CAL_TELL'][tstar]['PROFILE'] = profile
#			spflx.applyTelluricCal(redux['CAL_TELL'][tstar])
			redux[tstar] = spflx
# save this to pickle file
			f = open('{}/cal_tell_{}_{}.pkl'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']),'wb')
			pickle.dump(redux['CAL_TELL'][tstar],f)
			f.close()
# plot and save
			fig = spflx.plot(ylim=[0,1.2*numpy.quantile(spflx.flux.value,0.95)])
			fig.figure.savefig('{}/kast{}_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],tstar,redux['PARAMETERS']['DATE']))
			plt.close()
#			plt.clf()
			spflx.toFile('{}/kast{}_{}_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],tstar,redux['PARAMETERS']['DATE']))
			spflx.toFile('{}/kast{}_{}_{}.txt'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],tstar,redux['PARAMETERS']['DATE']))

# determine second-order flux corrections from telluric standards
	if ('CAL_REFLUX' not in list(redux.keys()) or reset==True or kwargs.get('reset_fluxrecal',False) == True) and 'TELLURIC' in list(redux['PARAMETERS'].keys()) and redux['PARAMETERS']['MODE'] != 'BLUE':
		if verbose==True: print('\nComputing second order flux corrections from telluric standards')
		redux['CAL_REFLUX'] = {}
		for tstar in list(redux['PARAMETERS']['TELLURIC'].keys()):
			redux['CAL_REFLUX'][tstar] = fluxReCalibrate(redux[tstar],spt=redux['PARAMETERS']['TELLURIC'][tstar]['SPT'],plot_file='{}/diagnostic_reflux_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
# save this to pickle file
			f = open('{}/cal_reflux_{}_{}.pkl'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']),'wb')
			pickle.dump(redux['CAL_REFLUX'][tstar],f)
			f.close()

# now analyze science targets - ONLY DOING FIRST FILE
	for src in list(redux['PARAMETERS']['SOURCE'].keys()):
		if src not in list(redux.keys()) or reset==True or kwargs.get('reset_source',False) == True:
			if verbose==True: print('\nExtracting {} spectrum of {}'.format(redux['PARAMETERS']['MODE'],src))
#			print(redux['PARAMETERS']['SOURCE'][src])
			ims,hds = readFiles(redux['PARAMETERS']['SOURCE'][src]['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'],verbose=verbose,)
			if len(redux['PARAMETERS']['SOURCE'][src]['FILES']) > 1:
				im = crRejectCombine(ims,verbose=verbose)
				hd = hds[0]
			else:
				im = ims
				hd = hds
# analyze image data
			imr,var = reduceScienceImage(im,redux['BIAS'],redux['FLAT'],redux['MASK'],hd=hd,verbose=verbose)
#			imr,var = reduceScienceImage(im,redux['BIAS'],redux['FLAT'],redux['MASK'],hd['EXPTIME'],verbose=verbose,)
# set extraction parameters
			src_wnd = copy.deepcopy(src_wnd0)
			if 'WINDOW' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				src_wnd = redux['PARAMETERS']['SOURCE'][src]['WINDOW']
			bck_wnd = copy.deepcopy(bck_wnd0)
			if 'BACK' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				bck_wnd = redux['PARAMETERS']['SOURCE'][src]['BACK']
			if 'CENTER' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				cntr = int(redux['PARAMETERS']['SOURCE'][src]['CENTER'])
				print('using input center {}'.format(cntr))
			else: 
				cntr = findPeak(imrect,plot_file='{}/diagnostic_profile_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],src))
				print('found new center {}'.format(cntr))
			tcorr,ttrace = False,False
			if 'TELLURIC' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				tcorr,ttrace,reflux = True,True,True
				tstar = redux['PARAMETERS']['SOURCE'][src]['TELLURIC']
			if 'APPLY_TRACE' in list(redux['PARAMETERS']['SOURCE'][src].keys()) and tcorr==True: 
				ttrace = redux['PARAMETERS']['SOURCE'][src]['APPLY_TRACE'].upper()=='TRUE'
			cflag = True
			if 'RECENTER' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				cflag = redux['PARAMETERS']['SOURCE'][src]['RECENTER'].upper()=='TRUE'
			# pflag = 'SOURCE'
			# if 'PROFILE' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
			# 	pflag = redux['PARAMETERS']['SOURCE'][src]['PROFILE'].upper()
# trace source or apply input trace
			if ttrace==False: 
				if verbose==True: print('Tracing science target'.format(src))
				trace = traceDispersion(imr,cntr=cntr,window=src_wnd,method='maximum',verbose=verbose,plot_file='{}/diagnostic_trace_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],src))
			else: 
				if verbose==True: print('Using trace from telluric standard {}'.format(tstar))
				trace=redux['CAL_TELL'][tstar]['TRACE']
			imrect = rectify(imr,trace)
			varrect = rectify(var,trace)
			maskrect = rectify(redux['MASK'],trace)
# recenter?
			if cflag == True or cntr<0: 
				cntr = findPeak(imrect,cntr=cntr,window=src_wnd,plot_file='{}/diagnostic_profile_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],src,redux['PARAMETERS']['MODE']))
				if verbose==True: print('Recentered source to {}'.format(cntr))
# choose profile and extract
			profile = []
			# if pflag.upper()=='SOURCE': 
			# 	if verbose==True: print('Using spatial profile measured from source')
			# 	profile = spatialProfile(imrect,cntr=cntr,window=src_wnd)
			if 'PROFILE' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				if redux['PARAMETERS']['SOURCE'][src]['PROFILE'].upper() == 'TELLURIC' and tcorr==True:
					profile = redux['CAL_TELL'][tstar]['PROFILE']
			# if pflag.upper()=='TELLURIC' and tcorr==True: 
			# 	if verbose==True: print('Using spatial profile measured from telluric standard {}'.format(tstar))
			# else: 
			# 	if verbose==True: print('Using flat spatial profile')
			# 	profile = numpy.ones(int(2*src_wnd+1))
			spflx = extractSpectrum(imrect,var=varrect,mask=maskrect,cntr=cntr,src_wnd=src_wnd,bck_wnd=bck_wnd,profile=profile,method=redux['PARAMETERS']['SOURCE'][src]['METHOD'],plot_file='{}/diagnostic_extraction_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],src,redux['PARAMETERS']['MODE']))
			spflx.name = src
			spflx.header = hd
# reapply arc solution
			arcdp,arcdphd = readFiles(redux['PARAMETERS']['ARC_DEEP']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
			arcrect = rectify(arcdp,trace)
			arcrecal = waveCalibrateArcs(arcrect,trace=trace,prior=redux['CAL_WAVE'],mode=redux['PARAMETERS']['MODE'],plot_file='{}/diagnostic_wavecal_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],src,redux['PARAMETERS']['MODE']))
			spflx.applyWaveCal(arcrecal)
# apply flux calibration
			if 'CAL_FLUX' in list(redux.keys()): 
				print('Applying flux calibration')
				spflx.applyFluxCal(redux['CAL_FLUX'],log=True)
			else: print('Warning: no flux calibration applied')
# apply second order flux calibration
# note: could also do this with wave and correction keywords with applyTelluricCal
			if 'CAL_REFLUX' in list(redux.keys()) and 'TELLURIC' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				print('Applying second-order flux calibration from {}'.format(tstar))
				spflx.applyFluxCal(redux['CAL_REFLUX'][tstar],log=False)
			else: print('Warning: no second order flux calibration applied')
# apply telluric correction
			if tcorr==True: 
				print('Applying telluric correction from {}'.format(tstar))				
				spflx.applyTelluricCal(redux['CAL_TELL'][tstar])
			else: print('Warning: no telluric correction applied')
			redux[src] = spflx
# plot and export spectrum
			fig = spflx.plot(ylim=[0,1.2*numpy.quantile(spflx.flux.value,0.95)])
			fig.figure.savefig('{}/kast{}_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],src,redux['PARAMETERS']['DATE']))
			plt.close()
#			plt.clf()
			spflx.toFile('{}/kast{}_{}_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],src,redux['PARAMETERS']['DATE']))
			spflx.toFile('{}/kast{}_{}_{}.txt'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],src,redux['PARAMETERS']['DATE']))

# save reduction information to a pickle file
	f = open('{}/reduction_structure_{}.pkl'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']),'wb')
	pickle.dump(redux,f)
	f.close()

	return




############################################################
# KAST SPECTRAL ANALYSIS FUNCTIONS
# These functions perform some basic analysis on spectra
############################################################

# PROGRAM CONSTANTS FOR ANALYSIS
INDEX_SETS = {
	'kirkpatrick1999': {'altname': ['kirkpatrick','kirkpatrick19','kir99'], 'bibcode': '1999ApJ...519..802K', 'indices': {\
		'Rb-a': {'ranges': ([.77752,.77852]*u.micron,[.78152,.78252]*u.micron,[.77952,.78052]*u.micron), 'method': 'line', 'sample': 'integrate'},\
		'Rb-b': {'ranges': ([.79226,.79326]*u.micron,[.79626,.79726]*u.micron,[.79426,.79526]*u.micron), 'method': 'line', 'sample': 'integrate'},\
		'Na-a': {'ranges': ([.81533,.81633]*u.micron,[.81783,.81883]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'Na-b': {'ranges': ([.81533,.81633]*u.micron,[.81898,.81998]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'Cs-a': {'ranges': ([.84961,.85061]*u.micron,[.85361,.85461]*u.micron,[.85161,.85261]*u.micron), 'method': 'line', 'sample': 'integrate'},\
		'Cs-b': {'ranges': ([.89185,.89285]*u.micron,[.89583,.89683]*u.micron,[.89385,.89485]*u.micron), 'method': 'line', 'sample': 'integrate'},\
		'TiO-a': {'ranges': ([.7033,.7048]*u.micron,[.7058,.7073]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'TiO-b': {'ranges': ([.8400,.8415]*u.micron,[.8435,.8470]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'VO-a': {'ranges': ([.7350,.7370]*u.micron,[.7550,.7570]*u.micron,[.7430,.7470]*u.micron), 'method': 'sumnum', 'sample': 'integrate'},\
		'VO-b': {'ranges': ([.7860,.7880]*u.micron,[.8080,.8100]*u.micron,[.7960,.8000]*u.micron), 'method': 'sumnum', 'sample': 'integrate'},\
		'CrH-a': {'ranges': ([.8580,.8600]*u.micron,[.8621,.8641]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'CrH-b': {'ranges': ([.9940,.9960]*u.micron,[.9970,.9990]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'FeH-a': {'ranges': ([.8660,.8680]*u.micron,[.8700,.8720]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'FeH-b': {'ranges': ([.9863,.9883]*u.micron,[.9908,.9928]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'Color-a': {'ranges': ([.9800,.9850]*u.micron,[.7300,.7350]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'Color-b': {'ranges': ([.9800,.9850]*u.micron,[.7000,.7050]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'Color-c': {'ranges': ([.9800,.9850]*u.micron,[.8100,.8150]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'Color-d': {'ranges': ([.9675,.9850]*u.micron,[.7350,.7550]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
	}},\
	'martin1999': {'altname': ['martin','martin99','mar99'], 'bibcode': '1999AJ....118.2466M', 'indices': {\
		'PC3': {'ranges': ([.823,.827]*u.micron,[.754,.758]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'PC6': {'ranges': ([.909,.913]*u.micron,[.650,.654]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'CrH1': {'ranges': ([.856,.860]*u.micron,[.861,.865]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'CrH2': {'ranges': ([.984,.988]*u.micron,[.997,1.001]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'FeH1': {'ranges': ([.856,.860]*u.micron,[.8685,.8725]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'FeH2': {'ranges': ([.984,.988]*u.micron,[.990,.994]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'H2O1': {'ranges': ([.919,.923]*u.micron,[.928,.932]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'TiO1': {'ranges': ([.700,.704]*u.micron,[.706,.710]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'TiO2': {'ranges': ([.838,.842]*u.micron,[.844,.848]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'VO1': {'ranges': ([.754,.758]*u.micron,[.742,.746]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
		'VO2': {'ranges': ([.799,.803]*u.micron,[.790,.794]*u.micron), 'method': 'ratio', 'sample': 'integrate'},\
	}},\
	'gizis1997': {'altname': ['gizis','gizis97','giz97'], 'bibcode': '', 'indices': {\
		'CaH1': {'ranges': [[6380,6390]*u.Angstrom,[6410,6420]*u.Angstrom,[6345,6355]*u.Angstrom],'method': 'avedenom','sample': 'average'},\
		'CaH2': {'ranges': [[6814, 6846]*u.Angstrom,[7042, 7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'CaH3': {'ranges': [[6960,6990]*u.Angstrom,[7042,7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'TiO5': {'ranges': [[7126,7135]*u.Angstrom,[7042,7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
	}},\
	'lepine2003': {'altname': ['lepine','lepine03','lep03'], 'bibcode': '', 'indices': {\
		'CaH1': {'ranges': [[6380,6390]*u.Angstrom,[6410,6420]*u.Angstrom,[6345,6355]*u.Angstrom],'method': 'avedenom','sample': 'average'},\
		'CaH2': {'ranges': [[6814, 6846]*u.Angstrom,[7042, 7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'CaH3': {'ranges': [[6960,6990]*u.Angstrom,[7042,7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'TiO5': {'ranges': [[7126,7135]*u.Angstrom,[7042,7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'VO1': {'ranges': [[7430, 7470]*u.Angstrom,[7550,7570]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'TiO6': {'ranges': [[7550,7570]*u.Angstrom,[7745,7765]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'VO2': {'ranges': [[7920, 7960]*u.Angstrom,[8440, 8470]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'TiO7': {'ranges': [[8440, 8470]*u.Angstrom,[8400, 8420]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'Color-M': {'ranges': [[8105, 8155]*u.Angstrom,[6510, 6560]*u.Angstrom],'method': 'ratio','sample': 'average'},\
	}},\
	'reid1995': {'altname': ['reid','reid95','rei95'], 'bibcode': '', 'indices': {\
		'TiO1': {'ranges': [[6718, 6723]*u.Angstrom,[6703, 6708]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'TiO2': {'ranges': [[7058, 7061]*u.Angstrom,[7043, 7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'TiO3': {'ranges': [[7092, 7097]*u.Angstrom,[7079, 7084]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'TiO4': {'ranges': [[7130, 7135]*u.Angstrom,[7115, 7120]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'TiO5': {'feature': [[7126, 7135]*u.Angstrom,[7042, 7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'CaH1': {'ranges': [[6380,6390]*u.Angstrom,[6410,6420]*u.Angstrom,[6345,6355]*u.Angstrom],'method': 'avdenom','sample': 'average'},\
		'CaH2': {'ranges': [[6814, 6846]*u.Angstrom,[7042, 7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'CaH3': {'ranges': [[6960,6990]*u.Angstrom,[7042,7046]*u.Angstrom],'method': 'ratio','sample': 'average'},\
		'CaOH': {'ranges': [[6230, 6240]*u.Angstrom,[6345, 6354]*u.Angstrom],'method': 'ratio','sample': 'average'},\
	}},\
	'burgasser2003': {'altname': ['burgasser','burgasser03','bur03'], 'bibcode': '', 'indices': {\
		'CsI-A': {'range': [[8496.1, 8506.1]*u.Angstrom, [8536.1, 8546.1]*u.Angstrom,[8516.1, 8626.1]*u.Angstrom],'method': 'sumnum_twicedenom','sample': 'average'},\
		'CsI-B': {'range': [[8918.5, 8928.5]*u.Angstrom, [8958.3, 8968.3]*u.Angstrom,[8938.5, 8948.3]*u.Angstrom],'method': 'sumnum_twicedenom','sample': 'average'},\
		'H2O': {'range': [[9220, 9240]*u.Angstrom,[9280, 9300]*u.Angstrom],'method': 'ratio','sample': 'integrate'},\
		'CrH-A': {'range': [[8560, 8600]*u.Angstrom,[8610, 8650]*u.Angstrom],'method': 'ratio','sample': 'integrate'},\
		'CrH-B': {'range': [[9855, 9885]*u.Angstrom,[9970, 10000]*u.Angstrom],'method': 'ratio','sample': 'integrate'},\
		'FeH-A': {'range': [[8560, 8600]*u.Angstrom,[8685, 8725]*u.Angstrom],'method': 'ratio','sample': 'integrate'},\
		'FeH-B': {'range': [[9855, 9885]*u.Angstrom,[9905, 9935]*u.Angstrom],'method': 'ratio','sample': 'integrate'},\
		'Color-e': {'range': [[9140, 9240]*u.Angstrom,[8400, 8500]*u.Angstrom],'method': 'ratio','sample': 'average'},\
	}},\
}

ZETA_RELATIONS = {
	'lepine2007': {'altname': ['lepine','lepine07','lep07'], 'reference': 'Lepine et al. (2007)','bibcode':'2007ApJ...669.1235L', 'coeff': [-0.164,0.67,-0.118,-0.05],'range': []},
	'dhital2011': {'altname': ['dhital','dhital11','dhi11'], 'reference': 'Dhital et al. (2011)','bibcode':'', 'coeff': [-0.005,-0.183,0.694,-0.127,-0.047],'range': []},
}

EW_LINES = {
	'H I': {'altname': ['H','H1','HI'], 'lines': [4861*u.Angstrom,6563*u.Angstrom]},\
	'Li I': {'altname': ['Li','Li1','LiI'], 'lines': [6708*u.Angstrom]},\
	'K I': {'altname': ['K','K1','KI'], 'lines': [7665*u.Angstrom,7699*u.Angstrom]},\
	'Na I': {'altname': ['Na','Na1','NaI'], 'lines': [8183*u.Angstrom,8195*u.Angstrom]},\
	'Cs I': {'altname': ['Cs','Cs1','CsI'], 'lines': [8521*u.Angstrom,8943*u.Angstrom]},\
	'Rb I': {'altname': ['Rb','Rb1','RbI'], 'lines': [7800*u.Angstrom,7948*u.Angstrom]},\
	'Ca I': {'altname': ['Ca','Ca1','CaI'], 'lines': [7209*u.Angstrom,7213*u.Angstrom,7326*u.Angstrom]},\
	'Ca II': {'altname': ['Ca2','CaII'], 'lines': [8436*u.Angstrom]},\
	'Ti I': {'altname': ['Ti','Ti1','Ti I'], 'lines': [6573*u.Angstrom,8542*u.Angstrom]},\
}

CHI_RELATIONS = {
    'schmidt2014': {'altname': ['schmidt','schmidt14'], 'reference': 'Schmidt et al. (2014)','bibcode':'2014PASP..126..642S', 'sptoffset': 0, 'method': 'interpolate', 'scale': 1.e-6,
        'spt': [17,18,19,20,21,22,23,24,25,26,27], \
        'values': [10.28,4.26,2.52,1.98,2.25,2.11,1.67,1.16,1.46,1.23,0.73],\
        'scatter': [3.13,1.18,0.58,0.27,0.11,0.36,0.22,0.3,0.28,0.3,0.3],\
        },
    'douglas2014': {'altname': ['douglas','douglas14'], 'reference': 'Douglas et al. (2014)','bibcode':'2014ApJ...795..161D', 'sptoffset': 0, 'method': 'interpolate', 'scale': 1.e-5,
        'spt': [10,11,12,13,14,15,16,17,18,19], \
        'values': [6.6453,6.0334,5.2658,4.4872,3.5926,2.4768,1.7363,1.2057,0.6122,0.3522],\
        'scatter': [0.6207,0.5326,0.5963,0.4967,0.5297,0.4860,0.3475,0.3267,0.2053,0.1432],\
        },
}

# HELPFUL FUNCTIONS
def padWhereArray(w,mx):
	'''
	Purpose:

		Pads the output of a numpy.where array to select (if available) one more index spot beyond limits

	'''
	if w[0][0] > 0: w = (numpy.insert(w[0],0,w[0][0]-1),)
	if w[0][-1] < mx: w = (numpy.append(w[0],w[0][-1]+1),)
	return w



def compareSpectra(sp1,sp2orig,fit_range=[],fitcycle=5,sclip=3.,plot=False,plot_file='',verbose=ERROR_CHECKING,**kwargs):
	'''
	Compares to spectra to each other and returns the best fit statistic and scale factor
	Input: spectra objects
	Output: fit stat and scale factor 

	NOTE
	provide a way of doing comparision without uncertainties
	'''
# make sure both spectra conform to same wave and flux units
	sp2 = copy.deepcopy(sp2orig)
	try:
		sp2.convertWave(sp1.wave.unit)
	except:
		raise ValueError('Cannot convert second spectrum with wave units {} to wave units {}'.format(sp2.wave.unit,sp1.wave.unit))
	try:
		sp2.convertFlux(sp1.flux.unit)
	except:
		raise ValueError('Cannot convert second spectrum with flux units {} to flux units {}'.format(sp2.flux.unit,sp1.flux.unit))

# interpolate onto a common wavelength scale
	wave = sp1.wave.value
	f1 = sp1.flux.value
	u1 = sp1.unc.value
	f2interp = interp1d(sp2.wave.value,sp2.flux.value,bounds_error=False,fill_value=0.)
	u2interp = interp1d(sp2.wave.value,sp2.unc.value,bounds_error=False,fill_value=0.)
	f2 = f2interp(wave)
	u2 = u2interp(wave)

# check uncertainties
	if numpy.isnan(numpy.median(u1))==True and numpy.isnan(numpy.median(u2))==False: vtot = u2**2
	elif numpy.isnan(numpy.median(u1))==False and numpy.isnan(numpy.median(u2))==True: vtot = u1**2
	elif numpy.isnan(numpy.median(u1))==True and numpy.isnan(numpy.median(u2))==True: vtot=numpy.ones(len(wave))
	else: vtot = u1**2+u2**2

# fit range
	weights = numpy.ones(len(wave))
	if len(fit_range)>1:
		weights[wave<numpy.nanmin(fit_range)]=0
		weights[wave>numpy.nanmax(fit_range)]=0

# using just chi2 for now
	scale_factor = numpy.nansum(weights*f1*f2/vtot)/numpy.nansum(weights*(f2**2)/vtot)
	stat = numpy.nansum(weights*(f1-f2*scale_factor)**2/vtot)

# add in a little iteration
	if fitcycle>1:
		for i in range(fitcycle):
			wdiff = numpy.absolute(weights*(f1-f2*scale_factor)**2/vtot)
			weights[wdiff>sclip]==0
			scale_factor = numpy.nansum(weights*f1*f2/vtot)/numpy.nansum(weights*(f2**2)/vtot)
			stat = numpy.nansum(weights*(f1-f2*scale_factor)**2/vtot)

	if plot==True or plot_file!='':
		uncplot = vtot**0.5
		if numpy.isnan(numpy.median(u1))==True and numpy.isnan(numpy.median(u2))==True: uncplot = numpy.zeros(len(wave))
		print(numpy.nanmedian(uncplot))
		xlim = kwargs.get('xlim',[numpy.nanmin(sp1.wave.value),numpy.nanmax(sp1.wave.value)])
		ylim = kwargs.get('ylim',[numpy.nanmin([0,numpy.quantile(f1,0.05)]),1.2*numpy.quantile(f1,0.95)])
		fig,(ax1,ax2) = plt.subplots(2,1,sharex='col',figsize=kwargs.get('figsize',[8,8]))
		ax1.plot(wave,f1,c=kwargs.get('color',PLOT_DEFAULTS['color']),ls=kwargs.get('ls',PLOT_DEFAULTS['ls']),alpha=kwargs.get('alpha',PLOT_DEFAULTS['alpha']))
		ax1.plot(wave,f2*scale_factor,c=kwargs.get('color',PLOT_DEFAULTS['comparison_color']),ls=kwargs.get('ls',PLOT_DEFAULTS['comparison_ls']),alpha=kwargs.get('alpha',PLOT_DEFAULTS['comparison_alpha']))
		ax1.legend(kwargs.get('legend',[sp1.name,sp2.name]),fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
		ax1.plot(wave,uncplot**0.5,c=kwargs.get('unc_color',PLOT_DEFAULTS['unc_color']),ls=kwargs.get('unc_ls',PLOT_DEFAULTS['unc_ls']),alpha=kwargs.get('unc_alpha',PLOT_DEFAULTS['unc_alpha']))
		ax1.plot(wave,numpy.zeros(len(wave)),c=kwargs.get('zero_color',PLOT_DEFAULTS['zero_color']),ls=kwargs.get('zero_ls',PLOT_DEFAULTS['zero_ls']),alpha=kwargs.get('zero_alpha',PLOT_DEFAULTS['zero_alpha']))
		ax1.fill_between(wave,numpy.zeros(len(wave))+ylim[0],(1-weights)*ylim[1],facecolor='grey',alpha=0.2)
#		ax1.set_xlim(kwargs.get('xlim',[numpy.nanmin(wave),numpy.nanmax(wave)]))
		ax1.set_ylim(ylim)
		ax1.set_ylabel(kwargs.get('ylabel',r'Flux ({})'.format(sp1.flux.unit)),fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
#		ax1.set_xticks(fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
#		ax1.set_yticks(fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
		diff = f1-f2*scale_factor
		ylim2 = kwargs.get('ylim2',[-3.*numpy.nanstd(diff),3.*numpy.nanstd(diff)])
		ax2.plot(wave,diff,c=kwargs.get('background_color',PLOT_DEFAULTS['background_color']),ls=kwargs.get('background_ls',PLOT_DEFAULTS['background_ls']),alpha=kwargs.get('background_alpha',PLOT_DEFAULTS['background_alpha']))
		ax2.plot(wave[weights==1],diff[weights==1],c=kwargs.get('color',PLOT_DEFAULTS['color']),ls=kwargs.get('ls',PLOT_DEFAULTS['ls']),alpha=kwargs.get('alpha',PLOT_DEFAULTS['alpha']))
		ax2.fill_between(wave,uncplot,-1.*uncplot,facecolor=kwargs.get('unc_color',PLOT_DEFAULTS['unc_color']),linestyle=kwargs.get('unc_ls',PLOT_DEFAULTS['unc_ls']),alpha=kwargs.get('unc_alpha',PLOT_DEFAULTS['unc_alpha']))
#		ax2.plot(wave,-1.*sp.unc.value,c=kwargs.get('unc_color',PLOT_DEFAULTS['unc_color']),ls=kwargs.get('unc_ls',PLOT_DEFAULTS['unc_ls']),alpha=kwargs.get('unc_alpha',PLOT_DEFAULTS['unc_alpha']))
		ax2.set_xlabel(kwargs.get('xlabel',r'Wavelength ({})'.format(sp1.wave.unit)),fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
		ax2.set_ylabel(kwargs.get('ylabel','O-C'),fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
		ax2.set_xlim(xlim)
		ax2.set_ylim(ylim2)
#		ax2.set_xticks(fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
#		ax2.set_yticks(fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
		ax2.plot(wave,numpy.zeros(len(wave)),c=kwargs.get('zero_color',PLOT_DEFAULTS['zero_color']),ls=kwargs.get('zero_ls',PLOT_DEFAULTS['zero_ls']),alpha=kwargs.get('zero_alpha',PLOT_DEFAULTS['zero_alpha']))
		if plot_file!='': 
			fig.savefig(plot_file)
			plt.close()
#		plt.clf()

	return stat, scale_factor



def compareSpectra_simple(sp1,sp2orig,fit_range=[],exclude_range=[],error_value=numpy.nan,plot=False,plot_file='',verbose=ERROR_CHECKING,**kwargs):
	'''
	A stripped down version of compareSpectra() to address errors 
	'''
	sp2 = copy.deepcopy(sp2orig)
	sp2.convertWave(sp1.wave.unit)
	sp2.convertFlux(sp1.flux.unit)
	wave = sp1.wave.value
	f1 = sp1.flux.value
	u1 = sp1.unc.value
	f2interp = interp1d(sp2.wave.value,sp2.flux.value,bounds_error=False,fill_value=0.)
	f2 = f2interp(wave)
	vtot = u1**2
	weights = numpy.ones(len(wave))
	if len(fit_range)>1:
		weights[wave<numpy.nanmin(fit_range)]=0
		weights[wave>numpy.nanmax(fit_range)]=0
	if len(exclude_range)>1:
		weights[numpy.logical_and(wave>numpy.nanmin(exclude_range),wave<numpy.nanmax(exclude_range))]=0
	if numpy.nansum(weights)>0:
		scale_factor = numpy.nansum(weights*f1*f2/vtot)/numpy.nansum(weights*f2*f2/vtot)
		stat = numpy.nansum(weights*(f1-f2*scale_factor)**2/vtot)
		if plot==True: 
			plt.clf()
			xlim = kwargs.get('xlim',[numpy.nanmin(wave),numpy.nanmax(wave)])
			ylim = kwargs.get('ylim',[-2.*numpy.nanmedian(u1),1.2*numpy.nanquantile(f1,0.95)])
			fig = plt.figure(figsize=[8,6])
			grid = plt.GridSpec(4, 1, hspace=0, wspace=0.2)
			ax_top = fig.add_subplot(grid[:-1,0])
			ax_btm = fig.add_subplot(grid[-1,0])
			ax_top.plot(wave,f1,'k-')
			ax_top.plot(wave,f2*scale_factor,'m-')
			ax_top.legend([sp1.name,sp2.name],fontsize=16)
			ax_top.plot(wave,u1,'k--')
			ax_top.plot(wave,numpy.zeros(len(wave)),'k:')
			ax_top.set_xlim(xlim)
			ax_top.set_ylim(ylim)
			ax_top.set_ylabel('Flux Density',fontsize=16)
			ax_btm.plot(wave,f1-f2*scale_factor,'k-')
			ax_btm.plot(wave,u1,'k--')
			ax_btm.plot(wave,-1.*u1,'k--')
			ax_btm.fill_between(wave,-1.*u1,u1,color='k',alpha=0.1)
			ax_btm.plot(wave,numpy.zeros(len(wave)),'k--')
			ax_btm.set_xlim(xlim)
			ax_btm.set_ylim([-5.*numpy.nanmedian(u1),5.*numpy.nanmedian(u1)])
			ax_btm.set_xlabel('Wavelength (Angstrom)',fontsize=16)
			ax_btm.set_ylabel('O-C',fontsize=16)
			wv = wave[weights==0]
			if len(wv)>0:
				for w in wv: 
					ax_top.plot(ylim,[w,w],c='grey',ls='-',alpha=0.3)
					ax_btm.plot([-3.*numpy.nanmedian(u1),3.*numpy.nanmedian(u1)],[w,w],c='grey',ls='-',alpha=0.3)
			if plot_file!='': 
				fig.savefig(plot_file)
				plt.close()

	else:
		scale_factor = error_value 
		stat = error_value
		if verbose==True: print('No data matched in fit range or was in non-excluded regions')


# 		fig,(ax1,ax2) = plt.subplots(2,1,sharex='col',figsize=kwargs.get('figsize',[8,8]))
# 		ax1.plot(wave,f1,c=kwargs.get('color',PLOT_DEFAULTS['color']),ls=kwargs.get('ls',PLOT_DEFAULTS['ls']),alpha=kwargs.get('alpha',PLOT_DEFAULTS['alpha']))
# 		ax1.plot(wave,f2*scale_factor,c=kwargs.get('color',PLOT_DEFAULTS['comparison_color']),ls=kwargs.get('ls',PLOT_DEFAULTS['comparison_ls']),alpha=kwargs.get('alpha',PLOT_DEFAULTS['comparison_alpha']))
# 		ax1.legend(kwargs.get('legend',[sp1.name,sp2.name]),fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
# 		ax1.plot(wave,u1,c=kwargs.get('unc_color',PLOT_DEFAULTS['unc_color']),ls=kwargs.get('unc_ls',PLOT_DEFAULTS['unc_ls']),alpha=kwargs.get('unc_alpha',PLOT_DEFAULTS['unc_alpha']))
# 		ax1.plot(wave,numpy.zeros(len(wave)),c=kwargs.get('zero_color',PLOT_DEFAULTS['zero_color']),ls=kwargs.get('zero_ls',PLOT_DEFAULTS['zero_ls']),alpha=kwargs.get('zero_alpha',PLOT_DEFAULTS['zero_alpha']))
# 		ax1.fill_between(wave,-1.*u1,u1,facecolor='k',alpha=0.1)
# #		ax1.set_xlim(kwargs.get('xlim',[numpy.nanmin(wave),numpy.nanmax(wave)]))
# 		ax1.set_xlim(xlim)
# 		ax1.set_ylim(ylim)
# 		ax1.set_ylabel(kwargs.get('ylabel',r'Flux ({})'.format(sp1.flux.unit)),fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
# #		ax1.set_xticks(fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
# #		ax1.set_yticks(fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
# 		diff = f1-f2*scale_factor
# 		ylim2 = kwargs.get('ylim2',[-3.*numpy.nanstd(diff),3.*numpy.nanstd(diff)])
# 		ax2.plot(wave,diff,c=kwargs.get('background_color',PLOT_DEFAULTS['background_color']),ls=kwargs.get('background_ls',PLOT_DEFAULTS['background_ls']),alpha=kwargs.get('background_alpha',PLOT_DEFAULTS['background_alpha']))
# 		ax2.plot(wave[weights==1],diff[weights==1],c=kwargs.get('color',PLOT_DEFAULTS['color']),ls=kwargs.get('ls',PLOT_DEFAULTS['ls']),alpha=kwargs.get('alpha',PLOT_DEFAULTS['alpha']))
# 		ax2.fill_between(wave,u1,-1.*u1,facecolor=kwargs.get('unc_color',PLOT_DEFAULTS['unc_color']),linestyle=kwargs.get('unc_ls',PLOT_DEFAULTS['unc_ls']),alpha=kwargs.get('unc_alpha',PLOT_DEFAULTS['unc_alpha']))
# #		ax2.plot(wave,-1.*sp.unc.value,c=kwargs.get('unc_color',PLOT_DEFAULTS['unc_color']),ls=kwargs.get('unc_ls',PLOT_DEFAULTS['unc_ls']),alpha=kwargs.get('unc_alpha',PLOT_DEFAULTS['unc_alpha']))
# 		ax2.set_xlabel(kwargs.get('xlabel',r'Wavelength ({})'.format(sp1.wave.unit)),fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
# 		ax2.set_ylabel(kwargs.get('ylabel','O-C'),fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
# 		ax2.set_xlim(xlim)
# 		ax2.set_ylim(ylim2)
# #		ax2.set_xticks(fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
# #		ax2.set_yticks(fontsize=kwargs.get('fontsize',PLOT_DEFAULTS['fontsize']))
# 		ax2.plot(wave,numpy.zeros(len(wave)),c=kwargs.get('zero_color',PLOT_DEFAULTS['zero_color']),ls=kwargs.get('zero_ls',PLOT_DEFAULTS['zero_ls']),alpha=kwargs.get('zero_alpha',PLOT_DEFAULTS['zero_alpha']))
#		plt.clf()

	return stat,scale_factor



def combineSpectra():
	'''
	Combines spectra together
	Input: spectra objects
	Output: combined spectrum object 
	'''
	pass

def stitchOrders():
	'''
	Stitches together multiple orders into one spectrum
	Input: spectra objects
	Output: combined spectrum object 
	'''
	pass

def initializeStandards(spt,sdss=True,folder=SPTSTDFOLDER,reset=False,sd=False,esd=False,usd=False,beta=False,giant=False,gamma=False,verbose=ERROR_CHECKING):
	'''
	Read spectral standards into standards structure
	'''
	if isinstance(spt,list) == False: spt = [spt]
	if len(spt) == 0: spt = [0,28]
	if isinstance(spt[0],float) == True or isinstance(spt[0],int) == True: 
		if len(spt) == 2: spt = numpy.arange(spt[0],spt[1]+1)
		spt = [typeToNum(s) for s in spt]

	if sd==True and 'sd' not in spt[0]: spt = ['sd'+s for s in spt]
	if esd==True and 'esd' not in spt[0]: spt = ['esd'+s for s in spt]
	if usd==True and 'usd' not in spt[0]: spt = ['usd'+s for s in spt]
	if beta==True and 'b' not in spt[0]: spt = [s+'b' for s in spt]
	if gamma==True and 'g' not in spt[0]: spt = [s+'g' for s in spt]
	if giant==True and 'I' not in spt[0]: spt = [s+'I' for s in spt]

	for s in spt:
		if s not in list(SPTSTDS.keys()) or reset==True:
			f = numpy.array(glob.glob('{}/{}*.txt'.format(folder,s.replace('.',''))))
			if len(f) > 1:
				if sdss==True and len(f[['SDSS' in x for x in f]])>0: f = f[['SDSS' in x for x in f]]
# for subdwarfs
				elif (sd==True or esd==True or usd==True) and len(f[['lepine2007' in x for x in f]])>0: 
#					print(len(f))
					f = f[['lepine2007' in x for x in f]]
				else: pass
			if len(f) > 0:
				SPTSTDS[s] = readSpectrum(f[0],name='{} STD'.format(s))
			else: 
				if verbose==True: print('Warning: cannot find a spectral standard for type {}'.format(s))
	return

def classifyTemplate(spec,spt=[],plot_file='',fit_range=[6500,8800],verbose=ERROR_CHECKING,**kwargs):
	'''
	Classification program using templates from resource library
	this is a simplified version of program
	'''
	if not isinstance(spec,Spectrum): raise('classifyTemplate requires a kastredux Spectrum object as an input')
	if not isinstance(spt,list): spt = [spt]
	if len(spt)==0: spt = [50,78]
	if len(spt)==1: spt = spt*2
	if len(spt)==2: initializeStandards(spt)

	sts = []
	stds = list(SPTSTDS.keys())
	for s in stds:
		st,scl = compareSpectra_simple(spec,SPTSTDS[s],fit_range=fit_range,error_value=numpy.nan,plot=False,**kwargs)
		if numpy.isfinite(st) == False: sts.append(1.e10)
		elif st <= 0: sts.append(1.e10)
		else: sts.append(st)
		if verbose==True: print('\t{}: {:.2f}'.format(SPTSTDS[s],st))
	spt = stds[numpy.argmin(sts)]
	if plot_file != '': st,scl = compareSpectra_simple(spec,SPTSTDS[spt],fit_range=fit_range,error_value=numpy.nan,plot=True,plot_file=plot_file,**kwargs)
	else: st,scl = compareSpectra_simple(spec,SPTSTDS[spt],fit_range=fit_range,error_value=numpy.nan,plot=False,**kwargs)

	return spt


def measureIndex(sp,ranges,sample='median',method='ratio',nsamples=100,noiseFlag=False,plot=False,verbose=ERROR_CHECKING,**pkwargs):
	"""
	Measures a spectral index on a spectrum based on defined methodology

	Parameters
	----------
	sp : Spectrum class object
		Single spectrum class object, should contain wave, flux and unc array elements

	ranges : array of astropy Quantites
		array of two-element lists that indicate the measurement regions in unitted wavelengths
		should be of the form [[3000,4000]*u.Angstrom,[4000,5000]*u.Angstrom,[5000,6000]*u.Angstrom]

	sample : string, default = 'integrate'
		defines the method by which flux in each region is combined; options are:

			* 'median': median value using numpy.nanmedian (default)
			* 'average': average value using numpy.nanmean
			* 'integrate': intergration using numpy.trapz
			* 'sum': total value using numpy.nansum
			* 'maximum' or 'max': maximum value using numpy.nanmax
			* 'minimum' or 'min': minimum value using numpy.nanmin

	method : string, default = 'ratio'
		defines the method by which regions are combined, proscribed as:

			* 'ratio': returns range[0]/range[1] (default)
			* 'single': returns range[0]
			* 'line' or 'avenum': returns 0.5 * (range[0]+range[1])/range[2]
			* 'inverse_line' or 'avedenum': returns 2 * range[0]/(range[1]+range[2])
			* 'change': returns 2 * (range[0]-range[1])/(range[0]+range[1])
			* 'sumnum': returns (range[0]+range[1])/range[2] = 2 * line
			* 'sumdenom': returns range[0]/(range[1]+range[2]) = 0.5 * inverse_line
			* 'doubleratio': returns (range[0]/range[1])/(range[1]/range[2])
			* 'allers': returns [range[2]*(lam0-lam1)/(lam2-lam1) + range[1]*(lam2-lam0)/(lam2-lam1)]/range[0]

	nsamples : int, default = 100
		number of samples for Monte Carlo uncertainty estimation

	noiseFlag : bool, default = False
		if True, do not conduct Monte Carlo uncertainty estimation 

	plot : bool, default = False
		if True, display plot of spectral region with index measurements indicated
		NOTE: this functionality is not yet implemented  

	verbose : bool, default = False
		if True, give extra feedback in computation process  

	Returns
	-------
	tuple of index measurement and uncertainty
		if noiseFlag = True, uncertainty will by numpy.nan

	Examples
	--------

	TBD

	See Also
	--------

	measureIndexSet : measures a predefined set of indices contained in INDICES dictionary

	""" 	


# error checking on number of arguments provided
	sample_types = ['integrate','sum','average','median','maximum','max','minimum','min']
	method_types = ['ratio','single','line','inverse_line','change','avenum','avedenom','sumnum','sumdenom','doubleratio','allers']
	if sample.lower() not in sample_types: raise ValueError('sample method {} not in allowed list: {}'.format(sample,sample_types))
	if verbose: print('Using measurement method {}'.format(sample))
	if method.lower() not in method_types: raise ValueError('index measurement method {} not in allowed list: {}'.format(method,method_types))
	if verbose: print('Using index combination method {}'.format(method))

	if len(ranges) == 0: raise ValueError('measureIndex needs at least 1 sample region; zero given')
	if len(ranges) == 1: method = 'single'
	if not isinstance(ranges[0],list) and not isinstance(ranges[0],numpy.ndarray) and not isinstance(ranges[0],tuple) and not isinstance(ranges[0],set):
		raise ValueError('measureIndex needs a list of wavelength ranges, you entered {} which has type {}'.format(ranges,type(ranges[0])))
	for r in ranges:
		if len(r) != 2: raise ValueError('Problem with range {} in input ranges {}: must be 2-element array'.format(r,ranges))
	if (len(ranges) < 2 and (method in ['ratio','change'])):
		raise ValueError('Index method {} needs at least 2 sample regions'.format(method))
	if (len(ranges) < 3 and (method in ['line','inverse_line','avenum','avedenom','sumnum','sumdenom','doubleratio','allers'])):
		raise ValueError('Index method {} needs at least 3 sample regions'.format(method))


# define the sample vectors
	value = numpy.zeros(len(ranges))
	value_sim = numpy.zeros((len(ranges),nsamples))

# loop over all sampling regions
	for i,waveRng in enumerate(ranges):

# convert units
		if isUnit(waveRng):
			waveRng = (waveRng.to(sp.wave.unit)).value
		elif isUnit(waveRng[0]):
			waveRng = [(w.to(sp.wave.unit)).value for w in waveRng]
		else:
			waveRng = ((waveRng*DEFAULT_WAVE_UNIT).to(sp.wave.unit)).value
		xNum = (numpy.arange(0,nsamples+1.0)/nsamples)* \
			(numpy.nanmax(waveRng)-numpy.nanmin(waveRng))+numpy.nanmin(waveRng)

# identify measureable regions
		w = numpy.where(numpy.logical_and(\
			numpy.logical_and(numpy.isnan(sp.flux.value) == False,numpy.isnan(sp.unc.value) == False),\
			numpy.logical_and(sp.wave.value >= numpy.nanmin(waveRng),sp.wave.value <= numpy.nanmax(waveRng))))
		if len(w[0]) == 0:
			noiseFlag = True
			w = numpy.where(numpy.logical_and(\
				numpy.isnan(sp.flux) == False,\
				numpy.logical_and(sp.wave.value >= numpy.nanmin(waveRng),sp.wave.value <= numpy.nanmax(waveRng))))
		if len(w[0]) == 0:
			if verbose: print('Warning: no data in the wavelength range {}'.format(waveRng))
			return numpy.nan,numpy.nan

# compute intepolated flux and noise
#		w = padWhereArray(w,len(sp.wave))
		f = interp1d(sp.wave.value[w],sp.flux.value[w],bounds_error=False,fill_value=numpy.nan)
		yNum = f(xNum)
		if noiseFlag == False:
			s = interp1d(sp.wave.value[w],sp.unc.value[w],bounds_error=False,fill_value=0.)
			yNum_e = s(xNum)

# first compute the actual value
		if (sample == 'integrate'): 
			w = numpy.where(numpy.logical_and(numpy.isfinite(xNum)==True,numpy.isfinite(yNum)==True))
			value[i] = trapz(yNum[w],xNum[w])
		elif (sample == 'average'): value[i] = numpy.nanmean(yNum)
		elif (sample == 'sum'): value[i] = numpy.nansum(yNum)
		elif (sample == 'median'): value[i] = numpy.nanmedian(yNum)
		elif (sample == 'maximum'): value[i] = numpy.nanmax(yNum)
		elif (sample == 'minimum'): value[i] = numpy.nanmin(yNum)
		else: value[i] = numpy.nanmean(yNum)
		value_sim[i,:] = value[i]

# now do Monte Carlo measurement of value and uncertainty
# THERE IS A PROBLEM HERE - ALL VALUE_SIM = NAN
		if noiseFlag == False: 
			for j in numpy.arange(nsamples):

				yVar = numpy.random.normal(yNum,yNum_e)
				if (sample == 'integrate'): 
					w = numpy.where(numpy.logical_and(numpy.isfinite(xNum)==True,numpy.isfinite(yVar)==True))
					value_sim[i,j] = trapz(yVar[w],xNum[w])
				elif (sample == 'average'): value_sim[i,j] = numpy.nanmean(yVar)
				elif (sample == 'sum'): value_sim[i,j] = numpy.nansum(yVar)
				elif (sample == 'median'): value_sim[i,j] = numpy.nanmedian(yVar)
				elif (sample == 'maximum'): value_sim[i,j] = numpy.nanmax(yVar)
				elif (sample == 'minimum'): value_sim[i,j] = numpy.nanmin(yVar)
				else: value_sim[i,j] = numpy.nanmean(yVar)

# compute index based on defined method
# default is a simple ratio
		if (method == 'single'): 
			val,vals = value[0],value_sim[0,:]
		elif (method == 'ratio'): 
			val,vals = value[0]/value[1],value_sim[0,:]/value_sim[1,:]
		elif (method == 'line' or method == 'avenum'): 
			val,vals = 0.5*(value[0]+value[1])/value[2],0.5*(value_sim[0,:]+value_sim[1,:])/value_sim[2,:]
		elif (method == 'inverse_line' or method == 'avedenom'):
			val,vals = 2.*value[0]/(value[1]+value[2]),2.*value_sim[0,:]/(value_sim[1,:]+value_sim[2,:])
		elif (method == 'change'): 
			val,vals = 2.*(value[0]-value[1])/(value[0]+value[1]),2.*(value_sim[0,:]-value_sim[1,:])/(value_sim[0,:]+value_sim[1,:])
		elif (method == 'sumnum'):
			val,vals = (value[0]+value[1])/value[2],(value_sim[0,:]+value_sim[1,:])/value_sim[2,:]
		elif (method == 'sumdenom'):
			val,vals = value[0]/(value[1]+value[2]),value_sim[0,:]/(value_sim[1,:]+value_sim[2,:])
		elif (method == 'doubleratio'):
			val,vals = (value[0]/value[1])/(value[1]/value[2]),(value_sim[0,:]/value_sim[1,:])/(value_sim[1,:]/value_sim[2,:])
		elif (method == 'allers'):
			val = (((numpy.mean(ranges[0])-numpy.mean(ranges[1]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value[2] \
				+ ((numpy.mean(ranges[2])-numpy.mean(ranges[0]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value[1]) \
				/value[0]
			vals = (((numpy.mean(ranges[0])-numpy.mean(ranges[1]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value_sim[2,:] \
				+ ((numpy.mean(ranges[2])-numpy.mean(ranges[0]))/(numpy.mean(ranges[2])-numpy.mean(ranges[1])))*value_sim[1,:]) \
				/value_sim[0,:]
		else:
			val,vals = value[0]/value[1],value_sim[0,:]/value_sim[1,:]

# PLOTTING/VISUALIZATION - not yet implemented
		# if plot == True:
		#	 bands = []
		#	 for r in ranges: bands.append(r)
		#	 inddict = {pkwargs.get('name','Index'): {'ranges': bands, 'value': value}}
		#	 plotSpectrum(sp,bands=bands,bandlabels=[pkwargs.get('name','') for b in bands],**pkwargs)		

# output mean, standard deviation

	if noiseFlag == True: 
		if verbose: print('Index value = {}'.format(val))
		return val, numpy.nan
	if verbose: print('Index value = {}+/-{}'.format(val,numpy.nanstd(vals)))
	return val, numpy.nanstd(vals)




def measureIndexSet(sp,ref='lepine2003',index_info={},range_keyword='ranges',method_keyword='method',sample_keyword='sample',info=False,verbose=ERROR_CHECKING,**kwargs):
	"""
	Measures a set of predefined or user-defined spectral indices

	Parameters
	----------
	sp : Spectrum class object
		Single spectrum class object, should contain wave, flux and unc array elements

	ref : string, default = 'lepine2003'
		named index set contained in the global variable INDEX_SETS
		alternate names: reference, set

	index_info : dict, default = {}
		optional dictionary to define index set information, which should be structured as follows:
		alternate names: indices, index, index_dict

			* series of index name strings, each pointing to a sub-dictionary
			* each sub-dictionary contains the following elements:
				* 'ranges' (or range_keyword): array of 2-element quantities indicating spectral measurement windows
				* 'sample' (or sample_keyword): method of sampling flux withn spectral windows
				* 'method' (or method_keyword): method of combining spectral measurements

	range_keyword : str, default = 'ranges'
		defines the dictionary name assigned to spectral measurement windows

	sample_keyword : str, default = 'sample'
		defines the dictionary name assigned to the method of sampling flux

	method_keyword : str, default = 'method'
		defines the dictionary name assigned to the method of combining spectral measurements

	info : bool, default = False
		If True, report out the pre-defined index sets
		alternate names: information, options, option, available

	verbose : bool, default = False
		if True, give extra feedback in computation process  

	Returns
	-------
	dict of index measurements
		keys correspond to the names of the indices, and point to 2-element tuple of (measurement, uncertainty) 

	Examples
	--------

	TBD

	See Also
	--------

	measureIndex : measures a single spectral index

	""" 	

# keyword parameters
	for k in ['reference','set']: ref = kwargs.get(k,ref)
	for k in ['indices','index','index_dict']: index_info = kwargs.get(k,index_info)
	for k in ['option','options','information','available']: info = kwargs.get(k,info)

# just return information on available sets
	if info==True:
		for k in list(INDEX_SETS.keys()):
			inds = list(INDEX_SETS[k]['indices'].keys())
			bibref = ''
			if 'bibcode' in list(INDEX_SETS[k].keys()): bibref = INDEX_SETS[k]['bibcode']
			r = inds[0]
			for i in inds[1:]: r=r+','+i
			print('{} (bibcode: {}): indices: {}'.format(k,bibref,r))
		return

# if index information is not passed, then check with INDEX_SET
	if len(index_info) == 0:
		tmp = checkDict(ref,INDEX_SETS)
		if tmp==False: raise ValueError('Index set {} is not one of the predefined index sets: {}'.format(ref,list(INDEX_SETS.keys())))
		index_info = copy.deepcopy(INDEX_SETS[tmp]['indices'])
		if verbose: print('Measuring indices from {} (bibcode: {})'.format(tmp,INDEX_SETS[tmp]['bibcode']))

# check that the relevant keywords are in the index dictionary
	names = list(index_info.keys())
	for n in names:
		for k in [range_keyword,method_keyword,sample_keyword]:
			if k not in list(index_info[n].keys()):
				raise ValueError('Keyword {} must be present for each defined index in input index dictionary, but is missing for index {}'.format(k,n))

# measure indices
	result = {}
	for n in names:
		ind,err = measureIndex(sp,index_info[n][range_keyword],method=index_info[n][method_keyword],sample=index_info[n][sample_keyword],verbose=verbose,**kwargs)
		result[n] = (ind,err)
		if verbose == True: print('Index {}: {:.3f}+/-{:.3f}'.format(n,ind,err))

	return result



def classifyIndices(spec,set='',verbose=ERROR_CHECKING):
	'''
	Classification program using spectral indices
	TBD
	'''
	pass






def measureEW(sp,lc,width=0.,recenter=True,absorption=True,continuum=0.,continuum_fit_order=1,npixline=2,output_unit=u.Angstrom,output='ew',plot=False,plot_file='',label='',nsamples=100,verbose=ERROR_CHECKING):
	'''
	Program to measure spectral equivalent width
	TBD
	'''

	nloop = 5
# input checks
	if isinstance(sp,Spectrum)==False: raise ValueError('First input should be a Spectrum object not type {}'.format(type(sp)))
	if isUnit(lc): line_center = lc.to(sp.wave.unit).value
	else: line_center = copy.deepcopy(lc)
	if not isinstance(line_center,int) and not isinstance(line_center,float):
		raise ValueError('Second input value should be a single float number; you entered {}'.format(lc))		
	if numpy.nanmin(sp.wave.value) > line_center or numpy.nanmax(sp.wave.value) < line_center:
		print('Wavelength {} is outside spectral data limits of {} to {}'.format(line_center*sp.wave.unit,numpy.nanmin(sp.wave),numpy.nanmax(sp.wave)))
		if output=='all': return {}
		else: return numpy.nan,numpy.nan

	if isUnit(width): line_width = width.to(sp.wave.unit).value
	else: line_width = copy.deepcopy(width)
	if not isinstance(line_width,int) and not isinstance(line_width,float):
		raise ValueError('Line width value should be a single float number; you entered {}'.format(width))		
# estimate line width if not provided
	if float(line_width) == 0.:
		ic = numpy.nanargmin(numpy.array([numpy.abs(a-line_center) for a in sp.wave.value]))
		line_width = numpy.absolute(sp.wave.value[ic+int(npixline)]-sp.wave.value[ic-int(npixline)])

# set up continuum		
	if isUnit(continuum): cont = continuum.to(sp.wave.unit).value
	else: cont = copy.deepcopy(continuum)
	if isinstance(cont,int) or isinstance(cont,float): cont = [line_width,line_width+numpy.abs(cont)]
	if cont[0] == cont[1]: cont = [line_width,line_width*2.]
	if numpy.nanmax(cont) < line_center: cont = [c+line_center for c in cont]
	if len(cont) < 4: cont = [2*line_center-cont[-1],2*line_center-cont[-2],cont[-1],cont[-2]]
	if verbose==True: print('Line center = {}, Line width = {}, Continuum = {}'.format(line_center,line_width,cont))

# preset fail condition
	ew = numpy.nan
	ew_unc = numpy.nan
	line_center_measure = numpy.nan
	line_center_measure_unc = numpy.nan
	rv = numpy.nan
	rv_unc = numpy.nan
	
# refine line centering
	line_center_measure = line_center
	if recenter == True:	
		for i in range(nloop):
			wc = numpy.where(numpy.logical_and(sp.wave.value >= line_center_measure-line_width,sp.wave.value <= line_center_measure+line_width))
			wv = numpy.array(sp.wave.value[wc])
			fl = numpy.array(sp.flux.value[wc])
			if absorption == True: line_center_measure = wv[numpy.nanargmin(fl)]
			else: line_center_measure = wv[numpy.nanargmax(fl)]
	rv = ((line_center_measure-line_center)/line_center)*const.c.to(u.km/u.s)
	cont = [c+line_center_measure-line_center for c in cont]

# sample range
	samplerng = [numpy.nanmin(cont)-0.1*(numpy.nanmax(cont)-numpy.nanmin(cont)),numpy.nanmax(cont)+0.1*(numpy.nanmax(cont)-numpy.nanmin(cont))]
	if numpy.nanmax(sp.wave.value) < samplerng[0] or numpy.nanmin(sp.wave.value) > samplerng[1]:
		print('Sample range {} to {} is outside the wavelength range of the spectrum {} to {}'.format(samplerng[0],samplerng[1],numpy.nanmin(sp.wave.value)*sp.wave.unit,numpy.nanmin(sp.wave.value)*sp.wave.unit))
		return numpy.nan, numpy.nan
		if output=='all': return {}
		else: return numpy.nan,numpy.nan

	w = numpy.where(numpy.logical_and(sp.wave.value >= samplerng[0],sp.wave.value <= samplerng[1]))
	if len(w[0]) > 0:
		f = interp1d(sp.wave.value[w],sp.flux.value[w],bounds_error=False,fill_value=0.)
		wline = numpy.linspace(line_center_measure-line_width,line_center_measure+line_width,nsamples)
		wcont = numpy.append(numpy.linspace(cont[0],cont[1],nsamples),numpy.linspace(cont[-2],cont[-1],nsamples))		
		fline = f(wline)
		fcont = f(wcont)
		pcont = numpy.poly1d(numpy.polyfit(wcont,fcont,continuum_fit_order))
		fcontfit = pcont(wline)
		# print(wline,fline)
		# print(wcont,fcont)
		ew = (trapz((numpy.ones(len(wline))-(fline/fcontfit)), wline)*sp.wave.unit).to(output_unit)

# plot this
		if plot == True:
			plt.clf()
			plt.plot(sp.wave.value,sp.flux.value,'k-')
			plt.fill_between(sp.wave.value,sp.flux.value-sp.unc.value,sp.flux.value+sp.unc.value,color='grey',alpha=0.2)
			plt.plot(wline,fline,'m-')
			plt.plot(wcont[wcont<line_center_measure-line_width],fcont[wcont<line_center_measure-line_width],'b-')
			plt.plot(wcont[wcont>line_center_measure+line_width],fcont[wcont>line_center_measure+line_width],'b-')
			plt.plot(wline,fcontfit,'b-')
			plt.plot([line_center_measure,line_center_measure],[numpy.nanmin(sp.flux.value[w]),numpy.nanmax(sp.flux.value[w])],'k--')
			plt.xlim(samplerng)
			plt.ylim([numpy.nanmin(sp.flux.value[w]),numpy.nanmax(sp.flux.value[w])])
			plt.title(label+' {:.4f}'.format(line_center_measure))
			plt.xlabel('Wavelength ({})'.format(sp.wave.unit))
			plt.ylabel('Flux Density ({})'.format(sp.flux.unit))
			if plot_file != '': plt.savefig(plot_file)


# MC for errors			
		if numpy.median(sp.unc.value[w]) != numpy.nan and numpy.median(sp.unc.value[w]) != 0.:
# default uncertainty = 2 x median uncertainty --> NOT USING
			ew_unc = (2.*numpy.nanmedian(sp.unc.value[w])/pcont(line_center_measure)*(numpy.nanmax(wline)-numpy.nanmin(wline))*sp.wave.unit).to(output_unit)
			ews,lns,rvs = [],[],[]
			for i in range(nsamples):
				wvv = sp.wave.value[w]
				flxv = numpy.random.normal(sp.flux.value[w],sp.unc.value[w])
#					spvar.flux = numpy.random.normal(sp.flux.value,sp.noise.value)*sp.flux.unit
				line_center_measure_var = line_center
				if recenter == True:	
					for i in range(5):
						wc = numpy.where(numpy.logical_and(wvv >= line_center_measure_var-line_width,wvv <= line_center_measure_var+line_width))
						wc = numpy.where(numpy.logical_and(wvv >= line_center_measure_var-line_width,wvv <= line_center_measure_var+line_width))
						wv = numpy.array(wvv[wc])
						fl = numpy.array(flxv[wc])
						if len(fl) > 0:
							if absorption == True: line_center_measure_var = wv[numpy.nanargmin(fl)]
							else: line_center_measure_var = wv[numpy.nanargmax(fl)]
				lns.append(line_center_measure_var)
				rvs.append((((line_center_measure_var-line_center)/line_center)*const.c.to(u.km/u.s)).value)

				f = interp1d(wvv,flxv,bounds_error=False,fill_value=0.)
				wline = numpy.linspace(line_center_measure_var-line_width,line_center_measure_var+line_width,nsamples)
				wcont = numpy.append(numpy.linspace(cont[0],cont[1],nsamples),numpy.linspace(cont[-2],cont[-1],nsamples))		
				fline = f(wline)
				fcont = f(wcont)
				pcont = numpy.poly1d(numpy.polyfit(wcont,fcont,continuum_fit_order))
				fcontfit = pcont(wline)
				ews.append(((trapz((numpy.ones(len(wline))-(fline/fcontfit)), wline)*sp.wave.unit).to(output_unit)).value)
			rv_unc = numpy.std(rvs)*u.km/u.s
#				ew_unc = numpy.sqrt(ew_unc**2+(numpy.std(ews)*output_unit)**2)
			ew_unc = (numpy.std(ews)*output_unit)
			line_center_measure_unc = numpy.std(lns)
							
	if output=='ew': return ew, ew_unc
	elif output=='rv': return rv, rv_unc
	elif output=='center': return line_center_measure*sp.wave.unit, line_center_measure_unc*sp.wave.unit
	else:
		return {'ew': ew, 
			'ew_unc': ew_unc,
			'line_center': line_center_measure*sp.wave.unit,
			'line_center_unc': line_center_measure_unc*sp.wave.unit,
			'rv': rv,
			'rv_unc': rv_unc}



def measureEWElement(sp,elem,verbose=ERROR_CHECKING,**kwargs):
	'''
	Program to measure spectral equivalent widths from a pre-defined set
	TBD
	'''
# input checks
	if isinstance(sp,Spectrum)==False: raise ValueError('First input should be a Spectrum object not type {}'.format(type(sp)))
	tmp = checkDict(elem,EW_LINES)
	if tmp==False: raise ValueError('Element {} is not one of the predefined line sets: {}'.format(elem,list(EW_LINES.keys())))
	lines = EW_LINES[tmp]['lines']
	if verbose: print('Measuring EWs of {} lines at {}'.format(tmp,lines))

# run through lines
	output = {}
	for l in lines:
		res = measureEW(sp,l,**kwargs)
		output['{}-{:.0f}'.format(tmp.replace(' ',''),l.value)] = res
	return output


def zeta(inp,ref='lepine2007',nsamples=100,noiseFlag=False,cah2_name='CaH2',cah3_name='CaH3',tio5_name='TiO5',verbose=ERROR_CHECKING,**kwargs):
	"""
	Measures the zeta paramter based on CaH2, CaH3, and TiO5 indices

	Parameters
	----------
	inp : Spectrum class object OR dict
		Can be either a single spectrum class object for which indices are measured, or a dictionary 
		containing index measurements, with same structure as output of measureIndexSet; i.e.,
		a series of index name strings, each pointing to 2-element array containing measurement and uncertainty
		dictionary keys should include the names assigned to cah2_name, cah3_name, and tio5_name should be presend

	ref : string, default = 'lepine2007'
		named index set contained in the global variable INDEX_SETS
		alternate names: reference, set


	cah2_name : str, default = 'CaH2'
		defines the dictionary name assigned to CaH2 index measurement

	cah3_name : str, default = 'CaH3'
		defines the dictionary name assigned to CaH3 index measurement

	tio5_name : str, default = 'TiO5'
		defines the dictionary name assigned to TiO5 index measurement

	verbose : bool, default = False
		if True, give extra feedback in computation process  

	nsamples : int, default = 100
		Number of Monte Carlo samples for uncertainty determination

	noiseFlag : bool, default = False
		if True, do not determine uncertainty through Monte Carlo sampling  

	Returns
	-------
	tuple of floats
		zeta value and its uncertainty

	Examples
	--------

	TBD

	See Also
	--------

	measureIndexSet : measures a set of spectral indices

	""" 	

# keyword parameters
	for k in ['reference','set']: ref = kwargs.get(k,ref)

# check ref is in ZETA_RELATIONS
	tmp = checkDict(ref,ZETA_RELATIONS)
	if tmp==False: raise ValueError('Reference {} is not one of the predefined calibrations: {}'.format(ref,list(ZETA_RELATIONS.keys())))
	coeff = ZETA_RELATIONS[tmp]['coeff']
	if verbose: print('Measuring zeta using the {} relation (bibcode: {})'.format(tmp,ZETA_RELATIONS[tmp]['bibcode']))

# measure indices if this is a spectrum
	if isinstance(inp,Spectrum):
		indices = measureIndexSet(inp,ref='gizis1997',noiseFlag=noiseFlag,verbose=verbose)
	elif isinstance(inp,dict):
		indices = copy.deepcopy(inp)
	else: raise ValueError('Input parameter should be a Spectrum object of dicitionary of index measurements; you passed a {}'.format(type(inp)))

# make sure relevant indices are present and check for uncertainties
	for k in [cah2_name,cah3_name,tio5_name]:
		if k not in list(indices.keys()): raise ValueError('Index {} is not among input indices: {}'.format(k,list(indices.keys())))
		if numpy.isfinite(indices[k][1])==False: noiseFlag=False

# compute zeta
	tio_ref = numpy.polyval(coeff,indices[cah2_name][0]+indices[cah3_name][0])
	zeta = (1.-indices[tio5_name][0])/(1.-tio_ref)
	e_zeta = numpy.nan
	if noiseFlag == False:
		tio_ref_samp = numpy.polyval(coeff,numpy.random.normal(indices[cah2_name][0],indices[cah2_name][1],nsamples)+numpy.random.normal(indices[cah3_name][0],indices[cah3_name][1],nsamples))
		zeta_samp = (numpy.ones(nsamples)-numpy.random.normal(indices[tio5_name][0],indices[tio5_name][1],nsamples))/(numpy.ones(nsamples)-tio_ref_samp)
		e_zeta = numpy.nanstd(zeta_samp)

	return zeta, e_zeta


def chiFactor(inp,ew=0.,e_ew=0.,ref='schmidt2014',output='loglhalbol',nsamples=100,noiseFlag=False,verbose=ERROR_CHECKING,**kwargs):
	"""
	Determines LHa/Lbol using equivalent width and chi factor

	Parameters
	----------

	Returns
	-------
	tuple of floats
		Lha/Lbol value and its uncertainty

	Examples
	--------

	TBD

	See Also
	--------

	measureIndexSet : measures a set of spectral indices

	""" 	

# keyword parameters
	for k in ['reference','set']: ref = kwargs.get(k,ref)


# measure everything if this is a spectrum
	if isinstance(inp,Spectrum):
		spt = typeToNum(classifyTemplate(inp))
		ew,e_ew = measureEW(inp,6563.,recenter=True,absorption=False)
	elif isinstance(inp,str):
		spt = typeToNum(inp)
	elif isinstance(inp,int) or isinstance(inp,float):
		spt = copy.deepcopy(inp)
	else: raise ValueError('Input parameter should be a Spectrum object or spectral type; you passed a {}'.format(type(inp)))

# check ew & ew_unc
	if ew == 0. or numpy.isfinite(ew)==False: output = 'chi'
	if ew < 0.: ew = numpy.absolute(ew)

# compute chi
# check ref is in CHI_RELATIONS
	tmp = checkDict(ref,CHI_RELATIONS)
	if tmp==False: raise ValueError('Reference {} is not one of the predefined calibrations: {}'.format(ref,list(CHI_RELATIONS.keys())))
	if verbose: print('Determining chi using the {} relation (bibcode: {})'.format(tmp,CHI_RELATIONS[tmp]['bibcode']))
# assuming interpolation; will need to generalize for fits
	chi_spt = CHI_RELATIONS[tmp]['spt']
	chi_val = CHI_RELATIONS[tmp]['values']
	chi_unc = CHI_RELATIONS[tmp]['scatter']
	chi_scale = CHI_RELATIONS[tmp]['scale']
	if spt-typeToNum('K0') < numpy.nanmin(chi_spt) or spt-typeToNum('K0') > numpy.nanmax(chi_spt): 
		print('Spectral type {} is outside range for {} relation: {} to {}'.format(typeToNum(spt),tmp,typeToNum(numpy.nanmin(chi_spt)+typeToNum('K0')),typeToNum(numpy.nanmax(chi_spt)+typeToNum('K0'))))
		return numpy.nan, numpy.nan
	f = interp1d(chi_spt,chi_val)
	chi = f(spt-typeToNum('K0'))*chi_scale
	f = interp1d(chi_spt,chi_unc)
	e_chi = f(spt-typeToNum('K0'))*chi_scale
	if verbose: print('chi = {}+/-{}'.format(chi,e_chi))
	if output=='chi': return chi, e_chi

# compute LHa/Lbol
	if isUnit(ew): ew = ew.to(u.Angstrom).value
	if isUnit(e_ew): e_ew = e_ew.to(u.Angstrom).value
	lhalbol = chi*ew
	e_lhalbol = numpy.nan
	if e_ew<= 0.or e_chi<=0 or numpy.isfinite(e_ew)==False or numpy.isfinite(e_chi)==False: noiseFlag=True
	if noiseFlag==False:
		samp = numpy.random.normal(chi,e_chi,nsamples)*numpy.random.normal(ew,e_ew,nsamples)
		e_lhalbol = numpy.nanstd(samp)

	if verbose: print('LHa/Lbol = {}+/-{}'.format(lhalbol,e_lhalbol))
	if output=='lhalbol': return lhalbol, e_lhalbol
	else: return numpy.log10(lhalbol), (e_lhalbol/lhalbol)/numpy.log(10)





