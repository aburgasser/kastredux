# -*- coding: utf-8 -*-
from __future__ import print_function

# WORKING COPY OF KAST REDUCTION CODE
# THINGS THAT NEED FIXING
# - documentation
# - function docstrings
# - if output dictionary exists, won't do anything! -> can't used saved files
# - add more dispersion options (esp. blue)
# - error on writing fits files - INSTRUMENT keyword (shorten to INSTRUM
#

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
VERSION = '2020.02.17'
__version__ = VERSION

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
# KAST INSTRUMENT CONSTANTS AND FUNCTIONS
############################################################

DEFAULT_WAVE_UNIT = u.Angstrom
DEFAULT_FLUX_UNIT = u.erg/u.cm/u.cm/u.Angstrom/u.s

MODES = {
	'RED': {'PREFIX': 'r', 'ALTNAME': ['r','rd','long','ir','nir'], 'NAME': 'KAST red', 'VERSION': 'kastr'},
	'BLUE': {'PREFIX': 'b', 'ALTNAME': ['b','bl','short','uv','vis'], 'NAME': 'KAST blue', 'VERSION': 'kastb'},
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
}
# assumes slow read
CCD_HEADER_KEYWORDS = {
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
CCD_PARAMETERS = {
	'RED-FAST': {'GAIN': 0.55, 'RN': 4.3},
	'RED-SLOW': {'GAIN': 1.9, 'RN': 3.7},
	'BLUE-FAST': {'GAIN': 1.3, 'RN': 6.5},
	'BLUE-SLOW': {'GAIN': 1.2, 'RN': 3.8},
}

# extract image mode from image header
def kastRBMode(hdr,keyword='VERSION'):
	if keyword not in list(hdr.keys()):
		raise ValueError('Header does not contain red/blue mode keyword {}'.format(keyword))
	md = hdr[keyword].strip()
	mode = ''
	for r in list(MODES.keys()):
		if MODES[r]['VERSION'] == md: mode=r
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
	elif keyword in list(CCD_HEADER_KEYWORDS.keys()): return hdr[CCD_HEADER_KEYWORDS[keyword]]
	elif keyword.upper() in list(CCD_HEADER_KEYWORDS.keys()): return hdr[CCD_HEADER_KEYWORDS[keyword.upper()]]
	else:
		raise ValueError('Cannot find keyword {} in header or header lookup table'.format(keyword))



# RESOURCE DATA
FLUXCALFOLDER = CODE_PATH+'/resources/flux_standards/'
FLUXCALS = {
	'FEIGE110': {'FILE' : 'ffeige110.dat'},
	'HILTNER600': {'FILE' : 'fhilt600.dat'},
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
ERROR_CHECKING = True


# SPECTRUM CLASS
class Spectrum(object):
	'''
	Container for spectrum object
	Includes wavelength, flux, variance, background, mask vectors; trace; and header
	Includes methods for combining spectrum objects together, reading/writing, conversion
	'''
	def __init__(self,**kwargs):
		core_attributes = {'instrument': 'KAST red','name': 'Unknown source','wave': [],'flux': [],'unc': [],'variance':[],'background': [],'mask': [],'header': {},}
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
		self.history = ['{} spectrum of {} successfully loaded'.format(self.instrument,self.name)]
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
		return '{} spectrum of {}'.format(self.instrument,self.name)

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
			out.variance = numpy.multiply(out.flux.value**2,((numpy.divide(flxs(out.wave.value),fself(out.wave.value))**2)+(numpy.divide(flxo(out.wave.value),fother(out.wave.value))**2)))*self.variance.unit*other.variance.unit
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
			setattr(out,k,numpy.divide(fself(out.wave.value),fother(out.wave.value))*(getattr(self,k).unit)/(getattr(other,k).unit))
			if k=='flux': flxs,flxo = fself,fother
# special for variance
		fself = interp1d(self.wave.value,self.variance.value,bounds_error=False,fill_value=0.)
		fother = interp1d(other.wave.value,other.variance.value,bounds_error=False,fill_value=0.)
		if numpy.random.choice(numpy.isfinite(self.variance.value))==True and numpy.random.choice(numpy.isfinite(other.variance.value))==True: 
			out.variance = numpy.multiply(out.flux.value**2,((numpy.divide(flxs(out.wave.value),fself(out.wave.value))**2)+(numpy.divide(flxo(out.wave.value),fother(out.wave.value))**2)))*self.variance.unit/other.variance.unit
		elif numpy.random.choice(numpy.isfinite(self.variance.value))==True:
			out.variance = (numpy.divide(fself(out.wave.value),flxo(out.wave.value))**2)*self.variance.unit/(other.flux.unit**2)
		elif numpy.random.choice(numpy.isfinite(other.variance.value))==True:
			out.variance = (numpy.divide(fother(out.wave.value),flxs(out.wave.value))**2)*self.variance.unit/(other.flux.unit**2)
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
		Do a boxcar smooth
		'''
		xsamp = numpy.arange(0,len(self.wave)-scale,scale)
		self.wave = numpy.array([self.wave.value[x+int(0.5*scale)] for x in xsamp])*self.wave.unit
		for k in ['flux','unc','background']:
			setattr(self,k,numpy.array([numpy.nanmedian(getattr(self,k).value[x:x+scale]) for x in xsamp])*getattr(self,k).unit)
		self.variance = self.unc**2
		self.history.append('Smoothed by a scale of {} pixels'.format(scale))
	
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

	def applyFluxCal(self,cal_flux):
		'''
		Apply flux calibration
		'''
# check for relevant parameters
		required = ['COEFF','UNIT']
		for r in required:
			if r not in list(cal_flux.keys()): raise ValueError('Required parameter {} not in flux calibration structure'.format(r))

		for k in ['flux','unc','background']:
			if k in list(self.__dict__.keys()):
				setattr(self,k,(getattr(self,k).value)*(10.**(numpy.polyval(cal_flux['COEFF'],self.wave.value)))*cal_flux['UNIT'])
		self.variance = self.unc**2
		self.history.append('Flux calibration applied')
		if 'NAME' in list(cal_flux.keys()): self.history.append('Flux calibrator {} used'.format(cal_flux['NAME']))
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
				if k.upper() not in ['HISTORY','COMMENT','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND'] and k.replace('#','') != '': # and k not in list(hdu.header.keys()):
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



# program to generate a log sheet
def makeLog(folder,output='',mode='RED',verbose=ERROR_CHECKING):
	'''
	:Purpose:

		Reads in the headers of files in a data folder and generates a log in pandas format

	:Required Inputs: 

		:param folder: folder for data

	:Optional Inputs: 

		:param output='': output file for log; can be .xlsx, .xls, .txt, .csv, .htm, .html, or .tex
		:param mode='RED': data mode (RED or BLUE)

	:Output: 

		If one file number provided, a 2D numpy array containing image and dictionary containing header
		If more than one file number provided, a 3D numpy array of [n files x n rows x n columns] and array of header dictionaries

	:Usage: 

		TBD

	'''
	gzflag = False
	hdkeys = ['OBSNUM','OBJECT','OBSTYPE','EXPTIME','DATE-OBS','RA','DEC','HA','AIRMASS','GRISM_N','GRATNG_N','SLIT_N','OBSERVER','DATASEC']
# check folder & mode	
	if os.path.exists(folder) == False: raise ValueError('Could not find data folder {}'.format(folder))

# get list of files	
	for mode in list(MODES.keys()):
		if checkDict(mode,MODES) == False: raise ValueError('Mode {} not defined'.format(mode))
		else: mode=checkDict(mode,MODES)
		files = glob.glob(folder+MODES[mode]['PREFIX']+'*.fits')
		if len(files)==0: 
			files = glob.glob(folder+MODES[mode]['PREFIX']+'*.gz')
			gzflag = True
		if len(files)==0: raise ValueError('Could not find data files {}*.fits or {}*.gz in {}'.format(MODES[mode]['PREFIX'],MODES[mode]['PREFIX'],folder))

# build up header information for files
		hd = {}
		hd['FILE'] = [os.path.basename(f) for f in files]
		for k in hdkeys: hd[k] = []
# NOTE: may need to utilize gzflag here
		for f in files:
			hdu = fits.open(f)
			for k in hdkeys: hd[k].append(hdu[0].header[k])
			hdu.close()
		dp = pandas.DataFrame(data=hd)
		dp.sort_values('FILE',inplace=True)
		dp.reset_index(inplace=True,drop=True)
		if output!='':
			ospl = output.split('.')
			outfile = output.replace('.{}'.format(ospl[-1]),'_{}.{}'.format(mode,ospl[-1]))
			if ospl[-1] == 'xls' or ospl[-1] == 'xlsx': dp.to_excel(outfile,index=False)
			elif ospl[-1] == 'csv': dp.to_csv(outfile,index=False,sep=',')
			elif ospl[-1] == 'txt': dp.to_csv(outfile,index=False,sep='\t')
			elif ospl[-1] == 'tex': dp.to_tex(outfile,index=False)
			elif ospl[-1] == 'html' or ospl[-1] == 'htm': dp.to_html(outfile,index=False)
			else: print('Cannot output to file {}; skipping'.format(outfile))
	return dp
		

def readKastFiles(num,folder='./',mode='',prefix='',rotate=True,verbose=ERROR_CHECKING):
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
	if len(files) == 0:
		if prefix=='': 
			if mode=='': raise ValueError('Mode (red/blue) or file prefix must be provided')
			if checkDict(mode,MODES) == False: raise ValueError('Mode {} not defined'.format(mode))
			else: mode=checkDict(mode,MODES)	
			prefix=MODES[mode]['PREFIX']
		files = ['{}/{}{}.fits'.format(folder,prefix,str(n).zfill(4)) for n in nlist]
# read in files
	for f in files:
		if not os.path.exists(f): print('Warning: cannot find data file {}; skipping'.format(f))
		else:
			hdulist = fits.open(f)
			im = hdulist[0].data
			hdr = hdulist[0].header
			hdulist.close()
			if mode=='':
				try: mode = kastRBMode(hdr)
				except: mode='UNKNOWN'
			if mode=='RED' and rotate==True: im = numpy.rot90(im,k=1)
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
		with images read in using readKastFiles()

	:Optional Inputs: 

		:param method='median': method for combining images (more than one of these can be included)
			* 'median': median combine
			* 'add': add all pixels
			* 'average: average value of all pixels
			* 'sigclip': sigma clipping, rejecting outliers N x std deviation, where N is specified by input ``sclip''
		:param axis=0.: axis upon which images are combined
		:param sclip=5.: number of standard deviations to apply signma clipping

		If reading in files, additional parameters for readKastFiles() can be included though **kwargs

	:Output: 

		Combined 2D image

	:Usage: 

		TBD

	'''
	images = copy.deepcopy(imarr)
	if isinstance(imarr,str) or len(numpy.shape(imarr))==1:
		try: 
			images,hds = readKastFiles(imarr,**kwargs)
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
				if checkDict(parameters[ref],MODES) == False: raise ValueError('Mode {} not defined'.format(parameters[ref]))
				else: parameters[ref]=checkDict(parameters[ref],MODES)
			elif ref in ['BIAS','FLAT','ARC','ARC_SHALLOW','ARC_DEEP','FLUXCAL']:
				dt = {}
				for p in parts[1:]:
					kv = p.split('=')
					dt[kv[0].upper().strip()] = kv[1].upper().strip()
				if 'FILES' not in list(dt.keys()): print('Warning: FILES keyword must be included in instruction file line for {}'.format(ref))
				else: dt['FILES'] = numberList(dt['FILES'])
				if 'BACK' in list(dt.keys()): 
					dt['BACK'] = dt['BACK'].replace('[','').replace(']','').replace('(','').replace(')','').split(',')
					dt['BACK'] = [int(x) for x in dt['BACK']]
				if 'WINDOW' in list(dt.keys()): dt['WINDOW'] = int(dt['WINDOW'])
				parameters[ref] = dt
			elif ref in ['SOURCE','TELLURIC']:
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



def makeBias(files,method='median',folder='./',mode='',prefix='',overwrite=True,verbose=ERROR_CHECKING,output=''):
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
	images,hds = readKastFiles(files,folder=folder,mode=mode,prefix=prefix)
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


def makeFlat(files,bias,method='median',quantile=0.9,folder='./',mode='',prefix='',verbose=ERROR_CHECKING,overwrite=True,output=''):
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
	images,hds = readKastFiles(files,folder=folder,mode=mode,prefix=prefix)
	if len(numpy.shape(images))==2:
		if verbose==True: print('Warning: image array is just a single image of dimensions {}; returning'.format(numpy.shape(images)))
		return images,hds

# read in bias if necessary
	if isinstance(bias,str): 
		if os.path.exists(bias): bias,bhd = readKastFiles(bias)
		elif os.path.exists(folder+'/'+bias): bias,bhd = readKastFiles(folder+'/'+bias)
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


def makeMask(bias,flat,sclip=3.,mode='',verbose=ERROR_CHECKING,overwrite=True,output=''):
	'''
	Creates a mask file: 0 = OK, 1 = BAD
	Input: flat, bias
	Output: mask image
	'''

# read in bias and flat if necessary
	if isinstance(bias,str): 
		if os.path.exists(bias): bias,bhd = readKastFiles(bias)
		elif os.path.exists(folder+'/'+bias): bias,bhd = readKastFiles(folder+'/'+bias)
		else: raise ValueError('Cannot find bias file {}'.format(bias))
		if mode=='': mode = kastRBMode(bhd)
	if isinstance(flat,str): 
		if os.path.exists(flat): flat,fhd = readKastFiles(flat)
		elif os.path.exists(folder+'/'+flat): flat,fhd = readKastFiles(folder+'/'+flat)
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
	mask[flat<(0.01*numpy.quantile(flat[mask==0],0.9))] = 1
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
		im,hd = readKastFiles(image,folder=folder,mode=mode)
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
		if os.path.exists(bias): bias,bhd = readKastFiles(bias)
		elif os.path.exists(folder+'/'+bias): bias,bhd = readKastFiles(folder+'/'+bias)
		else: raise ValueError('Cannot find bias file {}'.format(bias))
		if kastRBMode(bhd) != kastRBMode(hd): raise ValueError('Red/blue mode of bias frame is {} while that of science frames is {}'.format(kastRBMode(bhd),kastRBMode(hd)))
	if isinstance(flat,str): 
		if os.path.exists(flat): flat,fhd = readKastFiles(flat)
		elif os.path.exists(folder+'/'+flat): flat,fhd = readKastFiles(folder+'/'+flat)
		else: raise ValueError('Cannot find flat file {}'.format(flat))
		if kastRBMode(fhd) != kastRBMode(hd): raise ValueError('Red/blue mode of flat frame is {} while that of science frames is {}'.format(kastRBMode(bhd),kastRBMode(hd)))

# process image
	image_b = im-bias
	image_e = image_b*gain
# mask flat before dividing
	if len(mask)==0: mask = image_e*0
	fm = maskClean(flat,mask,replace=1.)
	image_f = image_e/(fm*exposure)
	variance = (image_e/fm+rn**2)/(exposure**2)
# mask image
	if mask_image==True:
		image_f = maskClean(image_f,mask,replace=numpy.nanmedian(image_f))	
		variance = maskClean(variance,mask,replace=numpy.nanmax(variance))	
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
	if cntr<0.: cntr=findPeak(im,verbose=verbose,**center_kwargs)   
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
	if method=='arc':
		subtract_background = False
	if method=='optimal' and len(profile)==0:
		profile = spatialProfile(im,cntr=cntr,window=src_wnd)
	if method=='boxcar' or len(profile)==0:
		profile = numpy.ones(2*src_wnd+1)
# force agreement between source window and extraction profile so profile can be used as input
	if len(profile)>(2*src_wnd+1): 
		profile = profile[int(0.5*(len(profile)-1)-src_wnd):int(0.5*(len(profile)-1)+src_wnd+1)]
	if len(profile)<(2*src_wnd+1): 
		src_wnd = int(0.5*(len(profile)-1))

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
		slc = slc-numpy.nanmin(slc)
		slc = slc/numpy.nanmax(slc)
		# xslc = numpy.arange(len(slc))+cntr-2*bck_wnd[1]
		# xref = int(numpy.nanmedian(x))
		# slc = numpy.nanmedian(im[int(trace[xref]-2*bck_wnd[1]):int(trace[xref]+2*bck_wnd[1]),int(xref-10):int(xref+10)],axis=1)

		xrng = [cntr-2*bck_wnd[1],cntr+2*bck_wnd[1]]
		yrng = [-0.05,1.2]
		plt.plot(numpy.arange(len(slc)),slc,'k-',alpha=0.5)
		plt.plot(numpy.arange(len(profile))+cntr-0.5*(len(profile)-1),profile,'b-')
		plt.legend(['Spatial Profile','Extraction Profile'])
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

		plt.tight_layout()
		try:
			plt.savefig(plot_file)
			plt.close()
		except:
			print('Warning: could not save flux calibration diagnostics to file {}'.format(plot_file))
#		plt.clf()
# generate a spectrum object
	return sp


def waveCalibrateArcs(arcim,deep=[],dispersion='',mode='',trace=[],prior={},fit_order=6,sfit_order=2,verbose=ERROR_CHECKING,middle=True,resolution=0.,lam0=[],pixel0=[],cntr=0.,fitcycle=5,sclip=2.,plot_file=''):
	'''
	Wavelength calibration from arc lamp
	Input: arc, list of lines
	Output: dictionary containing pixel->wave conversion, fit diagnostics
	THIS ASSUMES RED ONLY
	NEED TO INCLUDE HELIOCENTRIC CORRECTION EXPLICITLY
	'''
# check inputs
	if dispersion in list(DISPERSIONS.keys()):
		resolution = DISPERSIONS[dispersion]['RESOLUTION']
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
	elif mode=='BLUE':
		extwidth = 30
		swindow = 30
		dwindow = 25
		deep_threshold = 0.0015
		hehgcd = numpy.array([3261.05,3341.08,3403.65,3466.55,3610.51,3650.15,3662.88,3888.65,4046.56,4077.83,\
				4358.33,4471.50,4678.16,4799.92,4921.93,5015.68,5085.82,5460.74,5769.59,5790.65,5875.62])
		alines = hehgcd
		strong = numpy.array([3466.55,4046.56,4358.33,4678.16,4799.92,5085.82,5460.74,5875.62,5944.83])
		fit_order = 4
	else: raise ValueError('You must specify the red/blue mode (mode=RED or BLUE); mode={} was passed'.format(mode))

	cal_wave = {}

# extract arc traces from arc files
	if len(trace)==0:
		if cntr==0.: cntr = int(0.5*len(arcim[:,0]))
		trace = numpy.zeros(len(arcim[0,:]))+cntr
	tr = [int(t) for t in trace]
	ext = extractSpectrum(arcim,method='arc',subtract_background=False,trace=tr,src_wnd=extwidth,shift_trace=False,verbose=verbose)
	arctrace = ext.flux.value
	arctrace = arctrace-numpy.nanmin(arctrace)
	arctrace = arctrace/numpy.nanmax(arctrace)
	if len(deep)>0:
		ext = extractSpectrum(deep,method='arc',subtract_background=False,trace=tr,src_wnd=extwidth,shift_trace=False,verbose=verbose)
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
	if verbose==True: print('Final deep line pass: RMS={:.3f} Ang at {:.2f} Ang for {} lines and fit order {}, dRV = {:.1f} km/s'.format(rms,numpy.nanmedian(wave),len(lines_y),sfit_order,3.e5*rms/numpy.nanmedian(wave)))

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

def telluricCalibrate(tellsp,fitrange=[6200,8800],fit_order=5,fitcycle=10,sclip=3.,plot_file='',verbose=ERROR_CHECKING):
	'''
	Uses spectrum of telluric calibrator to determine corrections to telluric absorption
	Input: telluric cal spectrum, spectral type
	Output: dictionary containing telluric absorption correction, fit diagnostics
	'''
# regions of strong telluric absorption	
	tranges = [[5800,6000],[6250,6340],[6440,6620],[6850,6960],[7160,7350],[7580,7700],[8100,8350]]
# regions of G/A star line absorption	
	glines = [[3636,3656],[3960,3980],[4092,4112],[4330,4350],[4851,4871],[6553,6573],[8194,8214],[8480,8550],[8648,8668],[9536,9566]]

# mask out telluric and star lines to get fit of continuum
	wv = tellsp.wave.value
	flx = tellsp.flux.value
	mask = numpy.zeros(len(wv))
	mask[numpy.isnan(flx)==True] = 1
	mask[flx<=0] = 1
#	print(len(mask),numpy.nansum(mask))
	for t in tranges: mask[numpy.logical_and(wv>=t[0],wv<=t[1])]=1
#	print(len(mask),numpy.nansum(mask))
	for t in glines: mask[numpy.logical_and(wv>=t[0],wv<=t[1])]=1
#	print(len(mask),numpy.nansum(mask))
	mask[numpy.logical_or(wv<fitrange[0],wv>fitrange[1])]=1
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
	if numpy.nansum(tmask)>0:
		correction[numpy.logical_and(tmask==1,flx>0)] = continuum[numpy.logical_and(tmask==1,flx>0)]/flx[numpy.logical_and(tmask==1,flx>0)]

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
		plt.plot(wv[mask==0],flx[mask==0],'b-',)
		plt.plot(wv[mask==0],continuum[mask==0],'m-')
		plt.legend(['Observed','Masked','Continuum'])
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
	for src in list(parameters['SOURCE'].keys()):
		if verbose==True: print('Checking spatial profile of {}'.format(src))
		im,hd = readKastFiles(parameters['SOURCE'][src]['FILES'][0],folder=parameters['DATA_FOLDER'],mode=parameters['MODE'])
		if 'CENTER' in list(parameters['SOURCE'][src].keys()): cntr = int(parameters['SOURCE'][src]['CENTER'])
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
	if 'BIAS' not in list(redux.keys()) or reset==True:
		if bias_file!='':
			if os.path.exists(bias_file) == False: bias_file=redux['PARAMETERS']['REDUCTION_FOLDER']+'/'+bias_file
			if os.path.exists(bias_file) == False: bias_file=''
		if bias_file!='' and reset==False:
			if verbose==True: print('\nReading in bias frame from file {}'.format(bias_file))
			redux['BIAS'],redux['BIAS_HEADER'] = readKastFiles(bias_file,mode=redux['PARAMETERS']['MODE'],rotate=False)
		else:
			bias_file='{}/bias_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'])
			if verbose==True: print('\nReducing bias frames')
			redux['BIAS'],redux['BIAS_HEADER'] = makeBias(redux['PARAMETERS']['BIAS']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'],output=bias_file)

# read or create flat field frame		
	if 'FLAT' not in list(redux.keys()) or reset==True:
		if flat_file!='':
			if os.path.exists(flat_file) == False: flat_file=redux['PARAMETERS']['REDUCTION_FOLDER']+flat_file
			if os.path.exists(flat_file) == False: flat_file=''
		if flat_file!='' and reset==False:
			if verbose==True: print('\nReading in flat field frame from file {}'.format(flat_file))
			redux['FLAT'],redux['FLAT_HEADER'] = readKastFiles(flat_file,mode=redux['PARAMETERS']['MODE'],rotate=False)
		else:
			flat_file='{}/flat_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'])
			if verbose==True: print('\nReducing flat field frames')
			redux['FLAT'],redux['FLAT_HEADER'] = makeFlat(redux['PARAMETERS']['FLAT']['FILES'],redux['BIAS'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'],output=flat_file)

# read or create mask frame		
	if 'MASK' not in list(redux.keys()) or reset==True:
		if mask_file!='':
			if os.path.exists(mask_file) == False: mask_file=redux['PARAMETERS']['REDUCTION_FOLDER']+mask_file
			if os.path.exists(mask_file) == False: mask_file=''
		if mask_file!='' and reset==False:
			if verbose==True: print('\nReading in mask frame from file {}'.format(mask_file))
			redux['MASK'],redux['MASK_HEADER'] = readKastFiles(mask_file,mode=redux['PARAMETERS']['MODE'],rotate=False)
		else:
			mask_file='{}/mask_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'])
			if verbose==True: print('\nGenerating mask file')
			redux['MASK'] = makeMask(redux['BIAS'],redux['FLAT'],output=mask_file)
			if verbose==True: print('Masking {} bad pixels'.format(int(numpy.nansum(redux['MASK']))))

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
		arcsh,arcshhd = readKastFiles(redux['PARAMETERS']['ARC_SHALLOW']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
		arcdp,arcdphd = readKastFiles(redux['PARAMETERS']['ARC_DEEP']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
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
		if verbose==True: print('\nAnalyzing flux calibrator {}'.format(redux['PARAMETERS']['FLUXCAL']['NAME']))
		ims,hds = readKastFiles(redux['PARAMETERS']['FLUXCAL']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
		if len(redux['PARAMETERS']['FLUXCAL']['FILES']) > 1:
			im,hd = crRejectCombine(ims,verbose=verbose),hds[0]
		else: im,hd = ims,hds
# some parameters
		src_wnd = copy.deepcopy(src_wnd0)
		if 'WINDOW' in list(redux['PARAMETERS']['FLUXCAL'].keys()): src_wnd = redux['PARAMETERS']['FLUXCAL']['WINDOW']
		bck_wnd = copy.deepcopy(bck_wnd0)
		if 'BACK' in list(redux['PARAMETERS']['FLUXCAL'].keys()): bck_wnd = redux['PARAMETERS']['FLUXCAL']['BACK']
		fit_order = copy.deepcopy(fit_order0)
		if 'FIT_ORDER' in list(redux['PARAMETERS']['FLUXCAL'].keys()): fit_order = redux['PARAMETERS']['FLUXCAL']['FIT_ORDER']
# reduce data, determine trace and extract spectrum
#def reduceScienceImage(image, bias, flat, mask=[], hd={}, mode='', rdmode='', gain=0., rn=0., exposure='EXPTIME', mask_image=True, folder='./', verbose=ERROR_CHECKING):
		imr,var = reduceScienceImage(im,redux['BIAS'],redux['FLAT'],redux['MASK'],hd=hd)
		cntr = findPeak(imr)
		trace = traceDispersion(imr,cntr=cntr,window=src_wnd,method='maximum',plot_file='{}/diagnostic_trace_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['FLUXCAL']['NAME'],redux['PARAMETERS']['MODE']))
		imrect = rectify(imr,trace)
		varrect = rectify(var,trace)
		maskrect = rectify(redux['MASK'],trace)
		cntr = findPeak(imrect,plot_file='{}/diagnostic_profile_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['FLUXCAL']['NAME'],redux['PARAMETERS']['MODE']))
#		profile = spatialProfile(imrect,cntr=cntr,window=src_wnd)
		profile = numpy.ones(int(2*src_wnd+1))
		spflx = extractSpectrum(imrect,var=varrect,mask=maskrect,src_wnd=src_wnd,bck_wnd=bck_wnd,profile=profile,plot_file='{}/diagnostic_extraction_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['FLUXCAL']['NAME'],redux['PARAMETERS']['MODE']))
		spflx.name = redux['PARAMETERS']['FLUXCAL']['NAME']
		spflx.header = hd
#		spflx = extractSpectrum(imr,var,mask=redux['MASK'],trace=trace,src_wnd=10)
		arcdp,arcdphd = readKastFiles(redux['PARAMETERS']['ARC_DEEP']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
		arcrect = rectify(arcdp,trace)
		arcrecal = waveCalibrateArcs(arcrect,cntr=cntr,prior=redux['CAL_WAVE'],mode=redux['PARAMETERS']['MODE'],plot_file='{}/diagnostic_wavecal_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['FLUXCAL']['NAME'],redux['PARAMETERS']['MODE']))
		spflx.applyWaveCal(arcrecal)
		if redux['PARAMETERS']['MODE'] == 'BLUE':
			redux['CAL_FLUX'] = fluxCalibrate(spflx,redux['PARAMETERS']['FLUXCAL']['NAME'],fit_order=fit_order,fit_range=[3700.,5300.],plot_file='{}/diagnostic_fluxcal_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']))
		else:
			redux['CAL_FLUX'] = fluxCalibrate(spflx,redux['PARAMETERS']['FLUXCAL']['NAME'],fit_order=fit_order,plot_file='{}/diagnostic_fluxcal_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']))
		redux['CAL_FLUX']['TRACE'] = trace
		redux['CAL_FLUX']['PROFILE'] = profile
		spflx.applyFluxCal(redux['CAL_FLUX'])
		redux[redux['PARAMETERS']['FLUXCAL']['NAME']] = spflx
# save this to pickle file
		f = open('{}/cal_flux_{}.pkl'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE']),'wb')
		pickle.dump(redux['CAL_FLUX'],f)
		f.close()
# plot  and save
		fig = spflx.plot(ylim=[0,1.2*numpy.quantile(spflx.flux.value,0.95)])
		fig.figure.savefig('{}/kast{}_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],redux['PARAMETERS']['FLUXCAL']['NAME'],redux['PARAMETERS']['DATE']))
		plt.close()
#		plt.clf()
		spflx.toFile('{}/kast{}_{}_{}.fits'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],redux['PARAMETERS']['FLUXCAL']['NAME'],redux['PARAMETERS']['DATE']))
		spflx.toFile('{}/kast{}_{}_{}.txt'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],redux['PARAMETERS']['MODE'],redux['PARAMETERS']['FLUXCAL']['NAME'],redux['PARAMETERS']['DATE']))

# determine telluric calibrations
	if ('CAL_TELL' not in list(redux.keys()) or reset==True or kwargs.get('reset_tellcal',False) == True) and 'TELLURIC' in list(redux['PARAMETERS'].keys()):
		if verbose==True: print('\nComputing telluric corrections')
		redux['CAL_TELL'] = {}
		for tstar in list(redux['PARAMETERS']['TELLURIC'].keys()):
			ims,hds = readKastFiles(redux['PARAMETERS']['TELLURIC'][tstar]['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
			if len(redux['PARAMETERS']['TELLURIC'][tstar]['FILES']) > 1: im,hd = crRejectCombine(ims,verbose=verbose),hds[0]
			else: im,hd = ims,hds
# some parameters
			src_wnd = copy.deepcopy(src_wnd0)
			if 'WINDOW' in list(redux['PARAMETERS']['TELLURIC'][tstar].keys()): src_wnd = redux['PARAMETERS']['TELLURIC'][tstar]['WINDOW']
			bck_wnd = copy.deepcopy(bck_wnd0)
			if 'BACK' in list(redux['PARAMETERS']['TELLURIC'][tstar].keys()): bck_wnd = redux['PARAMETERS']['TELLURIC'][tstar]['BACK']
# reduce imaging data
			imr,var = reduceScienceImage(im,redux['BIAS'],redux['FLAT'],redux['MASK'],hd=hd)
#			imr,var = reduceScienceImage(im,redux['BIAS'],redux['FLAT'],redux['MASK'],hd['EXPTIME'])
# trace source and rectify
			cntr = findPeak(imr)
			trace = traceDispersion(imr,cntr=cntr,window=src_wnd,method='maximum',plot_file='{}/diagnostic_trace_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
			imrect = rectify(imr,trace)
			varrect = rectify(var,trace)
			maskrect = rectify(redux['MASK'],trace)
# find center
			cntr = findPeak(imrect,plot_file='{}/diagnostic_profile_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
# make profile
#			profile = spatialProfile(imrect,cntr=cntr,window=src_wnd)
			profile = numpy.ones(int(2*src_wnd+1))
# extract spectrum
			spflx = extractSpectrum(imrect,var=varrect,mask=maskrect,src_wnd=src_wnd,bck_wnd=bck_wnd,profile=profile,plot_file='{}/diagnostic_extraction_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
			spflx.name = tstar
			spflx.header = hd
# reidentify arc lines
			arcdp,arcdphd = readKastFiles(redux['PARAMETERS']['ARC_DEEP']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
			arcrect = rectify(arcdp,trace)
			arcrecal = waveCalibrateArcs(arcrect,trace=trace,prior=redux['CAL_WAVE'],mode=redux['PARAMETERS']['MODE'],plot_file='{}/diagnostic_wavecal_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
# apply wavelength calibration
			spflx.applyWaveCal(arcrecal)
#			spflx.applyWaveCal(redux['CAL_WAVE'])
# apply flux calibration
			spflx.applyFluxCal(redux['CAL_FLUX'])
# compute telluric corection
			redux['CAL_TELL'][tstar] = telluricCalibrate(spflx,plot_file='{}/diagnostic_telluic_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],tstar,redux['PARAMETERS']['MODE']))
			redux['CAL_TELL'][tstar]['NAME'] = tstar
			redux['CAL_TELL'][tstar]['TRACE'] = trace
			redux['CAL_TELL'][tstar]['PROFILE'] = profile
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

# now analyze science targets - ONLY DOING FIRST FILE
	for src in list(redux['PARAMETERS']['SOURCE'].keys()):
		if src not in list(redux.keys()) or reset==True or kwargs.get('reset_source',False) == True:
			if verbose==True: print('\nExtracting {} spectrum of {}'.format(redux['PARAMETERS']['MODE'],src))
#			print(redux['PARAMETERS']['SOURCE'][src])
			ims,hds = readKastFiles(redux['PARAMETERS']['SOURCE'][src]['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'],verbose=verbose,)
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
				tcorr,ttrace = True,True
				tstar = redux['PARAMETERS']['SOURCE'][src]['TELLURIC']
			if 'APPLY_TRACE' in list(redux['PARAMETERS']['SOURCE'][src].keys()) and tcorr==True: 
				ttrace = redux['PARAMETERS']['SOURCE'][src]['APPLY_TRACE'].upper()=='TRUE'
			cflag = True
			if 'RECENTER' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				cflag = redux['PARAMETERS']['SOURCE'][src]['RECENTER'].upper()=='TRUE'
			pflag = 'SOURCE'
			if 'PROFILE' in list(redux['PARAMETERS']['SOURCE'][src].keys()): 
				pflag = redux['PARAMETERS']['SOURCE'][src]['PROFILE'].upper()
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
			if pflag.upper()=='SOURCE': 
				if verbose==True: print('Using spatial profile measured from source')
				profile = spatialProfile(imrect,cntr=cntr,window=src_wnd)
			elif pflag.upper()=='TELLURIC' and tcorr==True: 
				if verbose==True: print('Using spatial profile measured from telluric standard {}'.format(tstar))
				profile = redux['CAL_TELL'][tstar]['PROFILE']
			else: 
				if verbose==True: print('Using flat spatial profile')
				profile = numpy.ones(int(2*src_wnd+1))
			spflx = extractSpectrum(imrect,var=varrect,mask=maskrect,cntr=cntr,src_wnd=src_wnd,bck_wnd=bck_wnd,profile=profile,plot_file='{}/diagnostic_extraction_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],src,redux['PARAMETERS']['MODE']))
			spflx.name = src
			spflx.header = hd
# reapply arc solution
			arcdp,arcdphd = readKastFiles(redux['PARAMETERS']['ARC_DEEP']['FILES'],folder=redux['PARAMETERS']['DATA_FOLDER'],mode=redux['PARAMETERS']['MODE'])
			arcrect = rectify(arcdp,trace)
			arcrecal = waveCalibrateArcs(arcrect,trace=trace,prior=redux['CAL_WAVE'],mode=redux['PARAMETERS']['MODE'],plot_file='{}/diagnostic_wavecal_{}_{}.pdf'.format(redux['PARAMETERS']['REDUCTION_FOLDER'],src,redux['PARAMETERS']['MODE']))
			spflx.applyWaveCal(arcrecal)
# apply flux calibration
			if 'CAL_FLUX' in list(redux.keys()): spflx.applyFluxCal(redux['CAL_FLUX'])
			else: print('Warning: no flux calibration applied')
# apply telluric correction
			if tcorr==True: spflx.applyTelluricCal(redux['CAL_TELL'][tstar])
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

	return redux


############################################################
# KAST SPECTRAL ANALYSIS FUNCTIONS
# These functions perform some basic analysis on spectra
############################################################

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

def compareSpectra_simple(sp1,sp2orig,fit_range=[],plot=False,plot_file='',**kwargs):
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
		ax_top.plot(wave,u1,'k--')
		ax_top.plot(wave,f1,'k-')
		ax_top.plot(wave,f2*scale_factor,'m-')
		ax_top.plot(wave,numpy.zeros(len(wave)),'k:')
		ax_top.legend([sp1.name,sp2.name],fontsize=16)
		ax_top.set_xlim(xlim)
		ax_top.set_ylim(ylim)
		ax_top.set_ylabel('Flux Density',fontsize=16)
		ax_btm.plot(wave,f1-f2*scale_factor,'k-')
		ax_btm.plot(wave,u1,'k--')
		ax_btm.plot(wave,-1.*u1,'k--')
		ax_btm.fill_between(wave,-1.*u1,u1,color='k',alpha=0.1)
		ax_btm.plot(wave,numpy.zeros(len(wave)),'k--')
		ax_btm.set_xlim(xlim)
		ax_btm.set_ylim([-3.*numpy.nanmedian(u1),3.*numpy.nanmedian(u1)])
		ax_btm.set_xlabel('Wavelength (Angstrom)',fontsize=16)
		ax_btm.set_ylabel('O-C',fontsize=16)

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
		if plot_file!='': 
			fig.savefig(plot_file)
			plt.close()
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

def classifyTemplate():
	pass

def classifyIndices():
	pass

def EW():
	pass





