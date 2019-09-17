#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Name: SEBCS_SA                                                               #
# Modul: SEBCS_SA.py                                                               #
# Author: Dr. Jakub Brom, University of South Bohemia, Faculty of Agriculture  #
# e-mail: jbrom@zf.jcu.cz                                                      #
# Date: 2019-03-14                                                             #
# License: (c) 2016 - 2019 Jakub Brom, University of South Bohemia,             #
#			Faculty of Agriculture                                             #
################################################################################

# Data import

import os
import glob
import string
import datetime
import sys
import time
import tempfile
import numpy as np

from osgeo import gdal

# Import SEBCSlib methods
from SEBCSlib import GeoIO, HeatFluxes, MeteoFeatures, SolarRadBalance, VegIndices, WindStability

gio = GeoIO()
ht = HeatFluxes()
mf = MeteoFeatures()
sr = SolarRadBalance()
vi = VegIndices()
ws = WindStability()


class SEBCS_SA:
	""" SEBCS_SA is class for calculation surface energy balance
	and crop water stress features.


	**Required parameters: Layers**

	:param band_red: Path to raster file in red region of spectra. 
		Spectral reflectance (rel.; 0-1).
	:type band_red: str
	:param band_nir: Path to raster file in NIR region of spectra.
		Spectral reflectance (rel.; 0-1).
	:type band_nir: str
	:param band_tir: Path to surface temperature layer (raster file).
		Temperature (grades of Celsius).
	:type band_tir: str
	:param DMT: Path to DMT layer (raster file). Metres.
	:type DMT: str
	:param ta: Path to air temperature layer (raster file). Temperature
		(grades of Celsius).
	:type ta: str


	**Required parameters: Constants**

	:param Rs_in: Value of actual incomming global radiation flux 
		(W.m-2) measured by pyranometer.
	:type Rs_in: float
	:param Rh: Value of actual air relative humidity (%).
	:type Rh: float
	:param latitude: Value of mean (center) latitude of area used
		for calculation	(dec. grades, e.g.: 49.47575471).
	:type latitude: float
	:param longitude: Value of mean (center) longitude of area used
		for calculation	(dec. grades, e.g.: 15.30221431).
	:type longitude: float
	:param date_acq: Date of data acquisition (ISO format,
		e.g.: '2016-05-09').
	:type date_acq: datetime.date
	:param timeGMT: Time of data acquisition in GMT (ISO
		format, e.g.: '14:46:00.0000').
	:type timeGMT: datetime.time


	**Required parameters: Outputs**

	:param output_list: List of output layers file names. Following
		names for the output layers should be used:
			Rs_in	- Incoming shortwave radiation (W.m-2)
			Rs_out	- Reflected shortwave radiation (W.m-2)
			albedo	- Albedo (rel.)
			RL_in	- Incoming longwave radiation (W.m-2)
			RL_emit - Reflected longwave radiation (W.m-2)
			Rn	- Total net radiation (W.m-2)
			Ts	- Surface temperature (Degree C)
			emis	- Surface emissivity (rel.)
			LE	- Latent heat flux (W.m-2)
			LEp	- Latent heat flux - Penman pot. (W.m-2)
			LE_PT	- Latent heat flux - Priestley-Taylor (W.m-2)
			H	- Sensible heat flux (W.m-2)
			G	- Ground heat flux (W.m-2)
			E_int	- Evapotranspiration intensity (mmol.m-2.s-1)
			EF	- Evaporative fraction (rel.)
			Bowen	- Bowen ratio (-)
			omega	- Decoupling coefficient (rel.)
			CWSI	- CWSI (-)
			frict_U	- Friction velocity (m.s-1)
			ra	- Aerodynamic surface resistance (s.m-1)
			rc	- Resistance for water vapour transfer (s.m-1)
			NDVI	- NDVI (-)
			MSAVI	- MSAVI (-)
			NDMI	- NDMI (-)
			slope	- Slope (Degree)
			aspect	- Aspect (Degree)
	:type output_list: list
	:param output_folder: Path to output folder where the outputs
		will be created.
	:type output_folder: str
	:param out_file_name: Name of exported multiband raster.
		Default is None.
	:type out_file_name: str


	**Default: Setting**

	:param method: Setting of calculation method. There are two
		calculation methods:
		- 'aero': aerodynamic method of calcultaion based
				on SEBAL method is used (default)
		- 'grad': Method grad is used for calculation based
				on gradient approach.
		See Manual for more details.
	:type method: int
	:param sat_type: Setting of satellite type which was used for
		calculation. There are three different satellite type:
		- 'L5': data from Landsat 5 TM (default);
		- 'L7': data from Landsta 7 ETM+;
		- 'other' - data from other satellite or aerial system.
		Satellite type setting is important for albedo calculation.
	:type sat_type: int
	:param emis_rule: Setting of emissivity calcultaion for correction
		of surface temperature layer. There are two variants:
		- 'No': No correction on emissivity (default)
		- 'Yes': Correction of the surface temperature layer
				on emissivity. In case the source data have
				emissivity equal to 1.
	:type emis_rule: str
	:param output_driver: GDAL driver for output files format. Different
		formats should be used. GTiff is default.
	:type output_driver: str
	:param multiband: Option of multiband raster creation. Default is False.
	:type multiband: bool


	**Optional: Spectral layers**

	:param band_blue: Path to raster file in blue region of spectra.
		Spectral reflectance (rel.; 0-1).
	:type band_blue: str
	:param band_green: Path to raster file in green region of spectra.
		Spectral reflectance (rel.; 0-1).
	:type band_green: str
	:param band_sw1: Path to raster file in SWIR region of spectra
		(ca 1.61 micro meters). Spectral reflectance (rel.; 0-1).
	:type band_sw1: str
	:param band_sw2: Path to raster file in SWIR region of spectra
		(ca 2.2 micro meters). Spectral reflectance (rel.; 0-1).
	:type band_sw2: str


	**Optional: Additional layers**

	:param U: Path to wind speed layer (raster file). Velocity (m.s-1).
		This layer should be set if method aero is used.
	:type U: str
	:param h_min: Path to minimal vegetation heigth layer (raster file).
		Heigth (m). This layer should be set if method aero is used.
	:type h_min: str
	:param h_max: Path to maximal vegetation heigth layer (raster file).
		Heigth (m). This layer should be set if method aero is used.
	:type h_max: str
	:param mask: Path to mask (raster layer) for the output. Mask should be
		real data type. Output layers are multiplied by mask where 0 are nodata 
		and 1 are data. Mask is used for data export only. Default is None.
	:type mask: str
	:param albedo: Path to albedo layer (raster file). This option
		is used prior to calculation of albedo by SEBCS SA.
		Deafult is None.
	:type albedo: str
	:param canopy: Path to canopy heigth layer. This option is used prior
		to calculation of canopy heigth from the min and max heigth values.
		Default is None.
	type canopy: str


	**Optional: Constants**

	:param Z_st: Heigth of wind speed heigth measurement (m).
		Default is 10 m.
	:type Z_st: float
	:param Z: Heigth of mixed layer (m). Default is 200 m
	:type Z: float
	:param cp: Specific heat capacity of dry air (J.kg-1.K-1).
		Default cp = 1012 J.kg-1.K-1 
	:type cp: float	


	**Returns**

	:returns: Raster file or files for calculated variables. Number
		of outputs and also type of outputs are set by output_list
		variable (list of names). Output file(s) is saved into
		output_folder.
	"""

	def __init__(self,
				band_red, band_nir, band_tir, DMT, ta,
				Rs_in = None, Rh = None, latitude = None, longitude = None, date_acq = None, timeGMT = None,
				output_list = None, output_folder = "/home", out_file_name = None,
				method = "aero", sat_type = "L5", emis_rule = "No", output_driver = "GTiff", multiband = False,
				band_blue = None, band_green = None, band_sw1 = None, band_sw2 = None,
				U = None, h_min = None, h_max = None, mask = None, albedo = None, canopy = None,
				Z_st=10.0, Z = 200, cp = 1012):
		
		# Variables declaration and import raster files:
		## Required layers	
		self.band_red = self.rasterToArray(band_red)
		if self.band_red is None:
			raise Exception("Could not open red band")
		self.band_nir = self.rasterToArray(band_nir)
		if self.band_nir is None:
			raise Exception("Could not open NIR band")
		self.band_tir = self.rasterToArray(band_tir)
		if self.band_tir is None:
			raise Exception("Could not open thermal band")
		self.DMT = self.rasterToArray(DMT)
		if self.DMT is None:
			raise Exception("Could not open DMT")
		self.ta = self.rasterToArray(ta)
		if self.ta is None:
			raise Exception("Could not open air temperature layer")
		
		## Required constants
		self.Rs_in = Rs_in
		if self.Rs_in is None:
			raise ValueError("Value of incomming shortwave radiation has not been set")
		self.Rh = Rh
		if self.Rh is None:
			raise ValueError("Value of air relative humidity has not been set")
		self.latitude = latitude
		if self.latitude is None:
			raise ValueError("Latitude has not been set")
		self.longitude = longitude
		if self.longitude is None:
			raise ValueError("Longitude has not been set")
		self.date_acq = date_acq
		if self.date_acq is None:
			raise ValueError("Date of data acquisition has not been set")
		self.timeGMT = timeGMT
		if self.timeGMT is None:
			raise ValueError("Time of data acquisition has not been set")
		
		## Required paths
		self.output_list = output_list
		if self.output_list is None:
			raise ValueError("List of output layers has not been set")
		self.output_folder = output_folder
		self.out_file_name = out_file_name
		
		# Default and otional setting
		self.method = method
		self.sat_type = sat_type
		self.emis_rule = emis_rule
		self.output_driver = output_driver
		self.multiband = multiband
		
		# Default and otional layers
		self.band_blue = self.rasterToArray(band_blue)
		self.band_green = self.rasterToArray(band_green)
		self.band_sw1 =	self.rasterToArray(band_sw1)
		self.band_sw2 =	self.rasterToArray(band_sw2)
		self.U = self.rasterToArray(U)
		self.h_min = self.rasterToArray(h_min)
		self.h_max = self.rasterToArray(h_max)
		self.canopy = self.rasterToArray(canopy)
		self.mask = self.rasterToArray(mask)
		self.albedo = self.rasterToArray(albedo)
		self.Z_st = Z_st
		
		# Constants
		self.Z = Z						# Heigth of mixed layer
		self.cp = cp					# Thermal capacity of dry air
		
		# Check for layers extension
		self.in_lyrs_list = [self.band_red, self.band_nir,
							self.band_tir, self.DMT, self.ta,
							self.band_blue, self.band_green,
							self.band_sw1, self.band_sw2,
							self.h_min, self.h_max, self.canopy, self.U,
							self.mask, self.albedo]

		
		# Input srs and geotransform
		self.gtransf, self.srs, self.xSize, self.ySize, self.EPSG = gio.readGeo(band_red)		# Reading geotransformation, SRS and pixel size

		# Vegetation indices
		self.ndvi = vi.vi_NDVI(self.band_red, self.band_nir)
		self.msavi = vi.vi_MSAVI(self.band_red, self.band_nir)
		self.ndmi = vi.vi_NDMI(self.band_nir, self.band_sw1)
		self.Fc = vi.fractVegCover(self.ndvi)
		
		# Slope and Aspect
		try:
			self.slope, self.aspect = sr.slopeAspect(self.DMT, self.xSize, self.ySize)
		except:
			self.slope = None
			self.aspect = None
			
		# Surface emissivity
		try:
			self.emis = mf.emissivity(self.Fc, self.band_red, self.ndvi)
		except:
			self.emis = None
			
		# Surface temperature
		self.ts = self.surfaceTemperature(self.band_tir, self.emis, self.emis_rule)
			
		# Auxiliary variables
		self.ta_Z = mf.airTemp(self.Z, self.Z_st, self.ta)
		self.P_Z = mf.airPress(self.Z, self.DMT)
		self.P_Zst = mf.airPress(self.Z_st, self.DMT)
		self.E_Z_sat = mf.satVapourPress(self.ta_Z)
		self.E_Zst_sat = mf.satVapourPress(self.ta)
		self.e_Z = mf.vapourPress(self.E_Z_sat, self.Rh)
		self.e_Zst = mf.vapourPress(self.E_Zst_sat, self.Rh)
		self.VPD = mf.vpd(self.E_Z_sat, self.e_Z)
		self.rho = mf.airDensity(self.ta_Z)
		self.latent = mf.latent(self.ta_Z)
		self.gamma = mf.gamma(self.P_Z, self.latent, self.cp)
		self.delta = mf.delta(self.ts, self.ta_Z)
		self.Es_sat = mf.satVapourPress(self.ts)
		
		# Radiation balance
		## Incomming shortwave radiation
		self.in_Rs = sr.solarInTopo(self.Rs_in, self.slope, self.aspect,
									self.latitude, self.longitude,
									self.date_acq, self.timeGMT)
		
		## Albedo
		if self.albedo is None:
			self.albedo = self.albedoCalc(self.sat_type, self.band_red, self.band_nir, 
									self.band_blue, self.band_green, self.band_sw1,
									self.band_sw2, self.ndvi, self.msavi)
		
		## Reflected radiation
		self.out_Rs = seb.reflectRs(self.in_Rs, self.albedo)
		
		## Longwave radiation balance
		self.RL_in, self.RL_out = self.longWaveRad(self.e_Z, self.ta_Z_K, self.ts_K, self.emis)
		
		## Net radiation
		self.Rn = self.in_Rs - self.out_Rs + self.RL_in - self.RL_out
		
		# Heat fluxes
		## Ground heat flux
		self.G = seb.groundFlux(self.ndvi, self.Rn, self.ts_C, self.albedo)
		self.Rn_G = seb.difRn_G(self.Rn, self.G)
		
		## H, LE, ra, friction velocity
		try:
			self.H, self.LE, self.EF, self.ra, self.frict = self.heatFluxes(self.ta_Z, self.ts_C,
										self.ts_K, self.Rn_G, self.albedo, self.ndvi, self.method,
										self.U, self.canopy, self.h_min, self.h_max, self.msavi, self.ro,
										self.Z, self.Z_st, self.cp)
		except:
			print("Fluxes no")
		
		
		## Evapotranspiration intensity
		self.E_int = seb.intensityE(self.LE, self.latent)
		
		## Potential LE
		self.LEp = seb.fluxLE_p(self.Rn_G, self.delta, self.VPD, self.ra, self.gamma, 
								self.ro, self.cp)
								
		## Equilibrium LE
		self.LE_EQ = seb.fluxLE_EQ(self.Rn_G, self.delta, self.gamma)
		
		## Priestley-Taylor LE (evaporation from wet surface)
		self.LE_PT = seb.fluxLE_PT(self.LE_EQ)
		
		# Heat flux indices
		## Bowen ratio
		self.bowen = seb.bowen(self.H, self.LE)
		
		## Omega factor (Decoupling coefficient)
		self.omega = seb.omega(self.LE, self.LEp)
		
		## rc - surface resistance for water vapour transfer
		self.rc = seb.rs(self.delta, self.gamma, self.omega, self.ra)
		
		## CWSI
		self.CWSI = self.cropWaterStressIndex(self.LEp, self.ra, self.rc,
									self.E_Z_sat, self.e_Z, self.ro, self.delta,
									self.gamma, self.cp)
		
		# Data export to raster
		self.dataExport(self.output_list, self.output_folder, self.output_driver, self.gtransf, self.srs, self.EPSG, self.out_file_name, self.multiband, self.mask)

	
	def rasterToArray(self, layer):
		"""Conversion of raster layer to numpy array.
		:param layer: Path to raster layer.
		:type layer: str
		
		:return in_layer: raster file converted to numpy array
		:type in_layer: numpy.ndarray
		"""
		try:
			if layer is not None:
				in_layer = gdal.Dataset.ReadAsArray(gdal.Open(layer)).astype(np.float32)
				in_layer = np.nan_to_num(in_layer)
			else:
				in_layer = None
			return in_layer
		except:
			in_layer = None
			return in_layer


	def lyrsExtent(self):
		"""Check differences between size of the input layers.
		"""
		
		in_lyrs_list_true = [i for i in self.in_lyrs_list if i != None]	# List of input layers without empty (None)
		
		len_list = []
		 
		for i in in_lyrs_list_true:
			len_list.append(len(i))
			
		if max(len_list) != min(len_list):
			raise Exception("Selected layers differ in spatial extent (number of columns or rows).")

	
	def vegIndices(self, bandRed, bandNIR, bandSWIR1 = None): 
		"""Vegetation indices calculation.
		:param bandRed: Reflectance in red spectral region.
		:type bandRed: numpy.ndarray
		:param bandNIR: Reflectance in NIR spectral region.
		:type bandNIR: numpy.ndarray
		:param bandSWIR1: Reflectance in SWIR1 spectral region
			(ca 1.6 micro meters).
		:type bandSWIR1: numpy.ndarray
		
		:return ndvi: NDVI according to Tucker (1979).
		:rtype ndvi: numpy.ndarray
		:return msavi: MSAVI2 according to Qi et al. (1994)
		:rtype msavi: numpy.ndarray
		:return ndmi: NDMI
		:rtype ndmi: numpy.ndarray
		:return Fc: Fractional vegetation cover.
		:rtype Fc: numpy.ndarray
		"""
		
		# NDVI, MSAVI
		try:
			if bandRed is not None and bandNIR is not None:
				ndvi = seb.vi_NDVI(bandRed, bandNIR)
				msavi = seb.vi_MSAVI(bandRed, bandNIR)
			else:
				if bandRed is not None:
					ndvi = np.zeros(bandRED.shape)
					msavi = np.zeros(bandRED.shape)								
				else:
					ndvi = np.zeros(bandNIR.shape)
					msavi = np.zeros(bandNIR.shape)
		except:
			ndvi = None
			msavi = None
					
		# NDMI
		try:
			if bandNIR is not None and bandSWIR1 is not None:
				ndmi = seb.vi_NDMI(bandNIR, bandSWIR1)
			else:
				if bandNIR is not None:
					ndmi = np.zeros(bandNIR.shape)
				else:
					ndmi = np.zeros(bandSWIR1.shape)
		except:
			ndmi = None
		
		# Fractional vegetation cover
		if ndvi is not None:
			Fc = seb.fractVegCover(ndvi)
		else:
			Fc = None
			
		return ndvi, msavi, ndmi, Fc
	
				
	def surfaceTemperature(self, tir_band, emissivity = None, emis_rule = "No"):
		"""Correction of surface temperature on emissivity.
		:param tir_band: Layer of surface temperature (grades of C).
		:type tir_band: numpy.ndarray
		:param emissivity: Layer of emissivity (rel.).
		:type emissivity: numpy.ndarray
		:param emis_rule: Setting if the correction will be done or not.
			No is default.
		:type emis_rule: str
		
		:return ts_C: Layer of corrected surface temperature (grades C).
		:rtype ts_C: numpy.ndarray
		:return ts_K: Layer of corrected surface temperature (grades K).
		:rtype ts_K: numpy.ndarray
		"""
		try:
			if tir_band is not None or emissivity is not None:
				if emis_rule != "No":
					ts = ((tir_band + 273.16)/emissivity**0.25) - 273.16
				else:
					ts = tir_band
			else:
				if tir_band is not None:
					ts = tir_band
				else:
					ts = None
		except:
			ts = None
		
		if ts is None:
			raise ValueError("Surface temperature has not been calculated.")
		
		return ts
				
			
	def albedoCalc(self, sat_type, band_red, band_nir, 
			band_blue = None, band_green = None, band_sw1 = None, band_sw2 = None,
			ndvi = None, msavi = None):
				
		"""Calculation of Albedo.
		"""
				
		try:
			if sat_type == "other" or band_blue is None or band_green is None or band_sw1 is None or band_sw2 is None:
				if ndvi is not None or msavi is not None: 
					albedo = sr.albedoBrom(ndvi, msavi)
				else:
					raise ValueError("Albedo has not been calculated.")
			else:
				albedo = sr.albedoLandsat(band_blue, band_green, band_red, band_nir, band_sw1, band_sw2, sat_type)
		except:
			raise ValueError("Albedo has not been calculated.")	
		
		return albedo


	def longWaveRad(self, e_Z, ta, ts, emissivity):
		"""Calculation of long wave radiation downward and upward
		fluxes.
		
		:param e_Z: Saturated water vapour pressure (kPa) at Z level
			(default at 200 m above surface).
		:type e_Z: numpy.ndarray
		:param ta_K: Air temperature at level Z (default at 200 m).
			in K.
		:type ta_K: numpy.ndarray
		:param ts_K: Surface temperature in K.
		:type ts_K: numpy.ndarray
		
		:return RL_in: Incomming longwave radiation (W.m-2).
		:rtype RL_in: numpy.ndarray
		:return RL_out: Emitted longwave radiation (W.m-2).
		:rtype RL_out: numpy.ndarray
		"""
		
		# Atmospheric emissivity
		emis_atm = sr.atmEmissivity(e_Z, ta)
		RL_in = sr.downRL(ta, emis_atm)
		RL_out = sr.outRL(ts, emissivity)
		
		return RL_in, RL_out
		

	def heatFluxes(self, ta_Z, ts_C, ts_K, Rn_G, albedo, ndvi, method = "aero",
					U = None, canopy = None, h_min = None, h_max = None, msavi = None, ro = None, 
					Z = 200.0, Z_st = 10.0, cp = 1012.0):
						
		"""Calculation of heat fluxes using aerodynamic or gradient
		method.
		"""

		if method == "aero":
			if U is None:
				raise ValueError("Missing data. SEBCS can not calculate heat fluxes.")
			else:
				if canopy is None:
					if h_min is None or h_max is None:
						raise ValueError("Missing data. SEBCS can not calculate heat fluxes.")
					else:
						h_eff = seb.vegHeight(h_min, h_max, msavi)
				else:
					h_eff = canopy
				disp = seb.zeroPlaneDis(h_eff)
				z0m = seb.z0m(ndvi, albedo)
				z0h = seb.z0h(z0m)
				Z_d = seb.z_d(Z, disp)
				Z_d_z0m = seb.ln_z_d_z0(Z_d, z0m)
				Z_d_z0h = seb.ln_z_d_z0(Z_d, z0h)
				U_Z = seb.windSpeedZ(U, Z, Z_st)
				psi_m, psi_h, frict = seb.stabCoef(U_Z, ta_Z, ts_C, Z_d, Z_d_z0m, Z_d_z0h)
				ra = seb.raThom(U_Z, psi_m, psi_h, Z_d_z0m, Z_d_z0h)
				Tmax = seb.maxT(Rn_G, ra, ta_Z, ro, cp)
				L_dry = seb.dryL(frict, ro, ta_Z, Rn_G, cp)
				FL = seb.functFL(L_dry, frict, ra)
				dT_dry = seb.dTdry(FL, ts_K)
				dT_b = seb.coef_b(Tmax, ta_Z, dT_dry)
				dT_a = seb.coef_a(ta_Z, dT_b)
				dT = seb.dT(ts_K, dT_a, dT_b)
			
				H = seb.fluxHAer(ra, ro, dT, cp)
				LE = seb.fluxLE(Rn_G, H)
				EF = seb.aeroEF(LE, Rn_G)
		else:
			EF = seb.gradEF(ts_C, ta_Z)
			LE = seb.gradLE(EF, Rn_G)
			H = seb.gradH(LE, Rn_G)
			ra = seb.raGrad(H, ro, ts_C, ta_Z, cp)
			frict = np.zeros(ra.shape)
# 	except:
# 			raise Exception("Heat fluxes has not been calculated.")
			
		return H, LE, EF, ra, frict
			
		
	def dataExport(self, output_list_names, out_folder, driver, gtransf, prj,
					EPSG, out_file_name, multiband, mask):
		"""Export output layers.
		"""
		
		out_names_all = ("Rs_in", "Rs_out", "albedo", "RL_in", "RL_emit", "Rn", 
					   "Ts", "emis", "LE", "LEp", "LE_PT", "H", "G", "E_int", 
					   "EF", "Bowen", "omega", "CWSI", "frict_U", "ra",
					   "rc", "NDVI", "MSAVI", "NDMI", "slope", "aspect")
		out_lyrs_all = (self.in_Rs, self.out_Rs, self.albedo, self.RL_in, self.RL_out,
						self.Rn, self.ts_C, self.emis, self.LE, self.LEp, self.LE_PT,
						self.H, self.G, self.E_int, self.EF, self.bowen, self.omega,
						self.CWSI, self.frict, self.ra, self.rc, self.ndvi,
						self.msavi, self.ndmi, self.slope, self.aspect)
		
		out_lyrs = []
		out_names = []
		
		# Create list of output layers indices
		index_list = []
		
		for i in out_names_all:
			rule = i in output_list_names
			if rule == True:
				index_list.append(out_names_all.index(i))
		
		# Export rasters:
		for i in index_list:
			if out_lyrs_all[i] is not None:
				# Chose correct layer and its name and path and create list
				# of out lyrs and its names
				
				out_lyrs.append(out_lyrs_all[i])
				out_names.append(out_names_all[i])
		
		# Mask application
		if mask is not None:
			for i in range(0, len(out_lyrs)):
				out_lyrs[i] = out_lyrs[i] * mask
		
		# Replace nan and inf values in output lyrs by zero.
		for i in range(0, len(out_lyrs)):
			out_lyrs[i] = np.nan_to_num(out_lyrs[i])
			out_lyrs[i] = np.where(out_lyrs[i] == np.isnan, 0.0, out_lyrs[i])
			out_lyrs[i] = np.where(out_lyrs[i] == -np.inf, 0.0, out_lyrs[i])
			out_lyrs[i] = np.where(out_lyrs[i] == np.inf, 0.0, out_lyrs[i])
				
		# Export data
		seb.arrayToRast(out_lyrs, out_names, prj, gtransf, EPSG, out_folder,
						driver, out_file_name, multiband)
			

# if __name__ == "__main__":
# 	band_red = "B3.tif"
# 	band_nir = "B4.tif"
# 	band_tir = "Ts.tif"
# 	DMT = "DMT.tif"
# 	ta = "Ta_2m.tif"
# 	Rs_in = 600
# 	Rh = 50
# 	latitude = 49.50833
# 	longitude = 15.50833
# 	date_acq = "2009-06-14"
# 	timeGMT = "9:53:0.000"
# 	output_list = ["Rs_in"]
# 	output_folder = os.path.join(os.path.dirname(__file__), "vystup")
# 	out_file_name = "pokus1"
# 	method="aero"
# 	sat_type="L5"
# 	emis_rule="No"
# 	output_driver="GTiff"
# 	multiband = False
# 	band_blue = "B1.tif"
# 	band_green="B2.tif"
# 	band_sw1="B5.tif"
# 	band_sw2="B7.tif"
# 	U= "U_2.tif"
# 	h_min = "LU_min.tif"
# 	h_max = "LU_max.tif"
# 	mask = "maska.tif"
# 	Z_st = 2
# 	albedo = None
# 	
# 	SEBCS_SA(band_red, band_nir, band_tir, DMT, ta,
# 				Rs_in, Rh, latitude, longitude, date_acq, timeGMT,
# 				output_list, output_folder, out_file_name,
# 				method, sat_type, emis_rule, output_driver, multiband,
# 				band_blue, band_green, band_sw1, band_sw2,
# 				U, h_min, h_max, mask, albedo,
# 				Z_st)
