# !/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
   /***************************************************************************
   SEBCSlib

   Collection of functions for energy balance and its components calculation.

                                -------------------
          begin                : 14-03-01
          date                 : 19-09-17
          git sha              : $Format:%H$
          copyright            : (C) 2019 Jakub Brom
          email                : jbrom@zf.jcu.cz

   ***************************************************************************/
   /***************************************************************************

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License  as published  by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License along
   with this program.  If not, see <https://www.gnu.org/licenses/>.

   ***************************************************************************/
"""

# Imports
import numpy as np
import math
import os
import warnings
from datetime import datetime
from scipy import ndimage
from osgeo import gdal
from osgeo import osr


# noinspection PyMethodMayBeStatic
class GeoIO:
	"""
	Class includes functions for reading of geotransformation features
	from raster and for writing ne rasters from numpy arrays. 
	"""

	def __init__(self):
		return

	def readGeo(self, rast):
		"""
		Reading geographical information from raster using
		GDAL.

		:param rast: Path to raster file in GDAL accepted format.
		:type rast: str

		:returns: The affine transformation coefficients.
		:rtype: tuple
		:returns: Projection information of the raster (dataset).
		:rtype: str
		:returns: Pixel width (m) on X axis.
		:rtype: float
		:returns: Pixel height (m) on Y axis.
		:rtype: float
		:returns: EPSG Geodetic Parameter Set code.
		:rtype: int
		"""

		try:
			ds = gdal.Open(rast)

			gtransf = ds.GetGeoTransform()
			prj = ds.GetProjection()
			x_size = gtransf[1]
			y_size = gtransf[5] * (-1)

			srs = osr.SpatialReference(wkt=prj)
			if srs.IsProjected:
				EPSG = int(srs.GetAttrValue("authority", 1))
			else:
				EPSG = None

			del ds

			return gtransf, prj, x_size, y_size, EPSG

		except IOError:
			warnings.warn("Geographical information has not been readed.", stacklevel=3)

			gtransf = None
			prj = None
			x_size = None
			y_size = None
			EPSG = None

			return gtransf, prj, x_size, y_size, EPSG

	def arrayToRast(self, arrays, names, prj, gtransf, EPSG, out_folder,
	                driver_name="GTiff", out_file_name=None, multiband=False):
		"""Export numpy 2D arrays to multiband or singleband raster
		files. Following common raster formats are accepted for export:\n
		
		- ENVI .hdr labeled raster format\n
		- Erdas Imagine (.img) raster format\n
		- Idrisi raster format (.rst)\n
		- TIFF / BigTIFF / GeoTIFF (.tif) raster format\n
		- PCI Geomatics Database File (.pix) raster format\n
		
		:param arrays: Numpy array or list of arrays for export
					   to raster.
		:type arrays: numpy.ndarray or list of numpy.ndarray
		:param names: Name or list of names of the exported bands 
					  (in case of multiband raster) or particular
					  rasters (in case of singleband rasters).
		:type names: str or list of str
		:param prj: Projection information of the exported raster
					(dataset).
		:type prj: str
		:param gtransf: The affine transformation coefficients.
		:type gtransf: tuple
		:param EPSG: EPSG Geodetic Parameter Set code.
		:type EPSG: int
		:param out_folder: Path to folder where the raster(s) will
						   be created.
		:type out_folder: str
		:param driver_name: GDAL driver. 'GTiff' is default.
		:type driver_name: str
		:param out_file_name: Name of exported multiband raster. Default
							  is None.
		:type out_file_name: str
		:param multiband: Option of multiband raster creation. Default 
						  is False.
		:type multiband: bool
		
		:returns: Raster singleband or multiband file(s)
		:rtype: raster
		"""

		# Convert arrays and names on list
		if type(arrays) is not list:
			arr_list = list()
			arr_list.append(arrays)
			arrays = arr_list
		if type(names) is not list:
			names_list = list()
			names_list.append(names)
			names = names_list

		if out_file_name is None:
			out_file_name = ""
			multiband = False

		# Drivers and suffixes
		driver_list = ["ENVI", "HFA", "RST", "GTiff", "PCIDSK"]     # GDAL driver for output files
		out_suffixes = ["", ".img", ".rst", ".tif", ".pix"]         # Suffixes of output names

		# Test driver
		if driver_name not in driver_list:
			raise ValueError("Unknown driver. Data could not be exported.")

		driver_index = driver_list.index(driver_name)
		suffix = out_suffixes[driver_index]

		if multiband is True and driver_name != "RST":
			out_file_name, ext = os.path.splitext(out_file_name)
			out_file = os.path.join(out_folder, out_file_name + suffix)

			try:
				driver = gdal.GetDriverByName(driver_name)
				ds = driver.Create(out_file, arrays[0].shape[1], arrays[0].shape[0], len(arrays), gdal.GDT_Float32)
				ds.SetProjection(prj)
				ds.SetGeoTransform(gtransf)
				if EPSG is not None:
					outRasterSRS = osr.SpatialReference()
					outRasterSRS.ImportFromEPSG(EPSG)
					ds.SetProjection(outRasterSRS.ExportToWkt())
				j = 1
				for i in arrays:
					ds.GetRasterBand(j).WriteArray(i)
					ds.GetRasterBand(j).SetDescription(names[j - 1])
					ds.GetRasterBand(j).SetMetadataItem("Band name", names[j - 1])
					ds.GetRasterBand(j).FlushCache()
					j = j + 1

				del ds

			except IOError:
				raise Exception("Raster file {} has not been created.".format(out_file_name + suffix))

		else:
			for i in range(0, len(arrays)):
				try:
					out_file_name, ext = os.path.splitext(names[i])
					out_file = os.path.join(out_folder, out_file_name + suffix)
					driver = gdal.GetDriverByName(driver_name)
					ds = driver.Create(out_file, arrays[i].shape[1], arrays[i].shape[0], 1, gdal.GDT_Float32)
					ds.SetProjection(prj)
					ds.SetGeoTransform(gtransf)
					if EPSG is not None:
						outRasterSRS = osr.SpatialReference()
						outRasterSRS.ImportFromEPSG(EPSG)
						ds.SetProjection(outRasterSRS.ExportToWkt())
					ds.GetRasterBand(1).WriteArray(arrays[i])
					ds.GetRasterBand(1).SetDescription(names[i])
					ds.GetRasterBand(1).SetMetadataItem("Band name", names[i])
					ds.GetRasterBand(1).FlushCache()

					del ds

				except IOError:
					raise Exception("Raster file {} has not been created.".format(names[i] + suffix))


# noinspection PyMethodMayBeStatic
class SolarRadBalance:
	"""
	Class contains functions for calculation of solar radiation balance
	and topographic features: incident radiation in dependence on surface geometry,
	slope of terrain, aspect of terrain, albedo, longwave radiation fluxes
	and atmospheric emissivity, shortwave radiation reflectance
	and total net radiation
	"""

	def __init__(self):
		return

	def slopeAspect(self, DMT, x_size, y_size):
		"""Slope and aspect of terrain (DMT) in degrees.

		:param DMT: Digital model of terrain (m a.s.l.)
		:type DMT: numpy.ndarray
		:param x_size: Size of pixel in x axis (m)
		:type x_size: float
		:param y_size: Size of pixel in y axis (m)
		:type y_size: float

		:returns: Slope of the terrain :math:`(\SI{}\degree)`
		:rtype: numpy.ndarray
		:returns: Aspect of the terrain :math:`(\SI{}\degree)`
		:rtype: numpy.ndarray
		"""

		x, y = np.gradient(DMT)
		slope = np.arctan(np.sqrt((x / x_size) ** 2.0 + (y / y_size) ** 2.0)) * 180 / np.pi
		aspect = 270 + np.arctan(x / y) * 180 / np.pi
		aspect = np.where(y > 0, aspect, aspect - 180)

		# Replacing nan values to 0 and inf to value
		slope = np.nan_to_num(slope)
		aspect = np.nan_to_num(aspect)
		del x
		del y
		return slope, aspect

	def solarInTopo(self, Rs_in, slope, aspect, latitude, longitude, date_acq, time_acq):

		"""Calculation of incident shortwave solar radiation flux
		according to the solar geometry, position (latitude
		and longitude) and shape of surface (slope and orientation).
		Flux of the solar energy :math:`(W.m^{-2})` is calculated on basis
		of the measured global radiation using pyranometer (incomming
		global radiation). Diffuse part of radiation is not separated
		in calculation.
		
		:param Rs_in: Global radiation measured by pyranometer :math:`(W.m^{-2})`.
		:type Rs_in: float
		:param slope: Slope of the terrain :math:`(\SI{}\degree)`.
		:type slope: numpy.ndarray
		:param aspect: Orientation of the terrain :math:`(\SI{}\degree)`.
		:type aspect: numpy.ndarray
		:param latitude: Mean latitude of the data in decimal degrees
		:type latitude: float
		:param longitude: Mean longitude of the data in decimal degrees
		:type longitude: float
		:param date_acq: Date of data acquisition in iso format
						 ('YYYY-mm-dd')
		:type date_acq: datetime.date
		:param time_acq: Time in GMT in datetime.time format
						 ('HH:MM:SS.SS')
		:type time_acq: datetime.time
		
		:returns: Incident shortwave radiation :math:`(W.m^{-2})`
				corrected on the terrain and solar geometry.
		:rtype: numpy.ndarray
		"""

		dat = datetime.strptime(str(date_acq), "%Y-%m-%d")  # conversion of iso date to datetime.date
		t = datetime.strptime(str(time_acq), "%H:%M:%S.%f")  # conversion of iso time to datetime.time
		dec_t = float(t.hour) + float(t.minute) / 60.0 + float(t.second) / 3600  # decimal time

		N = dat.strftime("%j")  # day of the year
		solar_time = dec_t + longitude / 360.0 * 24.0  # solar time
		hs = (12.0 - float(solar_time)) * 15.0  # local solar hour angle (°)
		ds = 23.45 * math.sin(360.0 * (284.0 + float(N)) / 365.0 * math.pi / 180.0)  # solar declination (°)
		ds_rad = math.radians(ds)  # solar declination (rad.)
		L_rad = math.radians(latitude)  # latitude (rad.)
		hs_rad = math.radians(hs)  # local solar hour angle (rad.)

		sin_alpha = (math.sin(L_rad) * math.sin(ds_rad) + math.cos(L_rad) * math.cos(ds_rad) * math.cos(
			hs_rad))  # sin of solar height angle
		if sin_alpha < 0.0:
			sin_alpha = 0.0

		slope_rad = np.radians(slope)
		asp_rad = np.radians((aspect - 180) * (-1))  # aspect transformation

		cosi = (np.sin(ds_rad) * (np.sin(L_rad) * np.cos(slope_rad)
		                          - np.cos(L_rad) * np.sin(slope_rad) * np.cos(asp_rad))
		        + np.cos(ds_rad) * np.cos(hs_rad) * (np.cos(L_rad)
		                                             * np.cos(slope_rad) + np.sin(L_rad) * np.sin(slope_rad)
		                                             * np.cos(asp_rad)) + np.cos(ds_rad) * np.sin(slope_rad)
		        * np.sin(asp_rad) * np.sin(hs_rad))

		Is = float(Rs_in) / sin_alpha  # radiation perpendicular to solar beam angle
		Rs_in_corr = Is * cosi  # radiation corrected on solar and terrain geometry

		return Rs_in_corr

	def atmEmissivity(self, e_Z, ta):
		"""
		Atmospheric emissivity calculated according to Idso
		(see Brutsaert 1982).

		:param e_Z: Atmospheric water vapour pressure (kPa)
		:type e_Z: numpy.ndarray, float
		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float

		:returns: Air emissivity (rel.)
		:rtype: numpy.ndarray, float
		"""

		try:
			emis_a = 1.24 * (e_Z * 10.0 / (ta + 273.16)) ** (1.0 / 7.0)
		except ArithmeticError:
			raise ArithmeticError("Air emissivity has not been calculated.")

		return emis_a

	def downRL(self, ta, emis_a):
		"""
		Funtion calculates downward flux of longwave radiation :math:`(W.m^{-2})`

		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float
		:param emis_a: Air emissivity (rel.)
		:type emis_a: numpy.ndarray, float

		:return RL_in: Downward flux of longwave radiation :math:`(W.m^{-2})`
		:rtype RL_in: numpy.ndarray, float
		"""

		try:
			RL_in = emis_a * 5.6703 * 10.0 ** (-8.0) * (ta + 273.16) ** 4
		except ArithmeticError:
			raise ArithmeticError("Downward longwave radiation flux has not been calculated.")

		return RL_in

	def outRL(self, ts, emiss):
		"""
		Upward flux of longwave radiation :math:`(W.m^{-2})`

		:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
		:type ts: numpy.ndarray
		:param emiss: Surface emissivity (rel.)
		:type emiss: numpy.ndarray

		:returns: Upward flux of longwave radiation :math:`(W.m^{-2})`
		:rtype: numpy.ndarray
		"""

		RL_out = emiss * 5.6703 * 10.0 ** (-8.0) * (ts + 273.16) ** 4

		return RL_out

	def albedoBrom(self, ndvi, msavi, c_a=0.08611, c_b=0.894716, c_c=5.558657,
	               c_d=-0.11829, c_e=-1.9818, c_f=-4.50339, c_g=-11.4625,
	               c_h=7.461454, c_i=5.299396, c_j=4.76657, c_k=-2.3127, c_l=-3.42739):
		"""
		Albedo (rel.) calculated according to Duffková and Brom et al. (2012)

		:param ndvi: Normalized Difference Vegetation Index (-)
		:type ndvi: numpy.array
		:param msavi: Modified Soil Adjusted Vegetation Index (-) according to Gao et al. 1996.
		:param c_a: constant. Default a = 0.08611
		:type c_a: float
		:param c_b: constant. Default a = 0.894716
		:type c_b: float
		:param c_c: constant. Default a = 5.558657
		:type c_c: float
		:param c_d: constant. Default a = -0.11829
		:type c_d: float
		:param c_e: constant. Default a = -1.9818
		:type c_e: float
		:param c_f: constant. Default a = -4.50339
		:type c_f: float
		:param c_g: constant. Default a = -11.4625
		:type c_g: float
		:param c_h: constant. Default a = 7.461454
		:type c_h: float
		:param c_i: constant. Default a = 5.299396
		:type c_i: float
		:param c_j: constant. Default a = 4.76657
		:type c_j: float
		:param c_k: constant. Default a = -2.3127
		:type c_k: float
		:param c_l: constant. Default a = -3.42739
		:type c_l: float

		:returns: Albedo (rel.)
		:rtype: numpy.ndarray
		"""

		albedo = (c_a + c_b * msavi + c_c * msavi ** 2 + c_d * ndvi + c_e * msavi ** 3
		          + c_f * msavi * ndvi + c_g * msavi ** 2 * ndvi
		          + c_h * msavi * ndvi ** 2 + c_i * msavi ** 2 * ndvi ** 2
		          + c_j * msavi ** 3 * ndvi + c_k * msavi ** 3 * ndvi ** 2
		          + c_l * msavi * ndvi ** 3)

		return albedo

	def albedoLandsat(self, blue, green, red, nir, swir1, swir2, sat_type):
		"""
		Albedo (rel.) calculated for Landsat satellite sensors. Albedo
		for Landsat 4 TM, 5 TM and Landsat 7 ETM+ is calculated
		according to Tasumi et al. (2008). Albedo for Landsat 8 OLI/TIRS
		is calculated according to Olmeo et al. (2017). Albedo is computed
		with spectral reflectance bands on relative scale (0 to 1).\n
		Note: This algorithm might be used with another data from different
		devices, however a comparable spectral data (bands) should be used.

		:param blue: Blue band (rel.)
		:type blue: numpy.ndarray
		:param green: Green band (rel.)
		:type green: numpy.ndarray
		:param red: Red band (rel.)
		:type red: numpy.ndarray
		:param nir: NIR band (rel.)
		:type nir: numpy.ndarray
		:param swir1: SWIR1 band on ca 1.61 :math:`\mu m` (rel.)
		:type swir1: numpy.ndarray
		:param swir2: SWIR2 band on ca 2.2 :math:`\mu m` (rel.)
		:type swir2: numpy.ndarray
		:param sat_type: Type of Landsat satellite: \n
						tL5 - Landsat 4 TM, 5 TM or Landsat 7 ETM+\n
						t\tL8 - Landsat 8 OLI/TIRS
		:type sat_type: str

		:returns: Albedo (rel.)
		:rtype: numpy.ndarray
		
		\n\n
		**References:**\n
		*G.F. Olmedo, S. Ortega-Farias, D. Fonseca-Luengo,
		D. de la Fuente-Saiz, F.F. Peñailillo 2018: Water: actual
		evapotranspiration with energy balance models. R Package
		Version 0.6 (2017)*\n
		*Tasumi, M., Allen, R.G., Trezza, R., 2008. At-Surface Reflectance
		and Albedo from Satellite for Operational Calculation of Land
		Surface Energy Balance. Journal of Hydrologic Engineering 13, 51–63.*
		https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51).
		"""

		# Constants
		bands = [blue, green, red, nir, swir1, swir2]
		if sat_type == "L8":
			wb = (0.246, 0.146, 0.191, 0.304, 0.105, 0.008)  # Constants for L8 according to Olmeo et. al 2017:
		# (G.F. Olmedo, S. Ortega-Farias, D. Fonseca-Luengo, D. de la Fuente-Saiz, F.F. Peñailillo
		# Water: actual evapotranspiration with energy balance models
		# R Package Version 0.6 (2017))
		else:
			wb = (0.254, 0.149, 0.147, 0.311, 0.103, 0.036)  # Constants according to Tasumi et al. 2008:
		# Tasumi, M., Allen, R.G., Trezza, R., 2008. At-Surface Reflectance and Albedo from Satellite
		# for Operational Calculation of Land Surface Energy Balance.
		# Journal of Hydrologic Engineering 13, 51–63. https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51)

		# Computing of broadband albedo
		alfa_list = []
		for i in range(0, len(bands)):
			alfa_band = bands[i] * wb[i]
			alfa_list.append(alfa_band)
		del alfa_band

		albedo = np.zeros(alfa_list[0].shape)
		i = 0
		while i != len(alfa_list):
			albedo = albedo + alfa_list[i]
			i = i + 1
		del alfa_list

		return albedo

	def reflectRs(self, Rs_in_corr, albedo):
		"""
		Amount of shortwave radiation reflected from surface :math:`(W.m^{-2})`
		"""
		Rs_out = Rs_in_corr * albedo

		return Rs_out

	def netRad(self, Rs_in_corr, Rs_out, RL_in, RL_out):
		"""
		Total net radiation :math:`(W.m^{-2})`
		"""

		Rn = Rs_in_corr - Rs_out + RL_in - RL_out

		return Rn


# noinspection PyMethodMayBeStatic,PyUnusedLocal
class VegIndices:
	"""
	Calculation of vegetation indices from spectral data.
	"""

	def __init__(self):
		return

	def viNDVI(self, red, nir):
		"""
		Normalized Difference Vegetation Index - NDVI.
		"""
		ignore_zero = np.seterr(all="ignore")
		ndvi = (nir - red) / (nir + red)
		ndvi[ndvi == np.inf] = 0  # replacement inf values by 0
		ndvi[ndvi == -np.inf] = 0  # replacement -inf values by 0
		return ndvi

	def viSAVI(self, red, nir, L=0.5):
		"""
		Soil Adjusted Vegetation Index - SAVI.
		"""

		ignore_zero = np.seterr(all="ignore")
		savi = (1 + L) * (nir - red) / (L + nir + red)

		return savi

	def viNDMI(self, nir, swir1):
		"""
		Normalized Vegetation Moisture Index - NDMI.
		"""

		ignore_zero = np.seterr(all="ignore")  # ignoring exceptions with dividing by zero
		ndmi = (nir - swir1) / (nir + swir1)
		ndmi[ndmi == np.inf] = 0  # replacement inf values by 0
		ndmi[ndmi == -np.inf] = 0  # replacement -inf values by 0

		return ndmi

	def viMSAVI(self, red, nir):
		"""
		Modified Soil Adjusted Vegetation Index - SAVI.
		"""

		ignore_zero = np.seterr(all="ignore")  # ignoring exceptions with dividing by zero
		msavi = 0.5 * ((2 * nir + 1) - ((2 * nir + 1) ** 2.0 - 8 * (nir - red)) ** 0.5)
		msavi[msavi == np.inf] = 0  # replacement inf values by 0
		msavi[msavi == -np.inf] = 0  # replacement -inf values by 0

		return msavi

	def fractVegCover(self, ndvi):
		"""
		Function returns fractional vegetation cover layer.
		"""

		Fc = ((ndvi - 0.2) / 0.3) ** 2

		return Fc


# noinspection PyMethodMayBeStatic,PyUnusedLocal
class HeatFluxes:
	"""
	Calculation of heat fluxes and heat balance features from spectral
	and thermal spatial data and meteorological measurements. The class
	contains a set of methods for calculation heat balance, e.g. ground
	heat flux, sensible heat flux, latent heat flux, evaporative fraction,
	omega factor (decoupling coefficient), surface resistance for water
	vapour transfer etc. Calculation for both aerodynamic and gradient
	method is included. 
	"""

	def __init__(self):
		return

	def fluxLE_p(self, Rn, G, delta, VPD, ra, gamma, ro, cp):
		"""
		Latent heat flux for potential evapotranspiration according to Penman (1948) :math:`(W.m^{-2})`.
		"""

		LE_p = (delta * (Rn - G) + ro * cp * VPD / ra) / (delta + gamma)

		return LE_p

	def fluxLE_EQ(self, Rn, G, delta, gamma):
		"""
		Equilibrium evaporation rate :math:`(W.m^{-2})`.
		"""

		LE_eq = delta / (delta + gamma) * (Rn - G)

		return LE_eq

	def fluxLE_PT(self, LE_eq, alpha=1.26):
		"""
		Evaporation from wet surface according to Priestley-Taylor (1972).
		"""

		LE_PT = LE_eq * alpha

		return LE_PT

	def fluxHAer(self, ra, rho, dT, cp=1012):
		"""
		Sensible heat flux :math:`(W.m^{-2})` calculated using aerodynamic method.
		"""

		H = rho * cp * dT / ra

		return H

	def gradH(self, LE, Rn, G):
		"""
		Sensible heat flux :math:`(W.m^{-2})` calculated using gradient method.
		"""

		H = Rn - G - LE

		return H

	def fluxLE(self, Rn, G, H):
		"""
		Latent heat flux :math:`(W.m^{-2})`
		"""

		LE = Rn - G - H

		return LE

	def gradLE(self, EF, Rn, G):
		"""
		Latent heat flux calculated using gradient method :math:`(W.m^{-2})`.
		"""

		LE = EF * (Rn - G)

		return LE

	def aeroEF(self, LE, Rn, G):
		"""
		Evaporative fraction calculated using aerodynamic method (rel.).
		"""

		EF = LE / (Rn - G)

		return EF

	def gradEF(self, ts, ta):
		"""
		Evaporative fraction calculated from gradient method according to Suleiman and Crago (2004).
		"""

		try:
			filt_ts = ndimage.median_filter(ts, 5)
		except ArithmeticError:
			warnings.warn("Median filter has not been used for Tmax estimation", stacklevel=3)
			filt_ts = Ts

		filt_ts = filt_ts[~np.isnan(filt_ts)]
		t_max = np.max(filt_ts)
		EF = (t_max - ts) / (t_max - ta)

		return EF

	def intensityE(self, LE, L):
		"""Evaporation intensity in mmol.m-2.s-1.
		
		:param LE: Latent heat flux :math:`(W.m^{-2})`
		:type LE: numpy.ndarray
		:param L: Latent heat of water evaporation (J.g-1).
		:type L: numpy.ndarray
		
		:returns E_int: Intensity of water evaporation (mmol.m-2.s-1).
		:rtype E_int: numpy.ndarray
		"""

		E_int = LE / L / 18 * 1000

		return E_int

	def omega(self, LE, LE_p):
		"""
		Omega factor (Decoupling coefficient) according to Jarvis and McNaughton (1985)
		"""

		omega = LE / LE_p

		return omega

	def rs(self, delta, gamma, omega, ra):
		"""
		Surface resistance for water vapour transfer (s.m-1)
		"""

		rs = (((delta + gamma) / omega - delta) * 1.0 / gamma - 1.0) * ra

		return rs

	def bowen(self, H, LE):
		"""
		Bowen ratio according to Bowen (1926)
		"""

		ignore_zero = np.seterr(all="ignore")
		bowen = H / LE
		bowen = np.nan_to_num(bowen)

		return bowen

	def rcp(self, E_Z_sat, e_Z, rho, LE_p, ra, gamma, cp):
		"""
		Surface resistance for water vapour transfer for potential evapotranspiration.
		"""

		rcp = (E_Z_sat - e_Z) * rho * cp / (gamma * LE_p) - ra

		return rcp

	def gamma_x(self, rcp, ra, gamma):
		"""
		Psychrometric constant corrected on the rcp and ra according to Jackson et al. (1981, 1988).
		"""

		gamma_x = gamma * (1.0 + rcp / ra)

		return gamma_x

	def cwsi(self, LEp, ra, rc, E_Z_sat, e_Z, rho, delta, gamma, cp=1012):
		"""
		Crop Water Stress Index calculated according to Jackson et al. (1981, 1988).
		"""
		try:
			rcp = self.rcp(E_Z_sat, e_Z, rho, LE_p, ra, gamma, cp)
			g_x = self.gamma_x(rcp, ra, gamma)
			cwsi = 1.0 - (delta + g_x) / (delta + gamma * (1.0 + rc / ra))
		except ArithmeticError:
			if LEp is not None:
				cwsi = np.zeros(LEp.shape)
			else:
				cwsi = None

		return cwsi

	def groundFlux(self, ndvi, Rn, ts, albedo):
		"""
		Ground heat flux according to Bastiaanssen et al. (1998)
		
		"""

		G = ts / albedo * (0.0038 * albedo + 0.0074 * albedo ** 2) * (1 - 0.98 * ndvi ** 4) * Rn

		return G


# noinspection PyMethodMayBeStatic
class MeteoFeatures:
	"""
	Calculation of basic meteorological features and miscellaneous
	variables.
	"""

	def __init__(self):
		return

	def airTemp(self, Z, Z_st, ta_x):
		"""
		Air temperature, recalculated on height Z (m) from data measured in height Z_st (m).
		ta_x - air temperature measured in height Z_st (°C)
		ta - air temperature measured in height Z (°C)
		"""
		ta = ta_x - (float(Z) - float(Z_st)) * 0.0065  # air temperature in C
		return ta

	def airPress(self, Z, DMT):
		"""
		Function returns atmospheric air pressure at level Z in kPa
		"""
		airP = 101.3 * ((293.0 - 0.0065 * (DMT + Z)) / 293.0) ** 5.26
		return airP

	# noinspection PyMethodMayBeStatic
	def satVapourPress(self, ta):
		"""
		Function returns layer of saturated water vapour pressure (kPa) at certain level.
		"""

		E_sat = 0.61121 * np.exp(17.502 * ta / (240.97 + ta))

		return E_sat

	def vapourPress(self, E_sat, Rh):
		"""
		Function returns layer of water vapour pressure in air (kPa) at certain level.
		"""

		e_abs = E_sat * Rh / 100.0

		return e_abs

	def vpd(self, E_sat, e_abs):
		"""
		Function returns water vapour pressure deficit (kPa).
		"""

		VPD = E_sat - e_abs

		return VPD

	def airDensity(self, ta):
		"""
		Function returns volumetric air density (kg.m-3).
		"""

		ro = 353.4 / (ta + 273.0)

		return ro

	def latent(self, ta):
		"""
		Function return latent heat for the water vapour exchange (J.g-1)
		"""

		latent = 2501 - 2.3723 * ta

		return latent

	def gamma(self, airP, latent, cp):
		"""
		Function returns psychrometric constant (kPa.K-1).
		"""

		gamma = cp * airP / (latent * 0.622) * 0.001

		return gamma

	def delta(self, ts, ta):
		"""
		Function calculates the slope of the humidity gradient to temperature gradient - delta function.
		Calculation according to Jackson et al. (1998).
		"""

		T = (ts + ta) / 2
		delta = (45.03 + 3.014 * T + 0.05345 * T ** 2 + 0.00224 * T ** 3) * 0.001
		del T

		return delta

	def emissivity(self, Fc, red, ndvi):
		"""
		Function returns emissivity according to Sobrino et al. (2004) NDVI Treshold Method.
		"""

		emiss = 0.004 * Fc + 0.986
		emiss = np.where(ndvi < 0.2, 1 - red,
		                 emiss)  # replacement of emissivity values for NDVI < 0.2 by values from red band
		emiss = np.where(ndvi > 0.5, 0.99, emiss)  # replacement of emissivity values for NDVI > 0.5 by 0.99
		emiss[emiss > 1] = 0.99  # replacement of values > 1 by 0.99
		emiss[emiss < 0.8] = 0.8  # replacement values < 0.8 by 0.8

		return emiss


# noinspection PyMethodMayBeStatic,PyUnusedLocal
class WindStability(MeteoFeatures, HeatFluxes):
	"""
	Atmospheric stability calculation. Class includes methods for calculation
	of boundary layer stability, friction velocity and aerodynamic
	resistance of the surface. Methods for calculation of stability
	parameters for both aerodynamic and gradient methods are included.
	Methods for SEBAL procedure are also used. Some another features
	are included.
	"""

	def vegHeight(self, h_min, h_max, msavi):
		"""
		height of vegetation cover (m) derived from MSAVI index
		according to Gao et al. (2011).
		"""
		minmsavi = np.min(msavi[msavi != 0])
		maxmsavi = np.max(msavi)
		h_eff = h_min + (msavi - minmsavi) / (minmsavi - maxmsavi) * (h_min - h_max)
		return h_eff

	def zeroPlaneDis(self, h_eff):
		"""
		Zero plane displacement (m).
		"""
		disp = 2.0 / 3.0 * h_eff
		return disp

	def z0m_METRIC(self, ndvi, albedo):
		"""
		Aerodynamic roughness of the surface for momentum transfer (m). 
		Calculated according to pySEBAL (see Jaafar, H.H., Ahmad, F.A.,
		2019. Time series trends of Landsat-based ET using automated
		calibration in METRIC and SEBAL: The Bekaa Valley, Lebanon. 
		Remote Sensing of Environment S0034425718305947.
		https://doi.org/10.1016/j.rse.2018.12.033)
		"""

		n = -5.037
		m = 1.096

		z0m = np.exp(m * ndvi / albedo + n)

		z0m = np.where(z0m < 0.001, 0.001, z0m)

		return z0m

	def z0m(self, h_eff=None, savi=None, ca=0.003, cb=5.26):
		"""
		Aerodynamic roughness of the surface for momentum transfer (m). 
		Calculated according to Thom (1975).
		"""

		if savi is None:
			z0m = h_eff * 0.123
		else:
			z0m = ca * np.exp(cb * savi)

		return z0m

	def z0h(self, z0m):
		"""
		Aerodynamic roughness of the surface for heat transfer (m)
		"""
		z0h = 0.1 * z0m
		return z0h

	def z_d(self, Z, disp):
		"""
		Z - d difference. Z is height of meteorological measurement (m) - Z = 200 m is used.
		"""
		Z_d = Z - disp
		return Z_d

	def windSpeedZ(self, U, Z=200.0, Z_st=2.0, h_st=0.12):
		"""
		Wind speed recalculated to height Z according to logarithmic law
		(Gao et al. 2011).

		Inputs:
		:param U: Wind speed measured on meteostation at level Z_st (m/s).
		:type U: float
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param Z_st: height of wind speed measurement (m). Default 2 m.
		:type Z_st: float
		:param h_st: height of vegetation cover under meteostation (m).
					 Default value is 0.12 m which corresponds with
					 reference cover used for meteostations.
		:type h_st: float

		Returns:
		:return Uz: Wind speed at mixinf layer (m/s)
		:rtype Uz: float
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			Uz = U * np.log(Z / (0.123 * h_st)) / np.log(Z_st / (0.123 * h_st))
		except ArithmeticError:
			raise ArithmeticError("Friction velocity has not been calculated")

		return Uz

	def frictVelo(self, Uz, z0m, Z=200.0, psi_m=0, kappa=0.41):
		"""
		Friction velocity of wind speed (m/s) corrected on atmospheric
		stability.

		Inputs:
		:param Uz: Wind speed at Z level (m/s)
		:type Uz: numpy.ndarray
		:param z0m: Surface roughness for momentum transfer (m)
		:type z0m: numpy.ndarray
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param psi_m: Stability parameter for momentum transfer.
					  Defaul 0.
		:type psi_m: numpy.ndarray
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float
		

		Returns:
		:return frict: Friction velocity (m/s)
		:rtype frict: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			frict = kappa * Uz / (np.log(Z / z0m) - psi_m)
		except ArithmeticError:
			raise ArithmeticError("Friction velocity has not been calculated")

		return frict

	def virtTemp(self, ta, ts, z0h, Z=200, psi_h=0, kappa=0.41):
		"""
		Virtual temperature (K) corrected on atmospheric
		stability.

		Inputs:
		:param ta: Air temperature at Z level (K, C degrees)
		:type ta: numpy.ndarray
		:param ts: Surface temperature (K, C degrees)
		:type ts: numpy.ndarray
		:param z0h: Surface roughness for heat transfer (m)
		:type z0h: numpy.ndarray
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param psi_h: Stability parameter for heat transfer, Default 0.
		:type psi_h: numpy.ndarray
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float
		

		:returns: Virtual temperature
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			t_virt = kappa * (ta - ts) / (np.log(Z / z0h) - psi_h)
		except ArithmeticError:
			raise ArithmeticError("Virtual temperature has not been calculated")

		return t_virt

	def psiM(self, L, X, Z=200, a=1.0, b=0.667, c=5.0, d=0.35):
		"""
		Calculation of stability parameter for momentum transfer (-)
		according to Beljaars et Holstag (1991) for stable conditions
		and Liu et al. (2007) for unstable and neutral conditions.

		Inputs:
		:param L: Monin-Obukhov length (m)
		:type L: numpy.ndarray
		:param X: X coefficient for stability calculation
		:type X: numpy.ndarray
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param a: Coefficient. Default a = 1.0
		:type a: float
		:param b: Coefficient. Default b = 0.667
		:type b: float
		:param c: Coefficient. Defalt c = 5.0
		:type c: float
		:param d: Coefficient. Default d = 0.35
		:type d: float

		Returns:
		:return psi_m: Stability parameter for momentum transfer (-)
		:rtype psi_m: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		dzeta = Z / L

		try:
			psi_m_stab = -(a * dzeta + b * (dzeta - c / d) * np.exp((-d) * dzeta) + b * c / d)
		except ArithmeticError:
			raise ArithmeticError(
				"Stability coefficient for momentum transfer for stable conditions has not been calculated")

		try:
			psi_m_unstab = 2.0 * np.log((1.0 + X) / 2.0) + np.log((1.0 + X ** 2.0) / 2.0) - 2.0 * np.arctan(
				X) + np.pi / 2.0
		except ArithmeticError:
			raise ArithmeticError(
				"Stability coefficient for momentum transfer for unstable conditions has not been calculated")

		psi_m = np.where(dzeta < 0.0, psi_m_unstab, psi_m_stab)

		return psi_m

	def psiH(self, L, X, Z=200, a=1.0, b=0.667, c=5.0, d=0.35):
		"""
		Calculation of stability parameter for heat transfer (-)
		according to Beljaars et Holstag (1991) for stable conditions
		and Liu et al. (2007) for unstable and neutral conditions.

		Inputs:
		:param L: Monin-Obukhov length (m)
		:type L: numpy.ndarray
		:param X: X coefficient for stability calculation
		:type X: numpy.ndarray
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param a: Coefficient. Default a = 1.0
		:type a: float
		:param b: Coefficient. Default b = 0.667
		:type b: float
		:param c: Coefficient. Defalt c = 5.0
		:type c: float
		:param d: Coefficient. Default d = 0.35
		:type d: float

		Returns:
		:return psi_m: Stability parameter for momentum transfer (-)
		:rtype psi_m: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		dzeta = Z / L

		try:
			psi_h_stab = -(
					(1 + 2 * a * dzeta / 3) ** 1.5 + b * (dzeta - c / d) * np.exp((-d) * dzeta) + (b * c / d - 1))
		except ArithmeticError:
			raise ArithmeticError(
				"Stability coefficient for heat transfer for stable conditions has not been calculated")

		try:
			psi_h_unstab = 2.0 * np.log((1.0 + X ** 2.0) / 2.0)
		except ArithmeticError:
			raise ArithmeticError(
				"Stability coefficient for heat transfer for unstable conditions has not been calculated")

		psi_h = np.where(dzeta < 0.0, psi_h_unstab, psi_h_stab)

		return psi_h

	def lengthMO(self, frict, ts, H=None, rho=None, t_virt=None, cp=1012, kappa=0.41, gravit=9.81):
		"""
		Monin-Obukhov length (m)

		Inputs:
		:param frict: Friction velocity (m/s)
		:type frict: numpy.ndarray
		:param ts: Surface temperature (C degree)
		:type ts: numpy.ndarray
		:param H: Sensible heat flux (W.m2)
		:type H: numpy.ndarray
		:param rho: Specific air density (g/m3)
		:type rho: numpy.ndarray
		:param t_virt: Virtual temperature
		:type t_virt: numpy.ndarray
		:param cp: Thermal heat capacity of dry air (K.kg-1.K-1)
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float
		:param gravit: Gravitation forcing (m/s2). Default 9.81
		:type gravit: float

		Returns:
		:return L: Monin-Obukhov length (m)
		:rtype L: numpy.ndarray 
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			if H is None:
				L = frict ** 2.0 * (ts + 273.16) / (kappa * gravit * t_virt)
			else:
				L = -(frict ** 3 * (ts + 273.15) * rho * cp / (kappa * gravit * H))
		except ArithmeticError:
			raise ArithmeticError("Monin-Obukhov lenght has not been calculated")

		return L

	def coefX(self, Z, L):
		"""
		X coefficient for atmospheric stability calculation.

		Inputs:
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param L: Monin-Obukhov length (m)
		:type L: numpy.ndarray

		Returns:
		:return X: X coefficient for atmospheric stability calculation
		:rtype X: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			X = np.where((1.0 - 16.0 * Z / L) ** 0.25 < 0.0, 0.0, (1.0 - 16.0 * Z / L) ** 0.25)
			X = np.nan_to_num(X)
		except ArithmeticError:
			raise ArithmeticError("Coefficient X has not been calculated")

		return X

	def stabCoef(self, Uz, ta, ts, z0m, z0h, Z=200, L=-10000.0, n_iter=10, a=1.0, b=0.667, c=5.0, d=0.35, kappa=0.41,
	             gravit=9.81):
		"""
		Stability parameters calculation using iterative procedure
		described by Itier (1980).

		Inputs:
		:param Uz: Wind speed at level Z (m/s)
		:type Uz: numpy.ndarray
		:param ta: Air temperature at Z level (K, C degrees)
		:type ta: numpy.ndarray
		:param ts: Surface temperature (K, C degrees)
		:type ts: numpy.ndarray
		:param z0m: Surface roughness for momentum transfer (m)
		:type z0m: numpy.ndarray
		:param z0h: Surface roughness for heat transfer (m)
		:type z0h: numpy.ndarray
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float (Numpy array)
		:param L: Initial value of Monin-Obukhov length (m).
				  Default -10000.0
		:type L: float
		:param n_iter: Number of iteration
		:type n_iter: int
		:param a: Coefficient. Default a = 1.0
		:type a: float
		:param b: Coefficient. Default b = 0.667
		:type b: float
		:param c: Coefficient. Defalt c = 5.0
		:type c: float
		:param d: Coefficient. Default d = 0.35
		:type d: float
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float
		:param gravit: Gravitation forcing (m/s2). Default 9.81
		:type gravit: float

		Returns:
		:return psi_m: Stability parameter for momentum transfer (-)
		:rtype psi_m: numpy.ndarray
		:return psi_h: Stability parameter for heat transfer (-)
		:rtype psi_h: numpy.ndarray
		:return frict: Friction velocity (m/s).
		:rtype frict: numpy.ndarray
		:return L: Monin-Obukhov length (m)
		:rtype L: numpy.ndarray 
		"""

		for i in range(n_iter):
			# L outliers
			L = np.where(Z / L < -10000000000.0, ndimage.median_filter(L, 5), L)
			L = np.where(Z / L > 10000000000.0, ndimage.median_filter(L, 5), L)
			L = np.where(Z / L == np.inf, ndimage.median_filter(L, 5), L)

			# Calculation of X coefficient
			X = self.coefX(Z, L)

			# Stability parameters
			psi_m = self.psiM(L, X, Z, a, b, c, d)
			psi_h = self.psiH(L, X, Z, a, b, c, d)

			# Friction velocity
			frict = self.frictVelo(Uz, z0m, Z, psi_m, kappa)

			# Virtual tempetarure
			t_virt = self.virtTemp(ta, ts, z0h, Z, psi_h, kappa)

			# Monin-Obukhov length
			L = self.lengthMO(frict, ts, t_virt=t_virt, kappa=kappa, gravit=gravit)

		return psi_m, psi_h, frict, L

	def raThom(self, Uz, z0m, z0h, psi_m=0.0, psi_h=0.0, Z=200.0, kappa=0.41):
		"""
		Aerodynamic resistance for heat and momentum transfer (s.m-1)
		calculated according to Thom (1975).

		Inputs:
		:param Uz: Wind speed at level Z (m/s)
		:type Uz: numpy.ndarray
		:param z0m: Surface roughness for momentum transfer (m)
		:type z0m: numpy.ndarray
		:param z0h: Surface roughness for heat transfer (m)
		:type z0h: numpy.ndarray
		:param psi_m: Stability parameter for momentum transfer (-).
					  Default is 0.
		:type psi_m: numpy.ndarray
		:param psi_h: Stability parameter for heat transfer (-)
					  Default is 0.
		:type psi_h: numpy.ndarray
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float (Numpy array)
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float

		Returns:
		:return ra: Aerodynamic resistance for heat and momentum
					transfer (s.m-1) calculated according to Thom (1975)
		:rtype ra: numpy.ndarray
		"""

		try:
			ra = (np.log(Z / z0m) - psi_m) * (np.log(Z / z0h) - psi_h) / (kappa ** 2 * Uz)
		except ArithmeticError:
			raise ArithmeticError("Aerodynamic resistance has not been calculated")

		return ra

	def raSEBAL(self, frict, z1=0.1, z2=2.0, psiH_z1=0, psiH_z2=0, kappa=0.41):

		ra = (np.log(z2 / z1) - psiH_z2 + psiH_z1) / (frict * kappa)

		return ra

	def raGrad(self, H, ro, ts, ta, cp=1012):
		"""
		Aerodynamic resistance for heat and momentum transfer (s.m-1) calculated
		from conversion of sensible heat flux equation.
		"""

		ra = ro * cp * (ts - ta) / H
		return ra

	def maxT(self, Rn, G, ra, ta, rho, cp=1012):
		"""
		Maximal surface temperature calculated on physical basis (K).
		"""
		t_max = (Rn - G) * ra / (rho * cp) + ta
		return t_max

	def coef_b(self, t_dry, t_wet, t_max, ta):
		"""
		Coefficient b calculated from temperature gradient.
		"""
		cb = (t_max - ta) / (t_dry - t_wet)  # minimal temperature Tmin is equal to ta + 273.16 (in K)
		return cb

	def coef_a(self, t_wet, cb):
		"""
		Coefficient a calculated from temperature gradient.
		"""
		ca = -cb * (t_wet + 273.16)  # minimal temperature Tmin is equal to ta + 273.16 (in K)
		return ca

	def dT(self, ts, ta, t_wet, t_dry, Rn, G, ra, rho, cp=1012):
		"""
		Temperature gradient calculated according to SEBAL (Bastiaanssen et al. 1998)
		"""
		t_max = self.maxT(Rn, G, ra, ta, rho, cp)
		cb = self.coef_b(t_dry, t_wet, t_max, ta)
		ca = self.coef_a(t_wet, cb)
		dT = ca + cb * (ts + 273.16)

		return dT

	def aeroSEBAL(self, Uz, ta, ts, z0m, Rn, G, rho, t_dry, t_wet, Z=200, z2=2, z1=0.1, niter=10, cp=1012, kappa=0.41):
		"""

		"""

		ignore_zero = np.seterr(all="ignore")

		# friction velocity for meteostation
		z0m_init = np.array(([0.12 * 0.123]), dtype=np.float)
		frict = self.frictVelo(Uz, z0m_init, Z)
		# ra for meteostation (neutral stab)
		ra = self.raSEBAL(frict, z1, z2)

		for i in range(niter):
			diffT = self.dT(ts, ta, t_wet, t_dry, Rn, G, ra, rho, cp)
			H = self.fluxHAer(ra, rho, diffT, cp)
			L = self.lengthMO(frict, ts, H, rho)
			X_z = self.coefX(Z, L)
			X_z1 = self.coefX(z1, L)
			X_z2 = self.coefX(z2, L)
			psiM_z = self.psiM(L, X_z, Z)
			psiH_z1 = self.psiH(L, X_z1, z1)
			psiH_z2 = self.psiH(L, X_z2, z2)
			frict = self.frictVelo(Uz, z0m, Z, psiM_z, kappa)
			ra = self.raSEBAL(frict, z1, z2, psiH_z1, psiH_z2, kappa)  # ra_h - heat transfer

		return ra, H
