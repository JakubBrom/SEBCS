#  /***************************************************************************
#  SEBCS_lib.py
#
#  TODO: popis
#
#                                -------------------
#          begin                :
#          date                 :
#          git sha              : $Format:%H$
#          copyright            : (C) 2014-2019 Jakub Brom
#          email                : jbrom@zf.jcu.cz
#
#  ***************************************************************************/
#  /***************************************************************************
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License  as published  by
#  the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License along
#  with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#   ***************************************************************************/

#  /***************************************************************************
#  SEBCS_lib.py
#
#  Collection of functions for energy balance and their components calculation.
#
#
#  ***************************************************************************/
#  /***************************************************************************
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License  as published  by
#  the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License along
#  with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#   ***************************************************************************/


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

		:return: The affine transformation coefficients.
		:rtype: tuple
		:return: Projection information of the raster (dataset).
		:rtype: str
		:return: Pixel width (m) on X axis.
		:rtype: float
		:return: Pixel height (m) on Y axis.
		:rtype: float
		:return: EPSG Geodetic Parameter Set code.
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
		
		:param arrays: Numpy array or list of arrays for export to raster.
		:type arrays: numpy.ndarray, list
		:param names: Name or list of names of the exported bands (in case\
		of multiband raster) or particular rasters (in case of singleband\
		rasters).
		:type names: str, list
		:param prj: Projection information of the exported raster (dataset).
		:type prj: str
		:param gtransf: The affine transformation coefficients.
		:type gtransf: tuple
		:param EPSG: EPSG Geodetic Parameter Set code.
		:type EPSG: int
		:param out_folder: Path to folder where the raster(s) will be created.
		:type out_folder: str
		:param driver_name: GDAL driver. 'GTiff' is default.
		:type driver_name: str
		:param out_file_name: Name of exported multiband raster.
		:type out_file_name: str
		:param multiband: Option of multiband raster creation.
		:type multiband: bool

		:return: Raster singleband or multiband file(s)
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

	@staticmethod
	def rasterToArray(layer):
		"""Conversion of raster layer to numpy array.

		:param layer: Path to raster layer.
		:type layer: str

		:return: raster file converted to numpy array
		:rtype: numpy.ndarray
		"""

		lyr_name = os.path.split(layer)[1]

		try:
			if layer is not None or layer is not "":
				try:
					new_array = gdal.Dataset.ReadAsArray(gdal.Open(layer)).astype(
						np.float32)
					new_array = np.nan_to_num(new_array)
				except:
					new_array = None
			else:
				new_array = None

		except IOError:
			warnings.warn("Layer {lr} has not been readed. No data will be "
			              "used instead".format(lr=lyr_name), stacklevel=3)
			new_array = None

		return new_array

	def lyrsExtent(self, in_lyrs_list):
		"""
		Check differences between size of the input layers. Function compares
		rasters only on basis of number of columns and rows.

		:param in_lyrs_list: List of layers (Numpy 2D arrays).
		:type in_lyrs_list: list

		:return: True
		"""
		# TODO: Solution of different spatial extent of rasters...

		in_lyrs_list_true = [i for i in in_lyrs_list if
		                     i is not None]  # List of input layers without
		# empty (None)

		len_list = []

		for i in in_lyrs_list_true:
			len_list.append(len(i))

		if max(len_list) != min(len_list):
			raise ValueError(
				"Selected layers differ in spatial extent (number of columns "
				"or rows).")

		return


# noinspection PyMethodMayBeStatic,PyUnusedLocal,SpellCheckingInspection,SpellCheckingInspection
class VegIndices:
	"""
	Calculation of vegetation indices from spectral data.
	"""

	def __init__(self):
		return

	def viNDVI(self, red, nir):
		"""
		Normalized Difference Vegetation Index - NDVI.

		:param red: Spectral reflectance in RED region (rel.)
		:type red: numpy.ndarray
		:param nir: Spectral reflectance in NIR region (rel.)
		:type nir: numpy.ndarray

		:return: Normalized Difference Vegetation Index - NDVI (unitless)
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			ndvi = (nir - red) / (nir + red)
			ndvi = np.where(ndvi == np.inf, 0, ndvi)  # replacement inf values
			# by 0
			ndvi = np.where(ndvi == -np.inf, 0, ndvi)  # replacement -inf
		# values by 0
		except ArithmeticError:
			raise ArithmeticError("NDVI has not been calculated.")

		return ndvi

	# noinspection SpellCheckingInspection
	def viSAVI(self, red, nir, L=0.5):
		# noinspection SpellCheckingInspection
		"""
				Soil Adjusted Vegetation Index - SAVI (Huete, 1988).

				:param red: Spectral reflectance in RED region (rel.)
				:type red: numpy.ndarray
				:param nir: Spectral reflectance in NIR region (rel.)
				:type nir: numpy.ndarray
				:param L: Parameter L. Default L=0.5
				:type L: float

				:return: Soil Adjusted Vegetation Index - SAVI (unitless)
				:rtype: numpy.ndarray

				\n
				**References**\n
				*Huete A.R. (1988): A soil-adjusted vegetation index (SAVI) Remote
				Sensing of Environment 27, 47-57.*
				"""

		ignore_zero = np.seterr(all="ignore")
		
		try:
			savi = (1 + L) * (nir - red) / (L + nir + red)
		except ArithmeticError:
			raise ArithmeticError("SAVI has not been calculated.")

		return savi

	def viOSAVI(self, red, nir, L=0.16):
		"""
		Optimized Soil Adjusted Vegetation Index - OSAVI (Rondeaux et al. (
		1996)).

		:param red: Spectral reflectance in RED region (rel.)
		:type red: numpy.ndarray
		:param nir: Spectral reflectance in NIR region (rel.)
		:type nir: numpy.ndarray
		:param L: Parameter L. Default L=0.5
		:type L: float

		:return: Soil Adjusted Vegetation Index - OSAVI (unitless)
		:rtype: numpy.ndarray

		\n
		**References**\n
		*Rondeaux G., Steven M., Baret F. (1996): Optimisation of
		soil-adjusted vegetation indices Remote Sensing of Environment,
		55 (1996), pp. 95-107*
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			savi = (1 + L) * (nir - red) / (L + nir + red)
		except ArithmeticError:
			raise ArithmeticError("SAVI has not been calculated.")

		return savi

	def viNDMI(self, nir, swir1):
		"""
		Normalized Vegetation Moisture Index - NDMI.

		:param nir: Spectral reflectance in NIR region (rel.)
		:type nir: numpy.ndarray
		:param swir1: Spectral reflectance in SWIR region (approx. 1.61\
		:math:`\\mu{m}` (rel.)
		:type nir: numpy.ndarray

		:return: Normalized Vegetation Moisture Index - NDMI (unitless)
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")  # ignoring exceptions with dividing by zero

		if swir1 is not None:
			try:
				ndmi = (nir - swir1) / (nir + swir1)
				ndmi[ndmi == np.inf] = 0  # replacement inf values by 0
				ndmi[ndmi == -np.inf] = 0  # replacement -inf values by 0
			except ArithmeticError:
				raise ArithmeticError("NDMI has not been calculated.")
		else:
			ndmi = None

		return ndmi

	def viMSAVI(self, red, nir):
		"""
		Modified Soil Adjusted Vegetation Index - MSAVI (Qi et al. 1994).

		:param red: Spectral reflectance in RED region (rel.)
		:type red: numpy.ndarray
		:param nir: Spectral reflectance in NIR region (rel.)
		:type nir: numpy.ndarray

		:return: Modified Soil Adjusted Vegetation Index - SAVI (unitless)
		:rtype: numpy.ndarray

		\n
		**References**\n
		*Qi, J., Chehbouni, A., Huete, A.R., Kerr, Y.H., Sorooshian, S.,
		1994. A modified soil adjusted vegetation index. Remote Sensing of
		Environment 48, 119–126. https://doi.org/10.1016/0034-4257(94)90134-1*
		"""

		ignore_zero = np.seterr(all="ignore")  # ignoring exceptions with dividing by zero

		try:
			msavi = 0.5 * ((2 * nir + 1) - ((2 * nir + 1) ** 2.0 - 8 * (nir - red)) ** 0.5)
			msavi[msavi == np.inf] = 0  # replacement inf values by 0
			msavi[msavi == -np.inf] = 0  # replacement -inf values by 0
		except ArithmeticError:
			raise ArithmeticError("MSAVI has not been calculated.")

		return msavi

	def fractVegCover(self, ndvi):
		"""
		Fractional vegetation cover layer - Fc.

		:param ndvi: Normalized Difference Vegetation Index - NDVI (unitless)
		:type ndvi: numpy.ndarray

		:return: Fractional vegetation cover layer - Fc (unitless)
		:rtype: numpy.ndarray
		"""

		try:
			Fc = ((ndvi - 0.2) / 0.3) ** 2
		except ArithmeticError:
			raise ArithmeticError("Fractional cover index has not been "
			                      "calculated.")

		return Fc

	def LAI(self, red, nir, method=3):
		"""
		Leaf Area Index (LAI) calculated according to several methods.

		:param red: Spectral reflectance in RED region (rel.)
		:type red: numpy.ndarray
		:param nir: Spectral reflectance in NIR region (rel.)
		:type nir: numpy.ndarray
		:param method: Method of LAI calculation:\n
				* 1: Pôças
				* 2: Bastiaanssen
				* 3: Jafaar (default)
				* 4: Anderson
				* 5: vineyard
				* 6: Carrasco
				* 7: Turner
		:type method: int

		:return: Leaf Area Index (LAI) :math:`(m^2.m^{-2})`
		:rtype: numpy.ndarray

		\n
		**References**\n
		*Anderson, M., Neale, C., Li, F., Norman, J., Kustas, W., Jayanthi,
		H., Chavez, J., 2004. Upscaling ground observations of vegetation
		water content, canopy height, and leaf area index during SMEX02 using
		aircraft and Landsat imagery. Remote Sensing of Environment 92,
		447–464. https://doi.org/10.1016/j.rse.2004.03.019* \n
		*Bastiaanssen, W.G.M., Menenti, M., Feddes, R.A., Holtslag, A.A.M.,
		1998. A remote sensing surface energy balance algorithm for land (
		SEBAL). 1. Formulation. Journal of Hydrology 212–213, 198–212.
		https://doi.org/10.1016/S0022-1694(98)00253-4* \n
		*Carrasco-Benavides, M., Ortega-Farías, S., Lagos, L., Kleissl, J.,
		Morales-Salinas, L., Kilic, A., 2014. Parameterization of the
		Satellite-Based Model (METRIC) for the Estimation of Instantaneous
		Surface Energy Balance Components over a Drip-Irrigated Vineyard.
		Remote Sensing 6, 11342–11371. https://doi.org/10.3390/rs61111342* \n
		*Jaafar, H.H., Ahmad, F.A., 2019. Time series trends of Landsat-based
		ET using automated calibration in METRIC and SEBAL: The Bekaa
		Valley, Lebanon. Remote Sensing of Environment S0034425718305947.
		https://doi.org/10.1016/j.rse.2018.12.033* \n
		*Pôças, I., Paço, T.A., Cunha, M., Andrade, J.A., Silvestre, J.,
		Sousa, A., Santos, F.L., Pereira, L.S., Allen, R.G., 2014.
		Satellite-based evapotranspiration of a super-intensive olive
		orchard:  Application of METRIC algorithms. Biosystems Engineering
		128, 69–81. https://doi.org/10.1016/j.biosystemseng.2014.06.019* \n
		*Turner, D.P., Cohen, W.B., Kennedy, R.E., Fassnacht, K.S., Briggs,
		J.M., 1999. Relationships between Leaf Area Index and Landsat TM
		Spectral Vegetation Indices across Three Temperate Zone Sites.
		Remote Sensing of Environment 70, 52–68.
		https://doi.org/10.1016/S0034-4257(99)00057-7*
		"""

		savi = self.viSAVI(red, nir, 0.5)
		ndvi = self.viNDVI(red, nir)
		osavi = self.viOSAVI(red, nir)
		msavi = self.viMSAVI(red, nir)
		# methods = {1:'Pôças', 2:'Bastiaanssen', 3:'Jafaar', 4:'Anderson',
		# 5:'vineyard', 6:'Carrasco', 7:'Turner'}

		if method is 1:     # Pôças
			LAI = np.where(savi > 0, 11.0 * savi**3, 0)
			LAI = np.where(savi > 0.817, 6, LAI)

		elif method is 2:   # Bastiaanssen
			LAI = np.where(savi > 0, np.log((0.61 - savi)/0.51)/0.91 * (-1), 0)
			LAI = np.where(savi >= 0.61, 6, LAI)

		elif method is 3:   # Jafaar
			LAI_1 = np.where(savi > 0, 11.0 * savi ** 3, 0)
			LAI_1 = np.where(savi > 0.817, 6, LAI_1)

			LAI_2 = np.where(savi > 0, np.log((0.61 - savi) / 0.51) / 0.91 * (
				-1), 0)
			LAI_2 = np.where(savi >= 0.61, 6, LAI_2)

			LAI = (LAI_1 + LAI_2)/2

		elif method is 4:   # Anderson
			LAI = (4 * osavi - 0.8) * (1 + 4.73e-6 * np.exp(15.64 * osavi))

		elif method is 5:   # vineyard
			LAI = 4.9 * ndvi - 0.46

		elif method is 6:   # Carrasco
			LAI = 1.2 - 3.08 * np.exp(-2013.35 * ndvi**6.41)

		elif method is 7:   # Turner
			LAI = 0.5724 + 0.0989 * ndvi - 0.0114 * ndvi**2 + 0.0004 * ndvi**3

		# elif method is 8:
		# 	LAI = 6.0/(1 + np.exp(-(8 * savi - 5)))
		# 	# Method proposed by Brom (mathematical definition only,
		# 	not tested). Good approximation with another methods (1, 2, 3) but
		# 	this method	is very sensitive on parameters in exponent
		else:
			raise ArithmeticError("LAI has not been calculated.")

		LAI = np.where(LAI <= 0, 0, LAI)

		return LAI

	def vegHeight(self, h_min, h_max, msavi):
		"""
		Height of effective vegetation cover (m) derived from MSAVI index
		according to Gao et al. (2011).

		:param h_min: Maximal height of vegetation cover (m)
		:type h_max: numpy.ndarray
		:param h_max: Minimal height of vegetation cover (m)
		:type h_min: numpy.ndarray
		:param msavi: Modified Soil Adjusted Vegetation Index (MSAVI)
		:type msavi: numpy.ndarray

		:return: Effective vegetation cover height (m)
		:rtype: numpy.ndarray

		\n
		**References**\n
		*Gao, Z.Q., Liu, C.S., Gao, W., Chang, N.-B., 2011. A coupled remote
		sensing and the Surface Energy Balance with Topography Algorithm
		(SEBTA) to estimate actual evapotranspiration over heterogeneous
		terrain. Hydrol. Earth Syst. Sci. 15, 119–139.
		https://doi.org/10.5194/hess-15-119-2011*
		"""

		try:
			minmsavi = np.min(msavi[msavi != 0])
			maxmsavi = np.max(msavi)
			h_eff = h_min + (msavi - minmsavi) / (minmsavi - maxmsavi) * (h_min - h_max)

			return h_eff

		except ArithmeticError:
			raise ArithmeticError("Vegetation height has not been calculated.")


# noinspection PyMethodMayBeStatic
class SolarRadBalance(VegIndices):
	"""
	Class contains functions for calculation of solar radiation balance
	and topographic features: incident radiation in dependence on surface
	geometry, slope of terrain, aspect of terrain, albedo, longwave radiation
	fluxes and atmospheric emissivity, shortwave radiation reflectance
	and total net radiation
	"""

	def slopeAspect(self, DMT, x_size, y_size):
		"""
		Slope and aspect of terrain (DMT) in degrees.

		:param DMT: Digital model of terrain (m a.s.l.)
		:type DMT: numpy.ndarray
		:param x_size: Size of pixel in x axis (m)
		:type x_size: float
		:param y_size: Size of pixel in y axis (m)
		:type y_size: float

		:return: Slope of the terrain :math:`(\SI{}\degree)`
		:rtype: numpy.ndarray
		:return: Aspect of the terrain :math:`(\SI{}\degree)`
		:rtype: numpy.ndarray
		"""

		try:
			x, y = np.gradient(DMT)
			slope = np.arctan(np.sqrt((x / x_size) ** 2.0 + (y / y_size) **
			                          2.0)) * 180 / np.pi
			aspect = 270 + np.arctan(x / y) * 180 / np.pi
			aspect = np.where(y > 0, aspect, aspect - 180)

			# Replacing nan values to 0 and inf to value
			slope = np.nan_to_num(slope)
			aspect = np.nan_to_num(aspect)
			del x
			del y
		except ArithmeticError:
			raise ArithmeticError("Slope and aspect has not been calculated.")

		return slope, aspect

	def solarInTopo(self, Rs_in, slope, aspect, latitude, longitude,
	                date_acq, time_acq):

		"""Calculation of incident shortwave solar radiation flux
		according to the solar geometry, position (latitude
		and longitude) and shape of surface (slope and orientation).
		Flux of the solar energy :math:`(W.m^{-2})` is calculated on basis
		of the measured global radiation using pyranometer (incomming
		global radiation). Diffuse part of radiation is not separated
		in calculation.

		:param Rs_in: Global radiation measured by pyranometer\
		:math:`(W.m^{-2})`.
		:type Rs_in: float
		:param slope: Slope of the terrain :math:`(\SI{}\degree)`.
		:type slope: numpy.ndarray
		:param aspect: Orientation of the terrain :math:`(\SI{}\degree)`.
		:type aspect: numpy.ndarray
		:param latitude: Mean latitude of the data in decimal degrees
		:type latitude: float
		:param longitude: Mean longitude of the data in decimal degrees
		:type longitude: float
		:param date_acq: Date of data acquisition in iso format ('YYYY-mm-dd')
		:type date_acq: datetime.date
		:param time_acq: Time in GMT in datetime.time format ('HH:MM:SS.SS')
		:type time_acq: datetime.time

		:returns: Incident shortwave radiation :math:`(W.m^{-2})` corrected\
		on the terrain and solar geometry.
		:rtype: numpy.ndarray
		"""

		# Date conversions
		try:
			dat = datetime.strptime(str(date_acq),
			                        "%Y-%m-%d")  # conversion of iso date to datetime.date
			t = datetime.strptime(str(time_acq),
			                      "%H:%M:%S.%f")  # conversion of iso time to datetime.time
			dec_t = float(t.hour) + float(t.minute) / 60.0 + float(
				t.second) / 3600  # decimal time
			N = dat.strftime("%j")  # day of the year
		except ArithmeticError:
			raise ArithmeticError("Date transformation has not been done.")

		# Solar radiation geometry calculation
		try:
			solar_time = dec_t + longitude / 360.0 * 24.0  # solar time
			hs = (12.0 - float(solar_time)) * 15.0  # local solar hour angle (°)
			ds = 23.45 * math.sin(360.0 * (284.0 + float(
				N)) / 365.0 * math.pi / 180.0)  # solar declination (°)
			ds_rad = math.radians(ds)  # solar declination (rad.)
			L_rad = math.radians(latitude)  # latitude (rad.)
			hs_rad = math.radians(hs)  # local solar hour angle (rad.)

			sin_alpha = (math.sin(L_rad) * math.sin(ds_rad) + math.cos(L_rad)
			             * math.cos(ds_rad) * math.cos(
						hs_rad))  # sin of solar height angle

			if sin_alpha < 0.0:
				sin_alpha = 0.0

			slope_rad = np.radians(slope)
			asp_rad = np.radians((aspect - 180) * (-1))  # aspect transformation

			cosi = (np.sin(ds_rad) * (np.sin(L_rad) * np.cos(slope_rad)
			                          - np.cos(L_rad) * np.sin(
						slope_rad) * np.cos(asp_rad))
			        + np.cos(ds_rad) * np.cos(hs_rad) * (np.cos(L_rad)
			                                             * np.cos(
								slope_rad) + np.sin(L_rad) * np.sin(slope_rad)
			                                             * np.cos(
								asp_rad)) + np.cos(ds_rad) * np.sin(slope_rad)
			        * np.sin(asp_rad) * np.sin(hs_rad))
		except ArithmeticError:
			raise ArithmeticError("Solar geometry corrections have not been "
			                      "done.")

		# Direct radiation intensity correction on DEM
		try:
			Is = float(
				Rs_in) / sin_alpha  # radiation perpendicular to solar beam angle
			Rs_in_corr = Is * cosi  # radiation corrected on solar and terrain geometry
		except ArithmeticError:
			raise ArithmeticError("Radiation intensity correction on DEM has "
			                      "not been done.")

		return Rs_in_corr

	def atmEmissivity(self, e_Z, ta):
		"""
		Atmospheric emissivity calculated according to Idso (see
		Brutsaert 1982).

		:param e_Z: Atmospheric water vapour pressure (kPa)
		:type e_Z: numpy.ndarray, float
		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float

		:return: Air emissivity (rel.)
		:rtype: numpy.ndarray, float
		"""

		try:
			emis_a = 1.24 * (e_Z * 10.0 / (ta + 273.16)) ** (1.0 / 7.0)
		except ArithmeticError:
			raise ArithmeticError("Air emissivity has not been calculated.")

		return emis_a

	def downRL(self, ta, emis_a):
		"""
		Funtion calculates downward flux of longwave radiation :math:`(W.m^{
		-2})`

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
			raise ArithmeticError("Downward longwave radiation flux has not "
			                      "been calculated.")

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

		try:
			RL_out = emiss * 5.6703 * 10.0 ** (-8.0) * (ts + 273.16) ** 4

		except ArithmeticError:
			raise ArithmeticError("Upward longwave radiation flux has not "
			                      "been calculated.")

		return RL_out

	def albedoBrom(self, ndvi, msavi, c_a=0.08611, c_b=0.894716, c_c=5.558657,
	               c_d=-0.11829, c_e=-1.9818, c_f=-4.50339, c_g=-11.4625,
	               c_h=7.461454, c_i=5.299396, c_j=4.76657, c_k=-2.3127,
	               c_l=-3.42739):
		"""
		Albedo (rel.) calculated according to Duffková and Brom et al. (2012)

		:param ndvi: Normalized Difference Vegetation Index (-)
		:type ndvi: numpy.array
		:param msavi: Modified Soil Adjusted Vegetation Index (-) according\
		to Gao et al. 1996.
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

		\n
		**References:**\n
		*Duffková, R., Brom, J., Žížala, D., Zemek, F., Procházka, J.,
		Nováková, E., Zajíček, A., Kvítek, T., 2012. Určení infiltračních
		oblastí pomocí vodního stresu vegetace na základě dálkového průzkumu
		Země a pozemních měření. Certifikovaná metodika. VÚMOP, v.v.i., Praha.*
		"""

		try:
			albedo = (c_a + c_b * msavi + c_c * msavi ** 2 + c_d * ndvi + c_e
			          * msavi ** 3 + c_f * msavi * ndvi + c_g * msavi ** 2
			          * ndvi + c_h * msavi * ndvi ** 2 + c_i * msavi ** 2
			          * ndvi ** 2 + c_j * msavi ** 3 * ndvi + c_k * msavi ** 3
			          * ndvi ** 2 + c_l * msavi * ndvi ** 3)

		except ArithmeticError:
			raise ArithmeticError("Albedo has not been calculated.")

		return albedo

	def albedoLandsat(self, blue, green, red, nir, swir1, swir2, sat_type="L8"):
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
		:param sat_type: Type of Landsat satellite: \n\n
				- L5 - Landsat 4 TM, 5 TM or Landsat 7 ETM+\n
				- L8 - Landsat 8 OLI/TIRS
		:type sat_type: str

		:return: Albedo (rel.)
		:rtype: numpy.ndarray

		\n
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
			wb = (0.246, 0.146, 0.191, 0.304, 0.105, 0.008)  # Constants for
		# L8 according to Olmeo et. al 2017: (G.F. Olmedo, S. Ortega-Farias,
		# D. Fonseca-Luengo, D. de la Fuente-Saiz, F.F. Peñailillo
		# Water: actual evapotranspiration with energy balance models
		# R Package Version 0.6 (2017))

		else:
			wb = (0.254, 0.149, 0.147, 0.311, 0.103, 0.036)  # Constants
		# according to Tasumi et al. 2008: Tasumi, M., Allen, R.G., Trezza,
		# R.,  2008. At-Surface Reflectance and Albedo from Satellite
		# for Operational Calculation of Land Surface Energy Balance.
		# Journal of Hydrologic Engineering 13, 51–63.
		# https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51)

		# Computing of broadband albedo
		try:
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
		except ArithmeticError:
			raise ArithmeticError("Albedo has not been calculated.")

		return albedo

	def albedo(self, band_red, band_nir, sat_type="L8", band_blue=None,
	           band_green=None, band_sw1=None, band_sw2=None):
		"""
		Calculation of Albedo according to data type (satellite data type)
		or data availability. Albedo can be calculated using Landsat data 
		or using any data including RED and NIR band. For the Landsat 8 data 
		Olmedo method is used, for tle Landsat 4, 5 and 7 Tasumi approach is 
		used. If only RED and NIR bands are available, Brom method is used.
		
		:param band_red: Red band (rel.)
		:type band_red: numpy.ndarray
		:param band_nir: NIR band (rel.)
		:type band_nir: numpy.ndarray
		:param sat_type: Type of Landsat satellite: \n\n
				- L5 - Landsat 4 TM, 5 TM or Landsat 7 ETM+\n
				- L8 - Landsat 8 OLI/TIRS
		:type sat_type: str
		:param band_blue: Blue band (rel.)
		:type band_blue: numpy.ndarray
		:param band_green: Green band (rel.)
		:type band_green: numpy.ndarray
		:param band_sw1: SWIR1 band on ca 1.61 :math:`\mu m` (rel.)
		:type band_sw1: numpy.ndarray
		:param band_sw2: SWIR2 band on ca 2.2 :math:`\mu m` (rel.)
		:type band_sw2: numpy.ndarray

		:return: Albedo (rel.)
		:rtype: numpy.ndarray
		
		\n
		**References:**\n
		*G.F. Olmedo, S. Ortega-Farias, D. Fonseca-Luengo,
		D. de la Fuente-Saiz, F.F. Peñailillo 2018: Water: actual
		evapotranspiration with energy balance models. R Package
		Version 0.6 (2017)*\n
		*Tasumi, M., Allen, R.G., Trezza, R., 2008. At-Surface Reflectance
		and Albedo from Satellite for Operational Calculation of Land
		Surface Energy Balance. Journal of Hydrologic Engineering 13, 51–63.
		https://doi.org/10.1061/(ASCE)1084-0699(2008)13:2(51).*\n
		*Duffková, R., Brom, J., Žížala, D., Zemek, F., Procházka, J.,
		Nováková, E., Zajíček, A., Kvítek, T., 2012. Určení infiltračních
		oblastí pomocí vodního stresu vegetace na základě dálkového průzkumu
		Země a pozemních měření. Certifikovaná metodika. VÚMOP, v.v.i., Praha.*
		"""

		try:
			if sat_type == "other" or band_blue is None or band_green is None\
					or band_sw1 is None or band_sw2 is None:
				ndvi = self.viNDVI(band_red, band_nir)
				msavi = self.viMSAVI(band_red, band_nir)
				albedo = self.albedoBrom(ndvi, msavi)

			else:
				albedo = self.albedoLandsat(band_blue, band_green, band_red,
				                            band_nir, band_sw1, band_sw2,
				                            sat_type)
		except ArithmeticError:
			raise ArithmeticError("Albedo has not been calculated.")

		return albedo

	def reflectRs(self, Rs_in_corr, albedo):
		"""
		Amount of shortwave radiation reflected from surface :math:`(W.m^{-2})`

		:param Rs_in_corr: Incomming global radiation corrected on DEM\
		:math:`(W.m^{-2})`
		:type Rs_in_corr: numpy.ndarray
		:param albedo: Surface albedo (rel.)
		:type albedo: numpy.ndarray
		:return: Amount of reflected global radiation :math:`(W.m^{-2})`
		:rtype: numpy.ndarray
		"""
		try:
			Rs_out = Rs_in_corr * albedo
		except ArithmeticError:
			raise ArithmeticError("Amount of reflected shortwave radiation "
			                      "has not been calculated.")

		return Rs_out

	def netRad(self, Rs_in_corr, Rs_out, RL_in, RL_out):
		"""
		Total net radiation.

		:param Rs_in_corr: Incomming global (shortwave) radiation\
		:math:`(W.m^{-2})`
		:type Rs_in_corr: numpy.ndarray
		:param Rs_out: Outgoing (reflected) shortwave radiation\
		:math:`(W.m^{-2})`
		:type Rs_out: numpy.ndarray
		:param RL_in: Incomming (downward) longwave radiation :math:`(W.m^{-2})`
		:type RL_in: numpy.ndarray
		:param RL_out: Outgoing (upward) longwave radiation :math:`(W.m^{-2})`
		:type RL_out: numpy.ndarray

		:return: Total net radiation flux :math:`(W.m^{-2})`
		:rtype: numpy.ndarray
		"""

		try:
			Rn = Rs_in_corr - Rs_out + RL_in - RL_out
		except ArithmeticError:
			raise ArithmeticError("Total net radiation has not been "
			                      "calculated.")

		return Rn


# noinspection PyMethodMayBeStatic
class MeteoFeatures(VegIndices):
	"""
	Calculation of basic meteorological features and miscellaneous
	variables.
	"""

	def __init__(self):
		return

	def airTemp(self, ta_st, st_altitude, DMT, adiabatic=0.0065):
		"""
		Conversion of air temperature calculated from temperature data and
		digital	elevation model :math:`(\SI{}\degreeCelsius)` to spatial
		near to ground temperature or temperature at different altitude.

		:param ta_st: Air temperature measured on meteostation :math:`(\SI{
		}\degreeCelsius)`
		:type ta_st: float
		:param st_altitude: Meteostation altitude (m a.s.l.)
		:type st_altitude: float
		:param DMT: Digital elevation model (m) or altitude for calculation
		of ta value.
		:type DMT: numpy.ndarray, float
		:param adiabatic: Adiabatic lapse rate :math:`(\SI{}\degreeCelsius)`.
		Default = 0.0065 :math:`\SI{}\degreeCelsius.m^{-1}`
		:type adiabatic: float
		:return: Spatial air temperature :math:`(\SI{}\degreeCelsius)`
		:rtype: numpy.ndarray, float
		"""

		try:
			ta = ta_st - (DMT - float(st_altitude)) * float(adiabatic)
		except ArithmeticError:
			raise ArithmeticError("Spatial air temperature has not been "
			                      "calculated.")

		return ta

	def airTemperatureBlending(self, ta_surface, Z=200.0, Z_st=2.0,
	                           adiabatic=0.0065):
		"""
		Conversion of spatial air temperature measured for near surface
		level to level Z above the surface. This approach is usually used
		for calculation of air temperature at blending height (mixing layer).

		:param ta_surface: Spatial air temperature for near ground level (
		Z_st; :math:`\SI{}\degreeCelsius)`
		:type ta_surface:
		:param Z: Height of a level which is used for calculation. Default
		value corresponds with "blending" height or height of mixing layer.
		Default is 200 m above surface.
		:type Z: float
		:param Z_st: Height of air temperature measurement (m)
		:type Z_st: float
		:param adiabatic: Adiabatic lapse rate :math:`(\SI{}\degreeCelsius)`.
		Default = 0.0065 :math:`\SI{}\degreeCelsius.m^{-1}`
		:type adiabatic: float
		:return:
		"""

		try:
			ta = ta_surface - (Z - Z_st) * adiabatic
		except ArithmeticError:
			raise ArithmeticError("Air temperature for level Z has not been "
			                      "calculated.")

		return ta

	def surfaceTemperature(self, tir_band, emissivity=1.0, emis_rule="No"):
		"""Correction of surface temperature on emissivity.
		:param tir_band: Layer of surface temperature :math:`(\SI{
		}\degreeCelsius)`.
		:type tir_band: numpy.ndarray
		:param emissivity: Layer of emissivity (rel.).
		:type emissivity: numpy.ndarray
		:param emis_rule: Setting if the correction will be done or not.
			No is default.
		:type emis_rule: str

		:return: Layer of corrected surface temperature :math:`(\SI{
		}\degreeCelsius)`.
		:rtype: numpy.ndarray
		"""
		try:
			if emis_rule is not "No":
				ts = ((tir_band + 273.16) / emissivity ** 0.25) - 273.16
			else:
				ts = tir_band
		except ValueError:
			raise ValueError("Surface temperature has not been calculated.")

		return ts

	def airPress(self, ta, DMT, Z=200.0, P0=101.325, adiabatic=0.0065):
		"""
		Atmospheric air pressure at level Z (kPa)

		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float
		:param DMT: Digital elevation model or altitude (m)
		:type DMT: numpy.ndarray, float
		:param Z: Height of measurement above surface (m)
		:type Z: float
		:param P0: Sea level air pressure (kPa)
		:type P0: float
		:param adiabatic: Adiabatic lapse rate :math:`(\SI{}\degreeCelsius)`
		:type adiabatic: float

		:return: Air pressure at level Z (kPa)
		:rtype: numpy.ndarray, float
		"""

		try:
			airP = P0 * ((ta + 273.15) / (ta + 273.15 + adiabatic * (DMT +
			                                                         Z))) ** 5.257

		except ArithmeticError:
			raise ArithmeticError("Air pressure has not been calculated.")

		return airP

	def satVapourPress(self, ta):
		"""
		Saturated water vapour pressure (kPa) calculated using Magnus-Tetens
		equation.

		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float

		:return: Saturated water vapour pressure (kPa)
		:rtype: numpy.ndarray, float
		"""

		try:
			E_sat = 0.61121 * np.exp(17.502 * ta / (240.97 + ta))
		except ArithmeticError:
			raise ArithmeticError("Saturated water vapour pressure has not "
			                      "been calculated.")

		return E_sat

	def vapourPress(self, E_sat, Rh):
		"""
		Water vapour pressure in air (kPa).

		:param E_sat: Saturated water vapour pressure (kPa)
		:type E_sat: numpy.ndarray, float
		:param Rh: Relative humidity of air (%)
		:type Rh: numpy.ndarray, float

		:return: Water vapour pressure in air (kPa)
		:rtype: numpy.ndarray, float
		"""

		try:
			e_abs = E_sat * Rh / 100.0
		except ArithmeticError:
			raise ArithmeticError("Water vapour pressure has not been "
			                      "calculated")

		return e_abs

	def vpd(self, E_sat, e_abs):
		"""
		Water vapour pressure deficit (kPa).

		:param E_sat: Saturated water vapour pressure (kPa)
		:type E_sat: numpy.array, float
		:param e_abs: Water vapour pressure in air (kPa)
		:type e_abs: numpy.array, float

		:return: Water vapour pressure deficit (kPa)
		:rtype: numpy.ndarray, float
		"""

		try:
			VPD = E_sat - e_abs
		except ArithmeticError:
			raise ArithmeticError("Water vapour pressure deficit has not "
			                      "been calculated.")

		return VPD

	def airDensity(self, ta):
		"""
		Volumetric dry air density :math:`(kg.m^{-3})`.

		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float

		:return: Volumetric dry air density :math:`(kg.m^{-3})`
		:rtype: numpy.ndarray, float
		"""

		try:
			rho = 353.4 / (ta + 273.15)
		except ArithmeticError:
			raise ArithmeticError("Volumetric dry air density has not been "
			                      "calculated")

		return rho

	def latent(self, ta):
		"""
		Latent heat for the water vapour exchange :math:`(J.g^{-1})`.

		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float

		:return: Latent heat for the water vapour exchange :math:`(J.g^{-1})`
		:rtype: numpy.ndarray, float
		"""

		try:
			latent = 2501 - 2.3723 * ta
		except ArithmeticError:
			raise ArithmeticError("Latent heat has not been calculated.")

		return latent

	def gamma(self, airP, latent, cp=1012.0):
		"""
		Psychrometric constant :math:`(kPa.K^{-1})`.

		:param airP: Atmospheric pressure (kPa)
		:type airP: numpy.ndarray, float
		:param latent: Latent heat for water vapour exchange :math:`(J.g^{-1})`
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float

		:return: Psychrometric constant :math:`(kPa.K^{-1})`
		:rtype: numpy.ndarray, float
		"""

		try:
			gamma = cp * airP / (latent * 0.622) * 0.001
		except ArithmeticError:
			raise ArithmeticError("Psychrometric constant has not been "
			                      "calculated.")

		return gamma

	def delta(self, ts, ta):
		"""
		Slope of water vapour pressure gradient to temperature gradient - delta
		function :math:`(kPa.K^{-1})`. Calculation according to Jackson et
		al. (1998).

		:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
		:type ts: numpy.ndarray
		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float

		:return: Slope of water vapour pressure gradient to temperature\
		gradient - delta function :math:`(kPa.K^{-1})`
		:rtype: numpy.ndarray
		"""

		try:
			T = (ts + ta) / 2
			delta = (
						        45.03 + 3.014 * T + 0.05345 * T ** 2 + 0.00224 * T ** 3) * 0.001
			del T
		except ArithmeticError:
			raise ArithmeticError("Delta has not been calculated.")

		return delta

	def emissivity(self, red, ndvi):
		"""
		Surface emissivity calculated according to Sobrino et al. (2004) NDVI
		Treshold Method.

		:param red: Spectral reflectance in RED region (rel.)
		:type red: numpy.ndarray
		:param ndvi: Spectral vegetation index NDVI (unitless)
		:type ndvi: numpy.ndarray

		:return: Surface emissivity (rel.)
		:rtype: numpy.ndarray
		"""

		try:
			Fc = self.fractVegCover(ndvi)
			emiss = 0.004 * Fc + 0.986
			emiss = np.where(ndvi < 0.2, 1 - red,
			                 emiss)  # replacement of emissivity values for NDVI
			# < 0.2 by values from red band

			emiss = np.where(ndvi > 0.5, 0.99, emiss)  # replacement of
			# emissivity values for NDVI > 0.5 by 0.99

			emiss[emiss > 1] = 0.99  # replacement of values > 1 by 0.99
			emiss[emiss < 0.8] = 0.8  # replacement values < 0.8 by 0.8

		except ArithmeticError:
			raise ArithmeticError("Surface emissivity has not been calculated.")

		return emiss


# noinspection PyMethodMayBeStatic,PyUnusedLocal,PyShadowingNames
class HeatFluxes(MeteoFeatures):
	"""
	Calculation of heat fluxes and heat balance features from spectral
	and thermal spatial data and meteorological measurements. The class
	contains a set of methods for calculation heat balance, e.g. ground
	heat flux, sensible heat flux, latent heat flux, evaporative fraction,
	omega factor (decoupling coefficient), surface resistance for water
	vapour transfer etc. Calculation for both aerodynamic and gradient
	method is included. 
	"""

	def fluxLE_p(self, Rn, G, delta, VPD, ra, gamma, rho, cp=1012.0):
		"""
		Latent heat flux for potential evapotranspiration according to
		Penman (1948) :math:`(W.m^{-2})`.

		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray
		:param delta: Delta function :math:`(kPa.K^{-1})`
		:type delta: numpy.ndarray
		:param VPD: Water vapour pressure in air (kPa)
		:type VPD: numpy.ndarray
		:param ra: Aerodynamic resistance :math:`(s.m^{-1})`
		:type ra: numpy.ndarray
		:param gamma: Psychrometric constant :math:`(kPa.K^{-1})`.
		:type gamma: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float

		:return: Potential evaporation latent heat flux :math:`(W.m^{-2})`.
		:rtype: numpy.ndarray
		"""

		try:
			LE_p = (delta * (Rn - G) + rho * cp * VPD / ra) / (delta + gamma)
		except ArithmeticError:
			raise ArithmeticError("Latent heat flux for potential evaporation "
			                      "according to Penman (1948) has not been "
			                      "calculated")

		return LE_p

	def fluxLE_EQ(self, Rn, G, delta, gamma):
		"""
		Equilibrium evaporation rate :math:`(W.m^{-2})`.

		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray
		:param delta: Delta function :math:`(kPa.K^{-1})`
		:type delta: numpy.ndarray
		:param gamma: Psychrometric constant :math:`(kPa.K^{-1})`.
		:type gamma: numpy.ndarray

		:return: Equilibrium evaporation latent heat flux :math:`(W.m^{-2})`.
		:rtype: numpy.ndarray
		"""

		try:
			LE_eq = delta / (delta + gamma) * (Rn - G)
		except ArithmeticError:
			raise ArithmeticError("Latent heat flux for equilibrium "
			                      "evaporation has not been calculated")

		return LE_eq

	def fluxLE_PT(self, LE_eq, alpha=1.26):
		"""
		Evaporation from wet surface according to Priestley-Taylor (1972).

		:param LE_eq: Equilibrium evaporation latent heat flux\
		:math:`(W.m^{-2})`
		:type LE_eq: numpy.ndarray
		:param alpha: Priestley-Taylor alpha.
		:type alpha: float

		:return: Latent heat flux for wet surface calculated by\
		Priestley-Taylor :math:`(W.m^{-2})`
		:rtype: numpy.ndarray
		"""

		try:
			LE_PT = LE_eq * alpha
		except ArithmeticError:
			raise ArithmeticError("Latent heat flux for Priestley-Taylor "
			                      "evaporation has not been calculated")

		return LE_PT

	def fluxLE_ref(self, Rn, G, ta, ts, U, Rh, DMT, veg_type="short",
	               cp = 1012.0):
		"""
		Latent heat flux for reference evapotranspiration according to FAO56
		method (Allen et al. 1998, ASCE-ET 2000)

		:param Rn: Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray, float
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray, float
		:param ta: Air temperature measured at meteostation at approx. 2 m.
		:math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float
		:param ts: Surface temperature measured at meteostation at approx. 2 m.
		:math:`(\SI{}\degreeCelsius)`
		:param U: Wind speed measured at meteostation :math:`(m.s^{-1})`
		:type U: numpy.ndarray, float
		:param Rh: Relative humidity (%)
		:type Rh: numpy.ndarray, float
		:param DMT: Digital model of terrain or altitude (m a.s.l.)
		:type DMT: numpy.ndarray, float
		:param veg_type: Type of vegetation - "short" or "tall"
		:type veg_type: str
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float

		:return: Latent heat flux for reference evapotranspiration.
		:rtype: numpy.ndarray, float

		\n
		**References**\n
		*Allen, R. et al., 1998. Crop Evapotranspiration: Guidelines for
		Computing Crop Water Requirements. United Nations Food and
		Agriculture Organization, Irrigation and Drainage Paper 56, Rome,
		Italy. 300 pp.*
		*ASCE-ET 2000: ASCE’s Standardized Reference Evapotranspiration
		Equation. National Irrigation Symposium in Phoenix, Arizona, US.
		"""

		if veg_type is "short":
			if Rn > 0.0:
				cn = 37.0
				cd = 0.24
			else:
				cn = 37.0
				cd = 0.96
		else:
			if Rn > 0:
				cn = 66.0
				cd = 0.25
			else:
				cn = 66.0
				cd = 1.7

		delta = self.delta(ts, ta)
		air_pressure = self.airPress(ta, DMT, 2.0)
		latent_heat = self.latent(ta)
		gamma = self.gamma(air_pressure, latent_heat, cp)
		E_sat = self.satVapourPress(ta)
		e_abs = self.vapourPress(E_sat, Rh)
		VPD = self.vpd(E_sat, e_abs)

		try:
			LE_ref = (408.0 * delta * (Rn - G) * 0.0036) + (gamma * cn/(ta +
			          273.15) * U * VPD)/(delta + gamma * (1 + cd * U)) * \
			         latent_heat/3600.0
		except ArithmeticError:
			raise ArithmeticError("Reference evapotranspiration has not been "
			                      "calculated.")

		return LE_ref

	def fluxHAer(self, ra, rho, dT, cp=1012.0):
		"""
		Sensible heat flux :math:`(W.m^{-2})` calculated using aerodynamic 
		method.
		
		:param ra: Aerodynamic resistance for heat and momentum transfer\
		:math:`(s.m^{-1})` calculated according to Thom (1975)
		:type ra: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param dT: Temperature gradient calculated according to SEBAL\
		(Bastiaanssen et al. 1998) :math:`(\SI{}\degreeCelsius)`
		:type dT: numpy.ndarray
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float

		:return: Sensible heat flux :math:`(W.m^{-2})`
		:rtype: numpy.ndarray

		"""

		try:
			flux_H = rho * cp * dT / ra
		except ArithmeticError:
			raise ArithmeticError("Sensible heat flux has not been calculated.")

		return flux_H

	def gradH(self, LE, Rn, G):
		"""
		Sensible heat flux :math:`(W.m^{-2})` calculated using gradient method.

		:param LE: Latent heat flux :math:`(W.m^{-2})`
		:type LE: numpy.ndarray
		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray

		:return: Sensible heat flux :math:`(W.m^{-2})`
		:rtype: numpy.ndarray
		"""

		try:
			flux_H = Rn - G - LE
		except ArithmeticError:
			raise ArithmeticError("Sensible heat flux has not been calculated.")

		return flux_H

	def fluxLE(self, Rn, G, flux_H):
		"""
		Latent heat flux :math:`(W.m^{-2})`

		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray
		:param flux_H: Sensible heat flux :math:`(W.m^{-2})`
		:type flux_H: numpy.ndarray

		:return: Latent heat flux :math:`(W.m^{-2})`
		:rtype: numpy.ndarray
		"""

		try:
			LE = Rn - G - flux_H
		except ArithmeticError:
			raise ArithmeticError("Latent heat flux has not been calculated.")

		return LE

	def gradLE(self, EF, Rn, G):
		"""
		Latent heat flux calculated using gradient method :math:`(W.m^{-2})`.

		:param EF: Evaporative fraction (rel.)
		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray

		:return: Latent heat flux :math:`(W.m^{-2})`
		:rtype: numpy.ndarray
		"""

		try:
			LE = EF * (Rn - G)
		except ArithmeticError:
			raise ArithmeticError("Latent heat flux has not been calculated.")

		return LE

	def aeroEF(self, LE, Rn, G):
		"""
		Evaporative fraction calculated using aerodynamic method (rel.).

		:param LE: Latent heat flux :math:`(W.m^{-2})`
		:type LE: numpy.ndarray
		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray

		:return: Evaporative fraction (rel.)
		:rtype: numpy.ndarray
		"""

		try:
			EF = LE / (Rn - G)
		except ArithmeticError:
			raise ArithmeticError("Evaporative fraction has not been "
			                      "calculated.")

		return EF

	def gradEF(self, ts, ta, mask=None):
		"""
		Evaporative fraction calculated from gradient method according
		to Suleiman and Crago (2004).

		:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
		:type ts: numpy.ndarray
		:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray, float
		:param mask: Mask of the area of interest. Number of rows and columns
		should be the same.  Format (1, 0) or (1, nan).
		:type mask: numpy.ndarray

		:return: Evaporative fraction (rel.)
		:rtype: numpy.ndarray
		"""

		try:
			t_max = WindStability().dryT(ts, mask)
			EF = (t_max - ts) / (t_max - ta)
		except ArithmeticError:
			raise ArithmeticError("Evaporative fraction has not been "
			                      "calculated.")

		return EF

	def heatFluxes(self, Rn, G, ts, ta, method="aero", Uz=None, h_eff=None,
	               LAI=None, z0m=None, z0h=None, rho=None, disp=None, mask=None,
	               air_pressure=101.3, Z=200.0, cp=1012, L=-10000.0, n_iter=10,
	               a=1.0, b=0.667, c=5.0, d=0.35, kappa=0.41, gravit=9.81):
		# noinspection SpellCheckingInspection
		"""
		Function provides a calculation for heat fluxes, aerodynimc
		resistance of surface and friction velocity for three different
		methods: \
			1. "aero" - Aerodynamic method based on Monin-Obukhov theory \
			2. "sebal" - Method based on SEBAL approach provided by
			Bastiaanssen et al.

				:param disp:
				:param Rn: Total net radiation :math:`(W.m^{-2})`
				:type Rn: numpy.ndarray
				:param G: Ground heat flux :math:`(W.m^{-2})`
				:type G: numpy.ndarray
				:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
				:type ts: numpy.ndarray
				:param ta: Air temperature :math:`(\SI{}\degreeCelsius)`
				:type ta: numpy.ndarray, float
				:param method: Method of heat fluxes calculation: \n\n
						- aero - Aerodynamic method based on  calculation
						of ra using approach proposed by Thom (1975) \n
						- SEBAL - SEBAL method proposed by Bastiaanssen et al. (1998) \n
						- grad - Gradient method proposed by Suleiman and Crago (2004)
				:type method: str
				:param Uz: Wind speed measured on meteostation at level Z_st :math:`(m.s^{-1})`.
				:type Uz: float, numpy.ndarray
				:param h_eff: Effective height of vegetation cover
				:type h_eff: numpy.ndarray, float
				:param LAI: Leaf area index :math:`(m^{2}.m^{-2})`
				:param z0m: Aerodynamic roughness of the surface for momentum
				transfer (m)
				:type z0m: numpy.ndarray
				:param z0h: Aerodynamic roughness of the surface for heat
				transfer (m)
				:type z0h: numpy.ndarray
				:param rho: Volumetric dry air density :math:`(kg.m^{-3})`
				:type rho: numpy.ndarray, float
				:param disp: Zero plane displacement (m)
				:type disp: numpy.ndarray, float
				:param mask: Mask of the area of interest. Number of rows
				and columns should be the same.  Format (1, 0) or (1, nan).
				:type mask: numpy.ndarray
				:param air_pressure: Air pressure  at level Z (kPa)
				:type air_pressure: float, numpy.ndarray
				:param Z: Blending height (m)
				:type Z: float
				:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
				:type cp: float
				:param L: Initial Monin-Obukhov lenght (m).
				:type L: float
				:param n_iter: Number of iteration in stability coefficient calculation.
				:type n_iter: int
				:param a: Constant for stability parameters calculation (Beljaars and
				Holtslag, 1991)
				:type a: float
				:param b: Constant for stability parameters calculation (Beljaars and
				Holtslag, 1991)
				:type b: float
				:param c: Constant for stability parameters calculation (Beljaars and
				Holtslag, 1991)
				:type c: float
				:param d: Constant for stability parameters calculation (Beljaars and
				Holtslag, 1991)
				:type d: float
				:param kappa: von Karman constant. Default 0.41
				:type kappa: float
				:param gravit: Gravitation forcing (m/s2). Default 9.81
				:type gravit: float

				:returns: Sensible heat flux :math:`(W.m^{-2})`
				:rtype: numpy.ndarray
				:returns: Latent heat flux :math:`(W.m^{-2})`
				:rtype: numpy.ndarray
				:returns: Evaporative fraction (rel.)
				:rtype: numpy.ndarray
				:returns: Latent heat flux for equilibrium evaporation :math:`(W.m^{
				-2})`
				:rtype: numpy.ndarray
				:returns: Latent heat flux for Priestley-Taylor evaporation :math:`(
				W.m^{-2})`
				:rtype: numpy.ndarray
				:returns: Aerodynamic resistance for heat and momentum transfer :math:`(
				s.m^{-1})`
				:rtype: numpy.ndarray
				:returns: Friction velocity :math:`(m.s^{-1})`

				\n
				**References**\n
				*Bastiaanssen, W.G.M., Menenti, M., Feddes, R.A., Holtslag, A.A.M.,
				1998. A remote sensing surface energy balance algorithm for land (
				SEBAL). 1. Formulation. Journal of Hydrology 212–213, 198–212.
				https://doi.org/10.1016/S0022-1694(98)00253-4*\n
				*Beljaars, A.C.M., Holtslag, A.A.M., 1991. Flux Parametrization over
				Land Surfaces for Atmospheric Models. Journal of Applied Meteorology
				30, 327–341. https://doi.org/10.1175/1520-0450(
				1991)030<0327:FPOLSF>2.0.CO;2*\n
				*Suleiman, A., Crago, R., 2004. Hourly and Daytime Evapotranspiration
				from Grassland Using Radiometric Surface Temperatures. Agronomy
				Journal 96, 384–390. https://doi.org/10.2134/agronj2004.3840*\n
				*Thom, A.S., 1975. Momentum, mass and heat exchange of plant
				communities, in: Monteith, J.L. (Ed.), Vegetation and the Atmosphere,
				Vol. 1 Principles. Academic Press, London, pp. 57–110.*
				"""

		ignore_zero = np.seterr(all="ignore")

		# Handling inputs
		if method is "aero" or method is "SEBAL":
			if Uz is None:
				raise IOError("Heat fluxes have not been calculated - wind "
				              "speed has not been set up.")
			else:
				if z0m is None or z0h is None:
					if h_eff is not None and LAI is not None:
						try:
							z0m = WindStability().z0m(h_eff, LAI)
							z0h = WindStability().z0h(z0m)
						except ArithmeticError:
							raise ArithmeticError("Parameters of surface "
							                      "aerodynamic roughness"
							                      "have not been calculated.")
					else:
						raise IOError("Parameters of effective height of "
						              "vegetation or leaf area index has not "
						              "been set up correctly")

				if disp is None:
					if h_eff is not None:
						disp = WindStability().zeroPlaneDis(h_eff)
					else:
						disp = 0.0
		try:
			if rho is None or rho is "":
				rho = self.airDensity(ta)
		except ArithmeticError:
			raise ArithmeticError("Air density has not been "
			                      "calculated.")

		# Heat fluxes calculation
		if method is "aero":
			try:
				psi_m, psi_h, frict, L = WindStability().stabCoef(Uz, ta, ts,
				                                                  z0m, z0h,
				                                                  disp, Z, L,
				                                                  n_iter, a,
				                                                  b, c, d,
				                                                  kappa, gravit)
				ra = WindStability().raThom(Uz, z0m, z0h, disp=disp,
				                            psi_m=psi_m, psi_h=psi_h, Z=Z,
				                            kappa=kappa)
				dT = WindStability().dT(ts, ta, Rn, G, ra, rho, cp, mask)
				flux_H = self.fluxHAer(ra, rho, dT, cp)
				LE = self.fluxLE(Rn, G, flux_H)
				EF = self.aeroEF(LE, Rn, G)

			except ArithmeticError:
				raise ArithmeticError("Heat fluxes for aerodynamic method "
				                      "has not been calculated")

		elif method is "SEBAL":
			try:
				psi_m, psi_h, frict, L = WindStability().stabCoef(Uz, ta, ts,
				                                                  z0m,
				                                                  z0h, disp, Z,
				                                                  L, n_iter,
				                                                  a, b, c, d,
				                                                  kappa, gravit)
				ra = WindStability().raThom(Uz, z0m, z0h, disp=disp,
				                            psi_m=psi_m, psi_h=psi_h, Z=Z,
				                            kappa=kappa)
				ra_h, flux_H = WindStability().aeroSEBAL(Uz, ta, ts, z0m, Rn, G,
				                                         rho, n_iter, Z, cp=cp,
				                                         kappa=kappa, mask=mask)
				LE = self.fluxLE(Rn, G, flux_H)
				EF = self.aeroEF(LE, Rn, G)

			except ArithmeticError:
				raise ArithmeticError("Heat fluxes for SEBAL method "
				                      "has not been calculated")

		elif method is "grad":
			EF = self.gradEF(ts, ta, mask)
			LE = self.gradLE(EF, Rn, G)
			flux_H = self.gradH(LE, Rn, G)
			ra = WindStability().raGrad(flux_H, rho, ts, ta, cp)
			frict = np.zeros_like(ra)

		else:
			raise ArithmeticError("Heat fluxes have not been calculated")

		# Correction of heat fluxes for case the ts-ta < 0
		# <-- changing stability (temperature gradient inversion, stable
		# atmosphere)
		# <-- usually water bodies <-- values of fluxes can be recalculated
		# according to Priestley-Taylor formula

		try:
			delta = self.delta(ts, ta)
			latent_heat = self.latent(ta)
			gamma = self.gamma(air_pressure, latent_heat, cp)
			LE_eq = self.fluxLE_EQ(Rn, G, delta, gamma)
			LE_PT = self.fluxLE_PT(LE_eq, 1.26)

			LE = np.where(ts-ta < 0, LE_PT, LE)
			flux_H = Rn - G - LE
		except ArithmeticError:
			warnings.warn("Correction of heat fluxes has not been done",
			              stacklevel = 3)

		return flux_H, LE, EF, LE_eq, LE_PT, ra, frict

	def intensityE(self, LE, latent):
		"""Evaporation intensity in :math:`mmol.m^{-2}.s^{-1}`.
		
		:param LE: Latent heat flux :math:`(W.m^{-2})`
		:type LE: numpy.ndarray
		:param latent: Latent heat of water evaporation :math:`(J.g^{-1})`.
		:type latent: numpy.ndarray
		
		:returns: Intensity of water evaporation :math:`(mmol.m^{-2}.s^{-1})`.
		:rtype: numpy.ndarray
		"""

		E_int = LE / latent / 18 * 1000

		return E_int

	def omega(self, LE, LE_p):
		"""
		Omega factor (Decoupling coefficient) according to Jarvis and
		McNaughton (1985)

		:param LE: Latent heat flux :math:`(W.m^{-2})`
		:type LE: numpy.ndarray
		:param LE_p: Potential evaporation latent heat flux :math:`(W.m^{-2})`
		:type LE_p: numpy.ndarray

		:return: Omega factor (unitless)
		:rtype: numpy.ndarray

		\n
		**References:**\n
		*Jarvis, P., McNaughton, K., 1986. Stomatal Control of Transpiration:
		Scaling Up from Leaf to Region, in: Advances in Ecological Research.
		Elsevier, pp. 1–49.*
		"""

		try:
			omega = LE / LE_p
		except ArithmeticError:
			raise ArithmeticError("Omega factor has not been calculated")

		return omega

	def rs(self, delta, gamma, omega, ra):
		"""
		Surface resistance for water vapour transfer :math:`(s.m^{-1})`

		:param delta: Delta function :math:`(kPa.K^{-1})`
		:type delta: numpy.ndarray
		:param gamma: Psychrometric constant :math:`(kPa.K^{-1})`.
		:type gamma: numpy.ndarray
		:param omega: Decoupling coefficient (rel.)
		:type omega: numpy.ndarray
		:param ra: Aerodynamic resistance :math:`(s.m^{-1})`
		:type ra: numpy.ndarray

		:return: Surface resistance for water vapour transfer :math:`(s.m^{-1})`
		:rtype: numpy.ndarray
		"""

		try:
			rs = (((delta + gamma) / omega - delta) * 1.0 / gamma - 1.0) * ra
		except ArithmeticError:
			raise ArithmeticError("Surface resistance for water vapour "
			                      "transfer has not been calculated")

		return rs

	def bowen(self, flux_H, LE):
		"""
		Bowen ratio according to Bowen (1926)

		:param flux_H: Sensible heat flux :math:`(W.m^{-2})`
		:type flux_H: numpy.ndarray
		:param LE: Latent heat flux :math:`(W.m^{-2})`
		:type LE: numpy.ndarray

		:return: Bowen ratio (unitless)
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			bowen = flux_H / LE
			bowen = np.nan_to_num(bowen)
		except ArithmeticError:
			raise ArithmeticError("Bowen ratio has not been calculated")

		return bowen

	def rcp(self, E_Z_sat, e_Z, rho, LE_p, ra, gamma, cp=1012.0):
		"""
		Surface resistance for water vapour transfer for potential
		evapotranspiration :math:`(s.m^{-1})`.

		:param E_Z_sat: Saturated water vapour pressure (kPa).
		:type E_Z_sat: numpy.ndarray
		:param e_Z: Water vapour pressure (kPa).
		:type e_Z: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param LE_p: Potential evaporation latent heat flux :math:`(W.m^{-2})`
		:type LE_p: numpy.ndarray
		:param ra: Aerodynamic resistance :math:`(s.m^{-1})`
		:type ra: numpy.ndarray
		:param gamma: Psychrometric constant :math:`(kPa.K^{-1})`.
		:type gamma: numpy.ndarray
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float

		:return: Surface resistance for water vapour transfer for potential\
		evapotranspiration :math:`(s.m^{-1})`.
		:rtype: numpy.ndarray
		"""

		try:
			rcp = (E_Z_sat - e_Z) * rho * cp / (gamma * LE_p) - ra
		except ArithmeticError:
			raise ArithmeticError("Surface resistance for water vapour "
			                      "transfer for potential evaporation has not "
			                      "been calculated")

		return rcp

	def gamma_x(self, rcp, ra, gamma):
		"""
		Psychrometric constant corrected on the rcp and ra according to
		Jackson et al. (1981, 1988).

		:param rcp: Surface resistance for water vapour transfer for potential\
		evapotranspiration :math:`(s.m^{-1})`.
		:type rcp: numpy.ndarray
		:param ra: Aerodynamic resistance :math:`(s.m^{-1})`
		:type ra: numpy.ndarray
		:param gamma: Psychrometric constant :math:`(kPa.K^{-1})`.
		:type gamma: numpy.ndarray

		:return: Psychrometric constant corrected on the rcp and ra\
		according to Jackson et al. (1981, 1988) :math:`(kPa.K^{-1})`..
		:rtype: numpy.ndarray
		"""

		try:
			gamma_x = gamma * (1.0 + rcp / ra)
		except ArithmeticError:
			raise ArithmeticError("Corrected psychrometric constant has not "
			                      "been calculated")

		return gamma_x

	def cwsi(self, LEp, ra, rc, E_Z_sat, e_Z, rho, delta, gamma, cp=1012):
		"""
		Crop Water Stress Index calculated according to Jackson et al. (
		1981, 1988).

		:param LEp: Potential evaporation latent heat flux :math:`(W.m^{-2})`
		:type LEp: numpy.ndarray
		:param ra: Aerodynamic resistance :math:`(s.m^{-1})`
		:type ra: numpy.ndarray
		:param rc: Surface resistance for water vapour transfer\
		 :math:`(s.m^{-1})`
		:param E_Z_sat: Saturated water vapour pressure (kPa).
		:type E_Z_sat: numpy.ndarray
		:param e_Z: Water vapour pressure (kPa).
		:type e_Z: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param delta: Delta function :math:`(kPa.K^{-1})`
		:type delta: numpy.ndarray
		:param gamma: Psychrometric constant :math:`(kPa.K^{-1})`.
		:type gamma: numpy.ndarray
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float

		:return: Crop Water Stress Index (CWSI)
		:rtype: numpy.ndarray
		"""

		try:
			rcp = self.rcp(E_Z_sat, e_Z, rho, LEp, ra, gamma, cp)
			g_x = self.gamma_x(rcp, ra, gamma)
			cwsi = 1.0 - (delta + g_x) / (delta + gamma * (1.0 + rc / ra))

		except ArithmeticError:
			warnings.warn("CWSI has not been calculated", stacklevel=3)
			if LEp is not None:
				cwsi = np.zeros(LEp.shape)
			else:
				cwsi = None

		return cwsi

	def groundFlux(self, ndvi, Rn, ts, albedo):
		"""
		Ground heat flux according to Bastiaanssen et al. (1998)

		:param ndvi: NDVI spectral vegetation index (unitless)
		:type ndvi: numpy.ndarray
		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
		:type ts: numpy.ndarray
		:param albedo: Albedo (rel.)
		:type albedo: numpy.ndarray

		:return: Ground heat flux :math:`(W.m^{-2})`
		:rtype: numpy.ndarray
		"""

		G = ts / albedo * (0.0038 * albedo + 0.0074 * albedo ** 2) * (1 - 0.98 * ndvi ** 4) * Rn

		return G


# noinspection PyUnusedLocal,PyShadowingNames
class WindStability(HeatFluxes, VegIndices):
	"""
	Atmospheric stability calculation. Class includes methods for calculation
	of boundary layer stability, friction velocity and aerodynamic
	resistance of the surface. Methods for calculation of stability
	parameters for both aerodynamic and gradient methods are included.
	Methods for SEBAL procedure are also used. Some another features
	are included.
	"""

	@staticmethod
	def zeroPlaneDis(h_eff):
		"""
		Zero plane displacement (m) calculated according to Thom (1975).

		:param h_eff: Effective vegetation cover height (m)
		:type h_eff: numpy.ndarray

		:return: Zero plane displacement (m)
		:rtype: numpy.ndarray

		\n
		**References**\n
		*Thom, A.S., 1975. Momentum, mass and heat exchange of plant
		communities, in: Monteith, J.L. (Ed.), Vegetation and the Atmosphere,
		Vol. 1 Principles. Academic Press, London, pp. 57–110.*
		"""

		disp = 2.0 / 3.0 * h_eff
		return disp

	@staticmethod
	def z0m(h_eff, LAI):
		"""
		Aerodynamic roughness of the surface for momentum transfer (m).
		z0m is calculated according to Tasumi (2003) for vegetation cover
		lower tha 1 m and according to Thom (1975) for vegetation cover
		higher than 1 m. The lowest value of z0m is 0.005, which is typical
		value for agricultural bare soils (Allen et al. 2007)

		:param h_eff: Effective vegetation cover height (m)
		:type h_eff: numpy.ndarray
		:param LAI: Leaf Area Index
		:type LAI: numpy.ndarray

		:return: Aerodynamic roughness of the surface for momentum transfer (m)
		:rtype: numpy.ndarray

		\n
		**References**\n
		*Allen, R.G., Tasumi, M., Trezza, R., 2007. Satellite-Based Energy
		Balance for Mapping Evapotranspiration with Internalized Calibration
		(METRIC)—Model. J. Irrig. Drain Eng. 133, 380–394.
		https://doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380)* \n
		*Tasumi, M., 2003. Progress in Operational Estimation of Regional
		Evapotranspiration Using Satellite Imagery (Ph.D. Thesis).
		University of Idaho, Moscow, Idaho.*\n
		*Thom, A.S., 1975. Momentum, mass and heat exchange of plant
		communities, in: Monteith, J.L. (Ed.), Vegetation and the Atmosphere,
		Vol. 1 Principles. Academic Press, London, pp. 57–110.*
		"""

		try:
			z0m = np.where(h_eff < 1.0, 0.018 * LAI, 0.123 * h_eff)
			z0m = np.where(z0m < 0.005, 0.005, z0m)

		except ArithmeticError:
			raise ArithmeticError("Aerodynamic roughness of the surface for "
			                      "momentum transfer has not been calculated")
		return z0m

	@staticmethod
	def z0h(z0m):
		"""
		Aerodynamic roughness of the surface for heat transfer (m) calculated
		according to Thom (1975).

		:param z0m: Aerodynamic roughness of the surface for momentum\
		transfer (m)
		:type z0m: numpy.ndarray, float

		:return: Aerodynamic roughness of the surface for heat transfer (m)
		:rtype: numpy.ndarray, float

		\n
		**References**\n
		*Thom, A.S., 1975. Momentum, mass and heat exchange of plant
		communities, in: Monteith, J.L. (Ed.), Vegetation and the Atmosphere,
		Vol. 1 Principles. Academic Press, London, pp. 57–110.*
		"""
		z0h = 0.1 * z0m
		return z0h

	@staticmethod
	def windSpeedZ(U, Z=200.0, Z_st=2.0, h_st=0.12, ws_homog=1):
		"""
		Wind speed recalculated to height Z according to logarithmic law
		(Gao et al. 2011). The results can contain simple number, homogenous
		or heterogenous matrix of data (Numpy array).

		:param U: Wind speed measured on meteostation at level Z_st :math:`(m.s^{-1})`.
		:type U: float, numpy.ndarray
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param Z_st: Height of wind speed measurement (m). Default 2 m.
		:type Z_st: float
		:param h_st: Height of vegetation cover under meteostation (m).
					 Default value is 0.12 m which corresponds with
					 reference cover used for meteostations.
		:type h_st: float
		:param ws_homog: Indicates if the result matrix contains uniform
		result (mean value of wind speed for the area) or if the results
		includes heterogeneity of input data (each pixel is calculated
		separately). 0 - homogenous, 1 - heterogenous
		:type ws_homog: int
		:return: Wind speed at mixing layer :math:`(m.s^{-1})`
		:rtype: float, numpy.ndarray

		\n
		**References**\n
		*Gao, Z.Q., Liu, C.S., Gao, W., Chang, N.B., 2011. A coupled remote
		sensing and the Surface Energy Balance with Topography Algorithm
		(SEBTA) to estimate actual evapotranspiration over heterogeneous
		terrain. Hydrol. Earth Syst. Sci. 15, 119–139.
		https://doi.org/10.5194/hess-15-119-2011*
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			if ws_homog is 0:
				U_float = np.mean(U)
				Uz = U_float * np.log(Z / (0.123 * h_st)) / np.log(Z_st / (
						0.123 * h_st))
				if type(U) is not float:
					Uz = np.full_like(U, Uz, float)
			else:
				Uz = U * np.log(Z / (0.123 * h_st)) / np.log(Z_st / (
						0.123 * h_st))
		except ArithmeticError:
			raise ArithmeticError("Friction velocity has not been calculated")

		return Uz

	@staticmethod
	def frictVelo(Uz, z0m, disp=0.0, Z=200.0, psi_m=0, kappa=0.41):
		"""
		Friction velocity of wind speed :math:`(m.s^{-1})` corrected on atmospheric
		stability.

		:param disp:
		:param Uz: Wind speed at Z level :math:`(m.s^{-1})`
		:type Uz: numpy.ndarray
		:param z0m: Surface roughness for momentum transfer (m)
		:type z0m: numpy.ndarray, float
		:param disp: Zero plane displacement (m)
		:type disp: numpy.ndarray, float
		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param psi_m: Stability parameter for momentum transfer.
					  Defaul 0.
		:type psi_m: numpy.ndarray
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float

		:return: Friction velocity :math:`(m.s^{-1})`
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			frict = kappa * Uz / (np.log((Z - disp)/ z0m) - psi_m)
		except ArithmeticError:
			raise ArithmeticError("Friction velocity has not been calculated")

		return frict

	@staticmethod
	def virtTemp(ta, ts, z0h, disp=0.0, Z=200, psi_h=0, kappa=0.41):
		"""
		Virtual temperature (K) corrected on atmospheric stability.

		:param disp:
		:param ta: Air temperature at Z level (K, C degrees)
		:type ta: numpy.ndarray
		:param ts: Surface temperature (K, C degrees)
		:type ts: numpy.ndarray
		:param z0h: Surface roughness for heat transfer (m)
		:type z0h: numpy.ndarray
		:param disp: Zero plane displacement (m)
		:type disp: numpy.ndarray, float
		:param Z: Blending height (mixing layer height) (m). Default 200 m.
		:type Z: float
		:param psi_h: Stability parameter for heat transfer, Default 0.
		:type psi_h: numpy.ndarray
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float

		:returns: Virtual temperature (K)
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			t_virt = kappa * (ta - ts) / (np.log((Z - disp) / z0h) - psi_h)
		except ArithmeticError:
			raise ArithmeticError("Virtual temperature has not been calculated")

		return t_virt

	@staticmethod
	def psiM(L, X, Z=200, a=1.0, b=0.667, c=5.0, d=0.35):
		"""
		Calculation of stability parameter for momentum transfer (-)
		according to Beljaars et Holstag (1991) for stable conditions
		and Liu et al. (2007) for unstable and neutral conditions.

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

		:return: Stability parameter for momentum transfer (-)
		:rtype: numpy.ndarray
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

	@staticmethod
	def psiH(L, X, Z=200, a=1.0, b=0.667, c=5.0, d=0.35):
		"""
		Calculation of stability parameter for heat transfer (-)
		according to Beljaars et Holstag (1991) for stable conditions
		and Liu et al. (2007) for unstable and neutral conditions.

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

		:return: Stability parameter for momentum transfer (-)
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		dzeta = Z / L

		try:
			psi_h_stab = -((1 + 2 * a * dzeta / 3) ** 1.5 + b
			               * (dzeta - c / d) * np.exp((-d) * dzeta)
			               + (b * c/d - 1))
		except ArithmeticError:
			raise ArithmeticError("Stability coefficient for heat transfer "
			                      "for stable conditions has not been "
			                      "calculated")

		try:
			psi_h_unstab = 2.0 * np.log((1.0 + X ** 2.0) / 2.0)
		except ArithmeticError:
			raise ArithmeticError(
				"Stability coefficient for heat transfer for unstable conditions has not been calculated")

		psi_h = np.where(dzeta < 0.0, psi_h_unstab, psi_h_stab)

		return psi_h

	@staticmethod
	def lengthMO(frict, ts, flux_H=None, rho=None, t_virt=None, cp=1012, kappa=0.41, gravit=9.81):
		"""
		Monin-Obukhov length (m)

		:param frict: Friction velocity :math:`(m.s^{-1})`
		:type frict: numpy.ndarray
		:param ts: Surface temperature (C degree)
		:type ts: numpy.ndarray
		:param flux_H: Sensible heat flux (W.m2)
		:type flux_H: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param t_virt: Virtual temperature
		:type t_virt: numpy.ndarray
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float
		:param gravit: Gravitation forcing (m/s2). Default 9.81
		:type gravit: float

		:return: Monin-Obukhov length (m)
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			if flux_H is None:
				L = frict ** 2.0 * (ts + 273.16) / (kappa * gravit * t_virt)
			else:
				L = -(frict ** 3 * (ts + 273.15) * rho * cp / (kappa * gravit * flux_H))
		except ArithmeticError:
			raise ArithmeticError("Monin-Obukhov lenght has not been calculated")

		return L

	@staticmethod
	def coefX(Z, L):
		"""
		X coefficient for atmospheric stability calculation.

		:param Z: Blending height (mixing layer height) (m).
				  Default 200 m.
		:type Z: float
		:param L: Monin-Obukhov length (m)
		:type L: numpy.ndarray

		:return: X coefficient for atmospheric stability calculation
		:rtype: numpy.ndarray
		"""

		ignore_zero = np.seterr(all="ignore")

		try:
			X = np.where((1.0 - 16.0 * Z / L) ** 0.25 < 0.0, 0.0, (1.0 - 16.0 * Z / L) ** 0.25)
			X = np.nan_to_num(X)
		except ArithmeticError:
			raise ArithmeticError("Coefficient X has not been calculated")

		return X

	def stabCoef(self, Uz, ta, ts, z0m, z0h, disp=0, Z=200, L=-10000.0,
	             n_iter=10, a=1.0, b=0.667, c=5.0, d=0.35, kappa=0.41,
	             gravit=9.81):
		"""
		Stability parameters calculation using iterative procedure
		described by Itier (1980).

		:param disp:
		:param Uz: Wind speed at level Z :math:`(m.s^{-1})`
		:type Uz: numpy.ndarray
		:param ta: Air temperature at Z level (K, C degrees)
		:type ta: numpy.ndarray
		:param ts: Surface temperature (K, C degrees)
		:type ts: numpy.ndarray
		:param z0m: Surface roughness for momentum transfer (m)
		:type z0m: numpy.ndarray
		:param z0h: Surface roughness for heat transfer (m)
		:type z0h: numpy.ndarray
		:param disp: Zero plane displacement (m)
		:type disp: numpy.ndarray, float
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

		:return: Stability parameter for momentum transfer (-)
		:rtype: numpy.ndarray
		:return: Stability parameter for heat transfer (-)
		:rtype: numpy.ndarray
		:return: Friction velocity :math:`(m.s^{-1})`.
		:rtype: numpy.ndarray
		:return: Monin-Obukhov length (m)
		:rtype: numpy.ndarray
		"""

		# definition of returned variables
		psi_m = None
		psi_h = None
		frict = None

		# Stability parameters calculation
		try:
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
				frict = self.frictVelo(Uz, z0m, disp=disp, Z=Z, psi_m=psi_m,
				                       kappa=kappa)

				# Virtual tempetarure
				t_virt = self.virtTemp(ta, ts, z0h, disp=disp, Z=Z, psi_h=psi_h,
				                       kappa=kappa)

				# Monin-Obukhov length
				L = self.lengthMO(frict, ts, t_virt=t_virt, kappa=kappa, gravit=gravit)

		except ArithmeticError:
			raise ArithmeticError("Stability parameters has not been "
			                      "calculated.")

		return psi_m, psi_h, frict, L

	@staticmethod
	def raThom(Uz, z0m, z0h, disp=0.0, psi_m=0.0, psi_h=0.0, Z=200.0,
	           kappa=0.41):
		"""
		Aerodynamic resistance for heat and momentum transfer (s.m-1)
		calculated according to Thom (1975).

		:param disp:
		:param Uz: Wind speed at level Z :math:`(m.s^{-1})`
		:type Uz: numpy.ndarray
		:param z0m: Surface roughness for momentum transfer (m)
		:type z0m: numpy.ndarray
		:param z0h: Surface roughness for heat transfer (m)
		:type z0h: numpy.ndarray
		:param disp: Zero plane displacement (m)
		:type disp: numpy.ndarray, float
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

		:return ra: Aerodynamic resistance for heat and momentum
					transfer (s.m-1) calculated according to Thom (1975)
		:rtype ra: numpy.ndarray

		\n
		**References**\n
		*Thom, A.S., 1975. Momentum, mass and heat exchange of plant
		communities, in: Monteith, J.L. (Ed.), Vegetation and the Atmosphere,
		Vol. 1 Principles. Academic Press, London, pp. 57–110.*
		"""

		try:
			ra = (np.log((Z - disp)/ z0m) - psi_m) * (np.log((Z - disp)/z0h)
			                                          - psi_h) / (kappa ** 2 * Uz)
		except ArithmeticError:
			raise ArithmeticError("Aerodynamic resistance has not been "
			                      "calculated")

		return ra

	@staticmethod
	def raSEBAL(frict, psiH_z1, psiH_z2, z1=0.1, z2=2.0,
	            kappa=0.41):
		"""
		Calculation of surface aerodynamic resistance for heat transfer based on
		SEBAL approach (Bastiaanssen, 1998).

		:param frict: Friction velocity :math:`(m.s^{-1})`
		:type frict: numpy.ndarray
		:param psiH_z1: Stability parameter for momentum transfer (-) at\
		level z1
		:type psiH_z1: numpy.ndarray
		:param psiH_z2: Stability parameter for momentum transfer (-) at\
		level z2
		:type psiH_z2: numpy.ndarray
		:param z1: First height above zero plane displacement (m)
		:type z1: float
		:param z2: Second height above zero plane displacement (m)
		:type z2: float
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float

		:return: Aerodynamic resistance for heat transfer based on\
		SEBAL approach (Bastiaanssen, 1998).
		:rtype: numpy.ndarray
		"""

		# ra calculation
		try:
			ra = (np.log(z2 / z1) - psiH_z2 + psiH_z1) / (frict * kappa)
		except ArithmeticError:
			raise ArithmeticError("Aerodynamic resistance has not been "
			                      "calculated!")

		return ra

	@staticmethod
	def raGrad(flux_H, rho, ts, ta, cp=1012):
		"""
		Aerodynamic resistance for heat and momentum transfer `(s.m^{-1})`
		calculated from conversion of sensible heat flux equation.

		:param flux_H: Sensible heat flux :math:`(W.m^{-2})`
		:type flux_H: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray
		:param ta: Air temperature at level Z :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float

		:return: Aerodynamic resistance for gradient model :math:`(s.m^{-1})`
		:rtype: numpy.ndarray
		"""

		try:
			ra = rho * cp * (ts - ta) / flux_H
		except ArithmeticError:
			raise ArithmeticError("Aerodynamic resistance for gradient model "
			                      "has not been calculated")

		return ra

	@staticmethod
	def maxT(Rn, G, ra, ta, rho, cp=1012.0):
		"""
		Maximal surface temperature calculated on physical basis (K).
		
		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray
		:param ra: Aerodynamic resistance for heat transfer :math:`(s.m^{-1})`
		:type ra: numpy.ndarray
		:param ta: Air temperature at level Z :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float
		:return: Maximal surface temperature calculated from energy balance\
		equation :math:`(\SI{}\degreeCelsius)`
		:rtype: numpy.ndarray
		"""

		t_max = (Rn - G) * ra / (rho * cp) + ta
		
		return t_max

	@staticmethod
	def wetT(ta):
		"""
		Surface temperature for wet surface. In this case surface temperature
		for wet surface is equal to air temperature. This statement follows 
		from imagine that the LE = Rn - G and thus H is close to zero. \n

		:param ta: Air temperature at level Z :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray

		:return: Temperature of wet surface.
		:rtype: numpy.ndarray
		"""
		# TODO: definition of surface temperature on basis of ETr calculation

		t_wet = ta
		
		return t_wet
	
	@staticmethod
	def dryT(ts, mask=None):
		"""
		Extraction of temperature for dry surface with no evaporation.
		
		:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
		:type ts: numpy.ndarray
		:param mask: Mask of the area of interest. Number of rows and columns
		should be the same.  Format (1, 0) or (1, nan).
		:type mask: numpy.ndarray
		
		:return: Temperature of dry surface derived from surface temperature \
		layer.
		:rtype: float
		"""

		if mask is not None:
			try:
				ts = ts * mask
			except ArithmeticError:
				warnings.warn("Mask for dryT has not been used", stacklevel=3)

		try:
			filt_ts = ndimage.median_filter(ts, 5)
		except ArithmeticError:
			warnings.warn("Median filter has not been used for estimation of "
			              "dry surface temperature", stacklevel=3)
			filt_ts = ts

		filt_ts = filt_ts[~np.isnan(filt_ts)]
		t_dry = float(np.max(filt_ts))

		return t_dry

	@staticmethod
	def coef_b(t_dry, t_wet, t_max, ta):
		"""
		Coefficient b calculated from temperature gradient.

		:param t_dry: Temperature for dry surface derived from the surface\
		temperature layer :math:`(\SI{}\degreeCelsius)`
		:type t_dry: float, numpy.ndarray
		:param t_wet: Temperature for wet surface derived from the surface\
		temperature layer :math:`(\SI{}\degreeCelsius)`
		:type t_wet: float, numpy.ndarray
		:param t_max: Maximal surface temperature calculated from energy\
		balance :math:`(\SI{}\degreeCelsius)`
		:type t_max: numpy.ndarray
		:param ta: Air temperature at level Z :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray

		:return: Coefficient a calculated from temperature gradient
		:rtype: numpy.ndarray
		"""

		try:
			cb = (t_max - ta) / (t_dry - t_wet)  # minimal temperature Tmin is equal to ta + 273.16 (in K)
		except ArithmeticError:
			raise ArithmeticError("Calculation of coefficient b has not been "
			                      "done.")

		return cb

	@staticmethod
	def coef_a(t_wet, cb):
		"""
		Coefficient a calculated from temperature gradient.

		:param t_wet: Temperature for wet surface :math:`(\SI{}\degreeCelsius)`
		:type t_wet: numpy.ndarray
		:param cb: Coefficient b calculated from temperature gradient
		:type cb: numpy.ndarray

		:return: Coefficient a calculated from temperature gradient
		:rtype: numpy.ndarray
		"""

		try:
			ca = -cb * (t_wet + 273.16)  # minimal temperature Tmin is equal to ta + 273.16 (in K)
		except ArithmeticError:
			raise ArithmeticError("Calculation of coefficient a has not been "
			                      "done.")

		return ca

	def dT(self, ts, ta, Rn, G, ra, rho, cp=1012.0, mask=None):
		"""
		Temperature gradient calculated according to SEBAL
		(Bastiaanssen et al. 1998).
		
		:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
		:type ts: numpy.ndarray
		:param ta: Air temperature at level Z :math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray
		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray
		:param ra: Aerodynamic resistance :math:`(s.m^{-1})` 
		:type ra: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float
		:param mask: Mask of the area of interest. Number of rows and columns
		should be the same.  Format (1, 0) or (1, nan).
		:type mask: numpy.ndarray
		:return: Temperature difference dT calculated according to\
		Bastiaanssen (1998)
		:rtype: numpy.ndarray
		"""
				
		t_wet = self.wetT(ta)
		t_dry = self.dryT(ts, mask)
		t_max = self.maxT(Rn, G, ra, ta, rho, cp)
		cb = self.coef_b(t_dry, t_wet, t_max, ta)
		ca = self.coef_a(t_wet, cb)
		dT = ca + cb * (ts + 273.16)

		return dT

	def aeroSEBAL(self, Uz, ta, ts, z0m, Rn, G, rho, niter=10, Z=200, z1=0.1,
	              z2=2, cp=1012.0, kappa=0.41, mask=None):
		"""
		Calculation of sensible heat flux and surface aerodynamic resistance
		according to Bastiaanssen et al. (1998).

		:param Uz: Wind speed at level Z :math:`(m.s^{-1})`
		:type Uz: numpy.ndarray
		:param ta: Air temperature at blending height\
		:math:`(\SI{}\degreeCelsius)`
		:type ta: numpy.ndarray
		:param ts: Surface temperature :math:`(\SI{}\degreeCelsius)`
		:type ts: numpy.ndarray
		:param z0m: Surface roughness for momentum transfer (m)
		:type z0m: numpy.ndarray
		:param Rn: Total net radiation :math:`(W.m^{-2})`
		:type Rn: numpy.ndarray
		:param G: Ground heat flux :math:`(W.m^{-2})`
		:type G: numpy.ndarray
		:param rho: Specific air density :math:`(g.m^{-3})`
		:type rho: numpy.ndarray
		:param niter: Number of iteration
		:type niter: int
		:param Z: Blending height (mixing layer height) (m). Default 200 m.
		:type Z: float (Numpy array)
		:param z1: First height above zero plane displacement (m)
		:type z1: float
		:param z2: Second height above zero plane displacement (m)
		:type z2: float
		:param cp: Thermal heat capacity of dry air :math:`(K.kg^{-1}.K^{-1})`
		:type cp: float
		:param kappa: von Karman constant. Default 0.41
		:type kappa: float
		:param mask: Mask of the area of interest. Number of rows and columns
		should be the same.  Format (1, 0) or (1, nan).
		:type mask: numpy.ndarray

		:return: Aerodynamic resistance for heat transfer :math:`(s.m^{-1})`\
		calculated according to SEBAL (Bastiaanssen et al. 1998)
		:rtype: numpy.ndarray
		:return: Sensible heat flux :math:`(W.m^{-2})` calculated according\
		to SEBAL (Bastiaanssen et al. 1998)
		:rtype: numpy.ndarray
		"""

		# Settings and definition of returns
		ignore_zero = np.seterr(all="ignore")

		# initial part before iterative procedure:
		# psiH definition for initial part of calculation
		psiH_z1 = np.zeros_like(Uz)
		psiH_z2 = np.zeros_like(Uz)

		# friction velocity for meteostation
		try:
			z0m_init = 0.12 * 0.123
			frict = self.frictVelo(Uz, z0m_init, Z=Z)
			# ra for meteostation (neutral stab)
			ra = self.raSEBAL(frict, psiH_z1, psiH_z2, z1, z2, kappa)
		except ArithmeticError:
			raise ArithmeticError("Initial calculation of aerodynamic "
			                      "resistance and sensible heat flux using "
			                      "SEBAL has not finished.")

		# Iterative procedure:
		try:
			for i in range(niter):
				diffT = self.dT(ts, ta, Rn, G, ra, rho, cp, mask)
				H = self.fluxHAer(ra, rho, diffT, cp)
				L = self.lengthMO(frict, ts, H, rho)
				X_z = self.coefX(Z, L)
				X_z1 = self.coefX(z1, L)
				X_z2 = self.coefX(z2, L)
				psiM_z = self.psiM(L, X_z, Z)
				psiH_z1 = self.psiH(L, X_z1, z1)
				psiH_z2 = self.psiH(L, X_z2, z2)
				frict = self.frictVelo(Uz, z0m, Z=Z, psi_m=psiM_z, kappa=kappa)
				ra = self.raSEBAL(frict, psiH_z1, psiH_z2, z1, z2, kappa)
		except ArithmeticError:
			raise ArithmeticError("Aerodynamic resistance and sensible latent "
			                      "heat flux has not been calculated using "
			                      "SEBAL")

		return ra, H
