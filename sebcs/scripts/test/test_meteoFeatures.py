"""
   /***************************************************************************
   test_meteoFeatures.py

   Test for methods form MeteoFeatures class

                                -------------------
          begin                : 19-09-20
          date                 : 19-09-20
          git sha              : $Format:%H$
          copyright            : (C) 2014-2019 Jakub Brom
          email                : jbrom@zf.jcu.cz

   ***************************************************************************/
   /***************************************************************************

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License  as published  by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License along
   with this program.  If not, see <https://www.gnu.org/licenses/>.

   ***************************************************************************/
"""


from unittest import TestCase
import numpy as np
import matplotlib.pyplot as plt
from sebcs.scripts.SEBCS_lib import MeteoFeatures

met = MeteoFeatures()


class TestMeteoFeatures(TestCase):
	ta = 22
	ta_st = 25.0
	DMT = 350
	Z = 200
	Rh = 57

	# Results
	ta_res = 23.713
	pressure = 95.111
	e_sat = 2.642956
	e_air = e_sat * Rh/100

	def test_airTemp(self):
		self.T = met.airTemp(self.ta_st, self.Z)
		self.assertEqual(self.ta_res, self.T)

	def test_airPress(self):
		self.pres = met.airPress(self.ta, self.DMT, self.Z)
		self.assertAlmostEqual(self.pressure, self.pres, 3)

	def test_satVapourPress(self):
		self.e_s = met.satVapourPress(self.ta)
		self.assertAlmostEqual(self.e_sat, self.e_s, 5)

	def test_vapourPress(self):
		self.eair = met.vapourPress(self.e_sat, self.Rh)
		self.assertAlmostEqual(self.e_air, self.eair)

	def test_vpd(self):
		self.fail()

	def test_airDensity(self):
		self.fail()

	def test_latent(self):
		self.fail()

	def test_gamma(self):
		self.fail()

	def test_delta(self):
		self.fail()

	def test_emissivity(self):
		self.fail()
