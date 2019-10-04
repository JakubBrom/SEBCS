#  /***************************************************************************
#  test_windStability.py
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

from unittest import TestCase
import numpy as np
import matplotlib.pyplot as plt
from sebcs.scripts.SEBCSlib import WindStability

ws = WindStability()


class TestWindStability(TestCase):

	h_effect = np.array((0.18, 0.45, 1.06, 0.76, 11.52, 0.31, 24.56, 0.32))
	red_band = np.array((0.09, 0.01, 0.02, 0.07, 0.02, 0.08, 0.01, 0.06))
	nir_band = np.array((0.53, 0.15, 0.29, 0.53, 0.22, 0.52, 0.13, 0.53))
	alb = np.array((0.23, 0.06, 0.12, 0.22, 0.09, 0.22, 0.05, 0.22))
	savi = np.array((0.58, 0.31, 0.5, 0.63, 0.41, 0.6, 0.27, 0.64))
	ndvi = np.array((0.7, 0.85, 0.86, 0.78, 0.84, 0.73, 0.81, 0.79))

	z0m_z = 1.0
	z0m_m1 = np.array((0.02214, 0.05535, 0.13038, 0.09348, 1.41696, 0.03813,
	                   3.02088, 0.03936))
	z0m_m2 = np.array((0.06657007, 0.01599456, 0.04162131, 0.08129356,
	                   0.02530608, 0.07042951, 0.0131706, 0.09008218))
	z0m_m3 = np.array((1.4116612197468847, 419296.39410866355,
	                   136.7068305156706, 2.1867008338321656,
	                   1225.5538560014168, 1.8521244532279109,
	                   6931393.467897867, 2.538483933348297))

	def test_vegHeight(self):
		self.fail()

	def test_zeroPlaneDis(self):
		self.fail()

	def test_z0m_METRIC(self):
		self.fail()

	def test_z0m(self):
		self.fail()

	def test_z0h(self):
		self.fail()

	def test_windSpeedZ(self):
		self.fail()

	def test_frictVelo(self):
		self.fail()

	def test_virtTemp(self):
		self.fail()

	def test_psiM(self):
		self.fail()

	def test_psiH(self):
		self.fail()

	def test_lengthMO(self):
		self.fail()

	def test_coefX(self):
		self.fail()

	def test_stabCoef(self):
		self.fail()

	def test_raThom(self):
		self.fail()

	def test_raSEBAL(self):
		self.fail()

	def test_raGrad(self):
		self.fail()

	def test_maxT(self):
		self.fail()

	def test_wetT(self):
		self.fail()

	def test_dryT(self):
		self.fail()

	def test_coef_b(self):
		self.fail()

	def test_coef_a(self):
		self.fail()

	def test_dT(self):
		self.fail()

	def test_aeroSEBAL(self):
		self.fail()
