#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Imports
import numpy as np
import matplotlib.pyplot as plt
from sebcs.scripts.SEBCSlib import MeteoFeatures, HeatFluxes, WindStability, VegIndices


if __name__ == "__main__":
	met = MeteoFeatures()
	ht = HeatFluxes()
	win = WindStability()
	vi = VegIndices()

	# Data
	ta_2m = np.linspace(20.0, 25.0, 100)
	ts = np.linspace(20.0, 35.0, 100)
	Rn = np.linspace(500.0, 500.0, 100)
	G = np.linspace(50.0, 50.0, 100)
	U = np.linspace(2.0, 2.0, 100)
	DMT = np.linspace(500.0, 500.0, 1)
	albedo = np.linspace(0.18, 0.18, 1)
	ndvi = np.linspace(0.9, 0.9, 1)
	savi = np.linspace(0.5, 0.5, 1)
	h_eff = np.linspace(1.2, 1.2, 1)

	h_effect = np.array((0.18, 0.45, 1.06, 0.76, 11.52, 0.31, 24.56, 0.32))
	red_band = np.array((0.09, 0.01, 0.02, 0.07, 0.02, 0.08, 0.01, 0.06))
	nir_band = np.array((0.53, 0.15, 0.29, 0.53, 0.22, 0.52, 0.13, 0.53))
	alb = np.array((0.23, 0.06, 0.12, 0.22, 0.09, 0.22, 0.05, 0.22))
	savi = vi.viSAVI(red_band, nir_band, 0.5)
	ndvi = vi.viNDVI(red_band, nir_band)

	LAI = vi.LAI(red_band, nir_band, 1)
	LAI2 = vi.LAI(red_band, nir_band, 2)
	LAI3 = vi.LAI(red_band, nir_band, 3)
	LAI4 = vi.LAI(red_band, nir_band, 4)
	LAI5 = vi.LAI(red_band, nir_band, 5)
	LAI6 = vi.LAI(red_band, nir_band, 6)
	LAI7 = vi.LAI(red_band, nir_band, 7)


	# Konstanty
	# Z = 200.0
	# z1 = 0.1
	# z2 = 2.0
	# Z_st = 10.0
	# cp = 1012.0

	# zakladni meteo
	# ta = met.airTemp(Z, Z_st, ta_2m)
	# rho = met.airDensity(ta)
	# lat = met.latent(ta)
	# airP = met.airPress(Z, DMT)
	# delta = met.delta(ts, ta)
	# gama = met.gamma(airP, lat)
	#
	# # vitr
	# Uz = win.windSpeedZ(U, Z, Z_st)
	# disp = win.zeroPlaneDis(h_eff)
	# Z = Z - disp
	# z0m_vysl = win.z0m(method=3, h_eff=h_effect, red=red_band,
	#                   nir=nir_band, albedo=alb)
	# print(z0m_vysl)
	# z0h = win.z0h(z0m)
	# Z_d = win.z_d(Z, disp)
	# psi_m, psi_h, frict, L = win.stabCoef(Uz, ta, ts, z0m, z0h, Z)
	#
	# X_z1 = win.coefX(z1, L)
	# X_z2 = win.coefX(z2, L)
	# psiH_z1 = win.psiH(L, X_z2, z1)
	# psiH_z2 = win.psiH(L, X_z2, z2)
	#
	# ra = win.raThom(Uz, z0m, z0h, psi_m, psi_h)
	# ra2, H2 = win.aeroSEBAL(Uz, ta, ts, z0m, Rn, G, rho, Z, z2, z1, niter=10)
	# ra3 = win.raSEBAL(frict, z1, z2, psiH_z1, psiH_z2)
	#
	# # toky
	# Rn_G = Rn - G
	# dT = win.dT(ts, ta, Rn, G, ra, rho)
	# H = ht.fluxHAer(ra, rho, dT, cp)
	# dT3 = win.dT(ts, ta, Rn, G, ra, rho)
	# H3 = ht.fluxHAer(ra3, rho, dT3, cp)
	# LE_eq = ht.fluxLE_EQ(Rn, G, delta, gama)
	#
	# H4 = Rn_G - LE_eq * 1.26
	# t_wet = H4 * ra / (rho*cp) + ta
	#
	# graf
	# plt.plot(U, H3)

	# plt.scatter(savi, h_effect)

	plt.scatter(savi, LAI)
	plt.scatter(savi, LAI2)
	# plt.scatter(savi, LAI3)
	# plt.scatter(savi, LAI4)
	# plt.scatter(savi, LAI5)
	# plt.scatter(savi, LAI6)
	plt.scatter(savi, LAI7)
	# plt.plot(ts, H2)
	# plt.plot(ts, ra3)
	plt.show()
