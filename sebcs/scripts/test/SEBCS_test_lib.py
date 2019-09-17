#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Imports
import numpy as np
import matplotlib.pyplot as plt
from sebcs.scripts.SEBCSlib import MeteoFeatures, HeatFluxes, WindStability


if __name__ == "__main__":
	met = MeteoFeatures()
	ht = HeatFluxes()
	win = WindStability()

	# Data
	ta_2m = np.linspace(20.0, 20.0, 1)
	ts = np.linspace(20.0, 35.0, 100)
	Rn = np.linspace(500.0, 500.0, 1)
	G = np.linspace(50.0, 50.0, 1)
	U = np.linspace(2.0, 2.0, 100)
	DMT = np.linspace(500.0, 500.0, 1)
	albedo = np.linspace(0.18, 0.18, 1)
	ndvi = np.linspace(0.9, 0.9, 1)
	savi = np.linspace(0.5, 0.5, 1)
	h_eff = np.linspace(0.6, 0.6, 1)
	Tdry = 45.0
	Twet = 20.0

	# Konstanty
	Z = 200.0
	z1 = 0.1
	z2 = 2.0
	Z_st = 10.0
	cp = 1012.0

	# zakladni meteo
	ta = met.airTemp(Z, Z_st, ta_2m)
	rho = met.airDensity(ta)
	lat = met.latent(ta)
	airP = met.airPress(Z, DMT)

	# vitr
	Uz = win.windSpeedZ(U, Z, Z_st)
	disp = win.zeroPlaneDis(h_eff)
	Z = Z - disp
	z0m = win.z0m(savi=savi)
	z0h = win.z0h(z0m)
	Z_d = win.z_d(Z, disp)
	psi_m, psi_h, frict, L = win.stabCoef(Uz, ta, ts, z0m, z0h, Z)

	X_z1 = win.coefX(z1, L)
	X_z2 = win.coefX(z2, L)
	psiH_z1 = win.psiH(L, X_z2, z1)
	psiH_z2 = win.psiH(L, X_z2, z2)

	ra = win.raThom(Uz, z0m, z0h, psi_m, psi_h)
	ra2, H2 = win.aeroSEBAL(Uz, ta, ts, z0m, Rn, G, rho, Tdry, Twet, Z, z2, z1, niter=10)
	ra3 = win.raSEBAL(frict, z1, z2, psiH_z1, psiH_z2)

	# toky
	Rn_G = Rn - G
	dT = win.dT(ts, ta, Twet, Tdry, Rn, G, ra, rho)
	H = ht.fluxHAer(ra, rho, dT, cp)
	dT3 = win.dT(ts, ta, Twet, Tdry, Rn, G, ra, rho)
	H3 = ht.fluxHAer(ra3, rho, dT3, cp)

	# graf
	# plt.plot(U, H3)
	plt.plot(ts, ra)
	plt.plot(ts, ra2)
	plt.plot(ts, ra3)
	plt.show()
