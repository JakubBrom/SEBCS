#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Imports
import numpy as np
import matplotlib.pyplot as plt
from sebcs.scripts.SEBCS_lib import MeteoFeatures, HeatFluxes, WindStability, VegIndices


if __name__ == "__main__":
	met = MeteoFeatures()
	vi = VegIndices()
	ht = HeatFluxes()
	win = WindStability()

	# Data
	ta_2m = np.full((20,20), 23.0)
	ts = np.linspace(20,35,400).reshape(20,20)
	Rn = np.full((20,20), 650.0)
	G = np.full((20,20), 50)
	U = np.random.uniform(0.1,8, (20,20))
	# DMT = np.linspace(420.0, 580.0, 100)
	# albedo = np.linspace(0.18, 0.18, 100)
	# ndvi = np.linspace(0.9, 0.9, 1)
	# savi = np.linspace(0.5, 0.5, 1)
	h_eff = np.full((20,20), 0.6)
	LAI = np.full((20,20), 5.0)

	z0m = win.z0m(h_eff, LAI)
	z0h = win.z0h(z0m)

	H_aer, LE_aer, EF_aer, eq_aer, PT_aer, ra_aer, frict_aer = ht.heatFluxes(
		Rn, G, ts, ta_2m, "aero", Uz=U, h_eff=h_eff, LAI=LAI, z0m=z0m, z0h=z0h)

	H_sebal, LE_sebal, EF_sebal, eq_sebal, PT_sebal, ra_sebal, frict_sebal = \
		ht.heatFluxes(Rn, G, ts, ta_2m, "SEBAL", Uz=U, h_eff=h_eff, LAI=LAI)

	H_grad, LE_grad, EF_grad, eq_grad, PT_grad, ra_grad, frict_grad = \
		ht.heatFluxes(Rn, G, ts, ta_2m, "grad")


	f1 = plt.figure()
	f2 = plt.figure()
	f3 = plt.figure()
	ax1 = f1.add_subplot(111)
	ax2 = f2.add_subplot(111)
	ax3 = f3.add_subplot(111)
	ax1.scatter(U, PT_aer)
	# ax1.set_ylim([0,100])
	# ax2.imshow(PT_aer)
	# ax3.imshow(H_grad)
	plt.show()
