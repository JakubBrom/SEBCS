#  /***************************************************************************
#  SEBCS_main.py
#
#  Collection of functions and preparation of outputs for SEBCS
#  calculation procedure, e.g. data imports, miscellaneous data calculation,
#  preparation of outputs etc. Module SEBCS_main.py is new version of
#  previous SEBCS_SA.py module which will be (has been) removed.
#
#                                -------------------
#          begin                : 19-10-04
#          date                 : 19-10-04
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

# TODO vstupy:

self.red_path, self.nir_path, self.thermal_path, self.dmt_path, self.ta_path,
				self.glob_rad, self.humid, self.latitude, self.longitude, self.acq_date, self.acq_time,
				self.out_lyrs_fnames, self.out_folder_path, self.out_file_name,
				self.rb_method, self.rb_sat, self.emis_rule, self.out_driver, self.multiband,
				self.blue_path, self.green_path, self.swir1_path, self.swir2_path,
				self.wind_path, self.hmin_path, self.hmax_path, self.mask_path, self.albedo_path, self.canopy_path,
				self.hwind, Z = 200, cp = 1012

# TODO výstupy:
out_names_all = ("Rs_in", "Rs_out", "albedo", "RL_in", "RL_emit", "Rn",
                 "Ts", "emis", "LE", "LEp", "LE_PT", "H", "G", "E_int",
                 "EF", "Bowen", "omega", "CWSI", "frict_U", "ra",
                 "rc", "NDVI", "MSAVI", "NDMI", "slope", "aspect")
out_lyrs_all = (self.in_Rs, self.out_Rs, self.albedo, self.RL_in, self.RL_out,
                self.Rn, self.ts_C, self.emis, self.LE, self.LEp, self.LE_PT,
                self.H, self.G, self.E_int, self.EF, self.bowen, self.omega,
                self.CWSI, self.frict, self.ra, self.rc, self.ndvi,
                self.msavi, self.ndmi, self.slope, self.aspect)
