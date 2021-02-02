#  /***************************************************************************
#  readGeo_t.py
#
#  TODO: popis
#
#                                -------------------
#          begin                :
#          date                 :
#          git sha              : $Format:%H$
#          copyright            : (C) 2014-2020 Jakub Brom
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

from SEBCS_lib import GeoIO
geo = GeoIO()



if __name__ == '__main__':
	rast = "b5.rst"
	geoattr = geo.readGeo(rast)
	rarray = geo.rasterToArray(rast)
	print(geoattr)
	print(rarray)