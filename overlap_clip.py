#!/usr/bin/env python
# -*- coding: cp1250 -*-

# ------------------------------------------------------------------------------
#  Author: Jakub Brom
#  Date: 2020 - 9 - 21
#
#  Copyright (c)  Jakub Brom, 2020 - 2020.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------------

import shutil
import tempfile
import os
import warnings
import copy
import numpy as np
from osgeo import gdal, osr


def readGeo(rast):
    """
    Reading geographical information from raster using GDAL.

    :param rast: Path to raster file in GDAL accepted format.
    :type rast: str

    :returns: List of geotransformation parameters and features of
        an input raster:
            0: The affine transformation coefficients; tuple
            1: Projection information of the raster (dataset); str
            2: Pixel width (m) on X axis; float
            3: Pixel height (m) on Y axis; float
            4: EPSG Geodetic Parameter Set code; str
            5: Coordinate system; str
            6: Number of columns in the raster; int
            7: Number of rows in the raster; int
            8: SRS - OSR spatial reference; object
    :rtype: list
    """

    try:
        ds = gdal.Open(rast)

        gtransf = ds.GetGeoTransform()
        prj = ds.GetProjection()
        x_size = gtransf[1]
        y_size = gtransf[5] * (-1)
        cols = ds.RasterXSize
        rows = ds.RasterYSize

        srs = osr.SpatialReference(wkt=prj)
        if srs.IsProjected:
            EPSG = srs.GetAttrValue("authority", 1)
            geogcs = srs.GetAttrValue("geogcs", 0)
            loccs = srs.GetAttrValue("local_cs", 0)
        else:
            EPSG = None
            geogcs = None
            loccs = None

        del ds

    except IOError:
        raise IOError("Geographical information have not been readed.")

    return gtransf, prj, x_size, y_size, EPSG, geogcs, loccs, cols, \
           rows, srs


def fixSRS(rasters_in, tmp_out=True):
    """
    Checking and fixing the SRSs of input rasters. Function looks
    for possible definition of SRS in group of rasters --> local SRS
    is possibly replaced by global SRS if it is available in the set
    of rasters. The empty SRS is replaced by most frequent SRS
    in the set of rasters.

    :param rasters_in: List of input rasters paths.
    :type rasters_in: list
    :param tmp_out: Selection if the output layers will be created
                as temporal scratch layers or the original raster
                layers will be replaced by newly defined rasters.
    :type tmp_out: bool
    :return: List of output layers paths
    :rtype: list
    """

    # Manage output list
    rasters_out = copy.copy(rasters_in)

    # 1. Get rasters info
    EPSG_list = []
    LLC_X_list = []
    geoGCS_list = []
    localCS_list = []

    for i in rasters_in:
        ginfo = readGeo(i)
        cols = ginfo[7]
        EPSG_list.append(ginfo[4])
        geoGCS_list.append(ginfo[5])
        localCS_list.append(ginfo[6])
        LLC_X_list.append(ginfo[0][0] + ginfo[2] * cols)

    # 2. Check EPSGs
    # Convert EPSG to int
    for i in range(len(EPSG_list)):
        if EPSG_list[i] is None:
            EPSG_list[i] = 0
    EPSGs_int = np.array(EPSG_list).astype(int)

    # Check differences in the EPSGs
    if np.min(EPSGs_int) is not np.max(EPSGs_int):
        # Make temporary folder
        try:
            tmp_folder = tempfile.mkdtemp()
        except OSError as err:
            raise err

        # 3. Test of the layers with local CS - replace localCS with
        # geoCS
        localCS_index = [i for i, x in enumerate(localCS_list) if x is
                         not None]

        for i in localCS_index:
            # Name of local CS which can cover geo CS
            localCS_name = localCS_list[i][0:5]

            try:
                # Search for the projection in geoCS which
                # corresponds with localCS of the tested image
                geoCS_notNone = [(i, x) for i, x in
                                 enumerate(geoGCS_list) if x is
                                 not None]
                geo_index = [i for i, x in geoCS_notNone if x[0:5] ==
                             localCS_name][0]
                # Fill EPSG for a layer
                EPSG_list[i] = EPSG_list[geo_index]
                # create temporary output file
                tmp_file = tempfile.NamedTemporaryFile(
                    suffix=".tif", dir=tmp_folder, delete=False)
                tmp_file_name = tmp_file.name
                # convert SRS from local to global --> create
                # temporary file
                gdal.Warp(tmp_file_name, rasters_out[i], srcSRS=str(
                    "EPSG:" + EPSG_list[i]), warpOptions=[
                    "NUM_THREADS=ALL_CPUS"], multithread=True)
                # replace path of file with local projection by the
                # temporary file with correct projection
                print("SRS for layer " + rasters_out[i] + " has been "
                                                          "fixed.")
                rasters_out[i] = tmp_file_name
                geoGCS_list[i] = geoGCS_list[geo_index]
                localCS_list[i] = None

            except Warning:
                warnings.warn("The SRS for layer " + rasters_out[i] +
                              " has not been found. We will try to "
                              "find out another way how to fix SRS.",
                              stacklevel=3)

        # 4. Testing layers with missing EPSG on basis of coordinates
        noEPSG_index = [i for i, x in enumerate(geoGCS_list) if
                        x is None]

        # Searching the same coordinates and replacing SRS
        for i in noEPSG_index:
            # search for the same coordinates
            coord = LLC_X_list[i]

            try:
                geoCS_notNone = [i for i, x in enumerate(geoGCS_list) if
                                 x is
                                 not None]
                # Testing coordinates with layers with existing CS
                for j in geoCS_notNone:

                    if LLC_X_list[j] == coord:
                        EPSG_list[i] = EPSG_list[j]

                        tmp_file = tempfile.NamedTemporaryFile(
                            suffix=".tif", dir=tmp_folder, delete=False)
                        tmp_file_name = tmp_file.name
                        # convert SRS --> create temporary file
                        gdal.Warp(tmp_file_name, rasters_out[i],
                                  srcSRS=str("EPSG:" + EPSG_list[i]),
                                  warpOptions=["NUM_THREADS=ALL_CPUS"],
                                  multithread=True)
                        # replace path of file with local projection
                        # by the temporary file with correct projection
                        print("SRS for layer " + rasters_out[i] +
                              " has been fixed.")
                        rasters_out[i] = tmp_file_name
                        geoGCS_list[i] = geoGCS_list[j]

                        break

            except RuntimeError:
                raise RuntimeError("Ooops...! Something wrong "
                                   "happened.")

        # 5. Converting remaining files to prevailing EPSG
        noEPSG_index = [i for i, x in enumerate(geoGCS_list) if
                        x is None]
        EPSG_maxfreq = max(set(EPSG_list), key=EPSG_list.count)

        for i in noEPSG_index:
            warnings.warn("Layer " + rasters_out[i] + " has still no "
                          "SRS. The EPSG: " + EPSG_maxfreq + " will "
                                                             "be used.",
                          stacklevel=3)
            try:
                tmp_file = tempfile.NamedTemporaryFile(
                    suffix=".tif", dir=tmp_folder, delete=False)
                tmp_file_name = tmp_file.name
                # convert SRS --> create temporary file
                gdal.Warp(tmp_file_name, rasters_out[i],
                          srcSRS=str("EPSG:" + EPSG_maxfreq),
                          warpOptions=["NUM_THREADS=ALL_CPUS"],
                          multithread=True)
                # replace path of file with local projection by the
                # temporary file with correct projection
                print("SRS for layer " + rasters_out[i] + " has "
                                                          "been fixed.")
                rasters_out[i] = tmp_file_name
                EPSG_list[i] = EPSG_maxfreq

            except RuntimeError:
                raise RuntimeError("Ooops...! Something wrong "
                                   "happened.")

        if tmp_out is False:
            warnings.warn("All the files with incorrect SRS have been "
                          "replaced by files with fixed SRS. Some "
                          "errors may occur!", stacklevel=3)
            for i in range(len(rasters_out)):
                if rasters_out[i] is not rasters_in[i]:
                    shutil.copy(rasters_out[i], rasters_in[i])
            rasters_out = rasters_in
            # remove temporary files
            shutil.rmtree(tmp_folder)

    else:
        print("No SRS has been changed.")

    return rasters_out

def uniformSRS(rasters_in, epsg=None, tmp_out=True):
    """
    The function unifies SRSs of a group of input raster layers. If
    output EPSG is not defined, the most frequent EPSG is defined as
    default. If SRSs are not included in the raster datasets either
    set EPSG or WGS 84 (EPSG:4326; i.e. if epsg=None) is used for
    transformation of SRS.

    :param rasters_in: List of input rasters paths.
    :type rasters_in: list
    :param epsg: EPSG definition for output rasters. If epsg=None
                the most frequent EPSG in the layers group will be set.
    :type epsg: int
    :param tmp_out: Selection if the output layers will be created
                as temporal scratch layers or the original raster
                layers will be replaced by newly defined rasters.
    :type tmp_out: bool

    :return: List of output layers paths
    :rtype: list
    """

    # TODO: poresit vyjimky - melo by smysl, aby bylo namisto vyjimky
    #  jen varovani, s tim, ze by z procedury vypadly puvodni
    #  neupraveny data?
    # TODO: pokud GDAL neumi pracovat s nekterymi epsg (treba
    #  S-JTSK), tak by mozna bylo dobry zkusit QGIS QgsGeometryEngine
    #  (https://qgis.org/api/classQgsGeometryEngine.html)
    # TODO: poresit mazani dat z tmp -

    # GDAL exceptions
    gdal_ver = gdal.__version__
    gdal.UseExceptions()

    # Manage output list
    rasters_orig = copy.copy(rasters_in)
    rasters_out = copy.copy(rasters_in)

    # 1. Get rasters info
    EPSG_list = []
    geoGCS_list = []

    for i in rasters_in:
        ginfo = readGeo(i)
        EPSG_list.append(ginfo[4])
        geoGCS_list.append(ginfo[5])

    # Make temporary folder
    try:
        tmp_folder = tempfile.mkdtemp()
    except OSError as err:
        raise err

    # 2. Check SRSs
    # Convert EPSG to int
    for i in range(len(EPSG_list)):
        if EPSG_list[i] is None:
            EPSG_list[i] = 0
    EPSGs_int = np.array(EPSG_list).astype(int)

    # list of layers without SRS
    noEPSG_index = [i for i, x in enumerate(geoGCS_list) if
                    x is None]

    # 3. Check differences in the EPSGs
    if np.min(EPSGs_int) is not np.max(EPSGs_int):  # --> different SRS

        # 4. testing if there are no EPSGs in the raster list
        if len(noEPSG_index) > 0:
            # Fix SRS
            rasters_in = fixSRS(rasters_in, tmp_out=True)
            for i in rasters_in:
                ginfo = readGeo(i)
                EPSG_list[rasters_in.index(i)] = ginfo[4]
                geoGCS_list[rasters_in.index(i)] = ginfo[5]
        else:
            pass

    # 5. Testing if no projection is set up
    else:
        if len(noEPSG_index) > 0:
            warnings.warn("The projection is not defined in any layer "
                          "or only local projection is defined in all"
                          "cases. The continuing process may cause "
                          "an error in projection settings!",
                          stacklevel=3)
        else:
            pass

    # 6. Is the input EPSG set up?
    # EPSG is not set
    if epsg is None:
        # list of layers without SRS
        noEPSG_index = [i for i, x in enumerate(geoGCS_list) if
                        x is None]

        # Empty SRS for all the rasters --> WGS 84:
        if len(noEPSG_index) == len(EPSG_list):
            warnings.warn("The WGS 84 projection (EPSG:4326) will be "
                          "set for all layers.", stacklevel=3)
            try:
                for i in range(len(EPSG_list)):
                    tmp_file = tempfile.NamedTemporaryFile(
                        suffix=".tif", dir=tmp_folder, delete=False)
                    tmp_file_name = tmp_file.name
                    # convert SRS --> create temporary file
                    gdal.Warp(tmp_file_name, rasters_in[i],
                              srcSRS=str("EPSG:4326"),
                              warpOptions=["NUM_THREADS=ALL_CPUS"],
                              multithread=True)
                    # Write new files to the output layers path list
                    rasters_out[i] = tmp_file_name
            except RuntimeError:
                raise RuntimeError("The input SRSs of rasters are"
                                   "probably not supported by GDAL "
                                   "version {gdv}".format(gdv=gdal_ver))

        # EPSG with the highest frequency
        else:
            EPSG_maxfreq = max(set(EPSG_list), key=EPSG_list.count)
            warnings.warn("The projection EPSG:" + str(EPSG_maxfreq) +
                          " will be set for all layers.",
                          stacklevel=3)

            for i in range(len(EPSG_list)):
                try:
                    tmp_file = tempfile.NamedTemporaryFile(
                        suffix=".tif", dir=tmp_folder, delete=False)
                    tmp_file_name = tmp_file.name
                    # convert SRS --> create temporary file
                    gdal.Warp(tmp_file_name, rasters_in[i],
                              dstSRS=str("EPSG:" + str(EPSG_maxfreq)),
                              warpOptions=["NUM_THREADS=ALL_CPUS"],
                              multithread=True)
                    # Write new files to the output layers path list
                    rasters_out[i] = tmp_file_name

                except RuntimeError:
                    raise RuntimeError("The input SRSs of rasters are"
                                       "probably not supported by GDAL "
                                       "version {gdv}".format(
                                        gdv=gdal_ver))

    # EPSG is set
    else:
        # All layers have the same EPSG corresponding to input EPSG
        EPSG_cor = [x for x in EPSG_list if x is epsg]
        if len(EPSG_cor) is len(EPSG_list):
            rasters_out = rasters_in

        # Different or various SRS
        else:
            for i in range(len(EPSG_list)):
                try:
                    tmp_file = tempfile.NamedTemporaryFile(
                        suffix=".tif", dir=tmp_folder, delete=False)
                    tmp_file_name = tmp_file.name
                    # convert SRS --> create temporary file
                    gdal.Warp(tmp_file_name, rasters_in[i],
                              srcSRS=str("EPSG:" + str(EPSG_list[i])),
                              dstSRS=str("EPSG:" + str(epsg)),
                              warpOptions=["NUM_THREADS=ALL_CPUS"],
                              multithread=True)
                    # Write new files to the output layers path list
                    rasters_out[i] = tmp_file_name

                except RuntimeError:
                    raise RuntimeError("The input SRSs of rasters are"
                                       "probably not supported by GDAL "
                                       "version {gdv}".format(
                                        gdv=gdal_ver))

    # If temporary model is not used, replace original files by
    # results and remove temporary files
    if tmp_out is False:
        for i in range(len(rasters_out)):
            if rasters_out[i] is not rasters_orig[i]:
                shutil.copy(rasters_out[i], rasters_orig[i])
        rasters_out = rasters_orig
        # remove temporary files
        shutil.rmtree(tmp_folder)

    return rasters_out


# noinspection PyUnboundLocalVariable
def clipOverlappingArea(rasters_in, output_folder=None,
                        suffix='_clipped',
                        epsg=None, tmp_out=False):
    """
    Function for clipping rasters by their overlapping area.

    :param rasters_in: List of input rasters paths.
    :type rasters_in: list
    :param output_folder: Path to output folder defined by user
    :type output_folder: str
    :param suffix: Suffix of the output rasters names. The name of
    new raster is constructed from original name of raster and suffix
    :type suffix: str
    :param epsg: EPSG definition for output rasters. If epsg=None
                the most frequent EPSG in the layers group will be set.
    :type epsg: int
    :param tmp_out: If True the output layers will be created
                as temporal scratch layers.
    :type tmp_out: bool

    :return: Paths of rasters clipped by overlapping area of all the
    rasters
    :rtype: list
    """

    # 1. Unify EPSG of data
    rasters_unif = uniformSRS(rasters_in, epsg, tmp_out=True)

    # Make temporary folder for outputs
    if tmp_out is True:
        try:
            tmp_folder = tempfile.mkdtemp()
        except OSError as err:
            raise err

    # 2. Looking for coordinates of overlapping area
    # Lists of coordinates:
    EPSG_list = []
    ULC_X_list = []  # upper left corner Xmin
    ULC_Y_list = []  # upper left corner Ymax
    LRC_X_list = []  # lower right corner Xmax
    LRC_Y_list = []  # lower right corner Ymin
    X_size_list = []
    Y_size_list = []

    # extract geoinformation
    for i in rasters_unif:
        ginfo = readGeo(i)
        EPSG_list.append(ginfo[4])
        cols = ginfo[7]
        rows = ginfo[8]
        ULC_X_list.append(ginfo[0][0])
        LRC_X_list.append(ginfo[0][0] + ginfo[2] * cols)
        ULC_Y_list.append(ginfo[0][3])
        LRC_Y_list.append(ginfo[0][3] - ginfo[3] * rows)
        X_size_list.append(ginfo[2])
        Y_size_list.append(ginfo[3])

    # extract the coordinates and pixel size for overlapping area
    ULC_X_ovlp = np.max(ULC_X_list)
    LRC_X_ovlp = np.min(LRC_X_list)
    ULC_Y_ovlp = np.min(ULC_Y_list)
    LRC_Y_ovlp = np.max(LRC_Y_list)

    X_size_ovlp = np.min(X_size_list)
    Y_size_ovlp = np.min(Y_size_list)

    # 3. Rescale all the rasters to the same spatial extent (the same
    # pixel size) and clip all the layers by the overlapping area
    out_files_list = []
    for i in range(0, len(rasters_unif)):
        # Define paths of the output files
        if tmp_out is False:
            rast_in_name = os.path.splitext(rasters_in[i])[0]
            rast_name = rast_in_name + suffix + ".tif"
            out_file = os.path.join(output_folder, rast_name)
            out_files_list.append(out_file)

        else:
            tmp_file = tempfile.NamedTemporaryFile(
                suffix=".tif", dir=tmp_folder, delete=False)
            out_file = tmp_file.name
            out_files_list.append(out_file)

        # Clip files and change resolution
        gdal.Warp(out_file, rasters_unif[i],
                  outputBounds=(ULC_X_ovlp, LRC_Y_ovlp, LRC_X_ovlp,
                                ULC_Y_ovlp), xRes=X_size_ovlp,
                  yRes=Y_size_ovlp, warpOptions=[
                "NUM_THREADS=ALL_CPUS"], multithread=True)

    # TODO: je potreba otestovat, jestli nějakej raster není mimo
    #  rozsah ostatnich rastru
    # TODO: poresit nedokonale prekryvy, tzn, ze je prekryvna plocha
    #  napr. prazdna nebo tak podobne
    # TODO: zkontrolovat, jestli jsou souradnice stejne pro vsechny
    #  vrstvy. Pokud ano, tak uz dal neprovadet vypocet -->
    #  out_files_list = rasters_in
    # TODO: pridat vyjimky

    return out_files_list


if __name__ == '__main__':
    pass
    # TODO: add possible access from terminal
