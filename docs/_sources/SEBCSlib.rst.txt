Programmers documentation
==========================

The SEBCS plugin for QGIS contains the SEBCS_lib.py library, which may interest programmers and people who want to use methods for working with remote sensing data.

SEBCS_lib library
---------------------

SEBCS_lib is a library of SEBCS software (Brom 2013-2023) designed for
geographical and remote sensing purposes, mainly for the surface energy
balance analysis.

SEBCS_lib contains several tools for:

* raster manipulation (import, export, retrieval of coordinate system and metadata information)
* calculations of vegetation indices and aboveground biomass characteristics
* calculations of meteorological features of the boundary layer (physical properties)
* calculation of atmospheric stability characteristics
* calculation of energy fluxes and their characteristics, evapotranspiration intensity and vegetation water stress

Most of these methods work with Numpy arrays.

Access to SEBCS_lib library
---------------------------

The library (package) can be accessed in standard pythonic way using import in the Python console in QGIS. Classes can be imported separately.

.. code-block::

    from SEBCS import SEBCS_lib

The SEBCS_lib can be used within QGIS plugins. If there is an issue with importing SEBCS_lib please reinstall the plug-in.


SEBCS_lib documentation
-----------------------

.. automodule:: SEBCS_lib
    :members:



