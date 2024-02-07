SEBCS usage
=============

The SEBCS for QGIS plug-in allows spatial calculation of energy balance components and other surface characteristics. The plug-in uses three different approaches to calculate the actual latent heat flux. These are the aerodynamic method according to Penman-Monteith, the SEBAL method according to `Bastiaanssen et al. (1998) <https://www.sciencedirect.com/science/article/abs/pii/S0022169498002534>`_ and the gradient method based on the formulation of `Suleiman & Crago (2004) <https://acsess.onlinelibrary.wiley.com/doi/full/10.2134/agronj2004.3840>`_.

Inputs
------

The computation of the variables requires some data inputs. These are raster remote sensing data (Spectral data), raster data describing the surface (Additional data) and numerical data, including meteorological data and metadata. The input data varies depending on the calculation method used. An overview of the input data depending on the method used is given in the following table.

.. csv-table::
    :header: "Data","PM","SEBAL","Gradient","Optional","Comment"

        **Spectral raster data**
        Blue band,x,x,x,x
        Green band,x,x,x,x
        Red band,x,x,x
        NIR band,x,x,x
        SWIR1 band (~1.61μm),x,x,x,x
        SWIR2 band (~2.2μm),x,x,x,x
        Surface temperature (°C),x,x,x,x
        **Additional raster data**
        Digital model of terrain (m a.s.l),x,x,x
        Air temperature (°C),x,x,x
        Wind speed (m/s),x,x
        Canopy height(m),x,x,,,If known
        Max. canopy height (m),x,x,,,If canopy height layer is not available
        Min. canopy height(m),x,x,,,If canopy height layer is not available
        Mask,x,x,x,x
        Albedo (rel.),x,x,x,x
        **Numerical data**
        Global radiation (:math:`W.m^{-2}`),x,x,x
        Air relative humidity (%),x,x,x
        Height of wind speed measurement (m),x,x
        Date of data acquisition,x,x,x
        Time of data acquisition (GMT),x,x,x
        Mean latitude of data,x,x,x,,Can be calculated directly from data
        Mean longitude of data,x,x,x,,Can be calculated directly from data


All raster data must have the same size, i.e. the same number of columns and rows, and the same geographic coordination system!


Outputs
-------

The SEBCS for QGIS plug-in allows calculation of the following features described in the next table. The calculation of the features is described :doc:`here <calculation>`.

.. csv-table::
    :header: Feature, Unit

        Incomming shortwave radiation,:math:`W.m^{-2}`
        Reflected shortwave radiation,:math:`W.m^{-2}`
        Albedo,rel.
        Incomming longwave radiation,:math:`W.m^{-2}`
        Outgoing (emitted) longwave radiation,:math:`W.m^{-2}`
        Total net radiation,:math:`W.m^{-2}`
        Surface temperature,°C
        Surface emissivity,rel.
        Latent heat flux,:math:`W.m^{-2}`
        Latent heat flux - Penman potential,:math:`W.m^{-2}`
        Latent heat flux - Priestley-Taylor potential,:math:`W.m^{-2}`
        Sensible heat flux,:math:`W.m^{-2}`
        Ground heat flux,:math:`W.m^{-2}`
        Evapotranspiration intensity,:math:`mmol.m^{-2}.s^{-1}`
        Evaporative fraction,rel.
        Bowen ratio,unitless
        Decoupling coefficient,rel.
        Crop Water Stress Index (CWSI),unitless
        Friction velocity,:math:`m.s^{-1}`
        Aerodynamic surface resistance,:math:`s.m^{-1}`
        Surface resistance for water transfer,:math:`s.m^{-1}`
        Normalized Difference Vegetation Index (NDVI),unitless
        Modified Soil Adjusted Vegetation Index (MSAVI),unitless
        Normalized Difference Moisture Index (NDMI),unitless
        Soil Adjusted Vegetation Index (SAVI),unitless
        Leaf Area Index (LAI),:math:`m^{2}.m^{-2}`
        Slope,°
        Aspect,°


Input data preparation
----------------------

SEBCS for QGIS is primarily designed to use Landsat satellite data, but also allows the use of other spatial data sources. Regardless of the data source, the same spatial resolution (number of pixels) and spatial extent (same area) must be maintained for all layers used. The choice of the geographical reference system does not affect the function of the module.
Data preparation is specific for some data inputs. The procedure is as follows:


Optical data
............
The optical data must be used in the form of the surface spectral reflectance as relative values, in the range 0 to 1.


Air temperature
...............
The air temperature is the temperature measured by a weather station at the reference height :math:`z`. It is usually the temperature measured at a height of 2 m above the surface. Air temperature can be understood as either homogeneous or heterogeneous in a given space, depending on the measurement capabilities. For this reason, the air temperature is specified as a temperature layer for a given area.
Assuming a homogeneous temperature distribution in space, the layer has only one value, which can be created using the Create constant raster layer function in QGIS Processing module or the Raster Calculator.
In the case of the assumption of a heterogeneous air temperature distribution in the space, the air temperature needs to be modeled. If a sufficient number of measurements in a given area is available, various interpolation and geostatistical methods can be used to create an air temperature layer. A possible way to create an air temperature map is to assume an adiabatic change in air temperature with altitude (e.g. a decrease of 0.6 °C per 100 m altitude). In this case, it is necessary to have a digital model of terrain (elevation map converted to raster form, see below) and to know the altitude of the weather station location. The surface temperature is then calculated using the Raster Calculator according to the formula:

.. math::
    :label: eq:Ta_dmt

    T_a = T_{st} + \Gamma(Z_{st}-DMT)


Wind speed
..........
The wind speed is given as a layer of wind velocity at height :math:`Z_{st}` in :math:`m \cdot s^{-1}`. The situation is similar to that of surface temperature. The user can use a raster layer with constant values or a heterogeneous raster airflow map. The wind speed over a given area is, together with the air temperature, one of the major factors that greatly influence the result of the computation, so careful attention should be paid to the data preparation.


Digital model of terrain
........................
The digital terrain model is a raster layer that captures the elevation of the area of interest. It is usually created by interpolating the elevation given by contours or elevation points. It is possible to use elevation from LiDAR data or radar altimetry data, e.g. SRTM. In the case where an elevation map is not available and the terrain is essentially flat, a constant elevation layer can be used.


Vegetation cover height
........................
Information on the (effective) height of the vegetation cover is an important parameter for the calculation of aerodynamic parameters such as aerodynamic surface roughness or atmospheric boundary layer stability. The vegetation cover height can be obtained e.g. by scanning the surface using LiDAR as a digital surface model. In the case that the vegetation height layer is not available, it can be estimated by scaling the MSAVI index values between the minimum (:math:`h_{min}`) and maximum (:math:`h_{max}`) vegetation height (Gao et al. 2011):

.. math::
    :label: eq:veg_height

        h = h_{min} + \frac{MSAVI - MSAVI_{min}}{MSAVI_{min} - MSAVI_{max}} (h_{min} – h_{max})

To get an idea of the minimum and maximum vegetation cover height values to be used, either values in the table listed below or an estimated values can be used.



.. csv-table:: Summary of maximum and minimum values of vegetation cover effective height (m). Modified from Gao et al. (2011).
    :header: , Max. (m), Min. (m)

        Arid agricultural areas, 0.75, 0.01
        Forest areas, 15, 1.50
        Medium to tall graslands, 0.46, 0.31
        Low and sparse grasslands, 0.35, 0.20
        Water bodies, 0.001, –
        Urban area, 10, –
        Resident municipalities, 5, –
        Scattered municipalities 5, –
        Bareland, 0.001, –


For Central European agricultural landscapes, the maximum stand height is approx. 1-1.5 m for non forested areas. The minimum value can be set to 0.1 m. For forest areas, it is appropriate to use the value based on knowledge of the stand or forest unit.
The layers of minimum and maximum vegetation cover height are used as a raster layers.


