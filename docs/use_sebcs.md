#### [Go back to SEBCS](index.md)

# SEBCS usage

The SEBCS for QGIS plug-in allows spatial calculation of energy balance components and other surface characteristics. The plug-in uses three different approaches to calculate the actual latent heat flux. These are the aerodynamic method according to Penman-Monteith, the SEBAL method according to [Bastiaanssen et al.](https://www.sciencedirect.com/science/article/abs/pii/S0022169498002534) and the gradient method based on the formulation of [Suleiman & Crago](https://acsess.onlinelibrary.wiley.com/doi/full/10.2134/agronj2004.3840).

## Inputs

The computation of the variables requires some data inputs. These are raster remote sensing data (Spectral data), raster data describing the surface (Additional data) and numerical data, including meteorological data and metadata. The input data varies depending on the calculation method used. An overview of the input data depending on the method used is given in the following table.


| Data                                 | PM  | SEBAL | Gradient | Optional | Comment                                |
|--------------------------------------|-----|-------|----------|---------|----------------------------------------|
| **Spectral raster data**             |
| Blue band                            | x   | x     | x        | x       |
| Green band                           | x   | x     | x        | x       |
| Red band                             | x   | x     | x        |         |
| NIR band                             | x   | x     | x        |         |
| SWIR1 band (~1.61 μm)                | x   | x     | x        | x       |
| SWIR2 band (~2.2 μm)                 | x   | x     | x        | x       |
| Surface temperature (°C)             | x   | x     | x        | x       |
| **Additional raster data**           |
| Digital model of terrain (m a.s.l)   | x   | x     | x        |
| Air temperature (°C)                 | x   | x     | x        |
| Wind speed (m/s)                     | x   | x     |
| Canopy height (m)                    | x   | x     |          |         | If known                               |
| Max. canopy height (m)               | x   | x     |          |         | If canopy height layer is not available |
| Min. canopy height (m)               | x   | x     |          |         | If canopy height layer is not available |
| Mask                                 | x   | x     | x        | x       | 
| Albedo (rel.)                        | x   | x     | x        | x       |
| **Numerical data**                   |
| Global radiation (W.m<sup>-2</sup>)  | x   | x     | x        |
| Air relative humidity (%)            | x   | x     | x        |
| Height of wind speed measurement (m) | x   | x     |          |         |
| Date of data acquisition             | x   | x     | x        |
| Time of data acqusition (GMT)        | x   | x     | x        |
| Mean latitude of data                | x   | x     | x        |         | Can be calculated directly from data   |
| Mean longiitude of data              | x   | x     | x        |         | Can be calculated directly from data   |


<span style="color:red;"><b>All raster data must have the same size, i.e. the same number of columns and rows, and the same geographic coordination system!</b></span>


## Outputs

The SEBCS for QGIS plug-in allows calculation of the following features described in the next table. The calculation of the features is described [here](calc.md).

| Feature                                         | Unit                               |
|-------------------------------------------------|------------------------------------|
| Incomming shortwave radiation                   | W.m<sup>-2</sup>                   |
| Reflected shortwave radiation                   | W.m<sup>-2</sup>                   |
| Albedo                                          | rel.                               |
| Incomming longwave radiation                    | W.m<sup>-2</sup>                   |
| Outgoing (emitted) longwave radiation           | W.m<sup>-2</sup>                   |
| Total net radiation                             | W.m<sup>-2</sup>                   |
| Surface temperature                             | °C                                 |
| Surface emissivity                              | rel.                               |
| Latent heat flux                                | W.m<sup>-2</sup>                   |
| Latent heat flux - Penman pot.                  | W.m<sup>-2</sup>                   |
| Latent heat flux - Priestley-Taylor pot.        | W.m<sup>-2</sup>                   |
| Sensible heat flux                              | W.m<sup>-2</sup>                   |
| Ground heat flux                                | W.m<sup>-2</sup>                   |
| Evapotranspiration intensity                    | mmol.m<sup>-2</sup>.s<sup>-1</sup> |
| Evaporative fraction                            | rel.                               |
| Bowen ratio                                     | unitless                           |
| Decoupling coefficient                          | rel.                               |
| Crop Water Stress Index (CWSI)                  | unitless                           |
| Friction velocity                               | m.s<sup>-1</sup>                   |
| Aerodynamic surface resistance                  | s.m<sup>-1</sup>                   |
| Surface resistance for water transfer           | s.m<sup>-1</sup>                   |
| Normalized difference vegetation index (NDVI)   | unitless                           |
| Modified soil adjusted vegetation index (MSAVI) | unitless                           |
| Normalized difference moisture index (NDMI)     | unitless                           |
| Soil adjusted vegetation index (SAVI)           | unitless                           |
| Leaf area index (LAI)                           | m<sup>2</sup>.m<sup>-2</sup>       |
| Slope                                           | °                                  |
| Aspect                                          | °                                  |


## Input data preparation

TODO