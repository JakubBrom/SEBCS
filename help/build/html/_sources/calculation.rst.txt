Calculation
===========

Using the SEBCS for QGIS plug-in, it is possible to calculate several characteristics related to energy exchange at the surface and energy transformation into individual heat fluxes. It also calculates meteorological features and vegetation indices that have some connection with the land microclimate. The following features can be calculated:

* meteorological features
* vegetation cover characteristics
* radiation balance features
* surface aerodynamic features and boundary layer (atmospheric) stability
* heat balance features


Radiation balance features
---------------------------

The calculation of the radiation balance features includes the balance of short-wave and long-wave radiation on the active surface and the resulting total net radiation. The calculation procedure includes the calculation of the albedo, surface reflectance and atmospheric emissivity. The SEBCS for QGIS includes surface geometry issues in the radiation balance features calculation.
The relationship of the individual radiation balance components is summarized by the equation :ref:`(1) <rnflux>`:

.. _rnflux:
.. math::
    :label: eq:Rn_fluxes

    Rn=Rs_{\downarrow}-Rs_{\uparrow}+Rl_{\downarrow}-Rl_{\uparrow}

The individual radiation balance components calculation procedure is described in the following text.


Incoming short-wave radiation flux
..................................

The incoming (incident) short-wave radiation flux (irradiance) is the main source of energy for evaporation and energy transformations on the Earth's surface. The flux is usually measured at meteorological stations as the energy of radiation incident on a horizontal surface (global radiation), but the amount of energy reaching the surface is significantly related to its shape and the geometry of the incident radiation. The SEBCS for QGIS module calculates the flux of incident shortwave radiation as a function of surface shape and radiation geometry. The computational approach summarized by `Kumar (1997) <https://www.tandfonline.com/doi/abs/10.1080/136588197242266>`_ was used as follows.
The value of the incident short-wave radiation flux on the surface was calculated using equation:

.. math::
    :label: eq:RS_in

    Rs_{\downarrow} = I_s \cdot \cos i

where

.. math::
    :label: eq:cosi

    \cos i = \sin \delta_s (\sin Lat \cos \beta_s - \cos Lat \sin \beta_s \cos a_{w}) \\
        + \cos \delta_s \cos H_s (\cos Lat \cos \beta_s
        + \sin Lat \sin \beta_s \cos a_{w}) \\
        + \cos \delta_s \sin \beta_s \sin a_{w} \sin H_s

The value of :math:`\delta_s` can be calculated using equation:

.. math::
    :label: eq:sol_dek

    \delta_s = 23.45 \sin \left(\frac{360(284+N)}{365}\right)

The hour angle corresponds to the condition where :math:`H_s = 0` at noon (when the Sun is on the meridian). Each hour plus represents a change in angle of +15°; each hour minus represents a change in angle of -15°. By convention morning values are positive and afternoon values are negative.

.. math::
    :label: eq:hour_ang

    H_s = (12 - S_t)\cdot 15

Depending on the longitude, the solar time (:math:`S_t`) was determined according to the relation:

.. math::
    :label: eq:sol_time

    S_t = GMT + \frac{24 \cdot Long}{360}


Since the module uses the global radiation measured per horizontal surface as input, the value of :math:`I_s` is calculated based on the relationship:


.. math::
    :label: eq:irrad

    I_s = \frac{R_{s\downarrow konst}}{\sin \alpha_s}

where

.. math::
    :label: eq:sol_ang

    \sin \alpha_s = \sin \delta_s \sin Lat + \cos \delta_s \cos Lat \cos H_s


Surface reflectance (albedo)
............................

The calculation of the broadband reflectance or albedo is based on empirical approache. An empirical relation was used for the calculation, which calculates the broadband albedo based on the spectral reflectance of the surface for individual spectral bands according to the relation (Tasumi et al. 2008):

.. math::
    :label: eq:albedo

    \alpha = \displaystyle\sum_{b=1}^{7} (\rho_{s\_b} \cdot w_b)

The :math:`w_b` constants for each Landsat satellite are given in the table:

.. table:: Values of :math:`w_b` for particular spectral bands according to Liang et al. (2001 and 2003) and Tasumi et al. (2008) for Landsat 4, 5 and 7 and according to Olmeo et al (2017) for Landsat 8 and 9.

    +----------------+-------+-------+-------+-------+-------+-------+
    |Spectral band   | Blue  | Green | Red   | NIR   | SWIR1 | SWIR2 |
    +================+=======+=======+=======+=======+=======+=======+
    |Landsat 8, 9    | 0.246 | 0.146 | 0.191 | 0.304 | 0.105 | 0.008 |
    +----------------+-------+-------+-------+-------+-------+-------+
    |Landsat 4, 5, 7 | 0.254 | 0.149 | 0.147 | 0.311 | 0.103 | 0.036 |
    +----------------+-------+-------+-------+-------+-------+-------+


Alternatively, for other data sources, an empirical approach based on the use of vegetation indices NDVI and MSAVI can be used (Duffková et al. 2012):

.. math::
    :label: eq:albedo_Brom

    \alpha =  0.08611 + 0.89472 \cdot MSAVI + 5.55866 \cdot  MSAVI^2 -0.1183 \cdot NDVI\\
        - 1.9818 \cdot MSAVI^3 - 4.5034 \cdot MSAVI \cdot NDVI - 11.463 \cdot MSAVI^2 \cdot NDVI\\
        + 7.46145 \cdot MSAVI \cdot NDVI^2 + 5.2994 \cdot MSAVI^2 \cdot NDVI^2\\
        + 4.76657 \cdot MSAVI^3 \cdot NDVI - 2.3127 \cdot MSAVI^3 \cdot NDVI^2\\
        - 3.4274 \cdot MSAVI \cdot NDVI^3

The reflected shortwave radiation flux is calculated as follows:

.. math::
    :label: eq:Rs_out

    Rs_\uparrow = \alpha \cdot Rs_\downarrow


Incoming long-wave radiation flux
..................................

The incoming longwave radiation emitted by the atmosphere is calculated from the air temperature measured at :math:`z`. The calculation is based on the Stefan-Boltzmann law:

.. math::
    :label: eq:Rl_down

    Rl_{\downarrow} = \varepsilon_{ac} \sigma {T_{a\_K}}^4

The emissivity of the atmosphere is calculated based on the air temperature at the :math:`z` level and the amount of water vapour in the air. According to Brutsaert (1982) it is calculated:

.. math::
    :label: eq:emis_atm

    \varepsilon_{ac} = 1.24 \left( \frac{e_a \cdot 10}{T_a + 273.16} \right)^{\frac{1}{7}}

This approach is largely approximate and can be used for clear-sky weather conditions using the calculated :math:`\varepsilon_{ac}` value. A more appropriate approach is to measure the value of :math:`Rl_{\uparrow}` directly.


Emitted long-wave radiation flux
................................

The longwave radiation emitted by a surface is determined by the temperature of the surface. The calculation is based on the Stefan-Boltzmann law:

.. math::
    :label: eq:Rl_up

    Rl_{\uparrow} = \varepsilon \sigma {T_{s\_K}}^4

The surface emissivity is calculated based on an empirical approach of emissivity determination using the NDVI Treshold Method - :math:`NDVI^{THM}` (Sobrino et al. 2004). The method uses the NDVI index (Normalized Difference Vegetation Index, Tucker 1979). For the emissivity determination, the range of index values is divided into three categories:

* :math:`NDVI < 0.2` - In this case the surface is considered as bare ground and the emissivity values are derived from the reflectance values in the red region of the spectrum
* :math:`NDVI > 0.5` - In this case the surface is fully covered by vegetation and a typical emissivity value of :math:`\varepsilon = 0.99` is determined.
* :math:`0.2 ≤ NDVI ≤ 0.5` - In this case the surface can be considered as a mixture of bare soil and vegetation cover.

The relationship between emissivity and surface cover is shown by the following relationship:

.. math::
    :label: eq:emis_all

    \varepsilon = \varepsilon_v P_v + \varepsilon_s (1 - P_v) +d\varepsilon

where fraction of vegetation cover is calculated as follows:

.. math::
    :label: eq:frac

    P_v = \left(\frac{NDVI - NDVI_{min}}{NDVI_{max} - NDVI_{min}}\right)^2

where :math:`NDVI_{min} = 0.2` and :math:`NDVI_{max} = 0.5`.
As a result, the emissivity can be calculated based on the empirical relation as follows (Sobrino et al. 2004):

.. math::
    :label: eq:emis_surf

    \varepsilon = 0.004 P_v + 0.986

.. TODO
    Surface aerodynamic parameters and atmospheric stability
    ---------------------------------------------------------
    Heat balance and energy fluxes
    --------------------------------

Vegetation cover characteristics
---------------------------------

The calculation of vegetation cover characteristics includes the estimation of vegetation spectral indices and leaf area index.

**Normalized Difference Vegetation Index (NDVI)** is one of the most widely used spectral vegetation indices. NDVI provides information on vegetation cover, its quality and possibly also its quantity (biomass). NDVI can be calculated as follows:

.. math::
    :label: eq:ndvi

    NDVI=\frac{R_{NIR}-R_{RED}}{R_{NIR}+R_{RED}}


**Modified Soil Adjusted Vegetation Index (MSAVI)** has similar uses to the NDVI. Unlike the NDVI, there is no oversaturation of values with higher vegetation cover. The MSAVI can also be used to estimate vegetation height. MSAVI can be calculated as follows:

.. math::
    :label: eq:msavi

    MSAVI=0.5\cdot(2R_{NIR}+1-\sqrt{(2R_{NIR}+1)^{2}-8\cdot(R_{NIR}-R_{RED})}


**Normalized Difference Moisture Index (NDMI)** can be used for surface moisture estimation. NDMI can be calculated as follows:



.. math::
    :label: eq:ndmi

    NDVI=\frac{R_{NIR}-R_{SWIR1}}{R_{NIR}+R_{SWIR1}}


**Soil Adjusted Vegetation Index (SAVI)** is next widely used vegetation index. It is similar to MSAVI. SAVI can be calculated as follows:

.. math::
    :label: eq:savi

    SAVI = \frac{NIR - R}{NIR + R + L} (1 + L)

where :math:`L` is soil brightness correction factor defined as 0.5 to accommodate most land cover types.
SEBCS for QGIS calculates leaf area index (:math:`LAI`) according to `Jafaar & Ahmad (2016) <https://www.sciencedirect.com/science/article/pii/S0034425718305947?via%3Dihub>`_:

.. math::
    LAI_1=
        \begin{cases}
            11 \cdot SAVI^3         & SAVI > 0;\  SAVI \leq 0.817\\
            6                       & SAVI > 0.817
        \end{cases}

.. math::
    LAI_2 =
      \begin{cases}
        -\frac{\ln \frac{0.61-SAVI}{0.51}}{0.91}
            & SAVI >0;\  SAVI\leq 0.61\\
        6                       & SAVI > 0.61
      \end{cases}

.. math::
    LAI = \frac{LAI_1 + LAI_2}{2}

The SEBCS_lib library contains some additional methods for calculation. See :doc:`documentation <SEBCSlib>`.


.. TODO
    Meteorological features
    ------------------------

    .. math::
        :label: eq:Rn_bil

        Rn=G+H+LE



    .. math::
        :label: eq:Ta

        T_a = T_{st} + \Gamma(Z_{st}-DMT)

    .. math::
        :label: eq:U

        U = U_{st}\frac{\ln\left(\frac{z}{z_{0m\_st}}\right)}{\ln\left(\frac{z_{st}}{z_{0m\_st}}\right)}

    .. math::
        :label: eq:z0m_st

        z_{0m\_st} = 0.123 h_{st}

    .. math::
        :label: eq:press

        P = 101.3 \left( \frac{293-\Gamma \cdot (DMT + z)}{293} \right)^{5.26}

    .. math::
        :label: eq:Ea_sat

        E_a = 0.61121 \cdot \exp \left( \frac{17.502 \cdot T_a}{240.97 + T_a} \right)

    pro účely výpočtu albeda je hodnota vypočtena pro teplotu vzduchu ve výšce zst.

    .. math::
        :label: eq:ea

        e_a = \frac{E_a \cdot Rh}{100}

    pro účely výpočtu albeda je hodnota vypočtena pro teplotu vzduchu ve výšce zst.

    .. math::
        :label: eq:rho

        \rho = \frac{353.4}{T_a + 273}

    .. math::
        :label: eq:latent

        \lambda = 2501 - 2.3723 \cdot T_a

    .. math::
        :label: eq:gama

        \gamma = \frac{c_p \cdot P}{\lambda \cdot 0.622}




    .. math::
        :label: eq:tsk

        T_s = T_{s\_K} - 273.16


    .. math::
        :label: eq:delta

        \Delta = 45.03 + 3.014 T + 0.05345 T^2 + 0.00224 T^3

    where

    .. math::
        :label: eq:t_mean

        T = \frac{T_a + T_s}{2}


    .. math::
        :label: eq:Es

        E_s = 0.61121 \cdot \exp{\left(\frac{17.502 \cdot T_s}{240.97 + T_s}\right)}



    .. math::
        :label: eq:Rn

        Rn = Rs_\downarrow - Rs_\uparrow + Rl_{\downarrow} - Rl_{\uparrow}

    .. math::
        :label: eq:vegheight

        h = h_{min} + \frac{MSAVI - MSAVI_{min}}{MSAVI_{min} - MSAVI_{max}} (h_{min} - h_{max})

    .. math::
        :label: eq:zeroplane

        d = \frac{2}{3}h

    .. math::
        :label: eq:zom

        z_{0m} = 0.123 \cdot h

    .. math::
        :label: eq:zoh

        z_{0h} = 0.1 \cdot z_{0m}

    .. math::
        :label: eq:MO_length

        L = \frac{{u^*}^3 \rho c_p T_{a\_K}}{\kappa g H} = \frac{{u^*}^2 T_{a\_K}}{\kappa g T^*}

    .. math::
        :label: eq:MOS

        \varsigma = \frac{z}{L}

    .. math::
        :label: eq:psim_instab

        \Psi_m(\varsigma) = 2 \ln \left( \frac{1 + x}{2} \right) + \ln \left( \frac{1+ x^2}{2}\right) - 2 \arctan (x) + \frac{\pi}{2}

    .. math::
        :label: eq:psih_instab

        \Psi_h(\varsigma) = 2 \ln \left( \frac{1 + x}{2} \right)

    .. math::
        :label: eq:x

        x = (1-16\varsigma)^{0.25}

    .. math::
        :label: eq:psim_stab

        \Psi_m(\varsigma) = -\left[a\varsigma + b \left( \varsigma - \frac{c}{d} \right) \exp(-d\varsigma) + \frac{bc}{d} \right]

    .. math::
        :label: eq:psih_stab

        \Psi_h(\varsigma) = - \left[ \left(1 + \frac{2a}{3}\varsigma \right) + b \left(\varsigma - \frac{c}{d} \right) \exp(-d\varsigma) + \left(\frac{bc}{d} -1 \right) \right]

    .. math::
        :label: eq:frict

        u^* = \frac{\kappa U}{\ln \left(\frac{z-d}{z_{0m}} \right) \Psi_m(\varsigma)}

    .. math::
        :label: eq:virtT

        T^* = \frac{\kappa(T_a - T_s)}{\ln \left(\frac{z-d}{z_{0h}} \right) \Psi_h(\varsigma)}

    .. math::
        :label: eq:ra_Thom

        r_a = \frac{\left[ \ln \left(\frac{z-d}{z_{0m}} \right) \Psi_m(\varsigma) \right]\left[ \ln \left(\frac{z-d}{z_{0h}} \right) \Psi_h(\varsigma) \right]}{U \kappa^2}

    .. math::
        :label: eq:ra_SEBAL

        r_{ah} = \frac{\ln \left( \frac{z_2}{z_1}\right) - \Psi_{h\_z_2}(\varsigma) + \Psi_{h\_z_1}(\varsigma)}{u^*\kappa}

    .. math::
        :label: eq:G

        G=\frac{T_{s}}{\alpha}\left(0.0038\alpha+0.0074\alpha^{2}\right)\left(1-0.98NDVI^{4}\right)Rn

    .. math::
        :label: eq:H

        H=\frac{\rho c_{p}\delta T}{r_{a}}


    .. math::
        :label: eq:LE_bil

        LE=Rn-H-G

    .. math::
        :label: eq:cwsi

        CWSI=1-\frac{\Delta+\gamma^{*}}{\Delta+\gamma\left(1+\frac{r_{c}}{r_{a}}\right)}

    .. math::
        :label: eq:gamma_stair

        \gamma^{*}=\gamma+\left(1+\frac{r_{cp}}{r_{a}}\right)

    .. math::
        :label: eq:r_c

        r_{c}=\left[\left(\frac{\Delta+\gamma}{\Omega}-\Delta\right)\frac{1}{\gamma}-1\right]r_{a}

    .. math::
        :label: eq:r_cp

        r_{cp}=\frac{\left(E_{s}-e_{a}\right)\rho c_{p}}{\gamma\cdot LE_{p}}-r_{a}

    .. math::
        :label: eq:EF

        EF=\frac{LE}{Rn-G}