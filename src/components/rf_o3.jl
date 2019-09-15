# --------------------------------------------------------
# Tropospheric Ozone Radiative Forcing.
# --------------------------------------------------------

@defcomp rf_o3 begin

    CH₄_0             = Parameter()             # Pre-industrial atmospheric methane concentration (ppb).
    OZ00CH4           = Parameter()             # Year 2000 tropospheric ozone forcing due to methane (Wm⁻²).
    TROZSENS          = Parameter()             # Tropospheric ozone radiative efficiency factor.
    OZNOX             = Parameter()             # Sensitivity coefficient of tropospheric ozone to nitrogen oxides emissions.
    OZCO              = Parameter()             # Sensitivity coefficient of tropospheric ozone to carbon monoxide emissions.
    OZVOC             = Parameter()             # Sensitivity coefficient of tropospheric ozone to Non-methane volatile organic compound emissions.
    TROZSENS          = Parameter()             # Tropospheric ozone radiative efficiency factor.
    OZCH4             = Parameter()             # Sensitivity coefficient of tropospheric ozone to methane concentration.
    index_2000::Int64 = Parameter()             # Index to access values in the year 2000 (just for convenience).
    CH₄               = Parameter(index=[time]) # Atmospheric methane concentration (ppb).
    NOx_emissions     = Parameter(index=[time]) # Nitrogen oxides emissions (Mt yr⁻¹).
    CO_emissions      = Parameter(index=[time]) # Carbon monoxide emissions (Mt yr⁻¹).
    NMVOC_emissions   = Parameter(index=[time]) # Non-methane volatile organic compound emissions (Mt yr⁻¹).
    FOSSHIST          = Parameter(index=[time]) # Historical fossil carbon dioxide emissions values used to calculate pre-2000 non-methane ozone forcing.

    QOZ               = Variable(index=[time])  # Tropospheric ozone radiative forcing not attributed to methane (Wm⁻²).
    QCH4OZ            = Variable(index=[time])  # Indirect forcing effect from methane due to methane-induced ozone production (Wm⁻²).
    rf_O₃             = Variable(index=[time])  # Total radiative forcing from tropospheric ozone (Wm⁻²).


    function run_timestep(p, v, d, t)

        if is_first(t)
            # Set initial values for O₃ forcing to 0.0.
            v.QCH4OZ[t] = 0.0
            v.QOZ[t] = 0.0
            v.rf_O₃[t] = 0.0

        else

            # The CH₄ emission cycle was converted to run on emissions for all time periods (to be able to use CH₄ concentration observations during calibration).
            # The non-methane component of the O₃ cycle follows the original version of MAGICC by switching on the primary emissions equations after year 2000.
            if gettime(t) < 2000
                # Calculate historic non-CH₄ tropospheric O₃ forcing.
                v.QOZ[t] = (0.33 - p.OZ00CH4) * (p.FOSSHIST[t] - p.FOSSHIST[1]) / (p.FOSSHIST[p.index_2000-1] - p.FOSSHIST[1])

            else
                # Calculate change in emissions for NOX, CO, and NMVOC relative to 2000.
                DENOX = p.NOx_emissions[t]   - p.NOx_emissions[p.index_2000]
                DECO  = p.CO_emissions[t]    - p.CO_emissions[p.index_2000]
                DEVOC = p.NMVOC_emissions[t] - p.NMVOC_emissions[p.index_2000]

                # Calculate tropospheric O₃ radiative forcing not attributed to CH₄.
                v.QOZ[t] = (0.33 - p.OZ00CH4) + (v.QOZ[p.index_2000-1] - v.QOZ[p.index_2000-2]) + p.TROZSENS*(p.OZNOX*DENOX + p.OZCO*DECO + p.OZVOC*DEVOC)
            end

            # Calculate indirect CH₄ forcing due to CH₄-induced O₃ production (based on pure emissions version of CH₄ cycle).
            v.QCH4OZ[t] = p.TROZSENS * p.OZCH4 * log(p.CH₄[t]/p.CH₄_0)

            # Calculate total tropospheric O₃ radiative forcing (from direct O₃ and indirect CH₄ effects).
            v.rf_O₃[t] = v.QOZ[t] + v.QCH4OZ[t]
        end
    end
end
