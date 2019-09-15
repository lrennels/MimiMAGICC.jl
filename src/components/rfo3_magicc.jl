# --------------------------------------------------------
# Tropospheric Ozone Concentration and Radiative Forcing.
# --------------------------------------------------------

@defcomp rf_o3 begin

    OZ00CH4           = Parameter()             # Year 2000 tropospheric ozone forcing due to methane (Wm⁻²).
    TROZSENS          = Parameter()             # Tropospheric ozone radiative efficiency factor.
    OZNOX             = Parameter()             # Sensitivity coefficient of tropospheric ozone to nitrogen oxides emissions.
    OZCO              = Parameter()             # Sensitivity coefficient of tropospheric ozone to carbon monoxide emissions.
    OZVOC             = Parameter()             # Sensitivity coefficient of tropospheric ozone to Non-methane volatile organic compound emissions.
    index_2000::Int64 = Parameter()             # Index to access values in the year 2000 (just for convenience).
    NOx_emissions     = Parameter(index=[time]) # Nitrogen oxides emissions (Mt yr⁻¹).
    CO_emissions      = Parameter(index=[time]) # Carbon monoxide emissions (Mt yr⁻¹).
    NMVOC_emissions   = Parameter(index=[time]) # Non-methane volatile organic compound emissions (Mt yr⁻¹).
    FOSSHIST          = Parameter(index=[time]) # Historical fossil carbon dioxide emissions values used to calculate pre-2000 non-methane ozone forcing.
    QCH4OZ            = Parameter(index=[time]) # Indirect forcing effect from methane due to methane-induced ozone production (Wm⁻²).

    QOZ               = Variable(index=[time])  # Tropospheric ozone radiative forcing not attributed to methane (Wm⁻²).
    rf_O₃             = Variable(index=[time])  # Total radiative forcing from tropospheric ozone (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Set a reference value that subtracts out reference CH₄ contribution.
        QREF = 0.33 - p.OZ00CH4

        # The CH₄ emission cycle was converted to run on emissions for all time periods (to use concentration observations during calibration).
        # The O₃ cycle follows the original version of MAGICC by switching to emissions inputs starting in 2000.
        if gettime(t) < 2000
            FOSS0 = p.FOSSHIST[1]
            v.QOZ[t] = QREF *(p.FOSSHIST[t] - FOSS0) / (p.FOSSHIST[p.index_2000-1] - FOSS0)
        else

            # Calcualte forcing difference in final periods before switching to emissions version.
            DQOZ = v.QOZ[p.index_2000-1] - v.QOZ[p.index_2000-2]

            # Calculate change in emissions for NOX, CO, and NMVOC relative to 2000.
            DENOX = p.NOx_emissions[t]   - p.NOx_emissions[p.index_2000]
            DECO  = p.CO_emissions[t]    - p.CO_emissions[p.index_2000]
            DEVOC = p.NMVOC_emissions[t] - p.NMVOC_emissions[p.index_2000]

            # Calculate tropospheric O₃ radiative forcing not attributed to CH₄.
            v.QOZ[t] = QREF + DQOZ + p.TROZSENS*(p.OZNOX*DENOX + p.OZCO*DECO + p.OZVOC*DEVOC)
        end

        # Calculate total tropospheric O₃ radiative forcing (from direct O₃ and indirect CH₄ effects).
        v.rf_O₃[t] = v.QOZ[t] + p.QCH4OZ[t]
    end
end
