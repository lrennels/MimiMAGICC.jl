# -------------------------------------------------------------------
# Straospheric Water Vapor from Methane Oxidation Radiative Forcing.
# -------------------------------------------------------------------

@defcomp rf_ch4h2o begin

    CH₄_0     = Parameter()             # Pre-industrial atmospheric methane concentration (ppb).
    STRATH2O  = Parameter()             # Share of direct methane radiative forcing that represents stratospheric water vapor forcing from methane oxidation.
    CH₄       = Parameter(index=[time]) # Atmospheric methane concentration (ppb).

    QCH4H2O   = Variable(index=[time])  # Radiative forcing for stratoshperic water vapor from methane oxidation (Wm⁻²).    rf_O₃             = Variable(index=[time])  # Total radiative forcing from tropospheric ozone (Wm⁻²).


    function run_timestep(p, v, d, t)

        if is_first(t)
            # Set initial forcing to 0.
            v.QCH4H2O[t] = 0.0
        else
            # Calculate indirect CH₄ forcing due to stratoshperic water vapor production from CH₄ oxidation.
            v.QCH4H2O[t] = p.STRATH2O * (0.036 * (sqrt(p.CH₄[t]) - sqrt(p.CH₄_0)))
        end
    end
end
