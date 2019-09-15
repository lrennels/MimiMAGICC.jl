# --------------------------------------------------
# Radiative Forcing from Methane.
# --------------------------------------------------

# Create a function 'interact' to account for the overlap in methane and nitrous oxide in their radiative absoprtion bands.
function interact(M, N)
    d = 1.0 + 0.636 * ((M * N)/(10^6))^0.75 + 0.007 * (M/(10^3)) * ((M * N)/(10^6))^1.52
    return 0.47 * log(d)
end

@defcomp rf_ch4 begin
    CH₄_0     = Parameter()             # Pre-industrial atmospheric methane concentration (ppb).
    N₂O_0     = Parameter()             # Pre-industrial atmospheric nitrous oxide concentration (ppb).
    STRATH2O  = Parameter()             # Share of direct methane radiative forcing that represents stratospheric water vapor forcing from methane oxidation.
    TROZSENS  = Parameter()             # Tropospheric ozone radiative efficiency factor.
    OZCH4     = Parameter()             # Sensitivity coefficient of tropospheric ozone to methane concentration.
    scale_CH₄ = Parameter()             # Scaling factor for uncertainty in methane radiative forcing.
    CH₄       = Parameter(index=[time]) # Atmospheric methane concentration (ppb).

    QCH4OZ    = Variable(index=[time])  # Indirect forcing effect from methane due to methane-induced ozone production (Wm⁻²).
    QCH4H2O   = Variable(index=[time])  # Radiative forcing for stratoshperic water vapor from methane oxidation (Wm⁻²).
    QMeth     = Variable(index=[time])  # Direct forcing from atmospheric methane concentrations (Wm⁻²).


    function run_timestep(p, v, d, t)

        if is_first(t)
            # Set initial forcings to 0.
            v.QMeth[t]   = 0.0
            v.QCH4H2O[t] = 0.0
            v.QCH4OZ[t]  = 0.0
        else
            # Direct radiative forcing from atmospheric CH₄ concentrations.
            QCH4 = 0.036 * (sqrt(p.CH₄[t]) - sqrt(p.CH₄_0))
            v.QMeth[t] = (QCH4 - (interact(p.CH₄[t], p.N₂O_0) - interact(p.CH₄_0, p.N₂O_0))) * p.scale_CH₄

            # Calculate indirect CH₄ forcing due to stratoshperic water vapor production from CH₄ oxidation.
            v.QCH4H2O[t] = p.STRATH2O * QCH4

            # Calculate indirect CH₄ forcing due to CH₄-induced O₃ production.
            v.QCH4OZ[t] = p.TROZSENS * p.OZCH4 * log(p.CH₄[t]/p.CH₄_0)
        end
    end
end
