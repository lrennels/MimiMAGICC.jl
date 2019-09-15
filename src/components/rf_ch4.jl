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
    scale_CH₄ = Parameter()             # Scaling factor for uncertainty in methane radiative forcing.
    CH₄       = Parameter(index=[time]) # Atmospheric methane concentration (ppb).

    QMeth     = Variable(index=[time])  # Direct forcing from atmospheric methane concentrations (Wm⁻²).


    function run_timestep(p, v, d, t)

        if is_first(t)
            # Set initial forcings to 0.
            v.QMeth[t]   = 0.0
        else
            # Direct radiative forcing from atmospheric CH₄ concentrations.
            v.QMeth[t] = (0.036 * (sqrt(p.CH₄[t]) - sqrt(p.CH₄_0)) - (interact(p.CH₄[t], p.N₂O_0) - interact(p.CH₄_0, p.N₂O_0))) * p.scale_CH₄
        end
    end
end
