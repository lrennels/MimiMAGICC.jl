@defcomp radiativeforcing begin
    deltat = Parameter()
    rf_co2 = Parameter(index=[time])
    rf_aerosol = Parameter(index=[time])
    QMeth = Parameter(index=[time]) #direct radiative forcing from atmospheric ch4 concentrations only
    rf_O3 = Parameter(index=[time]) #radiative forcing from tropospheric ozone (direct O3 and ch4-induced ozone production)
    QCH4H2O = Parameter(index=[time]) #indirect CH4 forcing due to production of stratospheric H20 from CH4 oxidation
    rf_other = Parameter(index=[time])
    alpha = Parameter()
    rf = Variable(index=[time])

    function run_timestep(p, v, d, t)
        if is_first(t)
            v.rf[t] = 0.0
        else    
            v.rf[t] = p.rf_co2[t] + p.QMeth[t]+ p.QCH4H2O[t] + p.rf_other[t] + (p.alpha * p.rf_aerosol[t]) + p.rf_O3[t]
        end    
    end

end
