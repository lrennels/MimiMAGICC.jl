#The function 'interact' accounts for the overlap in CH4 and N2O in their absoprtion bands
function interact(M, N)
    d = 1.0 + 0.636 * ((M * N)/(10^6))^0.75 + 0.007 * (M/(10^3)) * ((M * N)/(10^6))^1.52
    return 0.47 * log(d)
end

@defcomp rfch4magicc begin
    N₂O_0              = Parameter() #preindustrial nitrous oxide
    STRATH2O        = Parameter() #Ratio DELQH2O/DELQCH4, for Strat H2O from CH4 Oxidation
    OZ00CH4         = Parameter() #2000 OZONE FORCING DUE TO CH4 (ORIG .168, NOW .161)
    TROZSENS        = Parameter() #TROP OZONE FORCING SENSITIVITY (TAR=0.042, EMPIRICAL=0.0389)
    OZCH4           = Parameter() #CH4 TERM FACTOR IN OZONE EQU
    CH4_0           = Parameter()
    scale_CH₄ = Parameter() #Scaling factor for uncertainty in RF (following FAIR 1.3 approach)
    CH4             = Parameter(index=[time]) #Natural CH4 emissions (Tgrams)
    QCH4OZ          = Variable(index=[time]) #Enhancement of QCH4 due to CH4-induced ozone production.
    QCH4H2O         = Variable(index=[time]) #Additional indirect CH4 forcing due to production of stratospheric H20 from CH4 oxidation
    QMeth           = Variable(index=[time]) #direct radiative forcing from atmospheric ch4 concentrations only

    function run_timestep(p, v, d, t)
        #direct radiative forcing from atmospheric ch4 concentrations only
        QCH4 = 0.036 * (sqrt(p.CH4[t]) - sqrt(p.CH4_0))
    
        v.QMeth[t] = (QCH4 - (interact(p.CH4[t], p.N₂O_0) - interact(p.CH4_0, p.N₂O_0))) * p.scale_CH₄
    
        #Calculate indirect CH4 forcing due to production of stratospheric H20 from CH4 oxidation
        v.QCH4H2O[t] = p.STRATH2O * QCH4
    
        #Calculate enhnacement of direct CH4 forcing due to ch4-induced ozone production
        #NOTE: Old concentration based equation added 2000 O₃ focring due to CH₄. This term was removed for emissions version.
        v.QCH4OZ[t] = p.TROZSENS * p.OZCH4 * log(p.CH4[t]/p.CH4_0)
    
    end

end
