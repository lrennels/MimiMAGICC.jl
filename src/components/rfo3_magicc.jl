@defcomp rfo3magicc begin
    deltat          = Parameter()
    OZ00CH4         = Parameter() # 2000 OZONE FORCING DUE TO CH4 (ORIG .168, NOW .161)
    TROZSENS        = Parameter() #TROP OZONE FORCING SENSITIVITY (TAR=0.042, EMPIRICAL=0.0389)
    OZNOX           = Parameter() #NOX TERM FACTOR IN OZONE EQU (ORIG 0.17, NOW 0.125)
    OZCO            = Parameter() #CO TERM FACTOR IN OZONE EQU (ORIG .0014, NOW .0011)
    OZVOC           = Parameter() #VOC TERM FACTOR IN OZONE EQU (ORIG .0042, NOW .0033)
    NOX_emissions   = Parameter(index=[time]) # Global NOx emissions in Mt/yr
    CO_emissions    = Parameter(index=[time]) # Global CO emissions in Mt/yr
    NMVOC_emissions = Parameter(index=[time]) # Global non-methane VOC emissions in Mt/yr
    FOSSHIST        = Parameter(index=[time])
    QCH4OZ          = Parameter(index=[time]) #QCH4OZ IS THE TROP OZONE FORCING FROM CH4 CHANGES.
    QOZ             = Variable(index=[time]) #TROPOSPHERIC OZONE FORCING NOT ASSOCIATED WITH CH4.
    rf_O3           = Variable(index=[time]) #Total radiative forcing tropospheric ozone (from direct O3 and Ch4 oxidation)

    function run_timestep(p, v, d, t)
        QREF = 0.33 - p.OZ00CH4
    
        #Use concentrations for 1765-1999.
        if t <= 235 #index for year 1999
            FOSS0 = p.FOSSHIST[1]
            v.QOZ[t] = QREF *(p.FOSSHIST[t] - FOSS0) / (p.FOSSHIST[235] - FOSS0)
        else
    
            DQOZ = v.QOZ[235] - v.QOZ[234]
    
            #Calculate change in emissions for NOX, CO, and NMVOC
            DENOX = p.NOX_emissions[t] - p.NOX_emissions[236]
            DECO = p.CO_emissions[t] - p.CO_emissions[236]
            DEVOC = p.NMVOC_emissions[t] - p.NMVOC_emissions[236]
    
            #Non methane OZ
            v.QOZ[t] = QREF + DQOZ + p.TROZSENS*(p.OZNOX*DENOX + p.OZCO*DECO + p.OZVOC*DEVOC)
    
        end
    
        #Total trop O3 Radiative forcing (from direct O3 and indirect CH4 effect)
        v.rf_O3[t] = v.QOZ[t] + p.QCH4OZ[t]
    
     end
    
end
