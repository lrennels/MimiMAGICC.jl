@defcomp ch4cyclemagicc begin
    deltat          = Parameter()
    BBCH4           = Parameter() #Tg(CH4)/ppb unit conversion between emissions and concentrations Note:(1 Teragram = 1 Megatonne)
    ANOX            = Parameter()
    ACO             = Parameter()
    AVOC            = Parameter()
    TAUINIT         = Parameter()
    SCH4            = Parameter()
    TAUSOIL         = Parameter()
    TAUSTRAT        = Parameter()
    GAM             = Parameter()
    CH4_natural     = Parameter()
    fffrac          = Parameter() #fossil fuel fraction of total CH4 emissions
    CH4_hist        = Parameter(index=[time])
    temp            = Parameter(index=[time])
    CH4_emissions   = Parameter(index=[time]) # Global anthropogenic CH4 emissions in Mt/yr
    NOX_emissions   = Parameter(index=[time]) # Global NOx emissions in Mt/yr
    CO_emissions    = Parameter(index=[time]) # Global CO emissions in Mt/yr
    NMVOC_emissions = Parameter(index=[time]) # Global non-methane VOC emissions in Mt/yr
    TAU_OH          = Variable(index=[time]) #Lifetime OH sink (years)
    CH4             = Variable(index=[time]) #Atmospheric concetration of CH4
    emeth           = Variable(index=[time]) #Additional CO2 emissions due to CH4 oxidation

    function run_timestep(p, v, d, t)
        #Use concentrations for 1765-1999
        if t <= 235 #index for year 1999
            v.CH4[t] = p.CH4_hist[t]
            v.TAU_OH[t] = 0.0
            v.emeth[t] = 0.0
        else
            DELT00 = 2.0 * p.temp[235] - p.temp[234]
            TX = 2.0 * p.temp[t-1] - p.temp[t-2]
            D_Temp = TX - DELT00
    
            #Calculate change in emissions for NOX, CO, and NMVOC
            DENOX = p.NOX_emissions[t] - p.NOX_emissions[236]
            DECO = p.CO_emissions[t] - p.CO_emissions[236]
            DEVOC = p.NMVOC_emissions[t] - p.NMVOC_emissions[236]
    
            #Add natural emissions to methane emissions
            Total_CH4emissions = p.CH4_emissions[t] + p.CH4_natural
    
            #Calcualte Tauother (based on strat and soil lifetimes)
            TAUOTHER = 1.0 / (1.0/p.TAUSOIL + 1.0/p.TAUSTRAT)
    
            #Calculate base CH4 concentration
            CM00 = v.CH4[235]
    
            #######################################################
            #Calculate OH lifetime and CH4 concentration projection
            #######################################################
    
            B = v.CH4[t-1] * p.BBCH4
    
            #Convert historical concentrations to mass
            B00 = CM00 * p.BBCH4
    
            #calculate relative oh lifetime
            AAA = exp(p.GAM * (p.ANOX*DENOX + p.ACO*DECO + p.AVOC*DEVOC))
    
            X = p.GAM * p.SCH4
    
            #Calculate Tau_oh
            U = p.TAUINIT * AAA
    
            ##########################################
            #FIRST ITERATION
            ##########################################
    
            BBAR = B
    
            TAUBAR = U * (BBAR/B00) ^ X
    
            TAUBAR = p.TAUINIT / (p.TAUINIT/TAUBAR + 0.0316 * D_Temp)
    
            DB1 = Total_CH4emissions - (BBAR/TAUBAR) - (BBAR/TAUOTHER)
    
            B1 = B + DB1
    
    
            ##########################################
            #SECOND ITERATION
            ##########################################
            BBAR = (B + B1) / 2.0
    
            TAUBAR = U * (BBAR/B00) ^ X
    
            TAUBAR = TAUBAR * (1.0 - 0.5 * X * DB1/B)
    
            TAUBAR = p.TAUINIT / (p.TAUINIT/TAUBAR + 0.0316 * D_Temp)
    
            DB2 = Total_CH4emissions - (BBAR/TAUBAR) - (BBAR/TAUOTHER)
    
            B2 = B + DB2
    
            ##########################################
            #THIRD ITERATION
            ##########################################
            BBAR = (B + B2) / 2.0
    
            TAUBAR = U * (BBAR/B00) ^ X
    
            TAUBAR = TAUBAR * (1.0 - 0.5 * X * DB2/B)
    
            TAUBAR = p.TAUINIT / (p.TAUINIT/TAUBAR + 0.0316 * D_Temp)
    
            DB3 = Total_CH4emissions - (BBAR/TAUBAR) - (BBAR/TAUOTHER)
    
            B3 = B + DB3
    
            ##########################################
            #FOURTH ITERATION
            ##########################################
            BBAR = (B + B3) / 2.0
    
            TAUBAR = U * (BBAR/B00) ^ X
    
            TAUBAR = TAUBAR * (1.0 - 0.5 * X * DB3/B)
    
            TAUBAR = p.TAUINIT / (p.TAUINIT/TAUBAR + 0.0316 * D_Temp)
    
            DB4 = Total_CH4emissions - (BBAR/TAUBAR) - (BBAR/TAUOTHER)
    
            B4 = B + DB4
    
            v.TAU_OH[t] = TAUBAR
    
            v.CH4[t] = B4 / p.BBCH4
    
            #Calculate additional CO2 emissions due to CH4 oxidation
            v.emeth[t] = p.fffrac * 0.0020625 * (v.CH4[t] - 700.0) / p.TAUINIT
        end
    end
    
end
