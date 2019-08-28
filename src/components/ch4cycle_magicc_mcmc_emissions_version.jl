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
    fffrac          = Parameter()               #fossil fuel fraction of total CH4 emissions
    CH4_0           = Parameter() #Initial CH4 concentration.
    #CH4_2           = Parameter() #Period 2 CH4 concentration (necessary because of how code is set up for t > 2).
    temp            = Parameter(index=[time])
    CH4_emissions   = Parameter(index=[time])   # Global anthropogenic CH4 emissions in Mt/yr
    NOX_emissions   = Parameter(index=[time])   # Global NOx emissions in Mt/yr
    CO_emissions    = Parameter(index=[time])   # Global CO emissions in Mt/yr
    NMVOC_emissions = Parameter(index=[time])   # Global non-methane VOC emissions in Mt/yr
    CH4_natural     = Parameter()
    TAU_OH          = Variable(index=[time])   #Lifetime OH sink (years)
    CH4             = Variable(index=[time]) #Atmospheric concetration of CH4
    emeth           = Variable(index=[time]) #Additional CO2 emissions due to CH4 oxidation

    function run_timestep(p, v, d, t)
        # Set initial conditions before switching CH₄ emissions cycle on.
        # Note, period1 CH4 doesn't affect CH4 cycle results (period 2 does), due to subsequent equations needing values from 2 periods back.
        # Radiative forcing equations use Period 1 CH4 (so set to be the same).
        if t <= 2
            # Initialize first two periods of model (period 1 temp fixed)
            #v.CH4[1] = 0.0
            v.CH4[2] = p.CH4_0 # DOES THIS MAKE SENSE? Could also just say you initialize model for first two periods, and say period 1 temp fixed, and just say if t == 2, CH4 = CH4[0] and then set cH4 rf to period 2
            v.TAU_OH[t] = 0.0
            v.emeth[t] = 0.0
        else
            #DELT00 = 2.0 * p.temp[2] - p.temp[1]
            #TX = 2.0 * p.temp[t-1] - p.temp[t-2]
            #D_Temp = TX - DELT00
    
            # Model code had temperature sensitiity relative to temp increase above 2000, so we keep that?
            if t <= 235
                D_Temp = 0.0
            else
                DELT00 = 2.0 * p.temp[235] - p.temp[234]
                TX = 2.0 * p.temp[t-1] - p.temp[t-2]
                D_Temp = TX - DELT00
            end
    
            #Calculate change in emissions for NOX, CO, and NMVOC
            DENOX = p.NOX_emissions[t] - p.NOX_emissions[3]
            DECO = p.CO_emissions[t] - p.CO_emissions[3]
            DEVOC = p.NMVOC_emissions[t] - p.NMVOC_emissions[3]
    
            #Add natural emissions (from paper) to methane emissions
            Total_CH4emissions = p.CH4_emissions[t] + p.CH4_natural
    
            #Calcualte Tauother (based on strat and soil lifetimes)
            TAUOTHER = 1.0 / (1.0/p.TAUSOIL + 1.0/p.TAUSTRAT)
    
            #Calculate base CH4 concentration
            CM00 = v.CH4[2]
    
            #######################################################
            #Calculate OH lifetime and CH4 concentration projection
            #######################################################
    
            B = v.CH4[t-1] * p.BBCH4
    
            #Convert historical concentrations to mass
            B00 = CM00 * p.BBCH4
    
            #calculate relative oh lifetime see equation in paper A terms are coefficients given in param list. D terms are differences in emissions
            AAA = exp(p.GAM * (p.ANOX*DENOX + p.ACO*DECO + p.AVOC*DEVOC))
    
            #Change sign of ch4 parameter (does it do something else?)
            X = p.GAM * p.SCH4
    
            #Calculate Tau_oh
            U = p.TAUINIT * AAA
    
            ##########################################
            #FIRST ITERATION
            ##########################################
    
            #B defined above
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
