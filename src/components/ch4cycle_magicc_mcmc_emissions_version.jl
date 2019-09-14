# --------------------------------------------------
# Methane Cycle (with OH calculations).
# --------------------------------------------------

@defcomp ch4_cycle begin

    BBCH4             = Parameter()             # Tg to ppb unit conversion between emissions and concentrations Note:(1 Teragram = 1 Megatonne).
    ANOX              = Parameter()             # Hydroxyl lifetime coefficient for nitrogen oxides.
    ACO               = Parameter()             # Hydroxyl lifetime coefficient for carbon monoxide.
    AVOC              = Parameter()             # Hydroxyl lifetime coefficient for non-methane volatile organic compounds.
    SCH4              = Parameter()             # Hydroxyl lifetime coefficient for methane.
    TAUINIT           = Parameter()             # Initial methane tropospheric lifetime due to hydroxyl abundance (years).
    TAUSOIL           = Parameter()             # Methane soil sink lifetime (years).
    TAUSTRAT          = Parameter()             # Methane stratospheric lifetime (years).
    GAM               = Parameter()             # Hydroxyl lifetime parameter.
    CH4_natural       = Parameter()             # Natural methane emissions (Mt yr⁻¹).
    fffrac            = Parameter()             # Fossil fuel fraction of total methane emissions.
    CH4_0             = Parameter()             # Atmospheric methane pre-industrial concentration (ppb).
    index_2000::Int64 = Parameter()             # Index to access values in the year 2000 (just for convenience).
    temperature       = Parameter(index=[time]) # Global surface temperature anomaly (°C).
    CH4_emissions     = Parameter(index=[time]) # Global anthropogenic methane emissions (Mt yr⁻¹).
    NOX_emissions     = Parameter(index=[time]) # Global nitrogen oxides emissions (Mt yr⁻¹).
    CO_emissions      = Parameter(index=[time]) # Global carbon monoxide emissions (Mt yr⁻¹).
    NMVOC_emissions   = Parameter(index=[time]) # Global non-methane volatile organic compound emissions (Mt yr⁻¹).

    TAU_OH            = Variable(index=[time])  # Methane tropospheric lifetime due to hydroxyl abundance (years).
    CH4               = Variable(index=[time])  # Atmospheric methane concetration for current period (ppb).
    emeth             = Variable(index=[time])  # Methane that has been oxidized to carbon dioxide.


    function run_timestep(p, v, d, t)

        # Set initial conditions. NOTE: The CH₄ cycle equations require temperature values from the previous two timesteps (t-1 and t-2), so the cycle switches
            # on starting in timestep 3. Further, timestep 1 CH₄ concentrations do not affect timestep 1 temperature or the results below (but timestep 2 CH₄
            # concentrations do). We therefore treat period 2 CH₄ concentrations as an initial (and uncertain) condition.
        if t.t <= 2
            v.CH4[2] = p.CH4_0
            v.emeth[t] = 0.0
        else

            # Calculate changes for temperature sensitivity equation (maintain MAGICC's original scaling relative to increase above 2000).
            if gettime(t) < 2000
                D_Temp = 0.0
            else
                DELT00 = 2.0 * p.temperature[p.index_2000-1] - p.temperature[p.index_2000-2]
                TX = 2.0 * p.temperature[t-1] - p.temperature[t-2]
                D_Temp = TX - DELT00
            end

            # Calculate change in emissions for NOx, CO, and NMVOC (relative to period 3 when CH₄ cycle model engages due to t-1 and t-2 terms above).
            DENOX = p.NOX_emissions[t] - p.NOX_emissions[3]
            DECO  = p.CO_emissions[t] - p.CO_emissions[3]
            DEVOC = p.NMVOC_emissions[t] - p.NMVOC_emissions[3]

            # Add natural CH₄ emissions to anthropogenic CH₄ emissions.
            Total_CH4emissions = p.CH4_emissions[t] + p.CH4_natural

            # Calcualte TAUOTHER term (based on stratospheric and soil sink lifetimes).
            TAUOTHER = 1.0 / (1.0/p.TAUSOIL + 1.0/p.TAUSTRAT)

            ##########################################################
            # Iterate to Calculate OH lifetime and CH₄ concentration.
            # See equations A30-A32 in Meinsshausen et al.(2011).
            ##########################################################

            # Convert previous timestep concentration to mass.
            B = v.CH4[t-1] * p.BBCH4

            # Convert historical concentration to mass.
            B00 = p.CH4_0 * p.BBCH4

            # Account for other gases on OH abundance and tropospheric methane lifetime.
            AAA = exp(p.GAM * (p.ANOX*DENOX + p.ACO*DECO + p.AVOC*DEVOC))

            # Multiply two parameters together.
            X = p.GAM * p.SCH4

            # Multiply initial tropospheric lifetime by AAA.
            U = p.TAUINIT * AAA

            # Iterate to calculate a final CH₄ tropospheric lifetime (accounting for a temperature feedback).

            #----------------------------------------
            # FIRST ITERATION
            #----------------------------------------
            BBAR = B
            TAUBAR = U * (BBAR/B00) ^ X
            TAUBAR = p.TAUINIT / (p.TAUINIT/TAUBAR + 0.0316 * D_Temp)
            DB1 = Total_CH4emissions - (BBAR/TAUBAR) - (BBAR/TAUOTHER)
            B1 = B + DB1

            #----------------------------------------
            # SECOND ITERATION
            #----------------------------------------
            BBAR = (B + B1) / 2.0
            TAUBAR = U * (BBAR/B00) ^ X
            TAUBAR = TAUBAR * (1.0 - 0.5 * X * DB1/B)
            TAUBAR = p.TAUINIT / (p.TAUINIT/TAUBAR + 0.0316 * D_Temp)
            DB2 = Total_CH4emissions - (BBAR/TAUBAR) - (BBAR/TAUOTHER)
            B2 = B + DB2

            #----------------------------------------
            # THIRD ITERATION
            #----------------------------------------
            BBAR = (B + B2) / 2.0
            TAUBAR = U * (BBAR/B00) ^ X
            TAUBAR = TAUBAR * (1.0 - 0.5 * X * DB2/B)
            TAUBAR = p.TAUINIT / (p.TAUINIT/TAUBAR + 0.0316 * D_Temp)
            DB3 = Total_CH4emissions - (BBAR/TAUBAR) - (BBAR/TAUOTHER)
            B3 = B + DB3

            #----------------------------------------
            # FOURTH ITERATION
            #----------------------------------------
            BBAR = (B + B3) / 2.0
            TAUBAR = U * (BBAR/B00) ^ X
            TAUBAR = TAUBAR * (1.0 - 0.5 * X * DB3/B)
            TAUBAR = p.TAUINIT / (p.TAUINIT/TAUBAR + 0.0316 * D_Temp)
            DB4 = Total_CH4emissions - (BBAR/TAUBAR) - (BBAR/TAUOTHER)
            B4 = B + DB4

            # Calculate CH₄ tropospheric lifetime.
            v.TAU_OH[t] = TAUBAR

            # Calculate CH₄ concentration.
            v.CH4[t] = B4 / p.BBCH4

            #Calculate additional CO₂ emissions due to CH₄ oxidation.
            v.emeth[t] = p.fffrac * 0.0020625 * (v.CH4[t] - 700.0) / p.TAUINIT
        end
    end
end
