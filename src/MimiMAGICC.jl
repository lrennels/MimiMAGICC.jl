module MimiMAGICC

using Mimi
using DataFrames
using CSVFiles

include("components/ch4_cycle.jl")
include("components/rf_ch4.jl")
include("components/rf_ch4h2o.jl")
include("components/rf_o3.jl")


export get_magicc_ch4

function get_magicc_ch4(;rcp_scenario::String="RCP85", start_year::Int64=1765, end_year::Int64=2300)

    # ---------------------------------------------
    # Load and clean up necessary data.
    # ---------------------------------------------

    # Load RCP emissions and concentration scenario values (RCP options = "RCP26", "RCP45", "RCP60", and "RCP85").
    rcp_emissions      = DataFrame(load(joinpath(@__DIR__, "..", "data", rcp_scenario*"_EMISSIONS.csv"), skiplines_begin=36))
    rcp_concentrations = DataFrame(load(joinpath(@__DIR__, "..", "data", rcp_scenario*"_CONCENTRATIONS.csv"), skiplines_begin=37))

    # Load temperature scenario (mean response from SNEASY across 100,000 posterior parameter samples).
    rcp_temperature    = DataFrame(load(joinpath(@__DIR__, "..", "data", rcp_scenario*"_temperature_sneasy_1765_2300.csv")))

    # Load fossil CO₂ emissions scenario to calculate historic tropospheric O₃ forcing (this scenario is specific to how MAGICC calculates this forcing).
    foss_hist_for_O₃   = DataFrame(load(joinpath(@__DIR__, "..", "data", "FOSS_HIST_magicc.csv")))

    # Find indices for start and end years to crop scenarios to correct time frame.
    start_index, end_index = findall((in)([start_year, end_year]), rcp_emissions.YEARS)
    rcp_emissions    = rcp_emissions[start_index:end_index, :]
    rcp_temperature  = rcp_temperature[start_index:end_index, :]
    foss_hist_for_O₃ = foss_hist_for_O₃[start_index:end_index, :]

    # Set initial CH₄ and N₂O concentrations to RCP 1765 values.
    CH₄_0 = rcp_concentrations[1, :CH4]
    N₂O_0 = rcp_concentrations[1, :N2O]

    # Get index for year 2000 given model time horizon.
    index_2000 = findall(x -> x == 2000, start_year:end_year)[1]

    # ---------------------------------------------
    # Initialize Mimi model.
    # ---------------------------------------------

    # Create a Mimi model.
    m = Model()

    # Set time index.
    set_dimension!(m, :time, start_year:end_year)

    # ---------------------------------------------
    # Add components to model.
    # ---------------------------------------------
    add_comp!(m, ch4_cycle)
    add_comp!(m, rf_ch4)
    add_comp!(m, rf_ch4h2o)
    add_comp!(m, rf_o3)

    # ---------------------------------------------
    # Set component parameters.
    # ---------------------------------------------

    # ---- Methane Cycle ---- #
    set_param!(m, :ch4_cycle, :index_2000, index_2000)
    set_param!(m, :ch4_cycle, :BBCH4, 2.78)
    set_param!(m, :ch4_cycle, :ANOX, 0.0042)
    set_param!(m, :ch4_cycle, :ACO, -0.000105)
    set_param!(m, :ch4_cycle, :AVOC, -0.000315)
    set_param!(m, :ch4_cycle, :TAUINIT, 7.73)
    set_param!(m, :ch4_cycle, :SCH4, -0.32)
    set_param!(m, :ch4_cycle, :TAUSOIL, 160.)
    set_param!(m, :ch4_cycle, :TAUSTRAT, 120.)
    set_param!(m, :ch4_cycle, :GAM, -1.0)
    set_param!(m, :ch4_cycle, :CH4_natural, 266.5)
    set_param!(m, :ch4_cycle, :fffrac, 0.18)
    set_param!(m, :ch4_cycle, :CH₄_0, CH₄_0)
    set_param!(m, :ch4_cycle, :temperature, rcp_temperature.Temperature)
    set_param!(m, :ch4_cycle, :CH4_emissions, rcp_emissions.CH4)
    set_param!(m, :ch4_cycle, :NOX_emissions, rcp_emissions.NOx)
    set_param!(m, :ch4_cycle, :CO_emissions, rcp_emissions.CO)
    set_param!(m, :ch4_cycle, :NMVOC_emissions, rcp_emissions.NMVOC)

    # ---- Methane Radiative Forcing ---- #
    set_param!(m, :rf_ch4, :N₂O_0, N₂O_0)
    set_param!(m, :rf_ch4, :CH₄_0, CH₄_0)
    set_param!(m, :rf_ch4, :scale_CH₄, 1.0)

    # ---- Straospheric Water Vapor From Methane Radiative Forcing ---- #
    set_param!(m, :rf_ch4h2o, :CH₄_0, CH₄_0)
    set_param!(m, :rf_ch4h2o, :STRATH2O, 0.15)

    # ---- Tropospheric Ozone Radiative Forcing ---- #
    set_param!(m, :rf_o3, :CH₄_0, CH₄_0)
    set_param!(m, :rf_o3, :OZ00CH4, 0.161)
    set_param!(m, :rf_o3, :TROZSENS, 0.042)
    set_param!(m, :rf_o3, :OZNOX, 0.125)
    set_param!(m, :rf_o3, :OZCO, 0.0011)
    set_param!(m, :rf_o3, :OZVOC, 0.0033)
    set_param!(m, :rf_o3, :index_2000, index_2000)
    set_param!(m, :rf_o3, :NOx_emissions, rcp_emissions.NOx)
    set_param!(m, :rf_o3, :CO_emissions, rcp_emissions.CO)
    set_param!(m, :rf_o3, :NMVOC_emissions, rcp_emissions.NMVOC)
    set_param!(m, :rf_o3, :FOSSHIST, foss_hist_for_O₃.EFOSS)
    set_param!(m, :rf_o3, :TROZSENS, 0.042)
    set_param!(m, :rf_o3, :OZCH4, 5.0)

    # ---------------------------------------------
    # Create connections between Mimi components.
    # ---------------------------------------------
    connect_param!(m, :rf_ch4,    :CH₄, :ch4_cycle, :CH₄)
    connect_param!(m, :rf_ch4h2o, :CH₄, :ch4_cycle, :CH₄)
    connect_param!(m, :rf_o3,     :CH₄, :ch4_cycle, :CH₄)

    return(m)
end

end #module
