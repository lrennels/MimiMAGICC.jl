using Test
using CSVFiles
using DataFrames
using Mimi
using MimiMAGICC

@testset "MAGICC" begin

#------------------------------------------------------------------------------
#   1. Carry out test to check that the model runs.
#------------------------------------------------------------------------------

    @testset "MAGICC-model" begin

    m = MimiMAGICC.get_magicc_ch4()
    run(m)

    end # MAGICC-CH4 model run test.
end # All MAGICC-CH4 tests.
