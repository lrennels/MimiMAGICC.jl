@defcomp rfco2 begin
    atmco2 = Parameter(index=[time])
    rf_co2 = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.rf_co2[t] = 3.7 * log(p.atmco2[t]/p.atmco2[1])/log(2.0)
    end
    
end
