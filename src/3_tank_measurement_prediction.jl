#This is the measurement function without noise, i.e. perfect measurement of state

function measurement_prediction_3_tank_system(x, params, t)
    y = x[3]                            #in the case of 3 tank system, the level in tank 3 is measured = X[3]    
    return [y]
end

