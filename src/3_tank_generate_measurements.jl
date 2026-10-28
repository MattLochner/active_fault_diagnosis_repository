function generate_measurements_3_tank_system(Real_Process, i, vk_std)
    yk = Real_Process[3, i] + vk_std[1]*randn()                            #measurement function for 3 tank system adds noise to third tank level measurement
    return [yk]
end

#Measurement function with noise