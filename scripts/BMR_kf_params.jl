module Kf_model_params

export Kf_data, constant_inputs, Data2

struct Mol_W
    #molar masses of species in g/mol
    Ni::Float64
    S::Float64
    Cu::Float64
    O::Float64
    H::Float64
    Fe::Float64
    Rh::Float64
    Ru::Float64
    Ir::Float64
end

mol_W = Mol_W(58.6934, 32.065, 63.546, 15.999, 1.00784, 55.845, 102.906, 101.07, 192.217) #in g/mol

struct Shape_factors
    NiS::Float64
    Ni3S4::Float64
    Cu9S5::Float64
    FeOHSO4::Float64
    Rh2S3::Float64
    Rh::Float64
    RhO2::Float64 #What is supposed to happen to shape factors >>3? Should I just assume them to be 3
    RuS2::Float64
    Ru::Float64
    RuO2::Float64
    Ir2S3::Float64
    Ir::Float64
    IrO2::Float64
end

shape_factors = Shape_factors(1.6, 1.0, 3.3, 3.4, 1.9, 2.6, 3, 2.7, 1.9, 3, 2.1, 3.4, 3)

struct Reaction_kinetics
    k0::Vector{Float64}                     #Pre-exponential factors for calculation of the reaction rate constants according to Arrhenius equation
    activation_energies::Vector{Float64}    #kJ/mol, Activation energies for calculation of the reaction rate constants according to Arrhenius equation
    R::Float64                              #J/mol, Ideal gas constant
    shape_factors::Shape_factors            #structure containing values of the shape factors for various solid species
end

k0 = [2.33*1e9, 4.76*1e8, 1.68*1e10, 2*7.52*1e-5, 2*5.45*10, 1.43*1e3, 4.66*1e9, 1.89*1e8, 7.34*1e8, 1.2*1e5, 2.18*1e17, 2.38*1e6, 7.1*1e8, 2.99*1e17, 1.25*1e12, 9.16*1e6, 1.44*1e16, 9.29*1e14, 1.2*1e10, 2.47*1e22, 3.09*1e24]
#k0 = [2.33*1e9, 4.76*1e8, 1.68*1e10, 7.52*1e-5, 5.45*10, 1.43*1e3] #value in thesis is differet to value in model -> reaction rate pre-exponential factor is double value in thesis

activation_energies = [78.5*1000, 41.1*1000, 87.4*1000, -26.2*1000, -5.9*1000, 34.3*1000, 59.5*1000, 64.6*1000, 56.1*1000, 52.1*1000, 120.9*1000, 48.5*1000, 64.2*1000, 138.5*1000, 116.2*1000, 52.8*1000, 124.3*1000, 134.1*1000, 72.4*1000, 175.1*1000, 209.2*1000]

reaction_kinetics = Reaction_kinetics(k0, activation_energies, 8.314, shape_factors)

struct Spent_electrolyte_data
    #Specify composition of the spent electrolyte 
    Cu::Float64 #g/L
    Ni::Float64 #g/L
    Fe::Float64 #mg/L
    Rh::Float64 #mg/L 
    Ru::Float64 #mg/L 
    Ir::Float64 #mg/L
    H2SO4::Float64 #g/l
    density::Float64 #kg/L
end

#spent_e = Spent_electrolyte_data(22.8, 25.3, 857, 56.1, 193, 33.6, 36.4, 1.18)
spent_e = Spent_electrolyte_data(22.8, 25.3, 857, 56.1, 193, 33.6, 36.4, 1.18)  

struct Residue
    Cu::Float64 #Mass percentage Cu in first stage solids residue
    Ni::Float64 #Mass percentage Ni in first stage solids residue
    Fe::Float64 #Mass percentage Fe in first stage solids residue
    S::Float64 #Mass percentage S in first stage solids residue
    Rh::Float64 #ppm Rh in first stage solids residue
    Ru::Float64 #ppm Ru in first stage solids residue
    Ir::Float64 #ppm Ir in first stage solids residue
    density::Float64 #kg/L
    Cu_as_Cu9S5::Float64 #Fraction of copper in first stage residue present as Cu9S5, remainder as CuS
    Ni_as_NiS::Float64 #Fraction of nickel in first stage residue present as NiA, remainder as Ni3S4
    Rh_as_Rh2S3::Float64 #Fraction of rhodium in first stage residue present as Rh2S3, remainder as Rh
    Ru_as_RuS2::Float64 #Fraction of ruthenium in first stage residue present as RuS2, remainder as Ru
    Ir_as_Ir2S3::Float64 #Fraction of iridium in first stage residue present as Ir2S3, remainder as Ir
    Fe_as_FeOHSO4::Float64 #Fraction of iron in the first stage residue present as leachable Fe(OH)SO4, remainder not leachable
    CuS_leachable::Float64 #%, fraction of CuS fed to the autoclave that can be leached 
    #CuS_leachable is ignored in Dorfling (2012) and Miskin (2016)
end

residue = Residue(50, 13, 0.42, 26, 300, 880, 200, 4.45, 0.706, 0.85, 0.50, 0.7, 0.43, 1, 0.75)

mutable struct Stream
    massflow ::Float64 #kg/h
    T::Float64 #degC
    x_S ::Float64      #kg solid/kg total massflow          
    NiS ::Float64      #kg NiS/kg solid
    Ni3S4 ::Float64
    Cu9S5 ::Float64
    CuS  ::Float64
    FeOHSO4  ::Float64
    Rh2S3  ::Float64
    Rh_s  ::Float64
    RhO2  ::Float64
    RuS2  ::Float64
    Ru_s  ::Float64
    RuO2  ::Float64
    Ir2S3  ::Float64
    Ir_s  ::Float64
    IrO2  ::Float64    
    x_L  ::Float64     #kg liquid/kg total massflow
    Cu ::Float64       #kg solute / kg solvent
    Ni  ::Float64
    Fe  ::Float64
    Rh_l  ::Float64
    Ru_l  ::Float64
    Ir_l ::Float64
    H2SO4  ::Float64
    x_G :: Float64      #kg gas/kg total
    O2 ::Float64         #kg 02/kg gas
end

struct Autoclave
    comp_1_vol ::Float64 #m^3
    comp_2_vol ::Float64 #m^3
    comp_3_vol ::Float64 #m^3
    comp_4_vol ::Float64 #m^3
    comp_1_max_vol ::Float64 #m^3
    comp_2_max_vol ::Float64 #m^3
    comp_3_max_vol ::Float64 #m^3
    comp_4_max_vol ::Float64 #m^3
    total_volume ::Float64 #m^3
    pressure ::Float64 #bar
    outer_diameter ::Float64 #m
    outer_temperature ::Float64 # °C
    emissivity ::Float64 
    comp_1_length ::Float64 #m
    comp_2_length ::Float64 #m
    comp_3_length ::Float64 #m
    comp_4_length ::Float64 #m
    agitator_power ::Float64 #kW
    comp_width ::Float64 #m
    weir_height_AC1 ::Float64 #m
    weir_height_AC2 ::Float64 #m
    atmospheric_pressure ::Float64 #Pa
    tank_10_max_vol ::Float64 #m^3
    tank_20_max_vol ::Float64 #m^3
    tank_40_max_vol ::Float64 #m^3
    tank_50_max_vol ::Float64 #m^3
    tank_10_area ::Float64 #m^2
    tank_20_area ::Float64 #m^2
    tank_40_area ::Float64 #m^2
    tank_50_area ::Float64 #m^2
    AC3_area ::Float64 #m^2
    AC4_area ::Float64 #m^2
end

autoclave = Autoclave(14, 4.67, 4.67, 4.67, 18.67, 6.23, 6.23, 6.23, 46.7, 6.5, 2.5, 50, 0.15, 5.5, 1.8, 1.8, 1.7, 10, 1.8, 1.4, 1.4, 101325, 15.0, 35.0, 15.7, 8.3, 4.28, 7.53, 4.41, 2.89, 3.6, 3.4)

mutable struct Flow
    s1_solids_mass::Float64                #mass flow rate of first stage leach residue in stream 1
    s1_spent_volume::Float64               #volumetric flow rate of spent electrolyte in stream 1, l/h
    s1_solid_fraction_mean::Float64      #spent volume fraction in stream 1 (spent e volume/total volume)
    s2_volume::Float64                     #volumetric flow rate of spent electrolyte in stream 2, l/h
    s3_volume::Float64                     #volumetric flow rate of water in stream 3, l/h
    s3_mass_fraction_H2SO4::Float64        #mass fraction H2SO4 in formic filtrate stream (stream 3), l/h
    s23_volume::Float64                    #volumetric flow rate of sulphuric acid in stream 23, l/h, entering flash recycle tank
    s24_volume::Float64                    #volumetric flow rate of water in stream 23, l/h, entering flash recycle tank
    s5_temperature::Float64                #temperature of the stream 5, dC
    s9_mass::Float64                       #total mass flow rate of the flash recycle stream, kg/h
    oxygen_total_excess::Float64           #factor by which the theoretical amount of oxygen (required for complete leaching of NiS and Cu9S5) is multiplied
    oxygen_fraction_10::Float64            #fraction of fresh oxygen to compartment 1 (stream 10)
    oxygen_fraction_11::Float64            #fraction of fresh oxygen to compartment 2 (stream 11)
    oxygen_fraction_12::Float64            #fraction of fresh oxygen to compartment 3 (stream 12)
    s18_volume::Float64                    #volumetric flow rate of spent electrolyte in stream 18, l/h
    s19_volume::Float64                    #volumetric flow rate of water in stream 19, l/h
    s20_volume::Float64                    #volumetric flow rate of acid in stream 20, l/h
    oxygen_temp::Float64                   #oxygen temperature, dC
    steam_temp::Float64                    #Steam temperature, dC
    stream_21_temp::Float64                #temperature of stream 21, dC
    acid_temp::Float64                     #temperature of acid streams, dC
    formic_temp::Float64                   #temperature of formic filtrate stream, dC
    water_temp::Float64                    #temperature of water stream, dC
    spent_temp::Float64                    #temperature of spent electrolyte streams, dC
    slurry_temp::Float64                   #temperature of slurry stream, dC
    AC2_cooling_water_max::Float64         #max flow rate of cooling water
end

#flow = Flow(800, 680, 1500, 2800, 0.015, 10.0, 10.0, 50, 8000.0, 1.2, 0.55, 0.2, 1-0.55-0.2, 400, 800, 7.0, 25, 180, 35.0, 35.0, 35.0, 35.0, 50.0, 50.0)
flow = Flow(800, 1315, 0.34, 1650, 2010, 0.015, 40.0, 0.0, 50, 9000.0, 1.2, 0.55, 0.2, 1-0.55-0.2, 1100, 100, 22.0, 25, 180, 35.0, 35.0, 35.0, 35.0, 50.0, 50.0, 10000)

struct Acid
    concentration::Float64 #mass percentage H2SO4
    density::Float64 #density of acid, kg/l
end

acid = Acid(96, 1.8355)

struct Air
    temp::Float64                      #Ambient temperature, dC
    k::Float64                         #Thermal conductivity of air, W/m.K
    viscosity::Float64                 #Viscosity of air, m^2/s
    thermal_diff::Float64              #Thermal diffusivity of air, m^2/s
    Pr::Float64                        #Prandtl number of air at 300K
    thermal_expansion::Float64         #Thermal expansion coefficient, K^-1
end

air = Air(25, 0.0263, 15.9/1000000, 22.5/1000000, 0.707, 2.725/1000)

struct Heat_capacities
    H2SO4_a::Float64    #Heat capacity constant of sulphyric acid, Cp = a +bT, T in dC, Cp in kJ/mol.dC
    H2SO4_b::Float64    #Heat capacity constant of sulphyric acid, Cp = a +bT, T in dC, Cp in kJ/mol.dC
    H2O::Float64        #Heat capacity of water, kJ/mol.dC
    Steam_a::Float64    #Heat capacity of steam, Cp = a + bT +cT^2 + dT^3, T in dC, Cp i kJ/mol.dC
    Steam_b::Float64    #Heat capacity of steam, Cp = a + bT +cT^2 + dT^3, T in dC, Cp i kJ/mol.dC
    Steam_c::Float64    #Heat capacity of steam, Cp = a + bT +cT^2 + dT^3, T in dC, Cp i kJ/mol.dC
    Steam_d::Float64    #Heat capacity of steam, Cp = a + bT +cT^2 + dT^3, T in dC, Cp i kJ/mol.dC
    Oxygen_a::Float64   #Heat capacity of oxygen, Cp = a + bT +cT^2 + dT^3, T in dC, Cp i kJ/mol.dC
    Oxygen_b::Float64   #Heat capacity of oxygen, Cp = a + bT +cT^2 + dT^3, T in dC, Cp i kJ/mol.dC
    Oxygen_c::Float64   #Heat capacity of oxygen, Cp = a + bT +cT^2 + dT^3, T in dC, Cp i kJ/mol.dC
    Oxygen_d::Float64   #Heat capacity of oxygen, Cp = a + bT +cT^2 + dT^3, T in dC, Cp i kJ/mol.dC
    NiS::Float64        #Heat capacity of NiS, J/mol.K
    CuS::Float64        #Heat capacity of CuS, J/mol.K
    Cu9S5::Float64      #Heat capacity of Cu9S5, J/mol.K
    Ni3S4_a::Float64    #Heat capacity of Ni3S4, Cp = a + bT + cT^2 + dT^3 + e/(T^2), T in K, Cp in J/mol.K
    Ni3S4_b::Float64    #Heat capacity of Ni3S4, Cp = a + bT + cT^2 + dT^3 + e/(T^2), T in K, Cp in J/mol.K
    Ni3S4_c::Float64    #Heat capacity of Ni3S4, Cp = a + bT + cT^2 + dT^3 + e/(T^2), T in K, Cp in J/mol.K
    Ni3S4_d::Float64    #Heat capacity of Ni3S4, Cp = a + bT + cT^2 + dT^3 + e/(T^2), T in K, Cp in J/mol.K
    Ni3S4_e::Float64    #Heat capacity of Ni3S4, Cp = a + bT + cT^2 + dT^3 + e/(T^2), T in K, Cp in J/mol.K
end

cP = Heat_capacities(139.1/1000, 15.59/100000, 75.4/1000, 33.46/1000, 0.688/100000, 0.7604/(10^8), -3.593/(10^12), 29.1/1000, 1.158/100000, -0.6076/(10^8), 1.311/(10^12), 52, 43.7, 158, 121.9, 143.66/1000, 0.02537/(10^6), -0.01788/(10^9), 0.000434*(10^6))

struct Gen_constants
    H_Evap_H2O::Float64             #Heat of vaporisation for water at 100dC
    R::Float64                      #Ideal gas constant, J/mol.K
    atmospheric_pressure::Float64   #Atmospheric pressure, bar
end

gen_constants = Gen_constants(2257, 8.314, 1)

struct H_reaction
    h_reaction::Vector{Float64}
end

h_reaction = H_reaction([-283.8, -852.3, -6251.2, -756.1, -765.3, -83.1, -13685.6, -3110.2, -13814.1, -3142.3, -13800.0, -3140.0, -2037.6, -963.0, -15.4, -5117.3, -1139.8, 16.7, -2000.0, -1150.0, 10.0])

struct Q_removal
    delta_T_AC2::Float64
    delta_T_AC3::Float64
    delta_T_TK40::Float64
    m_13::Float64
end

q_removal = Q_removal(20.0, 20.0, 20.0, 50.0)

struct Flow_constants
    prep_tank2::Float64
    recycle_tank::Float64
    discharge_tank::Float64
    thickener::Float64
    prep_tank_3::Float64
    AC1 ::Float64
    AC2 ::Float64
    AC3 ::Float64
    AC4 ::Float64
    O2out ::Float64
    overflow ::Float64
end

flow_constants = Flow_constants(69.14, 207.85, 79.53, 79.53, 18.57, 0.062, 0.062, 91.0, 26, 0.00581, 25.0)

struct Control
    setpoint::Vector{Float64}
    bias::Vector{Float64}
    Kc::Vector{Float64}
    Ti::Vector{Float64}
end

setpointvector = [1.49490900588678, #1 - Tank 10 level
                    5*0.539208941805876, #2 - Tank 20 level
                    1.07417580142121, #3 - AC3 level
                    1.14804813024663, #4 - Tank 40 level
                    2.15338181620428, #5 - Tank 50 level
                    1.04843173186243, #6 - AC4 level
                    130.730542551484, #7 - AC1 temp
                    129.971628898129, #8 - AC2 temp
                    125.755146989207, #9 - AC3 temp
                    140.291446131407, #10 - AC4 temp
                    6.5, #11 - Pressure
                    4848.16511140736, #12 - Stream 5 flow rate
                    40.0, #13 - Stream 23 flow rate (acid addition to flash recycle tank)
                    0.0, #14 - Stream 24 flow rate
                    1100.0, #15 - Stream 18 flow rate
                    22.0, #16 - Stream 20 flow rate (acid addition to 3rd stage slurry prep tank)
                    1274.21320937669, #17 - Stream 21 flow rate
                    12.5, #18 - percentage solids in stream 5
                    42.5, #19 - spent volume percentage (input for volume ratio control between stream 2 and 3)
                    75.0 #20 - TK40 temperature
                    ]
#biasvector = [1494.77528089887, 11767.9759357334, 5133.43292214034, 5133.43292214032, 100.0, 1378.26367945271, 7261.39524426505, 500000, 500000, 50, 458.341013548438]
biasvector = [1472.235143, 11654.24591, 5143.601488, 5143.599924, 100.4388658, 1394.923108, 7138.315376, 5977.051433, 5982.974719, 69.50727363, 458.6299304, 15184.86177]
gainsvector = [784.0, -18642.0, -5457.0, -2139.4, 812/5, -1384, -471.25, -336.144947, -381.0643508, 10*1.13, 72.7, -389.2]
integraltimevector = [21.8, 1.62, 2.6385, 8.254, 14.22*5, 9.828, 2.336, 0.847, 0.527, 12.05/5, 1.312, 0.701]
control = Control(setpointvector, biasvector, gainsvector, integraltimevector)

#Generate autoregressive random function inputs and add them into the main parameters folder
using DataInterpolations
using Random
using Plots

Random.seed!(42)

    # ARG function
    """
        generate_ar1(t_span, dt; phi, sigma, mu) -> (t_ar, Fin2_ar)

    Generates a stationary AR(1) process:
        x[k] = mu + phi*(x[k-1] - mu) + sigma*randn()

    - phi:  AR coefficient (|phi| < 1 for stationarity)
    - sigma: noise standard deviation
    - mu:   mean (nominal Fin2 = 0.2)
    """
    function generate_ar1(t_span, dt; phi, sigma, mu)
        t_ar = collect(t_span[1]:dt:t_span[2])
        n    = length(t_ar)
        x    = zeros(n)
        x[1] = mu                              # start at mean
        for k in 2:n
            x[k] = mu + phi * (x[k-1] - mu) + sigma * randn()
            x[k] = max(0.0, x[k])             # flow can't be negative
        end
        return t_ar, x
    end

    t_span = [0.0, 6000]
    dt = 0.1 #h

    Random.seed!(42)
    t_array_stream_1, random_stream_1 = generate_ar1(t_span, dt; phi = 0.99, sigma = 1, mu = 44)
    t_array_stream_1, constant_stream_1 = generate_ar1(t_span, dt; phi = 0.0, sigma = 0.0, mu = 44)

    stream_1_linear_interp = LinearInterpolation(random_stream_1, t_array_stream_1)
    stream_1_constant = LinearInterpolation(constant_stream_1, t_array_stream_1)
    plot(stream_1_linear_interp, seriestype = :line)
    plot(t_array_stream_1, random_stream_1, seriestype = :line, xlab = "time (h)", ylab = "temperature (deg C)", lab = "stream 1")

    Random.seed!(43)
    t_array_stream_2, random_stream_2 = generate_ar1(t_span, dt; phi = 0.99, sigma = 1, mu = 44)
    t_array_stream_2, constant_stream_2 = generate_ar1(t_span, dt; phi = 0.0, sigma = 0.0, mu = 44)
    stream_2_linear_interp = LinearInterpolation(random_stream_2, t_array_stream_2)
    stream_2_constant = LinearInterpolation(constant_stream_2, t_array_stream_2)
    plot(stream_2_linear_interp)
    plot!(t_array_stream_2, random_stream_2, seriestype = :line, lab = "stream 2")

    Random.seed!(44)
    t_array_stream_3, random_stream_3 = generate_ar1(t_span, dt; phi = 0.99, sigma = 1, mu = 44)
    t_array_stream_3, constant_stream_3 = generate_ar1(t_span, dt; phi = 0.0, sigma = 0.0, mu = 44)    
    stream_3_linear_interp = LinearInterpolation(random_stream_3, t_array_stream_3)
    stream_3_constant = LinearInterpolation(constant_stream_3, t_array_stream_3)
    plot(stream_3_linear_interp)

    # Random.seed!(45)
    # #solid fraction in stream 1
    # t_array_solid_fraction_1, random_solid_fraction_1 = generate_ar1(t_span, dt; phi = 0.0, sigma = 0.0, mu = 0.45)    
    # stream_1_solid_fraction = LinearInterpolation(random_solid_fraction_1, t_array_solid_fraction_1)
    # plot(stream_1_solid_fraction)

    t_array_stream_23 = t_span[1]:dt:t_span[2]
    sinusoid_stream_23 = 20.5 .+ 1.5 .* sin.(2π .* t_array_stream_23 ./ 24)
    constant_stream_23 = 20.5 * fill(1, length(t_array_stream_23))
    stream_23_linear_interp = LinearInterpolation(sinusoid_stream_23, t_array_stream_23)
    stream_23_constant = LinearInterpolation(constant_stream_23, t_array_stream_23)
    plot(stream_23_linear_interp)

    Random.seed!(48)
    t_array_solid_frac_s1, solid_frac_s1_interp = generate_ar1(t_span, dt; phi = 0.99, sigma = 0.003, mu = flow.s1_solid_fraction_mean)
    t_array_solid_frac_s1, solid_frac_s1_constant = generate_ar1(t_span, dt; phi = 0.0, sigma = 0.0, mu = flow.s1_solid_fraction_mean)
    stream_1_solid_frac_interp = LinearInterpolation(solid_frac_s1_interp, t_array_solid_frac_s1)
    stream_1_solid_frac_constant = LinearInterpolation(solid_frac_s1_constant, t_array_solid_frac_s1)
    plot(t_array_solid_frac_s1, solid_frac_s1_interp, seriestype = :line, xlab = "time (h)", ylab = "stream 1 solid fraction")

    Random.seed!(45)
    t_array_copper, copper_linear_interp = generate_ar1(t_span, dt; phi = 0.999, sigma = 0.001, mu = 0.5)
    t_array_copper, copper_linear_constant = generate_ar1(t_span, dt; phi = 0.0, sigma = 0.0, mu = 0.5)
    copper_linear_interp = LinearInterpolation(copper_linear_interp, t_array_copper)
    copper_constant = LinearInterpolation(copper_linear_constant, t_array_copper)
    plot(copper_linear_interp)

    Random.seed!(47)
    t_array_ni, ni_linear_interp = generate_ar1(t_span, dt; phi = 0.999, sigma = 0.001, mu = 0.13)
    t_array_ni, ni_linear_constant = generate_ar1(t_span, dt; phi = 0.0, sigma = 0.0, mu = 0.13)    
    ni_linear_interp = LinearInterpolation(ni_linear_interp, t_array_ni)
    ni_constant = LinearInterpolation(ni_linear_constant, t_array_ni)

    plot(ni_linear_interp)

struct Auto_regressive_inputs{T1, T2, T3, T4, T5, T6, T7}
    stream_1_linear_interp::T1
    stream_1_solid_frac_interp::T2
    stream_2_linear_interp::T3
    stream_3_linear_interp::T4
    stream_23_linear_interp::T5
    copper_linear_interp::T6
    ni_linear_interp::T7
end

auto_regressive_inputs = Auto_regressive_inputs(stream_1_linear_interp, stream_1_solid_frac_interp, stream_2_linear_interp, stream_3_constant, stream_23_constant, copper_constant, ni_constant)
constant_inputs = Auto_regressive_inputs(stream_1_constant, stream_1_solid_frac_constant, stream_2_constant, stream_3_constant, stream_23_constant, copper_constant, ni_constant)

#Things to consider:
    #resampling linear interpolation every time could become a serious bottleneck. 
        #might need to
        #   pre-sample onto solver time grid
        #   reduce interpolation frequency
        #   cache values if reused within a step

mutable struct Data2{T<:Auto_regressive_inputs}
    mol_W::Mol_W
    reaction_kinetics::Reaction_kinetics
    cP::Heat_capacities
    gen_constants::Gen_constants
    air::Air
    h_reaction::H_reaction
    spent_e::Spent_electrolyte_data
    residue::Residue
    acid::Acid
    autoclave::Autoclave
    flow::Flow
    q_removal::Q_removal
    flow_constants::Flow_constants
    control::Control
    auto_regressive_inputs::T
    faultmode::Int
end

Kf_data = Data2(mol_W, reaction_kinetics, cP, gen_constants, air, h_reaction, spent_e, residue, acid, autoclave, flow, q_removal, flow_constants, control, auto_regressive_inputs,0)

end