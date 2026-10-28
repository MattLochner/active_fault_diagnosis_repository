function BMR_gt_measurement_prediction(X, data, t)
    
    faultmode = data.faultmode

    residue_solid_fraction = zeros(3)
            residue_solid_fraction[1] = data.auto_regressive_inputs.copper_linear_interp(t)
            residue_solid_fraction[2] = data.auto_regressive_inputs.ni_linear_interp(t)
    residue_solid_fraction[3] = data.residue.S/100

    Cu_as_Cu9S5 = data.residue.Cu_as_Cu9S5;
    Ni_as_NiS = data.residue.Ni_as_NiS;
    Rh_as_Rh2S3 = data.residue.Rh_as_Rh2S3;
    Ru_as_RuS2 = data.residue.Ru_as_RuS2;
    Ir_as_Ir2S3 = data.residue.Ir_as_Ir2S3;
    Fe_as_FeOHSO4 = data.residue.Fe_as_FeOHSO4;
    CuS_leachable = data.residue.CuS_leachable;

    spent_concentration = zeros(7)
    spent_concentration[1] = data.spent_e.Cu
    spent_concentration[2] = data.spent_e.Ni
    spent_concentration[3] = data.spent_e.H2SO4

    density_acid = data.acid.density
    density_spent = data.spent_e.density
    density_residue = data.residue.density
    acid_concentration = data.acid.concentration
    MW = data.mol_W
    Cp_H2SO4_a = 1000*data.cP.H2SO4_a / (2 * MW.H + MW.S + 4 * MW.O)
    Cp_H2SO4_b = 1000*data.cP.H2SO4_b / (2 * MW.H + MW.S + 4 * MW.O)
    Cp_H2O = 1000*data.cP.H2O / (2 * MW.H + MW.O)
    Cp_steam_a = 1000*data.cP.Steam_a/(2*MW.H+MW.O);
    Cp_steam_b = 1000*data.cP.Steam_b/(2*MW.H+MW.O);
    Cp_steam_c = 1000*data.cP.Steam_c/(2*MW.H+MW.O);
    Cp_steam_d = 1000*data.cP.Steam_d/(2*MW.H+MW.O);
    Cp_O2_a = 1000*data.cP.Oxygen_a/(2*MW.O);
    Cp_O2_b = 1000*data.cP.Oxygen_b/(2*MW.O);
    Cp_O2_c = 1000*data.cP.Oxygen_c/(2*MW.O);
    Cp_O2_d = 1000*data.cP.Oxygen_d/(2*MW.O);
    Cp_NiS = data.cP.NiS/(MW.Ni+MW.S);
    Cp_CuS = data.cP.CuS/(MW.Cu+MW.S);
    Cp_Cu9S5 = data.cP.Cu9S5/(9*MW.Cu+5*MW.S);
    Cp_Ni3S4_a = data.cP.Ni3S4_a/(3*MW.Ni+4*MW.S);
    Cp_Ni3S4_b = data.cP.Ni3S4_b/(3*MW.Ni+4*MW.S);
    Cp_Ni3S4_c = data.cP.Ni3S4_c/(3*MW.Ni+4*MW.S);
    Cp_Ni3S4_d = data.cP.Ni3S4_d/(3*MW.Ni+4*MW.S);
    Cp_Ni3S4_e = data.cP.Ni3S4_e/(3*MW.Ni+4*MW.S);

    H_Evap_H2O = data.gen_constants.H_Evap_H2O
    h_reaction = data.h_reaction

    s1_solid_fraction_mean = data.flow.s1_solid_fraction_mean
    s3_mass_fraction_H2SO4 = data.flow.s3_mass_fraction_H2SO4
    T_5 = data.flow.s5_temperature
    T_21 = data.flow.stream_21_temp
    fraction_oxygen_stream_10 = data.flow.oxygen_fraction_10
    fraction_oxygen_stream_11 = data.flow.oxygen_fraction_11

        oxygen_T = data.flow.oxygen_temp
        steam_T = data.flow.steam_temp
        formic_T = data.auto_regressive_inputs.stream_3_linear_interp(t)
        water_T = data.auto_regressive_inputs.stream_23_linear_interp(t)
        spent_T = data.auto_regressive_inputs.stream_2_linear_interp(t)
        slurry_T = data.auto_regressive_inputs.stream_1_linear_interp(t)
        acid_T = data.auto_regressive_inputs.stream_23_linear_interp(t)
        s1_solid_fraction = data.auto_regressive_inputs.stream_1_solid_frac_interp(t)

    agitator_power = 3600*data.autoclave.agitator_power
    autoclave_outer_diameter = data.autoclave.outer_diameter
    autoclave_surface_temperature = data.autoclave.outer_temperature
    autoclave_emissivity = data.autoclave.emissivity
    length_1 = data.autoclave.comp_1_length
    length_2 = data.autoclave.comp_2_length
    length_3 = data.autoclave.comp_3_length
    length_4 = data.autoclave.comp_4_length
    weir_height_AC1 = data.autoclave.weir_height_AC1
    weir_height_AC2 = data.autoclave.weir_height_AC2
    comp_width = data.autoclave.comp_width
    P_atm = data.autoclave.atmospheric_pressure
    volume_total = 1000*data.autoclave.total_volume

    tank_10_max_vol = data.autoclave.tank_10_max_vol * 1000
    tank_20_max_vol = data.autoclave.tank_20_max_vol * 1000
    tank_40_max_vol = data.autoclave.tank_40_max_vol * 1000
    tank_50_max_vol = data.autoclave.tank_50_max_vol * 1000

    AC1_max_vol = data.autoclave.comp_1_max_vol * 1000
    AC2_max_vol = data.autoclave.comp_2_max_vol * 1000
    AC3_max_vol = data.autoclave.comp_3_max_vol * 1000
    AC4_max_vol = data.autoclave.comp_4_max_vol * 1000

    tank_10_area = data.autoclave.tank_10_area
    tank_20_area = data.autoclave.tank_20_area
    tank_40_area = data.autoclave.tank_40_area
    tank_50_area = data.autoclave.tank_50_area
    AC3_area = data.autoclave.AC3_area
    AC4_area = data.autoclave.AC4_area

    volume_tank10 = X[2]*X[1]/density_residue + X[7]*X[1]/density_spent
    volume_tank20 = X[12]*X[11]/density_residue + X[17]*X[11]/density_spent
    volume_tank40 = X[53]*X[52]/density_residue + X[58]*X[52]/density_spent
    volume_tank50 = X[63]*X[62]/density_residue + X[68]*X[62]/density_spent

    volume_1 = X[22] * X[96] / density_residue + X[27] * X[96] / density_spent
    volume_2 = X[32] * X[97] / density_residue + X[37] * X[97] / density_spent
    volume_3 = X[42] * X[98] / density_residue + X[47] * X[98] / density_spent
    volume_4 = X[72] * X[99] / density_residue + X[77] * X[99] / density_spent
    
    slurry_height_1 = (volume_1/1000)/(length_1*comp_width)
    slurry_height_2 = (volume_2/1000)/(length_2*comp_width)

    tank_10_percent_full = volume_tank10 / tank_10_max_vol * 100
    tank_20_percent_full = volume_tank20 / tank_20_max_vol * 100
    tank_40_percent_full = volume_tank40 / tank_40_max_vol * 100
    tank_50_percent_full = volume_tank50 / tank_50_max_vol * 100
    comp_1_percent_full = volume_1 / AC1_max_vol * 100
    comp_2_percent_full = volume_2 / AC2_max_vol * 100
    comp_3_percent_full = volume_3 / AC3_max_vol * 100
    comp_3_percent_full_gt = comp_3_percent_full
    comp_4_percent_full = volume_4 / AC4_max_vol * 100

    level_tank10 = (volume_tank10/1000) / tank_10_area
    level_tank20 = (volume_tank20/1000) / tank_20_area
    level_tank40 = (volume_tank40/1000) / tank_40_area
    level_tank50 = (volume_tank50/1000) / tank_50_area
    level_AC3 = (volume_3/1000) / AC3_area
    level_AC4 = (volume_4/1000) / AC4_area

    V_vap = (volume_total - volume_1 - volume_2 - volume_3 - volume_4) / 1000
    T_vap = 1/(length_1+length_2+length_3+length_4)*(length_1*X[31]+length_2*X[41]+length_3*X[51]+length_4*X[81])

    if T_vap < 60
        P_vap_H2O = ((10 ^ (8.10765 - 1750.286 / (T_vap + 235))) / 760) * 101325
    else
        P_vap_H2O = ((10 ^ (7.96681 - 1668.210 / (T_vap + 228))) / 760) * 101325
    end

    P_O2_pascal = X[103] * 1000 * 8.314 * (T_vap + 273.15) / V_vap

    y_O2_vap    = P_O2_pascal / (P_O2_pascal + P_vap_H2O)
    P_autoclave = (P_vap_H2O + P_O2_pascal) / 100000

    air = data.air
    flow_constants = data.flow_constants

            setpoint = data.control.setpoint
            bias = data.control.bias
            Kc = data.control.Kc
            Ti = data.control.Ti

                V_5 = setpoint[12]
                V_23 = setpoint[13]
                V_24 = setpoint[14]
                V_18 = setpoint[15]
                V_20 = setpoint[16]
                V_21 = setpoint[17]

                level_error_tank_10 = level_tank10 - setpoint[1]
                V_1 = - Kc[1] * (level_error_tank_10 + 1/Ti[1] * X[104]) + bias[1]
                V_1 = max(0.0, V_1)
                if V_1 == 0.0
                    level_error_tank_10 = 0.0
                end

                level_error_tank_20 = level_tank20 - setpoint[2]
                V_7 = -Kc[2] * (level_error_tank_20 + 1/Ti[2] * X[105]) + bias[2]
                V_7 = max(0.0, V_7)
                if V_7 == 0.0
                    level_error_tank_20 = 0.0
                end

                if faultmode == 8
                    level_measured_AC3 = level_AC3 + 0.3
                    comp_3_percent_full = (level_measured_AC3 * 1000 * AC3_area) / AC3_max_vol * 100
                else
                    level_measured_AC3 = level_AC3
                end
                level_error_AC3 = level_measured_AC3 - setpoint[3]
                V_14 = -Kc[3] * (level_error_AC3 + 1/Ti[3] * X[106]) + bias[3]
                V_14 = max(0.0, V_14)
                if V_14 == 0.0
                    level_error_AC3 = 0.0
                end

                level_error_tank_40 = level_tank40 - setpoint[4]
                V_15 = -Kc[4] * (level_error_tank_40 + 1/Ti[4] * X[107]) + bias[4]
                V_15 = max(0.0, V_15)
                if V_15 == 0.0
                    level_error_tank_40 = 0.0
                end

                level_error_tank_50 = level_tank50 - setpoint[5]
                V_19 = -Kc[5] * (level_error_tank_50 + 1/Ti[5] * X[108]) + bias[5]
                V_19 = max(0.0, V_19)
                if V_19 == 0.0
                    level_error_tank_50 = 0.0
                end

                level_error_AC4 = level_AC4 - setpoint[6]
                V_22 = -Kc[6] * (level_error_AC4 + 1/Ti[6] * X[109]) + bias[6]
                V_22 = max(0.0, V_22)
                if V_22 == 0.0
                    level_error_AC4 = 0.0
                end

                temperature_error_AC1 = X[31] - setpoint[7]
                V_9 = -Kc[7] * (temperature_error_AC1 + 1/Ti[7] * X[110]) + bias[7]
                V_9 = max(0.0, V_9)
                if V_9 == 0.0
                    temperature_error_AC1 = 0.0
                end

                temperature_error_AC2 = X[41] - setpoint[8]
                V_H2O_AC2 = -Kc[8] * (temperature_error_AC2 + 1/Ti[8] * X[111]) + bias[8]
                V_H2O_AC2 = max(0.0, V_H2O_AC2)
                if V_H2O_AC2 == 0.0
                    temperature_error_AC2 = 0.0
                end
                V_H2O_AC2_max = data.flow.AC2_cooling_water_max
                if V_H2O_AC2 > V_H2O_AC2_max
                    V_H2O_AC2 = V_H2O_AC2_max
                    temperature_error_AC2 = 0.0
                end

                if faultmode == 5
                    scaling_factor = data.flow.scaling_factor_initial + data.flow.scaling_factor_rate_change * t
                else
                    scaling_factor = 1.0
                end

                Q_Removed_2 = scaling_factor * V_H2O_AC2 * 1.0 * Cp_H2O * data.q_removal.delta_T_AC2

                temperature_error_AC3 = X[51] - setpoint[9]
                V_H2O_AC3 = -Kc[9] * (temperature_error_AC3 + 1/Ti[9] * X[112]) + bias[9]
                V_H2O_AC3 = max(0.0, V_H2O_AC3)
                if V_H2O_AC3 == 0.0
                    temperature_error_AC3 = 0.0
                end

                Q_Removed_3 = V_H2O_AC3 * 1.0 * Cp_H2O * data.q_removal.delta_T_AC3
                
                temperature_error_TK40 = X[101] - setpoint[20]
                V_H2O_TK40 = - Kc[12] * (temperature_error_TK40 + 1/Ti[12] * X[115]) + bias[12]
                
                V_H2O_TK40 = max(0.0, V_H2O_TK40)
                if V_H2O_TK40 == 0.0
                    temperature_error_TK40 = 0.0
                end
                
                Q_Removed_TK40 = V_H2O_TK40 * 1.0 * Cp_H2O * data.q_removal.delta_T_TK40

                temperature_error_AC4 = X[81] - setpoint[10]
                m_13 = -Kc[10] * (temperature_error_AC4 + 1/Ti[10] * X[113]) + bias[10]
                m_13 = max(0.0, m_13)
                if m_13 == 0.0
                    temperature_error_AC4 = 0.0
                end

                pressure_error = P_autoclave - setpoint[11]
                m_O2_total = -Kc[11] * (pressure_error + 1/Ti[11] * X[114]) + bias[11]
                m_O2_total = max(0.0, m_O2_total)
                if m_O2_total == 0.0
                    pressure_error = 0.0
                end                

                s1_spent_volume_fraction = ((1-s1_solid_fraction)/data.spent_e.density)/((s1_solid_fraction/density_residue) + ((1-s1_solid_fraction)/data.spent_e.density))

                V_spent_1 = V_1*s1_spent_volume_fraction
                m_spent_1 = V_spent_1 * data.spent_e.density
                m_S_1 = (V_1 - V_spent_1)*density_residue

                density_stream_1 = 1 / (s1_solid_fraction_mean / density_residue + (1-s1_solid_fraction_mean)/data.spent_e.density)

                density_stream_2 = data.spent_e.density
                density_stream_3 = 1.0
                solid_percentage = setpoint[18]
                spent_percentage = setpoint[19]

                V_2 = ((density_stream_1*V_1*(1/density_stream_1 - 1/data.spent_e.density)/(1/density_residue - 1/data.spent_e.density))*100 - solid_percentage * density_stream_1*V_1)/(solid_percentage*(data.spent_e.density + density_stream_3 * (100/spent_percentage-1)))
                V_2 = max(0.0, V_2)

                V_3 = (100 / spent_percentage - 1) * V_2             
                V_3 = max(0.0, V_3)

        if faultmode == 1
           display("valve blockage not modelled yet")
        end
        if faultmode == 2
            display("valve wear fault not modelled")
        end
        if faultmode == 3
            display("valve stiction fault not modelled")
        end
        if faultmode == 4
            max_V_7 = 11000
            if V_7 > max_V_7
                V_7 = max_V_7
            end
        end
        # if faultmode == 5
        #     max_V_H2O_AC2 = 6500
        #     if V_H2O_AC2 > max_V_H2O_AC2
        #         V_H2O_AC2 = max_V_H2O_AC2
        #         Q_Removed_2 = V_H2O_AC2 * 1.0 * Cp_H2O * data.q_removal.delta_T_AC2
        #         temperature_error_AC2 = 0.0
        #     end
        # end
        if faultmode == 6
            m_S_1 = 0.0
            V_spent_1 = 0.0
            V_1 = 0.0
        end
        if faultmode == 7
            display("sulfuric acid controller misuse must be modelled in runODE's function")
        end
        if faultmode == 101
            V_3_leaksize = 100
            if V_3 > 100
                V_3 = V_3 - 100
            else
                V_3 = V_3
            end
        end
        if faultmode == 102
            V_2_maximum = 100
            if V_2 > V_2_maximum
                V_2 = V_2_maximum
            else
                V_2 = V_2
            end
        end
        if faultmode == 103
            m_S_1 = m_S_1 * 0.9
            V_spent_1 = V_spent_1 * 0.9
            V_1 = V_1 * 0.9
        end

        m_2 = V_2 * data.spent_e.density
        m_3 = V_3
        m_H2SO4_4 = m_3 * s3_mass_fraction_H2SO4

        m_5 = V_5 / (X[2] / density_residue + X[7] / density_spent)
        m_7 = V_7 / (X[12] / density_residue + X[17] / density_spent)
        m_9 = V_9 / (X[22] / density_residue + X[27] / density_spent)
        m_14 = V_14 / (X[42] / density_residue + X[47] / density_spent)
        m_15 = V_15 / (X[53] / density_residue + X[58] / density_spent)
        m_21 = V_21 / (X[63] / density_residue + X[68] / density_spent)
        m_22 = V_22 / (X[72] / density_residue + X[77] / density_spent)
       
        m_liquid_in = zeros(Float64, 4)
        m_liquid_in[1] = m_2 + m_3 + m_spent_1
        m_liquid_in[2:3] = (V_spent_1 + V_2) .* spent_concentration[1:2] / 1000
        m_liquid_in[4] = m_H2SO4_4 + (V_spent_1 + V_2) * spent_concentration[3]/1000

        m_solids_in = zeros(Float64, 4)
        m_solids_in[1] = m_S_1
        m_solids_in[2:4] = m_solids_in[1] .* residue_solid_fraction

        n_phases_in = zeros(Float64, 4)
        n_phases_in[1] = (m_solids_in[3]/MW.Ni)*Ni_as_NiS;
        n_phases_in[2] = ((m_solids_in[3]/MW.Ni)*(1-Ni_as_NiS))/3;
        n_phases_in[3] = ((m_solids_in[2]/MW.Cu)*Cu_as_Cu9S5)/9;
        n_phases_in[4] = (m_solids_in[2]/MW.Cu)*(1-Cu_as_Cu9S5);

        m_phases_in = zeros(Float64, 14)
        m_phases_in[1] = n_phases_in[1]*(MW.Ni+MW.S);
        m_phases_in[2] = n_phases_in[2]*(3*MW.Ni+4*MW.S);
        m_phases_in[3] = n_phases_in[3]*(9*MW.Cu+5*MW.S);
        m_phases_in[4] = n_phases_in[4]*(MW.Cu+MW.S);

        m_23 = V_23 * data.acid.density
        m_H2SO4_23 = data.acid.concentration/100 * m_23
        m_24 = V_24

        m_18 = V_18 * density_spent
        m_19 = V_19
        m_20 = V_20 * density_acid
        m_H2SO4_20 = (acid_concentration/100)*m_20

        m_liquid_TK50 = zeros(Float64, 8)
        m_liquid_TK50[1] = m_18 + m_19 + m_20
        m_liquid_TK50[2:3] = V_18 * spent_concentration[1:2]/1000
        m_liquid_TK50[4] = m_H2SO4_20 + V_18 * spent_concentration[3]/1000

        if T_vap < 60
            P_water_2nd_stage = ((10^(8.10765-1750.286/(T_vap+235)))/760)*1.01325
            else
            P_water_2nd_stage = ((10^(7.96681-1668.210/(T_vap+228)))/760)*1.01325
        end

    y_O2_2nd_stage = 1-P_water_2nd_stage/P_autoclave;

        Ra = (9.81*air.thermal_expansion*(autoclave_surface_temperature-air.temp)*(autoclave_outer_diameter^3))/(air.viscosity*air.thermal_diff)
        Nu = (0.6+(0.387*Ra^(1/6))/((1+(0.559/air.Pr)^(9/16))^(8/27)))^2
        h = air.k*Nu/autoclave_outer_diameter
        Q_loss = (h*22/7*autoclave_outer_diameter*(autoclave_surface_temperature-air.temp)+autoclave_emissivity*22/7*autoclave_outer_diameter*5.67*(10^-8)*((autoclave_surface_temperature+273.15)^4-(air.temp+273.15)^4))*3600/1000

        m_10 = m_O2_total * fraction_oxygen_stream_10;
        m_11 = m_O2_total * fraction_oxygen_stream_11;
        m_12 = m_O2_total - m_11 - m_10;

        spent_concentration[3]
        density_spent
        acid_concentration

        X_H2SO4_spent = spent_concentration[3] / (density_spent * 1000)
        X_H2SO4_acid = acid_concentration / 100

        Cp_L_1 = X_H2SO4_spent * Cp_H2SO4_a + (1-X_H2SO4_spent) * Cp_H2O
        Cp_L_2 = X_H2SO4_spent * Cp_H2SO4_a + (1-X_H2SO4_spent) * Cp_H2O
        Cp_L_3 = s3_mass_fraction_H2SO4 * Cp_H2SO4_a + (1-s3_mass_fraction_H2SO4) * Cp_H2O
        Cp_L_5 = X[10]*Cp_H2SO4_a + (1-X[10])*Cp_H2O;
        Cp_L_23 = X_H2SO4_acid * Cp_H2SO4_a + (1 - X_H2SO4_acid) * Cp_H2O
        Cp_L_24 = Cp_H2O
        Cp_L_7 = X[20]*Cp_H2SO4_a+(1-X[20])*Cp_H2O;
        Cp_L_9 = X[30]*Cp_H2SO4_a+(1-X[30])*Cp_H2O;
        Cp_L_AC1 = Cp_L_9;
        Cp_L_AC2 = X[40]*Cp_H2SO4_a+(1-X[40])*Cp_H2O;
        Cp_L_21 = X[71]*Cp_H2SO4_a+(1-X[71])*Cp_H2O;
        Cp_L_22 = X[80]*Cp_H2SO4_a+(1-X[80])*Cp_H2O;
        Cp_L_14 = X[50]*Cp_H2SO4_a+(1-X[50])*Cp_H2O;
        Cp_L_15 = X[61]*Cp_H2SO4_a+(1-X[61])*Cp_H2O;
        Cp_L_16 = Cp_L_15
        Cp_L_17 = Cp_L_15#Cp_H2O                                # Stream 17 (unused: multiplied by zero in dX_dt[102])
        Cp_L_18 = X_H2SO4_spent * Cp_H2SO4_a + (1-X_H2SO4_spent) * Cp_H2O
        Cp_L_19 = Cp_H2O
        Cp_L_20 = X_H2SO4_acid * Cp_H2SO4_a + (1 - X_H2SO4_acid) * Cp_H2O

    Normalx_NiS_S_1 = m_phases_in[1] / (m_phases_in[1]+m_phases_in[2]+m_phases_in[3]+m_phases_in[4])
    Normalx_Ni3S4_S_1 = m_phases_in[2] / (m_phases_in[1]+m_phases_in[2]+m_phases_in[3]+m_phases_in[4])
    Normalx_Cu9S5_S_1 = m_phases_in[3] / (m_phases_in[1]+m_phases_in[2]+m_phases_in[3]+m_phases_in[4])
    Normalx_CuS_S_1 = m_phases_in[4] / (m_phases_in[1]+m_phases_in[2]+m_phases_in[3]+m_phases_in[4])

    Normalx_NiS_S_5 = X[3] / (X[3] + X[4] + X[5] + X[6])
    Normalx_Ni3S4_S_5 = X[4] / (X[3] + X[4] + X[5] + X[6])
    Normalx_Cu9S5_S_5 = X[5] / (X[3] + X[4] + X[5] + X[6])
    Normalx_CuS_S_5 = X[6] / (X[3] + X[4] + X[5] + X[6])

    Normalx_NiS_S_7 = X[13] / (X[13] + X[14] + X[15] + X[16])
    Normalx_Ni3S4_S_7 = X[14] / (X[13] + X[14] + X[15] + X[16])
    Normalx_Cu9S5_S_7 = X[15] / (X[13] + X[14] + X[15] + X[16])
    Normalx_CuS_S_7 = X[16] / (X[13] + X[14] + X[15] + X[16])

    Normalx_NiS_S_9 = X[23] / (X[23] + X[24] + X[25] + X[26])
    Normalx_Ni3S4_S_9 = X[24] / (X[23] + X[24] + X[25] + X[26])
    Normalx_Cu9S5_S_9 = X[25] / (X[23] + X[24] + X[25] + X[26])
    Normalx_CuS_S_9 = X[26] / (X[23] + X[24] + X[25] + X[26])

    Normalx_NiS_S_AC2 = X[33] / (X[33] + X[34] + X[35] + X[36])
    Normalx_Ni3S4_S_AC2 = X[34] / (X[33] + X[34] + X[35] + X[36])
    Normalx_Cu9S5_S_AC2 = X[35] / (X[33] + X[34] + X[35] + X[36])
    Normalx_CuS_S_AC2 = X[36] / (X[33] + X[34] + X[35] + X[36])

    Normalx_NiS_S_21 = X[64] / (X[64] + X[65] + X[66] + X[67])
    Normalx_Ni3S4_S_21 = X[65] / (X[64] + X[65] + X[66] + X[67])
    Normalx_Cu9S5_S_21 = X[66] / (X[64] + X[65] + X[66] + X[67])
    Normalx_CuS_S_21 = X[67] / (X[64] + X[65] + X[66] + X[67])

    Normalx_NiS_S_22 = X[73] / (X[73] + X[74] + X[75] + X[76])
    Normalx_Ni3S4_S_22 = X[74] / (X[73] + X[74] + X[75] + X[76])
    Normalx_Cu9S5_S_22 = X[75] / (X[73] + X[74] + X[75] + X[76])
    Normalx_CuS_S_22 = X[76] / (X[73] + X[74] + X[75] + X[76])

    Normalx_NiS_S_14 = X[43] / (X[43] + X[44] + X[45] + X[46])
    Normalx_Ni3S4_S_14 = X[44] / (X[43] + X[44] + X[45] + X[46])
    Normalx_Cu9S5_S_14 = X[45] / (X[43] + X[44] + X[45] + X[46])
    Normalx_CuS_S_14 = X[46] / (X[43] + X[44] + X[45] + X[46])

    Normalx_NiS_S_15 = X[54] / (X[54] + X[55] + X[56] + X[57])
    Normalx_Ni3S4_S_15 = X[55] / (X[54] + X[55] + X[56] + X[57])
    Normalx_Cu9S5_S_15 = X[56] / (X[54] + X[55] + X[56] + X[57])
    Normalx_CuS_S_15 = X[57] / (X[54] + X[55] + X[56] + X[57])


    Normalx_NiS_S_21 = X[64] / (X[64] + X[65] + X[66] + X[67])
    Normalx_Ni3S4_S_21 = X[65] / (X[64] + X[65] + X[66] + X[67])
    Normalx_Cu9S5_S_21 = X[66] / (X[64] + X[65] + X[66] + X[67])
    Normalx_CuS_S_21 = X[67] / (X[64] + X[65] + X[66] + X[67])
    
    Cp_S_1 = Normalx_NiS_S_1*Cp_NiS + Normalx_Cu9S5_S_1*Cp_Cu9S5 + Normalx_CuS_S_1*Cp_CuS + Normalx_Ni3S4_S_1*Cp_Ni3S4_a;
    Cp_S_5 = Normalx_NiS_S_5*Cp_NiS + Normalx_Cu9S5_S_5*Cp_Cu9S5 + Normalx_CuS_S_5*Cp_CuS + Normalx_Ni3S4_S_5*Cp_Ni3S4_a;
    Cp_S_7 = Normalx_NiS_S_7*Cp_NiS + Normalx_Cu9S5_S_7*Cp_Cu9S5 + Normalx_CuS_S_7*Cp_CuS + Normalx_Ni3S4_S_7*Cp_Ni3S4_a;
    Cp_S_9 = Normalx_NiS_S_9*Cp_NiS + Normalx_Cu9S5_S_9*Cp_Cu9S5 + Normalx_CuS_S_9*Cp_CuS + Normalx_Ni3S4_S_9*Cp_Ni3S4_a;
    Cp_S_AC1 = Cp_S_9;
    Cp_S_AC2 = Normalx_NiS_S_AC2*Cp_NiS + Normalx_Cu9S5_S_AC2*Cp_Cu9S5 + Normalx_CuS_S_AC2*Cp_CuS + Normalx_Ni3S4_S_AC2*Cp_Ni3S4_a;
    Cp_S_21 = Normalx_NiS_S_21*Cp_NiS + Normalx_Cu9S5_S_21*Cp_Cu9S5 + Normalx_CuS_S_21*Cp_CuS + Normalx_Ni3S4_S_21*Cp_Ni3S4_a;
    Cp_S_22 = Normalx_NiS_S_22*Cp_NiS + Normalx_Cu9S5_S_22*Cp_Cu9S5 + Normalx_CuS_S_22*Cp_CuS + Normalx_Ni3S4_S_22*Cp_Ni3S4_a;
    Cp_S_14 = Normalx_NiS_S_14*Cp_NiS + Normalx_Cu9S5_S_14*Cp_Cu9S5 + Normalx_CuS_S_14*Cp_CuS + Normalx_Ni3S4_S_14*Cp_Ni3S4_a;
    Cp_S_15 = Normalx_NiS_S_15*Cp_NiS + Normalx_Cu9S5_S_15*Cp_Cu9S5 + Normalx_CuS_S_15*Cp_CuS + Normalx_Ni3S4_S_15*Cp_Ni3S4_a;
    Cp_S_17 = Cp_S_15
    Cp_S_21 = Normalx_NiS_S_21*Cp_NiS + Normalx_Cu9S5_S_21*Cp_Cu9S5 + Normalx_CuS_S_21*Cp_CuS + Normalx_Ni3S4_S_21*Cp_Ni3S4_a;

    Q_removed_L_9 = (m_9*X[27])*Cp_L_9*(X[31]-100);
    Q_removed_S_9 = (m_9*X[22])*Cp_S_9*(X[31]-100);
    Q_removal_9 = Q_removed_L_9+Q_removed_S_9;

    if (Q_removal_9>0)
        m_H2O_flash_evaporated = Q_removal_9/H_Evap_H2O;
    else
        m_H2O_flash_evaporated = 0;
    end;
    if (m_H2O_flash_evaporated>m_9*X[27]*(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]))
        m_H2O_flash_evaporated = m_9*X[27]*(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]);
        temp_flag = 1;
    else
        temp_flag = 0;
    end;

        liq_fractions_1_H2SO4 = X[30]
        liq_fractions_1 = liq_fractions(liq_fractions_1_H2SO4, 0.0, 0.0, 0.0)

        solids_in_1_NiS   = X[82] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent))
        solids_in_1_Ni3S4 = X[83] * X[22]  * (1/(X[22]/density_residue + X[27]/density_spent))
        solids_in_1_Cu9S5 = -10000.0

        solids_in_1 = solid_fractions(solids_in_1_NiS, solids_in_1_Ni3S4, solids_in_1_Cu9S5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        Cu9S5_fresh = X[2] * flow_constants.prep_tank2 * X[1]^0.5 * X[5] / V_5

        solids_out_1_NiS = X[22] * X[23] / (X[22] / density_residue + X[27] / density_spent)
        solids_out_1_Ni3S4 = X[22] * X[24] / (X[22] / density_residue + X[27] / density_spent)
        solids_out_1_Cu9S5 = X[22] * X[25] / (X[22] / density_residue + X[27] / density_spent)
        solids_out_1_CuS = X[22] * X[26] / (X[22] / density_residue + X[27] / density_spent)

        solids_out_1 = solid_fractions(solids_out_1_NiS, solids_out_1_Ni3S4, solids_out_1_Cu9S5, solids_out_1_CuS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        Cu_molal_9 = 1000*(X[28]/(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]))/MW.Cu;
        Ni_molal_9 = 1000*(X[29]/(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]))/MW.Ni;
        H2SO4_molal_9 = 1000*(X[30]/(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]))/(2*MW.H+MW.S+4*MW.O);
        O2_solubility_9 = (2*MW.O/1000)*(calculate_oxygen_solubility2(X[31],P_O2_pascal,Cu_molal_9,Ni_molal_9,H2SO4_molal_9));
        molal_O2_solubility_9 = (O2_solubility_9*1000)/(2*MW.O);

        r_1 = calculate_reaction_rates(X[31], liq_fractions_1, solids_in_1, solids_out_1, molal_O2_solubility_9, Cu9S5_fresh, data)

        extent_1 = zeros(Float64,5);
        extent_1[1] = volume_1*60*r_1[1]
        extent_1[2] = volume_1*60*r_1[2]
        extent_1[3] = volume_1*60*r_1[3]
        extent_1[4] = volume_1*60*r_1[4]
        extent_1[5] = volume_1*60*r_1[5]
        solid_consumption_1 = zeros(Float64,4);
        solid_consumption_1[1] = ((MW.Ni+MW.S)/1000)*(4*extent_1[1]+extent_1[2]);
        solid_consumption_1[2] = ((3*MW.Ni+4*MW.S)/1000)*(-extent_1[1]+2*extent_1[3]);
        solid_consumption_1[3] = ((9*MW.Cu+5*MW.S)/1000)*(extent_1[4]);
        solid_consumption_1[4] = ((MW.Cu+MW.S)/1000)*(-5*extent_1[4]+extent_1[5]);

        liquid_generation_1 = zeros(Float64,3);
        liquid_generation_1[1] = (4*extent_1[4]+extent_1[5])*(MW.Cu/1000);
        liquid_generation_1[2] = (extent_1[1]+extent_1[2]+6*extent_1[3])*(MW.Ni/1000);
        liquid_generation_1[3] = ((2*MW.H+MW.S+4*MW.O)/1000)*((-2*extent_1[1]+4*extent_1[3]-8*extent_1[4])/2);

        for i=(1:2)
            if (liquid_generation_1[i]<0.0)
                display("Problem: liquid metal species consumption!")
            end;
        end;

        O2_consumption_1 = (0.5*extent_1[1]+2*extent_1[2]+15*extent_1[3]+2*extent_1[4]+2*extent_1[5])*(2*MW.O/1000);
        n_O2_consumption_1 = O2_consumption_1/(2*MW.O);
        n_H2O_condensation_1 = (n_O2_consumption_1/y_O2_2nd_stage)*(1-y_O2_2nd_stage);
        H2O_condensation_1 = n_H2O_condensation_1*(2*MW.H+MW.O);

        V_out_AC1 = flow_constants.AC1*comp_width*real((complex(slurry_height_1-weir_height_AC1))^(3/2))
        m_AC1 = V_out_AC1*1000*(1/(X[22]/density_residue+X[27]/density_spent))*60*60

        Q_Reaction_1 = data.h_reaction.h_reaction[1:5].*extent_1; 
        Q_Cond_1 = H2O_condensation_1*(H_Evap_H2O+Cp_H2O*(100-25)+Cp_steam_a*(X[31]-100));        

        V_S_AC1 = X[22]*m_AC1/density_residue;
        V_L_AC1 = X[27]*m_AC1/density_spent;
        V_AC1 = V_S_AC1 + V_L_AC1;
        
        liq_fractions_2_H2SO4 = X[40];
        liq_fractions_2 = liq_fractions(liq_fractions_2_H2SO4, 0.0, 0.0, 0.0)

        solids_in_2_NiS = X[84] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent));
        solids_in_2_Ni3S4 = X[85] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent));
        solids_in_2_Cu9S5 = -10000.0

        solids_in_2 = solid_fractions(solids_in_2_NiS, solids_in_2_Ni3S4, solids_in_2_Cu9S5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        V_AC2 = zero(Float64)
        solid_consumption_2 = zeros(Float64,4)
        extent_2 = zeros(Float64,5)
        liquid_generation_2 = zeros(Float64,3);
        O2_consumption_2 = zero(Float64)
        H2O_condensation_2 = zero(Float64)
        m_AC2 = zero(Float64)

        solids_out_2_NiS = X[32] * X[33] / (X[32] / density_residue + X[37] / density_spent)
        solids_out_2_Ni3S4 = X[32] * X[34] / (X[32] / density_residue + X[37] / density_spent)
        solids_out_2_Cu9S5 = X[32] * X[35] / (X[32] / density_residue + X[37] / density_spent)
        solids_out_2_CuS = X[32] * X[36] / (X[32] / density_residue + X[37] / density_spent)

        solids_out_2 = solid_fractions(solids_out_2_NiS, solids_out_2_Ni3S4, solids_out_2_Cu9S5, solids_out_2_CuS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        Cu_molal_AC2 = 1000*(X[38]/(1-X[38]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[39]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[40]))/MW.Cu;
        Ni_molal_AC2 = 1000*(X[39]/(1-X[38]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[39]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[40]))/MW.Ni;
        H2SO4_molal_AC2 = 1000*(X[40]/(1-X[38]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[39]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[40]))/(2*MW.H+MW.S+4*MW.O);
        O2_solubility_AC2 = (2*MW.O/1000)*(calculate_oxygen_solubility2(X[41],P_O2_pascal,Cu_molal_AC2,Ni_molal_AC2,H2SO4_molal_AC2));
        molal_O2_solubility_AC2 = (O2_solubility_AC2*1000)/(2*MW.O);     
        
        r_2 = calculate_reaction_rates(X[41],liq_fractions_2,solids_in_2,solids_out_2,molal_O2_solubility_AC2,Cu9S5_fresh,data);

        extent_2 = zeros(Float64,5);
        extent_2[1] = volume_2*60*r_2[1]
        extent_2[2] = volume_2*60*r_2[2]
        extent_2[3] = volume_2*60*r_2[3]
        extent_2[4] = volume_2*60*r_2[4]
        extent_2[5] = volume_2*60*r_2[5]

        solid_consumption_2[1] = ((MW.Ni+MW.S)/1000)*(4*extent_2[1]+extent_2[2]);
        solid_consumption_2[2] = ((3*MW.Ni+4*MW.S)/1000)*(-extent_2[1]+2*extent_2[3]);
        solid_consumption_2[3] = ((9*MW.Cu+5*MW.S)/1000)*(extent_2[4]);
        solid_consumption_2[4] = ((MW.Cu+MW.S)/1000)*(-5*extent_2[4]+extent_2[5]);

        liquid_generation_2[1] = (4*extent_2[4]+extent_2[5])*(MW.Cu/1000);
        liquid_generation_2[2] = (extent_2[1]+extent_2[2]+6*extent_2[3])*(MW.Ni/1000);
        liquid_generation_2[3] = ((2*MW.H+MW.S+4*MW.O)/1000)*((-2*extent_2[1]+4*extent_2[3]-8*extent_2[4])/2);
    
        for i = 1:2
            if (liquid_generation_2[i]<0.0)           
                display("Problem: liquid metal species consumption!")
            end;
        end;

        O2_consumption_2 = (0.5*extent_2[1]+2*extent_2[2]+15*extent_2[3]+2*extent_2[4]+2*extent_2[5])*(2*MW.O/1000);
        O2_consumption_AC4 = O2_consumption_2-m_10;
        n_O2_consumption_AC4 = O2_consumption_AC4/(2*MW.O);
        n_H2O_condensation_2 = (n_O2_consumption_AC4/y_O2_2nd_stage)*(1-y_O2_2nd_stage);
        H2O_condensation_2 = n_H2O_condensation_2*(2*MW.H+MW.O);

        V_out_AC2 = flow_constants.AC2*comp_width*real((complex(slurry_height_2-weir_height_AC2))^(3/2))
        m_AC2 = V_out_AC2*1000*(1/(X[32]/density_residue+X[37]/density_spent))*60*60

        Q_Reaction_2 = data.h_reaction.h_reaction[1:5].*extent_2; 
        Q_Cond_2 = H2O_condensation_2*(H_Evap_H2O+Cp_H2O*(100-25)+Cp_steam_a*(X[41]-100));

        Cp_10 = Cp_O2_a;
        T_10 = oxygen_T;

        liq_fractions_3_H2SO4 = X[50];
        liq_fractions_3 = liq_fractions(liq_fractions_3_H2SO4, 0.0, 0.0, 0.0)

        solids_in_3_NiS = X[86] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent));
        solids_in_3_Ni3S4 = X[87] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent));
        solids_in_3_Cu9S5 = -10000.0

        solids_in_3 = solid_fractions(solids_in_3_NiS, solids_in_3_Ni3S4, solids_in_3_Cu9S5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        solid_consumption_3 = zeros(Float64,4)
        extent_3 = zeros(Float64,5)
        liquid_generation_3 = zeros(Float64,3);
        O2_consumption_3 = zero(Float64)
        H2O_condensation_3 = zero(Float64)

        solids_out_3_NiS = X[42] * X[43] / (X[42] / density_residue + X[47] / density_spent)
        solids_out_3_Ni3S4 = X[42] * X[44] / (X[42] / density_residue + X[47] / density_spent)
        solids_out_3_Cu9S5 = X[42] * X[45] / (X[42] / density_residue + X[47] / density_spent)
        solids_out_3_CuS = X[42] * X[46] / (X[42] / density_residue + X[47] / density_spent)
             
        solids_out_3 = solid_fractions(solids_out_3_NiS, solids_out_3_Ni3S4, solids_out_3_Cu9S5, solids_out_3_CuS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        Cu_molal_14 = 1000*(X[48]/(1-X[48]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[49]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[50]))/MW.Cu;
        Ni_molal_14 = 1000*(X[49]/(1-X[48]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[49]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[50]))/MW.Ni;
        H2SO4_molal_14 = 1000*(X[50]/(1-X[48]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[49]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[50]))/(2*MW.H+MW.S+4*MW.O);
        O2_solubility_14 = (2*MW.O/1000)*(calculate_oxygen_solubility2(X[51],P_O2_pascal,Cu_molal_14,Ni_molal_14,H2SO4_molal_14));
        molal_O2_solubility_14 = (O2_solubility_14*1000)/(2*MW.O);
        
        r_3 = calculate_reaction_rates(X[51],liq_fractions_3,solids_in_3,solids_out_3,molal_O2_solubility_14,Cu9S5_fresh,data);
        
            extent_3 = zeros(Float64,5);
            extent_3[1] = volume_3*60*r_3[1]
            extent_3[2] = volume_3*60*r_3[2]
            extent_3[3] = volume_3*60*r_3[3]
            extent_3[4] = volume_3*60*r_3[4]
            extent_3[5] = volume_3*60*r_3[5]        
        
            solid_consumption_3[1] = ((MW.Ni+MW.S)/1000)*(4*extent_3[1]+extent_3[2]);
            solid_consumption_3[2] = ((3*MW.Ni+4*MW.S)/1000)*(-extent_3[1]+2*extent_3[3]);
            solid_consumption_3[3] = ((9*MW.Cu+5*MW.S)/1000)*(extent_3[4]);
            solid_consumption_3[4] = ((MW.Cu+MW.S)/1000)*(-5*extent_3[4]+extent_3[5]);

            liquid_generation_3[1] = (4*extent_3[4]+extent_3[5])*(MW.Cu/1000);
            liquid_generation_3[2] = (extent_3[1]+extent_3[2]+6*extent_3[3])*(MW.Ni/1000);
            liquid_generation_3[3] = ((2*MW.H+MW.S+4*MW.O)/1000)*((-2*extent_3[1]+4*extent_3[3]-8*extent_3[4])/2);

            for i = 1:2
                if (liquid_generation_3[i]<0.0)           
                    display("Problem: liquid metal species consumption!")
                end;
            end;

        O2_consumption_3 = (0.5*extent_3[1]+2*extent_3[2]+15*extent_3[3]+2*extent_3[4]+2*extent_3[5])*(2*MW.O/1000);
        O2_consumption_AC5 = O2_consumption_3-m_11;
        n_O2_consumption_AC5 = O2_consumption_AC5/(2*MW.O);
        n_H2O_condensation_3 = (n_O2_consumption_AC5/y_O2_2nd_stage)*(1-y_O2_2nd_stage);
        H2O_condensation_3 = n_H2O_condensation_3*(2*MW.H+MW.O);

        Q_Reaction_3 = data.h_reaction.h_reaction[1:5].*extent_3; 

        Q_Cond_3 = H2O_condensation_3*(H_Evap_H2O+Cp_H2O*(100-25)+Cp_steam_a*(X[51]-100));

        Cp_11 = Cp_O2_a;
        T_11 = oxygen_T;

        liq_fractions_4_H2SO4 = X[80];
        liq_fractions_4 = liq_fractions(liq_fractions_4_H2SO4, 0.0, 0.0, 0.0)

        solids_in_4_NiS = X[92] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent));
        solids_in_4_Ni3S4 = X[93] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent));
        solids_in_4_Cu9S5 = -10000.0
    
        solids_in_4 = solid_fractions(solids_in_4_NiS, solids_in_4_Ni3S4, solids_in_4_Cu9S5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                     
        Cp_12 = Cp_O2_a;
        T_12 = oxygen_T;

        T_13 = steam_T;

        solid_consumption_4 = zeros(Float64,4)
        extent_4 = zeros(Float64,5)
        liquid_generation_4 = zeros(Float64,3);
        O2_consumption_4 = zero(Float64)
        H2O_condensation_4 = zero(Float64)

        solids_out_4_NiS = X[72] * X[73] / (X[72] / density_residue + X[77] / density_spent)
        solids_out_4_Ni3S4 = X[72] * X[74] / (X[72] / density_residue + X[77] / density_spent)
        solids_out_4_Cu9S5 = X[72] * X[75] / (X[72] / density_residue + X[77] / density_spent)
        solids_out_4_CuS = X[72] * X[76] / (X[72] / density_residue + X[77] / density_spent)

        solids_out_4 = solid_fractions(solids_out_4_NiS, solids_out_4_Ni3S4, solids_out_4_Cu9S5, solids_out_4_CuS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        
        Cu_molal_22 = 1000*(X[78]/(1-X[78]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[79]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[80]))/MW.Cu;
        Ni_molal_22 = 1000*(X[79]/(1-X[78]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[79]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[80]))/MW.Ni;
        H2SO4_molal_22 = 1000*(X[80]/(1-X[78]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[79]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[80]))/(2*MW.H+MW.S+4*MW.O);
        O2_solubility_22 = (2*MW.O/1000)*(calculate_oxygen_solubility2(X[81],P_O2_pascal,Cu_molal_22,Ni_molal_22,H2SO4_molal_22));
        molal_O2_solubility_22 = (O2_solubility_22*1000)/(2*MW.O);

        r_4 = calculate_reaction_rates(X[81],liq_fractions_4,solids_in_4,solids_out_4,molal_O2_solubility_22,Cu9S5_fresh,data);

            extent_4 = zeros(Float64,5);
            extent_4[1] = volume_4*60*r_4[1]
            extent_4[2] = volume_4*60*r_4[2]
            extent_4[3] = volume_4*60*r_4[3]
            extent_4[4] = volume_4*60*r_4[4]
            extent_4[5] = volume_4*60*r_4[5]       

            solid_consumption_4[1] = ((MW.Ni+MW.S)/1000)*(4*extent_4[1]+extent_4[2]);
            solid_consumption_4[2] = ((3*MW.Ni+4*MW.S)/1000)*(-extent_4[1]+2*extent_4[3]);
            solid_consumption_4[3] = ((9*MW.Cu+5*MW.S)/1000)*(extent_4[4]);
            solid_consumption_4[4] = ((MW.Cu+MW.S)/1000)*(-5*extent_4[4]+extent_4[5]);

            liquid_generation_4[1] = (4*extent_4[4]+extent_4[5])*(MW.Cu/1000);
            liquid_generation_4[2] = (extent_4[1]+extent_4[2]+6*extent_4[3])*(MW.Ni/1000);
            liquid_generation_4[3] = ((2*MW.H+MW.S+4*MW.O)/1000)*((-2*extent_4[1]+4*extent_4[3]-8*extent_4[4])/2);
            
        for i = 1:2
            if (liquid_generation_3[i]<0.0)           
                display("Problem: liquid metal species consumption!")
            end;
        end;

        O2_consumption_4 = (0.5*extent_4[1]+2*extent_4[2]+15*extent_4[3]+2*extent_4[4]+2*extent_4[5])*(2*MW.O/1000);
        m_O2_AC5 = m_12-O2_consumption_4;
        n_O2_AC5 = m_O2_AC5/(2*MW.O);
        n_H2O_AC5 = n_O2_AC5/y_O2_2nd_stage*(1-y_O2_2nd_stage);
        H2O_condensation_4 = -n_H2O_AC5*(2*MW.H+MW.O);

        Q_Reaction_4 = data.h_reaction.h_reaction[1:5].*extent_4;
        Q_Cond_4 = H2O_condensation_4*(H_Evap_H2O+Cp_H2O*(100-25)+Cp_steam_a*(X[81]-100));

    solid_fraction_stream_17 = 1

    m_17 = X[53] * m_15 * 1 / solid_fraction_stream_17
    m_16 = m_15 - m_17

    overflow_tank_10 = flow_constants.overflow*real(complex(volume_tank10 - tank_10_max_vol)^(3/2))
    overflow_tank_20 = flow_constants.overflow*real(complex(volume_tank20 - tank_20_max_vol)^(3/2))
    overflow_tank_40 = flow_constants.overflow*real(complex(volume_tank40 - tank_40_max_vol)^(3/2))

    overflow_tank_50 = flow_constants.overflow*real(complex(volume_tank50 - tank_50_max_vol)^(3/2))

    dX_dt = zeros(Float64, length(X))

    # dX_dt[1] = m_S_1+m_spent_1+m_2+m_3-(m_5+overflow_tank_10);
    # dX_dt[2] = (m_S_1 - X[2]*m_5 - X[2]*dX_dt[1])/X[1];
    # dX_dt[3:6] = (m_phases_in[1:4]-X[3:6].*X[2]*m_5-X[3:6].*(X[2]*dX_dt[1]+X[1]*dX_dt[2]))./(X[2]*X[1]);
    # dX_dt[7] = -1*dX_dt[2];
    # dX_dt[8:10] = (m_liquid_in[2:4] - X[7].*X[8:10]*m_5 - X[8:10].*(X[7]*dX_dt[1]+X[1]*dX_dt[7]))./(X[1]*X[7]);
      
    # dX_dt[100] = (m_S_1*Cp_S_1*(slurry_T-25)+m_spent_1*Cp_L_1*(slurry_T-25)+m_2*Cp_L_2*(spent_T-25)+m_3*Cp_L_3*(formic_T-25)- m_5*(X[7]*Cp_L_5+X[2]*Cp_S_5)*(X[100]-25) - (X[7]*Cp_L_5+X[2]*Cp_S_5)*(X[100]-25)*dX_dt[1])/(X[1]*(X[7]*Cp_L_5+X[2]*Cp_S_5))

    # dX_dt[11] = (m_5) + m_9 + m_23 + m_24 - m_H2O_flash_evaporated - (m_7+overflow_tank_20);
    # dX_dt[12] = (X[2]*(m_5)+m_9*X[22]-X[12]*(m_7+overflow_tank_20)-X[12]*dX_dt[11])/(X[11]);
    # dX_dt[13:16] = (((m_5))*X[2].*X[3:6]+m_9*X[22].*X[23:26]-((m_7+overflow_tank_20))*X[12].*X[13:16]-X[13:16].*(X[12]*dX_dt[11]+X[11]*dX_dt[12]))./(X[12]*X[11]);
    # dX_dt[94:95] = (((m_5))*X[2].*X[3:4]+m_9*X[22].*X[82:83]-((m_7+overflow_tank_20))*X[12].*X[94:95]-X[94:95].*(X[12]*dX_dt[11]+X[11]*dX_dt[12]))./(X[12]*X[11]);
    # dX_dt[17] = -1*dX_dt[12];
    # dX_dt[18:19] = (((m_5))*X[7].*X[8:9]+m_9*X[27].*X[28:29]-((m_7+overflow_tank_20))*X[17].*X[18:19]-X[18:19].*(X[17]*dX_dt[11]+X[11]*dX_dt[17]))./(X[17]*X[11]);
    # dX_dt[20] = (((m_5))*X[7].*X[10]+m_H2SO4_23+m_9*X[27].*X[30]-((m_7+overflow_tank_20))*X[17].*X[20]-X[20].*(X[17]*dX_dt[11]+X[11]*dX_dt[17]))./(X[17]*X[11]);
    # dX_dt[21] = (((m_5))*(X[7]*Cp_L_5+X[2]*Cp_S_5)*(X[100]-25)+m_23*Cp_L_23*(acid_T-25)+m_24*Cp_L_24*(water_T-25)+(m_9-m_H2O_flash_evaporated)*(X[27]*Cp_L_9+X[22]*Cp_S_9)*(100-25)-(m_7+overflow_tank_20)*(X[12]*Cp_S_7+X[17]*Cp_L_7)*(X[21]-25) - (X[12]*Cp_S_7+X[17]*Cp_L_7)*(X[21]-25)*dX_dt[11])    /   (X[11]*(X[12]*Cp_S_7+X[17]*Cp_L_7));
  
    # dX_dt[96] = m_7 - m_9 - m_AC1 + O2_consumption_1 + H2O_condensation_1;
    # dX_dt[22] = (X[12]*m_7-m_9*X[22]-X[22]*m_AC1-sum(solid_consumption_1)-X[22]*dX_dt[96])/(X[96]);
    # dX_dt[23:26] = (X[12]*m_7.*X[13:16]-m_9*X[22].*X[23:26]-m_AC1*X[22].*X[23:26]-solid_consumption_1-(X[96]*dX_dt[22]+X[22]*dX_dt[96]).*X[23:26])./(X[22]*X[96]);
    # dX_dt[82:83] = (X[12]*m_7.*X[94:95]-m_9*X[22].*X[82:83]-m_AC1*X[22].*X[82:83]-(X[96]*dX_dt[22]+X[22]*dX_dt[96]).*X[82:83])./(X[22]*X[96]);
    # dX_dt[27] = -1*dX_dt[22];
    # dX_dt[28:30] = (X[17]*m_7.*X[18:20]-m_9*X[27].*X[28:30]-m_AC1*X[27].*X[28:30]+liquid_generation_1-(X[96]*dX_dt[27]+X[27]*dX_dt[96]).*X[28:30])./(X[27]*X[96]);
    # dX_dt[31] = (3*agitator_power - Q_loss*length_1 + ((X[17]*Cp_L_7+X[12]*Cp_S_7)*(m_7)*(X[21]-25)) - ((X[31]-25)*(X[27]*Cp_L_9+X[22]*Cp_S_9)*(m_9+m_AC1)) + Q_Cond_1 + O2_consumption_1*Cp_O2_a*(X[31]-25) - sum(Q_Reaction_1)-(X[31]-25.0)*(X[27]*Cp_L_9+X[22]*Cp_S_9)*dX_dt[96])/(X[96]*(X[27]*Cp_L_9+X[22]*Cp_S_9));

    # dX_dt[97] = m_AC1 - m_AC2 + O2_consumption_2 + H2O_condensation_2
    # dX_dt[32] = (X[22]*m_AC1-X[32]*m_AC2-sum(solid_consumption_2)-X[32]*dX_dt[97])/(X[97]);
    # dX_dt[33:36] = (m_AC1*X[22].*X[23:26]-m_AC2*X[32].*X[33:36]-solid_consumption_2-(X[97]*dX_dt[32]+X[32]*dX_dt[97]).*X[33:36])./(X[97]*X[32]);
    # dX_dt[84:85] = (m_AC1*X[22].*X[82:83]-m_AC2*X[32].*X[84:85]-(X[97]*dX_dt[32]+X[32]*dX_dt[97]).*X[84:85])./(X[97]*X[32]);
    # dX_dt[37] = -1*dX_dt[32];
    # dX_dt[38:40] = (m_AC1*X[27].*X[28:30]-m_AC2*X[37].*X[38:40]+liquid_generation_2-(X[97]*dX_dt[37]+X[37]*dX_dt[97]).*(X[38:40]))./(X[37]*X[97]);
    # dX_dt[41] = (-1*Q_Removed_2-Q_loss*length_2+agitator_power+(m_AC1*(X[31]-25)*(X[22]*Cp_S_9+X[27]*Cp_L_9)+m_10*Cp_10*(T_10-25))-(m_AC2*(X[41]-25)*(X[37]*Cp_L_AC2+X[32]*Cp_S_AC2))-sum(Q_Reaction_2)+Q_Cond_2+(O2_consumption_2-m_10)*Cp_10*(X[41]-25)-(X[41]-25.0)*(X[37]*Cp_L_AC2+X[32]*Cp_S_AC2)*dX_dt[97])/(X[97]*(X[37]*Cp_L_AC2+X[32]*Cp_S_AC2));
    
    # dX_dt[98] = m_AC2 - m_14 + O2_consumption_3 + H2O_condensation_3
    # dX_dt[42] = (X[32]*m_AC2-X[42]*m_14-sum(solid_consumption_3)-X[42]*dX_dt[98])/(X[98]);
    # dX_dt[43:46] = (X[32]*m_AC2.*X[33:36]-m_14*X[42].*X[43:46]-solid_consumption_3-(X[98]*dX_dt[42]+X[42]*dX_dt[98]).*(X[43:46]))./(X[42]*X[98]);
    # dX_dt[86:87] = (X[32]*m_AC2.*X[84:85]-m_14*X[42].*X[86:87]-(X[98]*dX_dt[42]+X[42]*dX_dt[98]).*(X[86:87]))./(X[42]*X[98]);
    # dX_dt[47] = -1*dX_dt[42];
    # dX_dt[48:50] = (X[37]*m_AC2.*X[38:40]-m_14*X[47].*X[48:50]+liquid_generation_3-(X[98]*dX_dt[47]+X[47]*dX_dt[98]).*(X[48:50]))./(X[98]*X[47]);
    # dX_dt[51] = (-1*Q_Removed_3-Q_loss*length_3+agitator_power+(m_AC2*(X[41]-25)*(X[32]*Cp_S_AC2+X[37]*Cp_L_AC2)+m_11*Cp_11*(T_11-25))-(m_14*(X[51]-25)*(X[47]*Cp_L_14+X[42]*Cp_S_14))-sum(Q_Reaction_3)+Q_Cond_3+(O2_consumption_3-m_11)*Cp_11*(X[51]-25)-(X[51]-25.0)*(X[47]*Cp_L_14+X[42]*Cp_S_14)*dX_dt[98])/(X[98]*(X[47]*Cp_L_14+X[42]*Cp_S_14));
    
    # dX_dt[52] = m_14-(m_15+overflow_tank_40);
    # dX_dt[53] = (X[42]*m_14-X[53]*(m_15+overflow_tank_40)-X[53]*dX_dt[52])/(X[52]);
    # dX_dt[54:57] = (m_14*X[42].*X[43:46]-((m_15+overflow_tank_40))*X[53].*X[54:57]-(X[52]*dX_dt[53]+X[53]*dX_dt[52]).*X[54:57])./(X[53]*X[52]);
    # dX_dt[88:89] = (m_14*X[42].*X[86:87]-((m_15+overflow_tank_40))*X[53].*X[88:89]-(X[52]*dX_dt[53]+X[53]*dX_dt[52]).*X[88:89])./(X[53]*X[52]);
    # dX_dt[58] = -1*dX_dt[53];
    # dX_dt[59:61] = (m_14*X[47].*X[48:50]-((m_15+overflow_tank_40))*X[58].*X[59:61]-(X[52]*dX_dt[58]+X[58]*dX_dt[52]).*X[59:61])./(X[58]*X[52]);
    # dX_dt[101] = (-1*Q_Removed_TK40+m_14*(X[51]-25)*(X[47]*Cp_L_14+X[42]*Cp_S_14)-(m_15+overflow_tank_40)*(X[101]-25)*(X[53]*Cp_S_15+X[58]*Cp_L_15)-(X[53]*Cp_S_15+X[58]*Cp_L_15)*(X[101]-25)*dX_dt[52])/((X[53]*Cp_S_15+X[58]*Cp_L_15)*X[52])

    # dX_dt[62] = m_17+m_18+m_19+m_20-(m_21+overflow_tank_50);
    # dX_dt[63] = (solid_fraction_stream_17*m_17 - X[63]*(m_21+overflow_tank_50) - X[63]*dX_dt[62])/X[62];
    # dX_dt[64:67] = (solid_fraction_stream_17*m_17.*X[54:57]-X[64:67].*X[63]*(m_21+overflow_tank_50)-X[64:67].*(X[63]*dX_dt[62]+X[62]*dX_dt[63]))./(X[62]*X[63]);
    # dX_dt[90:91] = (solid_fraction_stream_17*m_17.*X[88:89]-X[90:91].*X[63]*(m_21+overflow_tank_50)-X[90:91].*(X[63]*dX_dt[62]+X[62]*dX_dt[63]))./(X[62]*X[63]);
    # dX_dt[68] = -1*dX_dt[63];
    # dX_dt[69:71] = ((1-solid_fraction_stream_17)*m_17*X[59:61] + m_liquid_TK50[2:4] - X[68].*X[69:71].*(m_21+overflow_tank_50) - X[69:71].*(X[68]*dX_dt[62]+X[62]*dX_dt[68]))./(X[62]*X[68]);
    # dX_dt[102] = (m_18*Cp_L_18*(spent_T-25)+m_19*Cp_L_19*(water_T-25)+m_20*Cp_L_20*(acid_T-25)+m_17*solid_fraction_stream_17*Cp_S_17*(X[101]-25)+m_17*(1-solid_fraction_stream_17)*Cp_L_17*(X[101]-25) - (m_21+overflow_tank_50)*(X[63]*Cp_S_21+X[68]*Cp_L_21)*(X[102]-25)-(X[63]*Cp_S_21+X[68]*Cp_L_21)*(X[102]-25)*dX_dt[62])/((X[63]*Cp_S_21+X[68]*Cp_L_21)*X[62])

    # dX_dt[99] = m_21 + m_13 + O2_consumption_4 + H2O_condensation_4 - m_22
    # dX_dt[72] = (m_21*X[63]-m_22*X[72]-sum(solid_consumption_4)-X[72]*dX_dt[99])/(X[99]);
    # dX_dt[73:76] = (m_21*X[63].*X[64:67]-m_22*X[72].*X[73:76]-solid_consumption_4-(X[99]*dX_dt[72]+X[72]*dX_dt[99]).*(X[73:76]))./(X[99]*X[72]);
    # dX_dt[92:93] = (m_21*X[63].*X[90:91]-m_22*X[72].*X[92:93]-(X[99]*dX_dt[72]+X[72]*dX_dt[99]).*(X[92:93]))./(X[99]*X[72]);
    # dX_dt[77] = -1*dX_dt[72];
    # dX_dt[78:80] = (m_21*X[68].*X[69:71]-m_22*X[77].*X[78:80]+liquid_generation_4-(X[99]*dX_dt[77]+X[77]*dX_dt[99]).*(X[78:80]))./(X[99]*X[77]);
    # dX_dt[81] = (agitator_power-Q_loss*length_4+m_21*(X[68]*Cp_L_21+X[63]*Cp_S_21)*(X[102]-25)+m_12*Cp_12*(T_12-25)+m_13*(Cp_H2O*(100-25)+H_Evap_H2O+Cp_steam_a*(T_13-100))-m_22*(X[77]*Cp_L_22+X[72]*Cp_S_22)*(X[81]-25)-(m_12-O2_consumption_4)*Cp_12*(X[81]-25)-sum(Q_Reaction_4)+Q_Cond_4-(X[81]-25.0)*(X[77]*Cp_L_22+X[72]*Cp_S_22)*dX_dt[99])/(X[99]*(X[77]*Cp_L_22+X[72]*Cp_S_22));

    # dX_dt[103] = m_O2_total / (2*MW.O) - (O2_consumption_1/(2*MW.O) + O2_consumption_2/(2*MW.O) + O2_consumption_3/(2*MW.O) + O2_consumption_4/(2*MW.O)) - y_O2_vap * flow_constants.O2out * sqrt(P_autoclave*100000-P_atm)
    
    # dX_dt[104] = level_error_tank_10
    # dX_dt[105] = level_error_tank_20
    # dX_dt[106] = level_error_AC3
    # dX_dt[107] = level_error_tank_40
    # dX_dt[108] = level_error_tank_50
    # dX_dt[109] = level_error_AC4
    # dX_dt[110] = temperature_error_AC1
    # dX_dt[111] = temperature_error_AC2
    # dX_dt[112] = temperature_error_AC3
    # dX_dt[113] = temperature_error_AC4
    # dX_dt[114] = pressure_error
    # dX_dt[115] = temperature_error_TK40
   
        mass_AC1 = volume_1/(X[22]/density_residue+X[27]/density_spent)
        mass_AC2 = volume_2/(X[32]/density_residue+X[37]/density_spent)
        mass_14 = volume_3/(X[42]/density_residue+X[47]/density_spent)
        mass_22 = volume_4/(X[72]/density_residue+X[77]/density_spent)

        H2SO4_concentration_AC1 = X[30] * density_spent * 1000  #g/l
        H_plus_concentration_AC1 = 2 * H2SO4_concentration_AC1 / (2*MW.H + MW.S + MW.O*4)   #mol/l
        H2SO4_concentration_AC2 = X[40] * density_spent * 1000  #g/l
        H_plus_concentration_AC2 = 2 * H2SO4_concentration_AC2 / (2*MW.H + MW.S + MW.O*4)   #mol/l
        H2SO4_concentration_AC3 = X[50] * density_spent * 1000  #g/l
        H_plus_concentration_AC3 = 2 * H2SO4_concentration_AC3 / (2*MW.H + MW.S + MW.O*4)   #mol/l
        H2SO4_concentration_AC4 = X[80] * density_spent * 1000  #g/l
        H_plus_concentration_AC4 = 2 * H2SO4_concentration_AC4 / (2*MW.H + MW.S + MW.O*4)   #mol/l

        Cu_concentration_AC1 = X[28] * density_spent * 1000 #g/l
        Cu_concentration_AC2 = X[38] * density_spent * 1000 #g/l
        Cu_concentration_AC3 = X[48] * density_spent * 1000 #g/l
        Cu_concentration_AC4 = X[78] * density_spent * 1000 #g/l
        
        NiS_init_concentration_TK10 = X[3] * X[2] * (1/(X[2]/density_residue + X[7]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_FRT = X[94] * X[12] * (1/(X[12]/density_residue + X[17]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_AC1 = X[82] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_AC2 = X[84] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_AC3 = X[86] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_AC4 = X[92] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent))

        NiS_conc_TK10 = X[3] * X[2] * (1/(X[2]/density_residue + X[7]/density_spent)) 
        NiS_conc_FRT = X[13] * X[12] * (1/(X[12]/density_residue + X[17]/density_spent))
        NiS_conc_AC1 = X[23] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent))
        NiS_conc_AC2 = X[33] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent))
        NiS_conc_AC3 = X[43] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent))
        NiS_conc_AC4 = X[73] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent))
        
        Ni3S4_init_concentration_TK10 = X[4] * X[2] * (1/(X[2]/density_residue + X[7]/density_spent)) #kg Ni3S4_init / l
        
        #solids_in_1_Ni3S4_check = X[12] * X[14] * (1/(X[12]/density_residue + X[17]/density_spent))       #Concentration of solids entering the first autoclave compartment, kg Ni3S4/l slurry
        #solids_in_1_Ni3S4 = X[12] * flow_constants.recycle_tank * X[11]^0.5 * X[14] / V_7          #Concentration of solids entering the first autoclave compartment, kg Ni3S4/l slurry
        
        Ni3S4_init_concentration_FRT = X[95] * X[12] * (1/(X[12]/density_residue + X[17]/density_spent)) #kg Ni3S4_init / l
        Ni3S4_init_concentration_AC1 = X[83] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent)) #kg Ni3S4_init / l
        Ni3S4_init_concentration_AC2 = X[85] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent)) #kg Ni3S4_init / l
        Ni3S4_init_concentration_AC3 = X[87] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent)) #kg Ni3S4_init / l
        Ni3S4_init_concentration_AC4 = X[93] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent)) #kg Ni3S4_init / l

        Ni3S4_conc_TK10 = X[4] * X[2] * (1/(X[2]/density_residue + X[7]/density_spent)) 
        Ni3S4_conc_FRT = X[14] * X[12] * (1/(X[12]/density_residue + X[17]/density_spent))
        Ni3S4_conc_AC1 = X[24] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent))
        Ni3S4_conc_AC2 = X[34] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent))
        Ni3S4_conc_AC3 = X[44] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent))
        Ni3S4_conc_AC4 = X[74] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent))  
        
        MM_NiS_1 = (1000*solids_out_1_NiS/(MW.Ni+MW.S))/(0.001+1000*solids_out_1_NiS/(MW.Ni+MW.S))
        MM_Ni3S4_1 = (1000*solids_out_1_Ni3S4/(3*MW.Ni+4*MW.S))/(0.0001+1000*solids_out_1_Ni3S4/(3*MW.Ni+4*MW.S))
        MM_Cu9S5_1 = (1000*solids_out_1_Cu9S5/(9*MW.Cu+5*MW.S))/(0.0001+1000*solids_out_1_Cu9S5/(9*MW.Cu+5*MW.S))
        MM_CuS_1 = (1000*solids_out_1_CuS/(MW.Cu+MW.S))/(0.0001+1000*solids_out_1_CuS/(MW.Cu+MW.S))

        MM_NiS_2 = (1000*solids_out_2_NiS/(MW.Ni+MW.S))/(0.001+1000*solids_out_2_NiS/(MW.Ni+MW.S))
        MM_Ni3S4_2 = (1000*solids_out_2_Ni3S4/(3*MW.Ni+4*MW.S))/(0.0001+1000*solids_out_2_Ni3S4/(3*MW.Ni+4*MW.S))
        MM_Cu9S5_2 = (1000*solids_out_2_Cu9S5/(9*MW.Cu+5*MW.S))/(0.0001+1000*solids_out_2_Cu9S5/(9*MW.Cu+5*MW.S))
        MM_CuS_2 = (1000*solids_out_2_CuS/(MW.Cu+MW.S))/(0.0001+1000*solids_out_2_CuS/(MW.Cu+MW.S))

        MM_NiS_3 = (1000*solids_out_3_NiS/(MW.Ni+MW.S))/(0.001+1000*solids_out_3_NiS/(MW.Ni+MW.S))
        MM_Ni3S4_3 = (1000*solids_out_3_Ni3S4/(3*MW.Ni+4*MW.S))/(0.0001+1000*solids_out_3_Ni3S4/(3*MW.Ni+4*MW.S))
        MM_Cu9S5_3 = (1000*solids_out_3_Cu9S5/(9*MW.Cu+5*MW.S))/(0.0001+1000*solids_out_3_Cu9S5/(9*MW.Cu+5*MW.S))
        MM_CuS_3 = (1000*solids_out_3_CuS/(MW.Cu+MW.S))/(0.0001+1000*solids_out_3_CuS/(MW.Cu+MW.S))

        MM_NiS_4 = (1000*solids_out_4_NiS/(MW.Ni+MW.S))/(0.001+1000*solids_out_4_NiS/(MW.Ni+MW.S))
        MM_Ni3S4_4 = (1000*solids_out_4_Ni3S4/(3*MW.Ni+4*MW.S))/(0.0001+1000*solids_out_4_Ni3S4/(3*MW.Ni+4*MW.S))
        MM_Cu9S5_4 = (1000*solids_out_4_Cu9S5/(9*MW.Cu+5*MW.S))/(0.0001+1000*solids_out_4_Cu9S5/(9*MW.Cu+5*MW.S))
        MM_CuS_4 = (1000*solids_out_4_CuS/(MW.Cu+MW.S))/(0.0001+1000*solids_out_4_CuS/(MW.Cu+MW.S))

        V_16 = m_16 / density_spent #it is assumed stream 16 contains no solids
        V_17 = solid_fraction_stream_17*m_17/density_residue + (1-solid_fraction_stream_17)*m_17/density_spent

        wt_percent_copper_initial = 100*(X[5]*MW.Cu*9/(MW.Cu*9+MW.S*5) + X[6]*MW.Cu/(MW.Cu+MW.S))
        wt_percent_Ni_initial = 100*(X[3]*MW.Ni/(MW.Ni+MW.S) + X[4]*MW.Ni*3/(MW.Ni*3+MW.S*4))
        wt_percent_S_initial = 100*(X[5]*MW.S*5/(MW.Cu*9+MW.S*5) + X[6]*MW.S/(MW.Cu+MW.S) + X[3]*MW.S/(MW.Ni+MW.S) + X[4]*MW.S*4/(MW.Ni*3+MW.S*4))

        wt_percent_copper_second_stage_leach = 100*(X[45]*MW.Cu*9/(MW.Cu*9+MW.S*5) + X[46]*MW.Cu/(MW.Cu+MW.S))
        wt_percent_Ni_second_stage_leach = 100*(X[43]*MW.Ni/(MW.Ni+MW.S) + X[44]*MW.Ni*3/(MW.Ni*3+MW.S*4))
        wt_percent_S_second_stage_leach = 100*(X[45]*MW.S*5/(MW.Cu*9+MW.S*5) + X[46]*MW.S/(MW.Cu+MW.S) + X[43]*MW.S/(MW.Ni+MW.S) + X[44]*MW.S*4/(MW.Ni*3+MW.S*4))

        wt_percent_copper_third_stage_leach = 100*(X[75]*MW.Cu*9/(MW.Cu*9+MW.S*5) + X[76]*MW.Cu/(MW.Cu+MW.S))
        wt_percent_Ni_third_stage_leach = 100*(X[73]*MW.Ni/(MW.Ni+MW.S) + X[74]*MW.Ni*3/(MW.Ni*3+MW.S*4))
        wt_percent_S_third_stage_leach = 100*(X[75]*MW.S*5/(MW.Cu*9+MW.S*5) + X[76]*MW.S/(MW.Cu+MW.S) + X[73]*MW.S/(MW.Ni+MW.S) + X[74]*MW.S*4/(MW.Ni*3+MW.S*4))
        
        density_stream_5 =  1 / (X[2] / density_residue + X[7] / density_spent)

        y_vect = zeros(54)

        y_vect[1] = V_1
        y_vect[2] = V_2
        y_vect[3] = V_3
        y_vect[4] = X[100]
        y_vect[5] = tank_10_percent_full
        y_vect[6] = V_5
        y_vect[7] = density_stream_5
        y_vect[8] = V_23
        y_vect[9] = V_7
        y_vect[10] = X[21]
        y_vect[11] = tank_20_percent_full
        y_vect[12] = V_9
        y_vect[13] = X[31]
        y_vect[14] = comp_1_percent_full
        y_vect[15] = X[41]
        y_vect[16] = comp_3_percent_full
        y_vect[17] = X[51]
        y_vect[18] = m_10
        y_vect[19] = m_11
        y_vect[20] = m_12
        y_vect[21] = P_autoclave
        y_vect[22] = T_vap
        y_vect[23] = V_14
        y_vect[24] = X[101]
        y_vect[25] = tank_40_percent_full
        y_vect[26] = V_15
        y_vect[27] = V_17
        y_vect[28] = V_18
        y_vect[29] = V_19
        y_vect[30] = V_20
        y_vect[31] = tank_50_percent_full
        y_vect[32] = X[102]
        y_vect[33] = V_21
        y_vect[34] = X[81]
        y_vect[35] = comp_4_percent_full
        y_vect[36] = V_22
        y_vect[37] = s1_solid_fraction
        y_vect[38] = V_H2O_AC2
        y_vect[39] = scaling_factor
        y_vect[40] = slurry_T
        y_vect[41] = spent_T
        y_vect[42] = comp_3_percent_full_gt
        y_vect[43] = level_AC3
        y_vect[44] = level_measured_AC3
        y_vect[45] = H2SO4_concentration_AC1
        y_vect[46] = H2SO4_concentration_AC2
        y_vect[47] = H2SO4_concentration_AC3
        y_vect[48] = H2SO4_concentration_AC4
        y_vect[49] = wt_percent_copper_initial
        y_vect[50] = wt_percent_copper_second_stage_leach
        y_vect[51] = wt_percent_copper_third_stage_leach
        y_vect[52] = density_stream_5
        y_vect[53] = V_H2O_AC3
        y_vect[54] = m_13

            # A1 = H2SO4_concentration_AC1,
            # A2 = H2SO4_concentration_AC2,
            # A3 = H2SO4_concentration_AC3,
            # A4 = H2SO4_concentration_AC4,
            # wt_percent_copper_initial = wt_percent_copper_initial,
            # wt_percent_Ni_initial = wt_percent_Ni_initial,
            # wt_percent_S_initial = wt_percent_S_initial,
            # wt_percent_copper_second_stage_leach = wt_percent_copper_second_stage_leach,
            # wt_percent_Ni_second_stage_leach = wt_percent_Ni_second_stage_leach, 
            # wt_percent_S_second_stage_leach = wt_percent_S_second_stage_leach,
            # wt_percent_copper_third_stage_leach = wt_percent_copper_third_stage_leach,
            # wt_percent_Ni_third_stage_leach = wt_percent_Ni_third_stage_leach,
            # wt_percent_S_third_stage_leach = wt_percent_S_third_stage_leach,
        
    return y_vect
end

function BMR_gt_measurement_prediction_with_samples(X, data, t)
    
    faultmode = data.faultmode

    residue_solid_fraction = zeros(3)
            residue_solid_fraction[1] = data.auto_regressive_inputs.copper_linear_interp(t)
            residue_solid_fraction[2] = data.auto_regressive_inputs.ni_linear_interp(t)
    residue_solid_fraction[3] = data.residue.S/100

    Cu_as_Cu9S5 = data.residue.Cu_as_Cu9S5;
    Ni_as_NiS = data.residue.Ni_as_NiS;
    Rh_as_Rh2S3 = data.residue.Rh_as_Rh2S3;
    Ru_as_RuS2 = data.residue.Ru_as_RuS2;
    Ir_as_Ir2S3 = data.residue.Ir_as_Ir2S3;
    Fe_as_FeOHSO4 = data.residue.Fe_as_FeOHSO4;
    CuS_leachable = data.residue.CuS_leachable;

    spent_concentration = zeros(7)
    spent_concentration[1] = data.spent_e.Cu
    spent_concentration[2] = data.spent_e.Ni
    spent_concentration[3] = data.spent_e.H2SO4

    density_acid = data.acid.density
    density_spent = data.spent_e.density
    density_residue = data.residue.density
    acid_concentration = data.acid.concentration
    MW = data.mol_W
    Cp_H2SO4_a = 1000*data.cP.H2SO4_a / (2 * MW.H + MW.S + 4 * MW.O)
    Cp_H2SO4_b = 1000*data.cP.H2SO4_b / (2 * MW.H + MW.S + 4 * MW.O)
    Cp_H2O = 1000*data.cP.H2O / (2 * MW.H + MW.O)
    Cp_steam_a = 1000*data.cP.Steam_a/(2*MW.H+MW.O);
    Cp_steam_b = 1000*data.cP.Steam_b/(2*MW.H+MW.O);
    Cp_steam_c = 1000*data.cP.Steam_c/(2*MW.H+MW.O);
    Cp_steam_d = 1000*data.cP.Steam_d/(2*MW.H+MW.O);
    Cp_O2_a = 1000*data.cP.Oxygen_a/(2*MW.O);
    Cp_O2_b = 1000*data.cP.Oxygen_b/(2*MW.O);
    Cp_O2_c = 1000*data.cP.Oxygen_c/(2*MW.O);
    Cp_O2_d = 1000*data.cP.Oxygen_d/(2*MW.O);
    Cp_NiS = data.cP.NiS/(MW.Ni+MW.S);
    Cp_CuS = data.cP.CuS/(MW.Cu+MW.S);
    Cp_Cu9S5 = data.cP.Cu9S5/(9*MW.Cu+5*MW.S);
    Cp_Ni3S4_a = data.cP.Ni3S4_a/(3*MW.Ni+4*MW.S);
    Cp_Ni3S4_b = data.cP.Ni3S4_b/(3*MW.Ni+4*MW.S);
    Cp_Ni3S4_c = data.cP.Ni3S4_c/(3*MW.Ni+4*MW.S);
    Cp_Ni3S4_d = data.cP.Ni3S4_d/(3*MW.Ni+4*MW.S);
    Cp_Ni3S4_e = data.cP.Ni3S4_e/(3*MW.Ni+4*MW.S);

    H_Evap_H2O = data.gen_constants.H_Evap_H2O
    h_reaction = data.h_reaction

    s1_solid_fraction_mean = data.flow.s1_solid_fraction_mean
    s3_mass_fraction_H2SO4 = data.flow.s3_mass_fraction_H2SO4
    T_5 = data.flow.s5_temperature
    T_21 = data.flow.stream_21_temp
    fraction_oxygen_stream_10 = data.flow.oxygen_fraction_10
    fraction_oxygen_stream_11 = data.flow.oxygen_fraction_11

        oxygen_T = data.flow.oxygen_temp
        steam_T = data.flow.steam_temp
        formic_T = data.auto_regressive_inputs.stream_3_linear_interp(t)
        water_T = data.auto_regressive_inputs.stream_23_linear_interp(t)
        spent_T = data.auto_regressive_inputs.stream_2_linear_interp(t)
        slurry_T = data.auto_regressive_inputs.stream_1_linear_interp(t)
        acid_T = data.auto_regressive_inputs.stream_23_linear_interp(t)
        s1_solid_fraction = data.auto_regressive_inputs.stream_1_solid_frac_interp(t)

    agitator_power = 3600*data.autoclave.agitator_power
    autoclave_outer_diameter = data.autoclave.outer_diameter
    autoclave_surface_temperature = data.autoclave.outer_temperature
    autoclave_emissivity = data.autoclave.emissivity
    length_1 = data.autoclave.comp_1_length
    length_2 = data.autoclave.comp_2_length
    length_3 = data.autoclave.comp_3_length
    length_4 = data.autoclave.comp_4_length
    weir_height_AC1 = data.autoclave.weir_height_AC1
    weir_height_AC2 = data.autoclave.weir_height_AC2
    comp_width = data.autoclave.comp_width
    P_atm = data.autoclave.atmospheric_pressure
    volume_total = 1000*data.autoclave.total_volume

    tank_10_max_vol = data.autoclave.tank_10_max_vol * 1000
    tank_20_max_vol = data.autoclave.tank_20_max_vol * 1000
    tank_40_max_vol = data.autoclave.tank_40_max_vol * 1000
    tank_50_max_vol = data.autoclave.tank_50_max_vol * 1000

    AC1_max_vol = data.autoclave.comp_1_max_vol * 1000
    AC2_max_vol = data.autoclave.comp_2_max_vol * 1000
    AC3_max_vol = data.autoclave.comp_3_max_vol * 1000
    AC4_max_vol = data.autoclave.comp_4_max_vol * 1000

    tank_10_area = data.autoclave.tank_10_area
    tank_20_area = data.autoclave.tank_20_area
    tank_40_area = data.autoclave.tank_40_area
    tank_50_area = data.autoclave.tank_50_area
    AC3_area = data.autoclave.AC3_area
    AC4_area = data.autoclave.AC4_area

    volume_tank10 = X[2]*X[1]/density_residue + X[7]*X[1]/density_spent
    volume_tank20 = X[12]*X[11]/density_residue + X[17]*X[11]/density_spent
    volume_tank40 = X[53]*X[52]/density_residue + X[58]*X[52]/density_spent
    volume_tank50 = X[63]*X[62]/density_residue + X[68]*X[62]/density_spent

    volume_1 = X[22] * X[96] / density_residue + X[27] * X[96] / density_spent
    volume_2 = X[32] * X[97] / density_residue + X[37] * X[97] / density_spent
    volume_3 = X[42] * X[98] / density_residue + X[47] * X[98] / density_spent
    volume_4 = X[72] * X[99] / density_residue + X[77] * X[99] / density_spent
    
    slurry_height_1 = (volume_1/1000)/(length_1*comp_width)
    slurry_height_2 = (volume_2/1000)/(length_2*comp_width)

    tank_10_percent_full = volume_tank10 / tank_10_max_vol * 100
    tank_20_percent_full = volume_tank20 / tank_20_max_vol * 100
    tank_40_percent_full = volume_tank40 / tank_40_max_vol * 100
    tank_50_percent_full = volume_tank50 / tank_50_max_vol * 100
    comp_1_percent_full = volume_1 / AC1_max_vol * 100
    comp_2_percent_full = volume_2 / AC2_max_vol * 100
    comp_3_percent_full = volume_3 / AC3_max_vol * 100
    comp_4_percent_full = volume_4 / AC4_max_vol * 100

    level_tank10 = (volume_tank10/1000) / tank_10_area
    level_tank20 = (volume_tank20/1000) / tank_20_area
    level_tank40 = (volume_tank40/1000) / tank_40_area
    level_tank50 = (volume_tank50/1000) / tank_50_area
    level_AC3 = (volume_3/1000) / AC3_area
    level_AC4 = (volume_4/1000) / AC4_area

    V_vap = (volume_total - volume_1 - volume_2 - volume_3 - volume_4) / 1000
    T_vap = 1/(length_1+length_2+length_3+length_4)*(length_1*X[31]+length_2*X[41]+length_3*X[51]+length_4*X[81])

    if T_vap < 60
        P_vap_H2O = ((10 ^ (8.10765 - 1750.286 / (T_vap + 235))) / 760) * 101325
    else
        P_vap_H2O = ((10 ^ (7.96681 - 1668.210 / (T_vap + 228))) / 760) * 101325
    end

    P_O2_pascal = X[103] * 1000 * 8.314 * (T_vap + 273.15) / V_vap

    y_O2_vap    = P_O2_pascal / (P_O2_pascal + P_vap_H2O)
    P_autoclave = (P_vap_H2O + P_O2_pascal) / 100000

    air = data.air
    flow_constants = data.flow_constants

            setpoint = data.control.setpoint
            bias = data.control.bias
            Kc = data.control.Kc
            Ti = data.control.Ti

                V_5 = setpoint[12]
                V_23 = setpoint[13]
                V_24 = setpoint[14]
                V_18 = setpoint[15]
                V_20 = setpoint[16]
                V_21 = setpoint[17]

                level_error_tank_10 = level_tank10 - setpoint[1]
                V_1 = - Kc[1] * (level_error_tank_10 + 1/Ti[1] * X[104]) + bias[1]
                V_1 = max(0.0, V_1)
                if V_1 == 0.0
                    level_error_tank_10 = 0.0
                end

                level_error_tank_20 = level_tank20 - setpoint[2]
                V_7 = -Kc[2] * (level_error_tank_20 + 1/Ti[2] * X[105]) + bias[2]
                V_7 = max(0.0, V_7)
                if V_7 == 0.0
                    level_error_tank_20 = 0.0
                end

                if faultmode == 8
                    level_measured_AC3 = level_AC3 + 0.3
                    comp_3_percent_full = (level_measured_AC3 * 1000 * AC3_area) / AC3_max_vol * 100
                else
                    level_measured_AC3 = level_AC3
                end
                level_error_AC3 = level_measured_AC3 - setpoint[3]
                V_14 = -Kc[3] * (level_error_AC3 + 1/Ti[3] * X[106]) + bias[3]
                V_14 = max(0.0, V_14)
                if V_14 == 0.0
                    level_error_AC3 = 0.0
                end

                level_error_tank_40 = level_tank40 - setpoint[4]
                V_15 = -Kc[4] * (level_error_tank_40 + 1/Ti[4] * X[107]) + bias[4]
                V_15 = max(0.0, V_15)
                if V_15 == 0.0
                    level_error_tank_40 = 0.0
                end

                level_error_tank_50 = level_tank50 - setpoint[5]
                V_19 = -Kc[5] * (level_error_tank_50 + 1/Ti[5] * X[108]) + bias[5]
                V_19 = max(0.0, V_19)
                if V_19 == 0.0
                    level_error_tank_50 = 0.0
                end

                level_error_AC4 = level_AC4 - setpoint[6]
                V_22 = -Kc[6] * (level_error_AC4 + 1/Ti[6] * X[109]) + bias[6]
                V_22 = max(0.0, V_22)
                if V_22 == 0.0
                    level_error_AC4 = 0.0
                end

                temperature_error_AC1 = X[31] - setpoint[7]
                V_9 = -Kc[7] * (temperature_error_AC1 + 1/Ti[7] * X[110]) + bias[7]
                V_9 = max(0.0, V_9)
                if V_9 == 0.0
                    temperature_error_AC1 = 0.0
                end

                temperature_error_AC2 = X[41] - setpoint[8]
                V_H2O_AC2 = -Kc[8] * (temperature_error_AC2 + 1/Ti[8] * X[111]) + bias[8]
                V_H2O_AC2 = max(0.0, V_H2O_AC2)
                if V_H2O_AC2 == 0.0
                    temperature_error_AC2 = 0.0
                end
                V_H2O_AC2_max = data.flow.AC2_cooling_water_max
                if V_H2O_AC2 > V_H2O_AC2_max
                    V_H2O_AC2 = V_H2O_AC2_max
                    temperature_error_AC2 = 0.0
                end

                if faultmode == 5
                    scaling_factor = data.flow.scaling_factor_initial + data.flow.scaling_factor_rate_change * t
                else
                    scaling_factor = 1.0
                end

                Q_Removed_2 = scaling_factor * V_H2O_AC2 * 1.0 * Cp_H2O * data.q_removal.delta_T_AC2

                temperature_error_AC3 = X[51] - setpoint[9]
                V_H2O_AC3 = -Kc[9] * (temperature_error_AC3 + 1/Ti[9] * X[112]) + bias[9]
                V_H2O_AC3 = max(0.0, V_H2O_AC3)
                if V_H2O_AC3 == 0.0
                    temperature_error_AC3 = 0.0
                end

                Q_Removed_3 = V_H2O_AC3 * 1.0 * Cp_H2O * data.q_removal.delta_T_AC3
                
                temperature_error_TK40 = X[101] - setpoint[20]
                V_H2O_TK40 = - Kc[12] * (temperature_error_TK40 + 1/Ti[12] * X[115]) + bias[12]
                
                V_H2O_TK40 = max(0.0, V_H2O_TK40)
                if V_H2O_TK40 == 0.0
                    temperature_error_TK40 = 0.0
                end
                
                Q_Removed_TK40 = V_H2O_TK40 * 1.0 * Cp_H2O * data.q_removal.delta_T_TK40

                temperature_error_AC4 = X[81] - setpoint[10]
                m_13 = -Kc[10] * (temperature_error_AC4 + 1/Ti[10] * X[113]) + bias[10]
                m_13 = max(0.0, m_13)
                if m_13 == 0.0
                    temperature_error_AC4 = 0.0
                end

                pressure_error = P_autoclave - setpoint[11]
                m_O2_total = -Kc[11] * (pressure_error + 1/Ti[11] * X[114]) + bias[11]
                m_O2_total = max(0.0, m_O2_total)
                if m_O2_total == 0.0
                    pressure_error = 0.0
                end                

                s1_spent_volume_fraction = ((1-s1_solid_fraction)/data.spent_e.density)/((s1_solid_fraction/density_residue) + ((1-s1_solid_fraction)/data.spent_e.density))

                V_spent_1 = V_1*s1_spent_volume_fraction
                m_spent_1 = V_spent_1 * data.spent_e.density
                m_S_1 = (V_1 - V_spent_1)*density_residue

                density_stream_1 = 1 / (s1_solid_fraction_mean / density_residue + (1-s1_solid_fraction_mean)/data.spent_e.density)

                density_stream_2 = data.spent_e.density
                density_stream_3 = 1.0
                solid_percentage = setpoint[18]
                spent_percentage = setpoint[19]

                V_2 = ((density_stream_1*V_1*(1/density_stream_1 - 1/data.spent_e.density)/(1/density_residue - 1/data.spent_e.density))*100 - solid_percentage * density_stream_1*V_1)/(solid_percentage*(data.spent_e.density + density_stream_3 * (100/spent_percentage-1)))
                V_2 = max(0.0, V_2)

                V_3 = (100 / spent_percentage - 1) * V_2             
                V_3 = max(0.0, V_3)

        if faultmode == 1
           display("valve blockage not modelled yet")
        end
        if faultmode == 2
            display("valve wear fault not modelled")
        end
        if faultmode == 3
            display("valve stiction fault not modelled")
        end
        if faultmode == 4
            max_V_7 = 11000
            if V_7 > max_V_7
                V_7 = max_V_7
            end
        end
        # if faultmode == 5
        #     max_V_H2O_AC2 = 6500
        #     if V_H2O_AC2 > max_V_H2O_AC2
        #         V_H2O_AC2 = max_V_H2O_AC2
        #         Q_Removed_2 = V_H2O_AC2 * 1.0 * Cp_H2O * data.q_removal.delta_T_AC2
        #         temperature_error_AC2 = 0.0
        #     end
        # end
        if faultmode == 6
            m_S_1 = 0.0
            V_spent_1 = 0.0
            V_1 = 0.0
        end
        if faultmode == 7
            display("sulfuric acid controller misuse must be modelled in runODE's function")
        end
        if faultmode == 101
            V_3_leaksize = 100
            if V_3 > 100
                V_3 = V_3 - 100
            else
                V_3 = V_3
            end
        end
        if faultmode == 102
            V_2_maximum = 100
            if V_2 > V_2_maximum
                V_2 = V_2_maximum
            else
                V_2 = V_2
            end
        end
        if faultmode == 103
            m_S_1 = m_S_1 * 0.9
            V_spent_1 = V_spent_1 * 0.9
            V_1 = V_1 * 0.9
        end

        m_2 = V_2 * data.spent_e.density
        m_3 = V_3
        m_H2SO4_4 = m_3 * s3_mass_fraction_H2SO4

        m_5 = V_5 / (X[2] / density_residue + X[7] / density_spent)
        m_7 = V_7 / (X[12] / density_residue + X[17] / density_spent)
        m_9 = V_9 / (X[22] / density_residue + X[27] / density_spent)
        m_14 = V_14 / (X[42] / density_residue + X[47] / density_spent)
        m_15 = V_15 / (X[53] / density_residue + X[58] / density_spent)
        m_21 = V_21 / (X[63] / density_residue + X[68] / density_spent)
        m_22 = V_22 / (X[72] / density_residue + X[77] / density_spent)
       
        m_liquid_in = zeros(Float64, 4)
        m_liquid_in[1] = m_2 + m_3 + m_spent_1
        m_liquid_in[2:3] = (V_spent_1 + V_2) .* spent_concentration[1:2] / 1000
        m_liquid_in[4] = m_H2SO4_4 + (V_spent_1 + V_2) * spent_concentration[3]/1000

        m_solids_in = zeros(Float64, 4)
        m_solids_in[1] = m_S_1
        m_solids_in[2:4] = m_solids_in[1] .* residue_solid_fraction

        n_phases_in = zeros(Float64, 4)
        n_phases_in[1] = (m_solids_in[3]/MW.Ni)*Ni_as_NiS;
        n_phases_in[2] = ((m_solids_in[3]/MW.Ni)*(1-Ni_as_NiS))/3;
        n_phases_in[3] = ((m_solids_in[2]/MW.Cu)*Cu_as_Cu9S5)/9;
        n_phases_in[4] = (m_solids_in[2]/MW.Cu)*(1-Cu_as_Cu9S5);

        m_phases_in = zeros(Float64, 14)
        m_phases_in[1] = n_phases_in[1]*(MW.Ni+MW.S);
        m_phases_in[2] = n_phases_in[2]*(3*MW.Ni+4*MW.S);
        m_phases_in[3] = n_phases_in[3]*(9*MW.Cu+5*MW.S);
        m_phases_in[4] = n_phases_in[4]*(MW.Cu+MW.S);

        m_23 = V_23 * data.acid.density
        m_H2SO4_23 = data.acid.concentration/100 * m_23
        m_24 = V_24

        m_18 = V_18 * density_spent
        m_19 = V_19
        m_20 = V_20 * density_acid
        m_H2SO4_20 = (acid_concentration/100)*m_20

        m_liquid_TK50 = zeros(Float64, 8)
        m_liquid_TK50[1] = m_18 + m_19 + m_20
        m_liquid_TK50[2:3] = V_18 * spent_concentration[1:2]/1000
        m_liquid_TK50[4] = m_H2SO4_20 + V_18 * spent_concentration[3]/1000

        if T_vap < 60
            P_water_2nd_stage = ((10^(8.10765-1750.286/(T_vap+235)))/760)*1.01325
            else
            P_water_2nd_stage = ((10^(7.96681-1668.210/(T_vap+228)))/760)*1.01325
        end

    y_O2_2nd_stage = 1-P_water_2nd_stage/P_autoclave;

        Ra = (9.81*air.thermal_expansion*(autoclave_surface_temperature-air.temp)*(autoclave_outer_diameter^3))/(air.viscosity*air.thermal_diff)
        Nu = (0.6+(0.387*Ra^(1/6))/((1+(0.559/air.Pr)^(9/16))^(8/27)))^2
        h = air.k*Nu/autoclave_outer_diameter
        Q_loss = (h*22/7*autoclave_outer_diameter*(autoclave_surface_temperature-air.temp)+autoclave_emissivity*22/7*autoclave_outer_diameter*5.67*(10^-8)*((autoclave_surface_temperature+273.15)^4-(air.temp+273.15)^4))*3600/1000

        m_10 = m_O2_total * fraction_oxygen_stream_10;
        m_11 = m_O2_total * fraction_oxygen_stream_11;
        m_12 = m_O2_total - m_11 - m_10;

        spent_concentration[3]
        density_spent
        acid_concentration

        X_H2SO4_spent = spent_concentration[3] / (density_spent * 1000)
        X_H2SO4_acid = acid_concentration / 100

        Cp_L_1 = X_H2SO4_spent * Cp_H2SO4_a + (1-X_H2SO4_spent) * Cp_H2O
        Cp_L_2 = X_H2SO4_spent * Cp_H2SO4_a + (1-X_H2SO4_spent) * Cp_H2O
        Cp_L_3 = s3_mass_fraction_H2SO4 * Cp_H2SO4_a + (1-s3_mass_fraction_H2SO4) * Cp_H2O
        Cp_L_5 = X[10]*Cp_H2SO4_a + (1-X[10])*Cp_H2O;
        Cp_L_23 = X_H2SO4_acid * Cp_H2SO4_a + (1 - X_H2SO4_acid) * Cp_H2O
        Cp_L_24 = Cp_H2O
        Cp_L_7 = X[20]*Cp_H2SO4_a+(1-X[20])*Cp_H2O;
        Cp_L_9 = X[30]*Cp_H2SO4_a+(1-X[30])*Cp_H2O;
        Cp_L_AC1 = Cp_L_9;
        Cp_L_AC2 = X[40]*Cp_H2SO4_a+(1-X[40])*Cp_H2O;
        Cp_L_21 = X[71]*Cp_H2SO4_a+(1-X[71])*Cp_H2O;
        Cp_L_22 = X[80]*Cp_H2SO4_a+(1-X[80])*Cp_H2O;
        Cp_L_14 = X[50]*Cp_H2SO4_a+(1-X[50])*Cp_H2O;
        Cp_L_15 = X[61]*Cp_H2SO4_a+(1-X[61])*Cp_H2O;
        Cp_L_16 = Cp_L_15
        Cp_L_17 = Cp_L_15#Cp_H2O                                # Stream 17 (unused: multiplied by zero in dX_dt[102])
        Cp_L_18 = X_H2SO4_spent * Cp_H2SO4_a + (1-X_H2SO4_spent) * Cp_H2O
        Cp_L_19 = Cp_H2O
        Cp_L_20 = X_H2SO4_acid * Cp_H2SO4_a + (1 - X_H2SO4_acid) * Cp_H2O

    Normalx_NiS_S_1 = m_phases_in[1] / (m_phases_in[1]+m_phases_in[2]+m_phases_in[3]+m_phases_in[4])
    Normalx_Ni3S4_S_1 = m_phases_in[2] / (m_phases_in[1]+m_phases_in[2]+m_phases_in[3]+m_phases_in[4])
    Normalx_Cu9S5_S_1 = m_phases_in[3] / (m_phases_in[1]+m_phases_in[2]+m_phases_in[3]+m_phases_in[4])
    Normalx_CuS_S_1 = m_phases_in[4] / (m_phases_in[1]+m_phases_in[2]+m_phases_in[3]+m_phases_in[4])

    Normalx_NiS_S_5 = X[3] / (X[3] + X[4] + X[5] + X[6])
    Normalx_Ni3S4_S_5 = X[4] / (X[3] + X[4] + X[5] + X[6])
    Normalx_Cu9S5_S_5 = X[5] / (X[3] + X[4] + X[5] + X[6])
    Normalx_CuS_S_5 = X[6] / (X[3] + X[4] + X[5] + X[6])

    Normalx_NiS_S_7 = X[13] / (X[13] + X[14] + X[15] + X[16])
    Normalx_Ni3S4_S_7 = X[14] / (X[13] + X[14] + X[15] + X[16])
    Normalx_Cu9S5_S_7 = X[15] / (X[13] + X[14] + X[15] + X[16])
    Normalx_CuS_S_7 = X[16] / (X[13] + X[14] + X[15] + X[16])

    Normalx_NiS_S_9 = X[23] / (X[23] + X[24] + X[25] + X[26])
    Normalx_Ni3S4_S_9 = X[24] / (X[23] + X[24] + X[25] + X[26])
    Normalx_Cu9S5_S_9 = X[25] / (X[23] + X[24] + X[25] + X[26])
    Normalx_CuS_S_9 = X[26] / (X[23] + X[24] + X[25] + X[26])

    Normalx_NiS_S_AC2 = X[33] / (X[33] + X[34] + X[35] + X[36])
    Normalx_Ni3S4_S_AC2 = X[34] / (X[33] + X[34] + X[35] + X[36])
    Normalx_Cu9S5_S_AC2 = X[35] / (X[33] + X[34] + X[35] + X[36])
    Normalx_CuS_S_AC2 = X[36] / (X[33] + X[34] + X[35] + X[36])

    Normalx_NiS_S_21 = X[64] / (X[64] + X[65] + X[66] + X[67])
    Normalx_Ni3S4_S_21 = X[65] / (X[64] + X[65] + X[66] + X[67])
    Normalx_Cu9S5_S_21 = X[66] / (X[64] + X[65] + X[66] + X[67])
    Normalx_CuS_S_21 = X[67] / (X[64] + X[65] + X[66] + X[67])

    Normalx_NiS_S_22 = X[73] / (X[73] + X[74] + X[75] + X[76])
    Normalx_Ni3S4_S_22 = X[74] / (X[73] + X[74] + X[75] + X[76])
    Normalx_Cu9S5_S_22 = X[75] / (X[73] + X[74] + X[75] + X[76])
    Normalx_CuS_S_22 = X[76] / (X[73] + X[74] + X[75] + X[76])

    Normalx_NiS_S_14 = X[43] / (X[43] + X[44] + X[45] + X[46])
    Normalx_Ni3S4_S_14 = X[44] / (X[43] + X[44] + X[45] + X[46])
    Normalx_Cu9S5_S_14 = X[45] / (X[43] + X[44] + X[45] + X[46])
    Normalx_CuS_S_14 = X[46] / (X[43] + X[44] + X[45] + X[46])

    Normalx_NiS_S_15 = X[54] / (X[54] + X[55] + X[56] + X[57])
    Normalx_Ni3S4_S_15 = X[55] / (X[54] + X[55] + X[56] + X[57])
    Normalx_Cu9S5_S_15 = X[56] / (X[54] + X[55] + X[56] + X[57])
    Normalx_CuS_S_15 = X[57] / (X[54] + X[55] + X[56] + X[57])


    Normalx_NiS_S_21 = X[64] / (X[64] + X[65] + X[66] + X[67])
    Normalx_Ni3S4_S_21 = X[65] / (X[64] + X[65] + X[66] + X[67])
    Normalx_Cu9S5_S_21 = X[66] / (X[64] + X[65] + X[66] + X[67])
    Normalx_CuS_S_21 = X[67] / (X[64] + X[65] + X[66] + X[67])
    
    Cp_S_1 = Normalx_NiS_S_1*Cp_NiS + Normalx_Cu9S5_S_1*Cp_Cu9S5 + Normalx_CuS_S_1*Cp_CuS + Normalx_Ni3S4_S_1*Cp_Ni3S4_a;
    Cp_S_5 = Normalx_NiS_S_5*Cp_NiS + Normalx_Cu9S5_S_5*Cp_Cu9S5 + Normalx_CuS_S_5*Cp_CuS + Normalx_Ni3S4_S_5*Cp_Ni3S4_a;
    Cp_S_7 = Normalx_NiS_S_7*Cp_NiS + Normalx_Cu9S5_S_7*Cp_Cu9S5 + Normalx_CuS_S_7*Cp_CuS + Normalx_Ni3S4_S_7*Cp_Ni3S4_a;
    Cp_S_9 = Normalx_NiS_S_9*Cp_NiS + Normalx_Cu9S5_S_9*Cp_Cu9S5 + Normalx_CuS_S_9*Cp_CuS + Normalx_Ni3S4_S_9*Cp_Ni3S4_a;
    Cp_S_AC1 = Cp_S_9;
    Cp_S_AC2 = Normalx_NiS_S_AC2*Cp_NiS + Normalx_Cu9S5_S_AC2*Cp_Cu9S5 + Normalx_CuS_S_AC2*Cp_CuS + Normalx_Ni3S4_S_AC2*Cp_Ni3S4_a;
    Cp_S_21 = Normalx_NiS_S_21*Cp_NiS + Normalx_Cu9S5_S_21*Cp_Cu9S5 + Normalx_CuS_S_21*Cp_CuS + Normalx_Ni3S4_S_21*Cp_Ni3S4_a;
    Cp_S_22 = Normalx_NiS_S_22*Cp_NiS + Normalx_Cu9S5_S_22*Cp_Cu9S5 + Normalx_CuS_S_22*Cp_CuS + Normalx_Ni3S4_S_22*Cp_Ni3S4_a;
    Cp_S_14 = Normalx_NiS_S_14*Cp_NiS + Normalx_Cu9S5_S_14*Cp_Cu9S5 + Normalx_CuS_S_14*Cp_CuS + Normalx_Ni3S4_S_14*Cp_Ni3S4_a;
    Cp_S_15 = Normalx_NiS_S_15*Cp_NiS + Normalx_Cu9S5_S_15*Cp_Cu9S5 + Normalx_CuS_S_15*Cp_CuS + Normalx_Ni3S4_S_15*Cp_Ni3S4_a;
    Cp_S_17 = Cp_S_15
    Cp_S_21 = Normalx_NiS_S_21*Cp_NiS + Normalx_Cu9S5_S_21*Cp_Cu9S5 + Normalx_CuS_S_21*Cp_CuS + Normalx_Ni3S4_S_21*Cp_Ni3S4_a;

    Q_removed_L_9 = (m_9*X[27])*Cp_L_9*(X[31]-100);
    Q_removed_S_9 = (m_9*X[22])*Cp_S_9*(X[31]-100);
    Q_removal_9 = Q_removed_L_9+Q_removed_S_9;

    if (Q_removal_9>0)
        m_H2O_flash_evaporated = Q_removal_9/H_Evap_H2O;
    else
        m_H2O_flash_evaporated = 0;
    end;
    if (m_H2O_flash_evaporated>m_9*X[27]*(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]))
        m_H2O_flash_evaporated = m_9*X[27]*(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]);
        temp_flag = 1;
    else
        temp_flag = 0;
    end;

        liq_fractions_1_H2SO4 = X[30]
        liq_fractions_1 = liq_fractions(liq_fractions_1_H2SO4, 0.0, 0.0, 0.0)

        solids_in_1_NiS   = X[82] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent))
        solids_in_1_Ni3S4 = X[83] * X[22]  * (1/(X[22]/density_residue + X[27]/density_spent))
        solids_in_1_Cu9S5 = -10000.0

        solids_in_1 = solid_fractions(solids_in_1_NiS, solids_in_1_Ni3S4, solids_in_1_Cu9S5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        Cu9S5_fresh = X[2] * flow_constants.prep_tank2 * X[1]^0.5 * X[5] / V_5

        solids_out_1_NiS = X[22] * X[23] / (X[22] / density_residue + X[27] / density_spent)
        solids_out_1_Ni3S4 = X[22] * X[24] / (X[22] / density_residue + X[27] / density_spent)
        solids_out_1_Cu9S5 = X[22] * X[25] / (X[22] / density_residue + X[27] / density_spent)
        solids_out_1_CuS = X[22] * X[26] / (X[22] / density_residue + X[27] / density_spent)

        solids_out_1 = solid_fractions(solids_out_1_NiS, solids_out_1_Ni3S4, solids_out_1_Cu9S5, solids_out_1_CuS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        Cu_molal_9 = 1000*(X[28]/(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]))/MW.Cu;
        Ni_molal_9 = 1000*(X[29]/(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]))/MW.Ni;
        H2SO4_molal_9 = 1000*(X[30]/(1-X[28]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[29]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[30]))/(2*MW.H+MW.S+4*MW.O);
        O2_solubility_9 = (2*MW.O/1000)*(calculate_oxygen_solubility2(X[31],P_O2_pascal,Cu_molal_9,Ni_molal_9,H2SO4_molal_9));
        molal_O2_solubility_9 = (O2_solubility_9*1000)/(2*MW.O);

        r_1 = calculate_reaction_rates(X[31], liq_fractions_1, solids_in_1, solids_out_1, molal_O2_solubility_9, Cu9S5_fresh, data)

        extent_1 = zeros(Float64,5);
        extent_1[1] = volume_1*60*r_1[1]
        extent_1[2] = volume_1*60*r_1[2]
        extent_1[3] = volume_1*60*r_1[3]
        extent_1[4] = volume_1*60*r_1[4]
        extent_1[5] = volume_1*60*r_1[5]
        solid_consumption_1 = zeros(Float64,4);
        solid_consumption_1[1] = ((MW.Ni+MW.S)/1000)*(4*extent_1[1]+extent_1[2]);
        solid_consumption_1[2] = ((3*MW.Ni+4*MW.S)/1000)*(-extent_1[1]+2*extent_1[3]);
        solid_consumption_1[3] = ((9*MW.Cu+5*MW.S)/1000)*(extent_1[4]);
        solid_consumption_1[4] = ((MW.Cu+MW.S)/1000)*(-5*extent_1[4]+extent_1[5]);

        liquid_generation_1 = zeros(Float64,3);
        liquid_generation_1[1] = (4*extent_1[4]+extent_1[5])*(MW.Cu/1000);
        liquid_generation_1[2] = (extent_1[1]+extent_1[2]+6*extent_1[3])*(MW.Ni/1000);
        liquid_generation_1[3] = ((2*MW.H+MW.S+4*MW.O)/1000)*((-2*extent_1[1]+4*extent_1[3]-8*extent_1[4])/2);

        for i=(1:2)
            if (liquid_generation_1[i]<0.0)
                display("Problem: liquid metal species consumption!")
            end;
        end;

        O2_consumption_1 = (0.5*extent_1[1]+2*extent_1[2]+15*extent_1[3]+2*extent_1[4]+2*extent_1[5])*(2*MW.O/1000);
        n_O2_consumption_1 = O2_consumption_1/(2*MW.O);
        n_H2O_condensation_1 = (n_O2_consumption_1/y_O2_2nd_stage)*(1-y_O2_2nd_stage);
        H2O_condensation_1 = n_H2O_condensation_1*(2*MW.H+MW.O);

        V_out_AC1 = flow_constants.AC1*comp_width*real((complex(slurry_height_1-weir_height_AC1))^(3/2))
        m_AC1 = V_out_AC1*1000*(1/(X[22]/density_residue+X[27]/density_spent))*60*60

        Q_Reaction_1 = data.h_reaction.h_reaction[1:5].*extent_1; 
        Q_Cond_1 = H2O_condensation_1*(H_Evap_H2O+Cp_H2O*(100-25)+Cp_steam_a*(X[31]-100));        

        V_S_AC1 = X[22]*m_AC1/density_residue;
        V_L_AC1 = X[27]*m_AC1/density_spent;
        V_AC1 = V_S_AC1 + V_L_AC1;
        
        liq_fractions_2_H2SO4 = X[40];
        liq_fractions_2 = liq_fractions(liq_fractions_2_H2SO4, 0.0, 0.0, 0.0)

        solids_in_2_NiS = X[84] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent));
        solids_in_2_Ni3S4 = X[85] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent));
        solids_in_2_Cu9S5 = -10000.0

        solids_in_2 = solid_fractions(solids_in_2_NiS, solids_in_2_Ni3S4, solids_in_2_Cu9S5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        V_AC2 = zero(Float64)
        solid_consumption_2 = zeros(Float64,4)
        extent_2 = zeros(Float64,5)
        liquid_generation_2 = zeros(Float64,3);
        O2_consumption_2 = zero(Float64)
        H2O_condensation_2 = zero(Float64)
        m_AC2 = zero(Float64)

        solids_out_2_NiS = X[32] * X[33] / (X[32] / density_residue + X[37] / density_spent)
        solids_out_2_Ni3S4 = X[32] * X[34] / (X[32] / density_residue + X[37] / density_spent)
        solids_out_2_Cu9S5 = X[32] * X[35] / (X[32] / density_residue + X[37] / density_spent)
        solids_out_2_CuS = X[32] * X[36] / (X[32] / density_residue + X[37] / density_spent)

        solids_out_2 = solid_fractions(solids_out_2_NiS, solids_out_2_Ni3S4, solids_out_2_Cu9S5, solids_out_2_CuS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        Cu_molal_AC2 = 1000*(X[38]/(1-X[38]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[39]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[40]))/MW.Cu;
        Ni_molal_AC2 = 1000*(X[39]/(1-X[38]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[39]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[40]))/MW.Ni;
        H2SO4_molal_AC2 = 1000*(X[40]/(1-X[38]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[39]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[40]))/(2*MW.H+MW.S+4*MW.O);
        O2_solubility_AC2 = (2*MW.O/1000)*(calculate_oxygen_solubility2(X[41],P_O2_pascal,Cu_molal_AC2,Ni_molal_AC2,H2SO4_molal_AC2));
        molal_O2_solubility_AC2 = (O2_solubility_AC2*1000)/(2*MW.O);     
        
        r_2 = calculate_reaction_rates(X[41],liq_fractions_2,solids_in_2,solids_out_2,molal_O2_solubility_AC2,Cu9S5_fresh,data);

        extent_2 = zeros(Float64,5);
        extent_2[1] = volume_2*60*r_2[1]
        extent_2[2] = volume_2*60*r_2[2]
        extent_2[3] = volume_2*60*r_2[3]
        extent_2[4] = volume_2*60*r_2[4]
        extent_2[5] = volume_2*60*r_2[5]

        solid_consumption_2[1] = ((MW.Ni+MW.S)/1000)*(4*extent_2[1]+extent_2[2]);
        solid_consumption_2[2] = ((3*MW.Ni+4*MW.S)/1000)*(-extent_2[1]+2*extent_2[3]);
        solid_consumption_2[3] = ((9*MW.Cu+5*MW.S)/1000)*(extent_2[4]);
        solid_consumption_2[4] = ((MW.Cu+MW.S)/1000)*(-5*extent_2[4]+extent_2[5]);

        liquid_generation_2[1] = (4*extent_2[4]+extent_2[5])*(MW.Cu/1000);
        liquid_generation_2[2] = (extent_2[1]+extent_2[2]+6*extent_2[3])*(MW.Ni/1000);
        liquid_generation_2[3] = ((2*MW.H+MW.S+4*MW.O)/1000)*((-2*extent_2[1]+4*extent_2[3]-8*extent_2[4])/2);
    
        for i = 1:2
            if (liquid_generation_2[i]<0.0)           
                display("Problem: liquid metal species consumption!")
            end;
        end;

        O2_consumption_2 = (0.5*extent_2[1]+2*extent_2[2]+15*extent_2[3]+2*extent_2[4]+2*extent_2[5])*(2*MW.O/1000);
        O2_consumption_AC4 = O2_consumption_2-m_10;
        n_O2_consumption_AC4 = O2_consumption_AC4/(2*MW.O);
        n_H2O_condensation_2 = (n_O2_consumption_AC4/y_O2_2nd_stage)*(1-y_O2_2nd_stage);
        H2O_condensation_2 = n_H2O_condensation_2*(2*MW.H+MW.O);

        V_out_AC2 = flow_constants.AC2*comp_width*real((complex(slurry_height_2-weir_height_AC2))^(3/2))
        m_AC2 = V_out_AC2*1000*(1/(X[32]/density_residue+X[37]/density_spent))*60*60

        Q_Reaction_2 = data.h_reaction.h_reaction[1:5].*extent_2; 
        Q_Cond_2 = H2O_condensation_2*(H_Evap_H2O+Cp_H2O*(100-25)+Cp_steam_a*(X[41]-100));

        Cp_10 = Cp_O2_a;
        T_10 = oxygen_T;

        liq_fractions_3_H2SO4 = X[50];
        liq_fractions_3 = liq_fractions(liq_fractions_3_H2SO4, 0.0, 0.0, 0.0)

        solids_in_3_NiS = X[86] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent));
        solids_in_3_Ni3S4 = X[87] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent));
        solids_in_3_Cu9S5 = -10000.0

        solids_in_3 = solid_fractions(solids_in_3_NiS, solids_in_3_Ni3S4, solids_in_3_Cu9S5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        solid_consumption_3 = zeros(Float64,4)
        extent_3 = zeros(Float64,5)
        liquid_generation_3 = zeros(Float64,3);
        O2_consumption_3 = zero(Float64)
        H2O_condensation_3 = zero(Float64)

        solids_out_3_NiS = X[42] * X[43] / (X[42] / density_residue + X[47] / density_spent)
        solids_out_3_Ni3S4 = X[42] * X[44] / (X[42] / density_residue + X[47] / density_spent)
        solids_out_3_Cu9S5 = X[42] * X[45] / (X[42] / density_residue + X[47] / density_spent)
        solids_out_3_CuS = X[42] * X[46] / (X[42] / density_residue + X[47] / density_spent)
             
        solids_out_3 = solid_fractions(solids_out_3_NiS, solids_out_3_Ni3S4, solids_out_3_Cu9S5, solids_out_3_CuS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        Cu_molal_14 = 1000*(X[48]/(1-X[48]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[49]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[50]))/MW.Cu;
        Ni_molal_14 = 1000*(X[49]/(1-X[48]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[49]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[50]))/MW.Ni;
        H2SO4_molal_14 = 1000*(X[50]/(1-X[48]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[49]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[50]))/(2*MW.H+MW.S+4*MW.O);
        O2_solubility_14 = (2*MW.O/1000)*(calculate_oxygen_solubility2(X[51],P_O2_pascal,Cu_molal_14,Ni_molal_14,H2SO4_molal_14));
        molal_O2_solubility_14 = (O2_solubility_14*1000)/(2*MW.O);
        
        r_3 = calculate_reaction_rates(X[51],liq_fractions_3,solids_in_3,solids_out_3,molal_O2_solubility_14,Cu9S5_fresh,data);
        
            extent_3 = zeros(Float64,5);
            extent_3[1] = volume_3*60*r_3[1]
            extent_3[2] = volume_3*60*r_3[2]
            extent_3[3] = volume_3*60*r_3[3]
            extent_3[4] = volume_3*60*r_3[4]
            extent_3[5] = volume_3*60*r_3[5]        
        
            solid_consumption_3[1] = ((MW.Ni+MW.S)/1000)*(4*extent_3[1]+extent_3[2]);
            solid_consumption_3[2] = ((3*MW.Ni+4*MW.S)/1000)*(-extent_3[1]+2*extent_3[3]);
            solid_consumption_3[3] = ((9*MW.Cu+5*MW.S)/1000)*(extent_3[4]);
            solid_consumption_3[4] = ((MW.Cu+MW.S)/1000)*(-5*extent_3[4]+extent_3[5]);

            liquid_generation_3[1] = (4*extent_3[4]+extent_3[5])*(MW.Cu/1000);
            liquid_generation_3[2] = (extent_3[1]+extent_3[2]+6*extent_3[3])*(MW.Ni/1000);
            liquid_generation_3[3] = ((2*MW.H+MW.S+4*MW.O)/1000)*((-2*extent_3[1]+4*extent_3[3]-8*extent_3[4])/2);

            for i = 1:2
                if (liquid_generation_3[i]<0.0)           
                    display("Problem: liquid metal species consumption!")
                end;
            end;

        O2_consumption_3 = (0.5*extent_3[1]+2*extent_3[2]+15*extent_3[3]+2*extent_3[4]+2*extent_3[5])*(2*MW.O/1000);
        O2_consumption_AC5 = O2_consumption_3-m_11;
        n_O2_consumption_AC5 = O2_consumption_AC5/(2*MW.O);
        n_H2O_condensation_3 = (n_O2_consumption_AC5/y_O2_2nd_stage)*(1-y_O2_2nd_stage);
        H2O_condensation_3 = n_H2O_condensation_3*(2*MW.H+MW.O);

        Q_Reaction_3 = data.h_reaction.h_reaction[1:5].*extent_3; 

        Q_Cond_3 = H2O_condensation_3*(H_Evap_H2O+Cp_H2O*(100-25)+Cp_steam_a*(X[51]-100));

        Cp_11 = Cp_O2_a;
        T_11 = oxygen_T;

        liq_fractions_4_H2SO4 = X[80];
        liq_fractions_4 = liq_fractions(liq_fractions_4_H2SO4, 0.0, 0.0, 0.0)

        solids_in_4_NiS = X[92] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent));
        solids_in_4_Ni3S4 = X[93] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent));
        solids_in_4_Cu9S5 = -10000.0
    
        solids_in_4 = solid_fractions(solids_in_4_NiS, solids_in_4_Ni3S4, solids_in_4_Cu9S5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                     
        Cp_12 = Cp_O2_a;
        T_12 = oxygen_T;

        T_13 = steam_T;

        solid_consumption_4 = zeros(Float64,4)
        extent_4 = zeros(Float64,5)
        liquid_generation_4 = zeros(Float64,3);
        O2_consumption_4 = zero(Float64)
        H2O_condensation_4 = zero(Float64)

        solids_out_4_NiS = X[72] * X[73] / (X[72] / density_residue + X[77] / density_spent)
        solids_out_4_Ni3S4 = X[72] * X[74] / (X[72] / density_residue + X[77] / density_spent)
        solids_out_4_Cu9S5 = X[72] * X[75] / (X[72] / density_residue + X[77] / density_spent)
        solids_out_4_CuS = X[72] * X[76] / (X[72] / density_residue + X[77] / density_spent)

        solids_out_4 = solid_fractions(solids_out_4_NiS, solids_out_4_Ni3S4, solids_out_4_Cu9S5, solids_out_4_CuS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        
        Cu_molal_22 = 1000*(X[78]/(1-X[78]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[79]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[80]))/MW.Cu;
        Ni_molal_22 = 1000*(X[79]/(1-X[78]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[79]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[80]))/MW.Ni;
        H2SO4_molal_22 = 1000*(X[80]/(1-X[78]*(MW.Cu+MW.S+4*MW.O)/MW.Cu-X[79]*(MW.Ni+MW.S+4*MW.O)/MW.Ni-X[80]))/(2*MW.H+MW.S+4*MW.O);
        O2_solubility_22 = (2*MW.O/1000)*(calculate_oxygen_solubility2(X[81],P_O2_pascal,Cu_molal_22,Ni_molal_22,H2SO4_molal_22));
        molal_O2_solubility_22 = (O2_solubility_22*1000)/(2*MW.O);

        r_4 = calculate_reaction_rates(X[81],liq_fractions_4,solids_in_4,solids_out_4,molal_O2_solubility_22,Cu9S5_fresh,data);

            extent_4 = zeros(Float64,5);
            extent_4[1] = volume_4*60*r_4[1]
            extent_4[2] = volume_4*60*r_4[2]
            extent_4[3] = volume_4*60*r_4[3]
            extent_4[4] = volume_4*60*r_4[4]
            extent_4[5] = volume_4*60*r_4[5]       

            solid_consumption_4[1] = ((MW.Ni+MW.S)/1000)*(4*extent_4[1]+extent_4[2]);
            solid_consumption_4[2] = ((3*MW.Ni+4*MW.S)/1000)*(-extent_4[1]+2*extent_4[3]);
            solid_consumption_4[3] = ((9*MW.Cu+5*MW.S)/1000)*(extent_4[4]);
            solid_consumption_4[4] = ((MW.Cu+MW.S)/1000)*(-5*extent_4[4]+extent_4[5]);

            liquid_generation_4[1] = (4*extent_4[4]+extent_4[5])*(MW.Cu/1000);
            liquid_generation_4[2] = (extent_4[1]+extent_4[2]+6*extent_4[3])*(MW.Ni/1000);
            liquid_generation_4[3] = ((2*MW.H+MW.S+4*MW.O)/1000)*((-2*extent_4[1]+4*extent_4[3]-8*extent_4[4])/2);
            
        for i = 1:2
            if (liquid_generation_3[i]<0.0)           
                display("Problem: liquid metal species consumption!")
            end;
        end;

        O2_consumption_4 = (0.5*extent_4[1]+2*extent_4[2]+15*extent_4[3]+2*extent_4[4]+2*extent_4[5])*(2*MW.O/1000);
        m_O2_AC5 = m_12-O2_consumption_4;
        n_O2_AC5 = m_O2_AC5/(2*MW.O);
        n_H2O_AC5 = n_O2_AC5/y_O2_2nd_stage*(1-y_O2_2nd_stage);
        H2O_condensation_4 = -n_H2O_AC5*(2*MW.H+MW.O);

        Q_Reaction_4 = data.h_reaction.h_reaction[1:5].*extent_4;
        Q_Cond_4 = H2O_condensation_4*(H_Evap_H2O+Cp_H2O*(100-25)+Cp_steam_a*(X[81]-100));

    solid_fraction_stream_17 = 1

    m_17 = X[53] * m_15 * 1 / solid_fraction_stream_17
    m_16 = m_15 - m_17

    overflow_tank_10 = flow_constants.overflow*real(complex(volume_tank10 - tank_10_max_vol)^(3/2))
    overflow_tank_20 = flow_constants.overflow*real(complex(volume_tank20 - tank_20_max_vol)^(3/2))
    overflow_tank_40 = flow_constants.overflow*real(complex(volume_tank40 - tank_40_max_vol)^(3/2))

    overflow_tank_50 = flow_constants.overflow*real(complex(volume_tank50 - tank_50_max_vol)^(3/2))

    # dX_dt = zeros(Float64, length(X))

    # dX_dt[1] = m_S_1+m_spent_1+m_2+m_3-(m_5+overflow_tank_10);
    # dX_dt[2] = (m_S_1 - X[2]*m_5 - X[2]*dX_dt[1])/X[1];
    # dX_dt[3:6] = (m_phases_in[1:4]-X[3:6].*X[2]*m_5-X[3:6].*(X[2]*dX_dt[1]+X[1]*dX_dt[2]))./(X[2]*X[1]);
    # dX_dt[7] = -1*dX_dt[2];
    # dX_dt[8:10] = (m_liquid_in[2:4] - X[7].*X[8:10]*m_5 - X[8:10].*(X[7]*dX_dt[1]+X[1]*dX_dt[7]))./(X[1]*X[7]);
      
    # dX_dt[100] = (m_S_1*Cp_S_1*(slurry_T-25)+m_spent_1*Cp_L_1*(slurry_T-25)+m_2*Cp_L_2*(spent_T-25)+m_3*Cp_L_3*(formic_T-25)- m_5*(X[7]*Cp_L_5+X[2]*Cp_S_5)*(X[100]-25) - (X[7]*Cp_L_5+X[2]*Cp_S_5)*(X[100]-25)*dX_dt[1])/(X[1]*(X[7]*Cp_L_5+X[2]*Cp_S_5))

    # dX_dt[11] = (m_5) + m_9 + m_23 + m_24 - m_H2O_flash_evaporated - (m_7+overflow_tank_20);
    # dX_dt[12] = (X[2]*(m_5)+m_9*X[22]-X[12]*(m_7+overflow_tank_20)-X[12]*dX_dt[11])/(X[11]);
    # dX_dt[13:16] = (((m_5))*X[2].*X[3:6]+m_9*X[22].*X[23:26]-((m_7+overflow_tank_20))*X[12].*X[13:16]-X[13:16].*(X[12]*dX_dt[11]+X[11]*dX_dt[12]))./(X[12]*X[11]);
    # dX_dt[94:95] = (((m_5))*X[2].*X[3:4]+m_9*X[22].*X[82:83]-((m_7+overflow_tank_20))*X[12].*X[94:95]-X[94:95].*(X[12]*dX_dt[11]+X[11]*dX_dt[12]))./(X[12]*X[11]);
    # dX_dt[17] = -1*dX_dt[12];
    # dX_dt[18:19] = (((m_5))*X[7].*X[8:9]+m_9*X[27].*X[28:29]-((m_7+overflow_tank_20))*X[17].*X[18:19]-X[18:19].*(X[17]*dX_dt[11]+X[11]*dX_dt[17]))./(X[17]*X[11]);
    # dX_dt[20] = (((m_5))*X[7].*X[10]+m_H2SO4_23+m_9*X[27].*X[30]-((m_7+overflow_tank_20))*X[17].*X[20]-X[20].*(X[17]*dX_dt[11]+X[11]*dX_dt[17]))./(X[17]*X[11]);
    # dX_dt[21] = (((m_5))*(X[7]*Cp_L_5+X[2]*Cp_S_5)*(X[100]-25)+m_23*Cp_L_23*(acid_T-25)+m_24*Cp_L_24*(water_T-25)+(m_9-m_H2O_flash_evaporated)*(X[27]*Cp_L_9+X[22]*Cp_S_9)*(100-25)-(m_7+overflow_tank_20)*(X[12]*Cp_S_7+X[17]*Cp_L_7)*(X[21]-25) - (X[12]*Cp_S_7+X[17]*Cp_L_7)*(X[21]-25)*dX_dt[11])    /   (X[11]*(X[12]*Cp_S_7+X[17]*Cp_L_7));
  
    # dX_dt[96] = m_7 - m_9 - m_AC1 + O2_consumption_1 + H2O_condensation_1;
    # dX_dt[22] = (X[12]*m_7-m_9*X[22]-X[22]*m_AC1-sum(solid_consumption_1)-X[22]*dX_dt[96])/(X[96]);
    # dX_dt[23:26] = (X[12]*m_7.*X[13:16]-m_9*X[22].*X[23:26]-m_AC1*X[22].*X[23:26]-solid_consumption_1-(X[96]*dX_dt[22]+X[22]*dX_dt[96]).*X[23:26])./(X[22]*X[96]);
    # dX_dt[82:83] = (X[12]*m_7.*X[94:95]-m_9*X[22].*X[82:83]-m_AC1*X[22].*X[82:83]-(X[96]*dX_dt[22]+X[22]*dX_dt[96]).*X[82:83])./(X[22]*X[96]);
    # dX_dt[27] = -1*dX_dt[22];
    # dX_dt[28:30] = (X[17]*m_7.*X[18:20]-m_9*X[27].*X[28:30]-m_AC1*X[27].*X[28:30]+liquid_generation_1-(X[96]*dX_dt[27]+X[27]*dX_dt[96]).*X[28:30])./(X[27]*X[96]);
    # dX_dt[31] = (3*agitator_power - Q_loss*length_1 + ((X[17]*Cp_L_7+X[12]*Cp_S_7)*(m_7)*(X[21]-25)) - ((X[31]-25)*(X[27]*Cp_L_9+X[22]*Cp_S_9)*(m_9+m_AC1)) + Q_Cond_1 + O2_consumption_1*Cp_O2_a*(X[31]-25) - sum(Q_Reaction_1)-(X[31]-25.0)*(X[27]*Cp_L_9+X[22]*Cp_S_9)*dX_dt[96])/(X[96]*(X[27]*Cp_L_9+X[22]*Cp_S_9));

    # dX_dt[97] = m_AC1 - m_AC2 + O2_consumption_2 + H2O_condensation_2
    # dX_dt[32] = (X[22]*m_AC1-X[32]*m_AC2-sum(solid_consumption_2)-X[32]*dX_dt[97])/(X[97]);
    # dX_dt[33:36] = (m_AC1*X[22].*X[23:26]-m_AC2*X[32].*X[33:36]-solid_consumption_2-(X[97]*dX_dt[32]+X[32]*dX_dt[97]).*X[33:36])./(X[97]*X[32]);
    # dX_dt[84:85] = (m_AC1*X[22].*X[82:83]-m_AC2*X[32].*X[84:85]-(X[97]*dX_dt[32]+X[32]*dX_dt[97]).*X[84:85])./(X[97]*X[32]);
    # dX_dt[37] = -1*dX_dt[32];
    # dX_dt[38:40] = (m_AC1*X[27].*X[28:30]-m_AC2*X[37].*X[38:40]+liquid_generation_2-(X[97]*dX_dt[37]+X[37]*dX_dt[97]).*(X[38:40]))./(X[37]*X[97]);
    # dX_dt[41] = (-1*Q_Removed_2-Q_loss*length_2+agitator_power+(m_AC1*(X[31]-25)*(X[22]*Cp_S_9+X[27]*Cp_L_9)+m_10*Cp_10*(T_10-25))-(m_AC2*(X[41]-25)*(X[37]*Cp_L_AC2+X[32]*Cp_S_AC2))-sum(Q_Reaction_2)+Q_Cond_2+(O2_consumption_2-m_10)*Cp_10*(X[41]-25)-(X[41]-25.0)*(X[37]*Cp_L_AC2+X[32]*Cp_S_AC2)*dX_dt[97])/(X[97]*(X[37]*Cp_L_AC2+X[32]*Cp_S_AC2));
    
    # dX_dt[98] = m_AC2 - m_14 + O2_consumption_3 + H2O_condensation_3
    # dX_dt[42] = (X[32]*m_AC2-X[42]*m_14-sum(solid_consumption_3)-X[42]*dX_dt[98])/(X[98]);
    # dX_dt[43:46] = (X[32]*m_AC2.*X[33:36]-m_14*X[42].*X[43:46]-solid_consumption_3-(X[98]*dX_dt[42]+X[42]*dX_dt[98]).*(X[43:46]))./(X[42]*X[98]);
    # dX_dt[86:87] = (X[32]*m_AC2.*X[84:85]-m_14*X[42].*X[86:87]-(X[98]*dX_dt[42]+X[42]*dX_dt[98]).*(X[86:87]))./(X[42]*X[98]);
    # dX_dt[47] = -1*dX_dt[42];
    # dX_dt[48:50] = (X[37]*m_AC2.*X[38:40]-m_14*X[47].*X[48:50]+liquid_generation_3-(X[98]*dX_dt[47]+X[47]*dX_dt[98]).*(X[48:50]))./(X[98]*X[47]);
    # dX_dt[51] = (-1*Q_Removed_3-Q_loss*length_3+agitator_power+(m_AC2*(X[41]-25)*(X[32]*Cp_S_AC2+X[37]*Cp_L_AC2)+m_11*Cp_11*(T_11-25))-(m_14*(X[51]-25)*(X[47]*Cp_L_14+X[42]*Cp_S_14))-sum(Q_Reaction_3)+Q_Cond_3+(O2_consumption_3-m_11)*Cp_11*(X[51]-25)-(X[51]-25.0)*(X[47]*Cp_L_14+X[42]*Cp_S_14)*dX_dt[98])/(X[98]*(X[47]*Cp_L_14+X[42]*Cp_S_14));
    
    # dX_dt[52] = m_14-(m_15+overflow_tank_40);
    # dX_dt[53] = (X[42]*m_14-X[53]*(m_15+overflow_tank_40)-X[53]*dX_dt[52])/(X[52]);
    # dX_dt[54:57] = (m_14*X[42].*X[43:46]-((m_15+overflow_tank_40))*X[53].*X[54:57]-(X[52]*dX_dt[53]+X[53]*dX_dt[52]).*X[54:57])./(X[53]*X[52]);
    # dX_dt[88:89] = (m_14*X[42].*X[86:87]-((m_15+overflow_tank_40))*X[53].*X[88:89]-(X[52]*dX_dt[53]+X[53]*dX_dt[52]).*X[88:89])./(X[53]*X[52]);
    # dX_dt[58] = -1*dX_dt[53];
    # dX_dt[59:61] = (m_14*X[47].*X[48:50]-((m_15+overflow_tank_40))*X[58].*X[59:61]-(X[52]*dX_dt[58]+X[58]*dX_dt[52]).*X[59:61])./(X[58]*X[52]);
    # dX_dt[101] = (-1*Q_Removed_TK40+m_14*(X[51]-25)*(X[47]*Cp_L_14+X[42]*Cp_S_14)-(m_15+overflow_tank_40)*(X[101]-25)*(X[53]*Cp_S_15+X[58]*Cp_L_15)-(X[53]*Cp_S_15+X[58]*Cp_L_15)*(X[101]-25)*dX_dt[52])/((X[53]*Cp_S_15+X[58]*Cp_L_15)*X[52])

    # dX_dt[62] = m_17+m_18+m_19+m_20-(m_21+overflow_tank_50);
    # dX_dt[63] = (solid_fraction_stream_17*m_17 - X[63]*(m_21+overflow_tank_50) - X[63]*dX_dt[62])/X[62];
    # dX_dt[64:67] = (solid_fraction_stream_17*m_17.*X[54:57]-X[64:67].*X[63]*(m_21+overflow_tank_50)-X[64:67].*(X[63]*dX_dt[62]+X[62]*dX_dt[63]))./(X[62]*X[63]);
    # dX_dt[90:91] = (solid_fraction_stream_17*m_17.*X[88:89]-X[90:91].*X[63]*(m_21+overflow_tank_50)-X[90:91].*(X[63]*dX_dt[62]+X[62]*dX_dt[63]))./(X[62]*X[63]);
    # dX_dt[68] = -1*dX_dt[63];
    # dX_dt[69:71] = ((1-solid_fraction_stream_17)*m_17*X[59:61] + m_liquid_TK50[2:4] - X[68].*X[69:71].*(m_21+overflow_tank_50) - X[69:71].*(X[68]*dX_dt[62]+X[62]*dX_dt[68]))./(X[62]*X[68]);
    # dX_dt[102] = (m_18*Cp_L_18*(spent_T-25)+m_19*Cp_L_19*(water_T-25)+m_20*Cp_L_20*(acid_T-25)+m_17*solid_fraction_stream_17*Cp_S_17*(X[101]-25)+m_17*(1-solid_fraction_stream_17)*Cp_L_17*(X[101]-25) - (m_21+overflow_tank_50)*(X[63]*Cp_S_21+X[68]*Cp_L_21)*(X[102]-25)-(X[63]*Cp_S_21+X[68]*Cp_L_21)*(X[102]-25)*dX_dt[62])/((X[63]*Cp_S_21+X[68]*Cp_L_21)*X[62])

    # dX_dt[99] = m_21 + m_13 + O2_consumption_4 + H2O_condensation_4 - m_22
    # dX_dt[72] = (m_21*X[63]-m_22*X[72]-sum(solid_consumption_4)-X[72]*dX_dt[99])/(X[99]);
    # dX_dt[73:76] = (m_21*X[63].*X[64:67]-m_22*X[72].*X[73:76]-solid_consumption_4-(X[99]*dX_dt[72]+X[72]*dX_dt[99]).*(X[73:76]))./(X[99]*X[72]);
    # dX_dt[92:93] = (m_21*X[63].*X[90:91]-m_22*X[72].*X[92:93]-(X[99]*dX_dt[72]+X[72]*dX_dt[99]).*(X[92:93]))./(X[99]*X[72]);
    # dX_dt[77] = -1*dX_dt[72];
    # dX_dt[78:80] = (m_21*X[68].*X[69:71]-m_22*X[77].*X[78:80]+liquid_generation_4-(X[99]*dX_dt[77]+X[77]*dX_dt[99]).*(X[78:80]))./(X[99]*X[77]);
    # dX_dt[81] = (agitator_power-Q_loss*length_4+m_21*(X[68]*Cp_L_21+X[63]*Cp_S_21)*(X[102]-25)+m_12*Cp_12*(T_12-25)+m_13*(Cp_H2O*(100-25)+H_Evap_H2O+Cp_steam_a*(T_13-100))-m_22*(X[77]*Cp_L_22+X[72]*Cp_S_22)*(X[81]-25)-(m_12-O2_consumption_4)*Cp_12*(X[81]-25)-sum(Q_Reaction_4)+Q_Cond_4-(X[81]-25.0)*(X[77]*Cp_L_22+X[72]*Cp_S_22)*dX_dt[99])/(X[99]*(X[77]*Cp_L_22+X[72]*Cp_S_22));

    # dX_dt[103] = m_O2_total / (2*MW.O) - (O2_consumption_1/(2*MW.O) + O2_consumption_2/(2*MW.O) + O2_consumption_3/(2*MW.O) + O2_consumption_4/(2*MW.O)) - y_O2_vap * flow_constants.O2out * sqrt(P_autoclave*100000-P_atm)
    
    # dX_dt[104] = level_error_tank_10
    # dX_dt[105] = level_error_tank_20
    # dX_dt[106] = level_error_AC3
    # dX_dt[107] = level_error_tank_40
    # dX_dt[108] = level_error_tank_50
    # dX_dt[109] = level_error_AC4
    # dX_dt[110] = temperature_error_AC1
    # dX_dt[111] = temperature_error_AC2
    # dX_dt[112] = temperature_error_AC3
    # dX_dt[113] = temperature_error_AC4
    # dX_dt[114] = pressure_error
    # dX_dt[115] = temperature_error_TK40
   
        mass_AC1 = volume_1/(X[22]/density_residue+X[27]/density_spent)
        mass_AC2 = volume_2/(X[32]/density_residue+X[37]/density_spent)
        mass_14 = volume_3/(X[42]/density_residue+X[47]/density_spent)
        mass_22 = volume_4/(X[72]/density_residue+X[77]/density_spent)

        H2SO4_concentration_AC1 = X[30] * density_spent * 1000  #g/l
        H_plus_concentration_AC1 = 2 * H2SO4_concentration_AC1 / (2*MW.H + MW.S + MW.O*4)   #mol/l
        H2SO4_concentration_AC2 = X[40] * density_spent * 1000  #g/l
        H_plus_concentration_AC2 = 2 * H2SO4_concentration_AC2 / (2*MW.H + MW.S + MW.O*4)   #mol/l
        H2SO4_concentration_AC3 = X[50] * density_spent * 1000  #g/l
        H_plus_concentration_AC3 = 2 * H2SO4_concentration_AC3 / (2*MW.H + MW.S + MW.O*4)   #mol/l
        H2SO4_concentration_AC4 = X[80] * density_spent * 1000  #g/l
        H_plus_concentration_AC4 = 2 * H2SO4_concentration_AC4 / (2*MW.H + MW.S + MW.O*4)   #mol/l

        Cu_concentration_AC1 = X[28] * density_spent * 1000 #g/l
        Cu_concentration_AC2 = X[38] * density_spent * 1000 #g/l
        Cu_concentration_AC3 = X[48] * density_spent * 1000 #g/l
        Cu_concentration_AC4 = X[78] * density_spent * 1000 #g/l
        
        NiS_init_concentration_TK10 = X[3] * X[2] * (1/(X[2]/density_residue + X[7]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_FRT = X[94] * X[12] * (1/(X[12]/density_residue + X[17]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_AC1 = X[82] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_AC2 = X[84] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_AC3 = X[86] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent)) #kg NiS_init / l
        NiS_init_concentration_AC4 = X[92] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent))

        NiS_conc_TK10 = X[3] * X[2] * (1/(X[2]/density_residue + X[7]/density_spent)) 
        NiS_conc_FRT = X[13] * X[12] * (1/(X[12]/density_residue + X[17]/density_spent))
        NiS_conc_AC1 = X[23] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent))
        NiS_conc_AC2 = X[33] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent))
        NiS_conc_AC3 = X[43] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent))
        NiS_conc_AC4 = X[73] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent))
        
        Ni3S4_init_concentration_TK10 = X[4] * X[2] * (1/(X[2]/density_residue + X[7]/density_spent)) #kg Ni3S4_init / l
        
        #solids_in_1_Ni3S4_check = X[12] * X[14] * (1/(X[12]/density_residue + X[17]/density_spent))       #Concentration of solids entering the first autoclave compartment, kg Ni3S4/l slurry
        #solids_in_1_Ni3S4 = X[12] * flow_constants.recycle_tank * X[11]^0.5 * X[14] / V_7          #Concentration of solids entering the first autoclave compartment, kg Ni3S4/l slurry
        
        Ni3S4_init_concentration_FRT = X[95] * X[12] * (1/(X[12]/density_residue + X[17]/density_spent)) #kg Ni3S4_init / l
        Ni3S4_init_concentration_AC1 = X[83] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent)) #kg Ni3S4_init / l
        Ni3S4_init_concentration_AC2 = X[85] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent)) #kg Ni3S4_init / l
        Ni3S4_init_concentration_AC3 = X[87] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent)) #kg Ni3S4_init / l
        Ni3S4_init_concentration_AC4 = X[93] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent)) #kg Ni3S4_init / l

        Ni3S4_conc_TK10 = X[4] * X[2] * (1/(X[2]/density_residue + X[7]/density_spent)) 
        Ni3S4_conc_FRT = X[14] * X[12] * (1/(X[12]/density_residue + X[17]/density_spent))
        Ni3S4_conc_AC1 = X[24] * X[22] * (1/(X[22]/density_residue + X[27]/density_spent))
        Ni3S4_conc_AC2 = X[34] * X[32] * (1/(X[32]/density_residue + X[37]/density_spent))
        Ni3S4_conc_AC3 = X[44] * X[42] * (1/(X[42]/density_residue + X[47]/density_spent))
        Ni3S4_conc_AC4 = X[74] * X[72] * (1/(X[72]/density_residue + X[77]/density_spent))  
        
        MM_NiS_1 = (1000*solids_out_1_NiS/(MW.Ni+MW.S))/(0.001+1000*solids_out_1_NiS/(MW.Ni+MW.S))
        MM_Ni3S4_1 = (1000*solids_out_1_Ni3S4/(3*MW.Ni+4*MW.S))/(0.0001+1000*solids_out_1_Ni3S4/(3*MW.Ni+4*MW.S))
        MM_Cu9S5_1 = (1000*solids_out_1_Cu9S5/(9*MW.Cu+5*MW.S))/(0.0001+1000*solids_out_1_Cu9S5/(9*MW.Cu+5*MW.S))
        MM_CuS_1 = (1000*solids_out_1_CuS/(MW.Cu+MW.S))/(0.0001+1000*solids_out_1_CuS/(MW.Cu+MW.S))

        MM_NiS_2 = (1000*solids_out_2_NiS/(MW.Ni+MW.S))/(0.001+1000*solids_out_2_NiS/(MW.Ni+MW.S))
        MM_Ni3S4_2 = (1000*solids_out_2_Ni3S4/(3*MW.Ni+4*MW.S))/(0.0001+1000*solids_out_2_Ni3S4/(3*MW.Ni+4*MW.S))
        MM_Cu9S5_2 = (1000*solids_out_2_Cu9S5/(9*MW.Cu+5*MW.S))/(0.0001+1000*solids_out_2_Cu9S5/(9*MW.Cu+5*MW.S))
        MM_CuS_2 = (1000*solids_out_2_CuS/(MW.Cu+MW.S))/(0.0001+1000*solids_out_2_CuS/(MW.Cu+MW.S))

        MM_NiS_3 = (1000*solids_out_3_NiS/(MW.Ni+MW.S))/(0.001+1000*solids_out_3_NiS/(MW.Ni+MW.S))
        MM_Ni3S4_3 = (1000*solids_out_3_Ni3S4/(3*MW.Ni+4*MW.S))/(0.0001+1000*solids_out_3_Ni3S4/(3*MW.Ni+4*MW.S))
        MM_Cu9S5_3 = (1000*solids_out_3_Cu9S5/(9*MW.Cu+5*MW.S))/(0.0001+1000*solids_out_3_Cu9S5/(9*MW.Cu+5*MW.S))
        MM_CuS_3 = (1000*solids_out_3_CuS/(MW.Cu+MW.S))/(0.0001+1000*solids_out_3_CuS/(MW.Cu+MW.S))

        MM_NiS_4 = (1000*solids_out_4_NiS/(MW.Ni+MW.S))/(0.001+1000*solids_out_4_NiS/(MW.Ni+MW.S))
        MM_Ni3S4_4 = (1000*solids_out_4_Ni3S4/(3*MW.Ni+4*MW.S))/(0.0001+1000*solids_out_4_Ni3S4/(3*MW.Ni+4*MW.S))
        MM_Cu9S5_4 = (1000*solids_out_4_Cu9S5/(9*MW.Cu+5*MW.S))/(0.0001+1000*solids_out_4_Cu9S5/(9*MW.Cu+5*MW.S))
        MM_CuS_4 = (1000*solids_out_4_CuS/(MW.Cu+MW.S))/(0.0001+1000*solids_out_4_CuS/(MW.Cu+MW.S))

        V_16 = m_16 / density_spent #it is assumed stream 16 contains no solids
        V_17 = solid_fraction_stream_17*m_17/density_residue + (1-solid_fraction_stream_17)*m_17/density_spent

        wt_percent_copper_initial = 100*(X[5]*MW.Cu*9/(MW.Cu*9+MW.S*5) + X[6]*MW.Cu/(MW.Cu+MW.S))
        wt_percent_Ni_initial = 100*(X[3]*MW.Ni/(MW.Ni+MW.S) + X[4]*MW.Ni*3/(MW.Ni*3+MW.S*4))
        wt_percent_S_initial = 100*(X[5]*MW.S*5/(MW.Cu*9+MW.S*5) + X[6]*MW.S/(MW.Cu+MW.S) + X[3]*MW.S/(MW.Ni+MW.S) + X[4]*MW.S*4/(MW.Ni*3+MW.S*4))

        wt_percent_copper_second_stage_leach = 100*(X[45]*MW.Cu*9/(MW.Cu*9+MW.S*5) + X[46]*MW.Cu/(MW.Cu+MW.S))
        wt_percent_Ni_second_stage_leach = 100*(X[43]*MW.Ni/(MW.Ni+MW.S) + X[44]*MW.Ni*3/(MW.Ni*3+MW.S*4))
        wt_percent_S_second_stage_leach = 100*(X[45]*MW.S*5/(MW.Cu*9+MW.S*5) + X[46]*MW.S/(MW.Cu+MW.S) + X[43]*MW.S/(MW.Ni+MW.S) + X[44]*MW.S*4/(MW.Ni*3+MW.S*4))

        wt_percent_copper_third_stage_leach = 100*(X[75]*MW.Cu*9/(MW.Cu*9+MW.S*5) + X[76]*MW.Cu/(MW.Cu+MW.S))
        wt_percent_Ni_third_stage_leach = 100*(X[73]*MW.Ni/(MW.Ni+MW.S) + X[74]*MW.Ni*3/(MW.Ni*3+MW.S*4))
        wt_percent_S_third_stage_leach = 100*(X[75]*MW.S*5/(MW.Cu*9+MW.S*5) + X[76]*MW.S/(MW.Cu+MW.S) + X[73]*MW.S/(MW.Ni+MW.S) + X[74]*MW.S*4/(MW.Ni*3+MW.S*4))
        
        density_stream_5 =  1 / (X[2] / density_residue + X[7] / density_spent)

        y_vect = zeros(49)

        y_vect[1] = V_1
        y_vect[2] = V_2
        y_vect[3] = V_3
        y_vect[4] = X[100]
        y_vect[5] = tank_10_percent_full
        y_vect[6] = V_5
        y_vect[7] = density_stream_5
        y_vect[8] = V_23
        y_vect[9] = V_7
        y_vect[10] = X[21]
        y_vect[11] = tank_20_percent_full
        y_vect[12] = V_9
        y_vect[13] = X[31]
        y_vect[14] = comp_1_percent_full
        y_vect[15] = X[41]
        y_vect[16] = comp_3_percent_full
        y_vect[17] = X[51]
        y_vect[18] = m_10
        y_vect[19] = m_11
        y_vect[20] = m_12
        y_vect[21] = P_autoclave
        y_vect[22] = T_vap
        y_vect[23] = V_14
        y_vect[24] = X[101]
        y_vect[25] = tank_40_percent_full
        y_vect[26] = V_15
        y_vect[27] = V_17
        y_vect[28] = V_18
        y_vect[29] = V_19
        y_vect[30] = V_20
        y_vect[31] = tank_50_percent_full
        y_vect[32] = X[102]
        y_vect[33] = V_21
        y_vect[34] = X[81]
        y_vect[35] = comp_4_percent_full
        y_vect[36] = V_22
        y_vect[37] = H2SO4_concentration_AC1
        y_vect[38] = H2SO4_concentration_AC2
        y_vect[39] = H2SO4_concentration_AC3
        y_vect[40] = H2SO4_concentration_AC4
        y_vect[41] = wt_percent_copper_initial = wt_percent_copper_initial
        y_vect[42] = wt_percent_Ni_initial = wt_percent_Ni_initial
        y_vect[43] = wt_percent_S_initial = wt_percent_S_initial
        y_vect[44] = wt_percent_copper_second_stage_leach = wt_percent_copper_second_stage_leach
        y_vect[45] = wt_percent_Ni_second_stage_leach = wt_percent_Ni_second_stage_leach
        y_vect[46] = wt_percent_S_second_stage_leach = wt_percent_S_second_stage_leach
        y_vect[47] = wt_percent_copper_third_stage_leach = wt_percent_copper_third_stage_leach
        y_vect[48] = wt_percent_Ni_third_stage_leach = wt_percent_Ni_third_stage_leach
        y_vect[49] = wt_percent_S_third_stage_leach = wt_percent_S_third_stage_leach
        
        
    return y_vect
end


function BMR_gt_generate_measurements(X, data, t, vkstd)
    #add noise sampled from normal distribution to y_prediction measurement
    y_prediction = BMR_gt_measurement_prediction(X, data, t)

    y_pred = y_prediction[1:length(vkstd)]       #this allows me to add variables of interest to gt_measurment_prediction and observe variables in the ODE, without changing the measurment function...
    yk = y_pred .+ vkstd .* randn(length(y_pred))
    return yk
end


function BMR_gt_generate_measurements_with_samples(X, data, t, vkstd)
    #add noise sampled from normal distribution to y_prediction measurement
    y_prediction = BMR_gt_measurement_prediction_with_samples(X, data, t)

    y_pred = y_prediction[1:length(vkstd)]       #this allows me to add variables of interest to gt_measurment_prediction and observe variables in the ODE, without changing the measurment function...
    yk = y_pred .+ vkstd .* randn(length(y_pred))
    return yk
end