mutable struct Params
    u1 #m^3/s
    u2 #m^3/s
    c1 #Normal(1.0, 0.0025)
    Sp #m^2
    g #m/s^2
    A #m^2
    c2 #Normal(0.8, 0.0025)
    c3 #Normal(1.0, 0.0025)
    alpha #Normal(0.6, 4e-4)
    r # Normal(0.002, 10e-6)
    faultmode #normal (0), fault 1 (1), fault 2 (2)
end

#RP stands for real process
#MP stands for model process

function Three_tank_GT_process!(x, p, t)
    u1, u2, c1, Sp, g, A, c2, c3, alpha, r, faultmode = p.u1, p.u2, p.c1, p.Sp, p.g, p.A, p.c2, p.c3, p.alpha, p.r, p.faultmode

    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    qf = 0

    if faultmode == 1
        u1 = u1 + (alpha - 1) * u1
    end
    if faultmode == 2
        qf = c2 * pi * r ^ 2 * sqrt(2 * g * x2)
    end

    dx = zeros(length(x))

    dx[1] = (-c1 * Sp * sign(x1 - x3) * sqrt(2 * g * abs(x1 - x3)) + u1)/(A)
    dx[2] = (-c3 * Sp * sign(x2 - x3) * sqrt(2 * g * abs(x2 - x3)) - c2 * Sp * sqrt(2 * g * x2) - qf + u2)/(A)
    dx[3] = (c1 * Sp * sign(x1 - x3) * sqrt(2 * g * abs(x1 - x3)) - c3 * Sp * sign(x3 - x2) * sqrt(2 * g * abs(x3 - x2)))/(A)

    return dx

end

function Three_tank_model_process!(x, p, t)
    u1, u2, c1, Sp, g, A, c2, c3, alpha, r, faultmode = p.u1, p.u2, p.c1, p.Sp, p.g, p.A, p.c2, p.c3, p.alpha, p.r, p.faultmode
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    qf = 0

    if faultmode == 1
        u1 = u1 + (alpha - 1) * u1
    end
    if faultmode == 2
        qf = c2 * pi * r ^ 2 * sqrt(2 * g * x2)
    end    

    dx = zeros(length(x))
    dx[1] = (-c1 * Sp * sign(x1 - x3) * sqrt(2 * g * abs(x1 - x3)) + u1)/(A)
    dx[2] = (-c3 * Sp * sign(x2 - x3) * sqrt(2 * g * abs(x2 - x3)) - c2 * Sp * sqrt(2 * g * x2) - qf + u2)/(A)
    dx[3] = (c1 * Sp * sign(x1 - x3) * sqrt(2 * g * abs(x1 - x3)) - c3 * Sp * sign(x3 - x2) * sqrt(2 * g * abs(x3 - x2)))/(A)
    return dx

end
