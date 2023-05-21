function f(X, Q, u, t)
    # Unpack state variables
    x, y, psi, v = X
    q = Q[1]

    # Define constants
    c = 1.5
    rho = 1.226
    m = 300
    Iz = 1000
    L = 6
    B = 2
    g = 9.81

    # Compute forces
    T = u[1]
    phi = u[2]
    beta = psi - phi
    Fa = 0.5 * rho * v^2 * c * L * B
    Fr = -0.5 * rho * v^2 * L * B * sin(beta)
    Fh = T
    Fx = Fa + Fr - Fh * sin(beta)
    Fy = Fh * cos(beta) - m * g
    Mz = Fh * L * cos(beta) - Fr * L * sin(beta)

    # Compute state derivatives
    dxdt = v * cos(psi)
    dydt = v * sin(psi)
    dpsidt = v * L / Iz * sin(beta)
    dvdt = 1 / m * Fy

    dXdt = [dxdt; dydt; dpsidt; dvdt]
    return dXdt
end
