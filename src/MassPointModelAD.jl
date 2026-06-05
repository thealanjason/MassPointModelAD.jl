module MassPointModelAD

using Interpolations: Interpolations
using OrdinaryDiffEq: OrdinaryDiffEq, ODEProblem, Tsit5, solve

export MPMParams, DEM, load_dem, dem_from_func, MPM_rhs, simulate

include("dem.jl")
include("utils.jl")

"""
    MPMParams(coulomb_friction, turbulent_friction, g)

Parameters for the mass point model.

- `coulomb_friction`: dry Coulomb friction coefficient (dimensionless)
- `turbulent_friction`: turbulent friction coefficient (m/s²); use `Inf` to disable
- `g`: gravitational acceleration magnitude (m/s²)
"""
struct MPMParams{T<:Real}
    coulomb_friction::T
    turbulent_friction::T
    g::T
end

MPMParams(coulomb_friction, turbulent_friction, g) =
    MPMParams(promote(coulomb_friction, turbulent_friction, g)...)

function _dem_derivatives(dem::DEM, x, y)
    return dem.fgrad_x(x, y),
    dem.fgrad_y(x, y),
    dem.fgrad_xx(x, y),
    dem.fgrad_xy(x, y),
    dem.fgrad_yy(x, y)
end

function _velocity_magnitude(ux, uy, gradx, grady)
    cross = 2 * ux * uy * gradx * grady / sqrt(1 + gradx^2) / sqrt(1 + grady^2)
    return sqrt(max(ux^2 + uy^2 + cross, zero(ux)))
end

"""
    MPM_rhs(alpha, params, dem; curvature=false)

Right-hand side of the mass point ODE system (eqs. 12–15).

`alpha = [x, y, ux, uy]`. Returns `dalphadt`.
"""
function MPM_rhs(alpha, params::MPMParams, dem::DEM; curvature::Bool = false)
    x, y, ux, uy = alpha[1], alpha[2], alpha[3], alpha[4]
    coulomb_friction = params.coulomb_friction
    turbulent_friction = params.turbulent_friction
    g = params.g

    gradx, grady, gradxx, gradxy, gradyy = _dem_derivatives(dem, x, y)

    sqrt_x = sqrt(1 + gradx^2)
    sqrt_y = sqrt(1 + grady^2)
    sqrt_xy = sqrt(1 + gradx^2 + grady^2)

    gx = gradx * (-g) / sqrt_x
    gy = grady * (-g) / sqrt_y
    gn = g / sqrt_xy

    if curvature
        kx = gradxx / (1 + gradx^2) / sqrt_xy
        kxy = gradxy / sqrt_x / sqrt_y / sqrt_xy
        ky = gradyy / (1 + grady^2) / sqrt_xy
        uku = kx * ux^2 + 2 * kxy * ux * uy + ky * uy^2
    else
        uku = zero(gn)
    end

    u = _velocity_magnitude(ux, uy, gradx, grady)

    dalphadt = similar(alpha)
    dalphadt[1] = ux / sqrt_x
    dalphadt[2] = uy / sqrt_y

    if u != 0
        friction = coulomb_friction * (gn + uku) + g * u^2 / turbulent_friction
        dalphadt[3] = gx - (ux / u) * friction
        dalphadt[4] = gy - (uy / u) * friction
    else
        friction = coulomb_friction * gn
        dalphadt[3] =
            gx +
            friction * gradx * (1 + gradx^2 + grady^2) / sqrt_x /
            sqrt(gradx^2 + grady^2 + (gradx^2 + grady^2)^2)
        dalphadt[4] =
            gy +
            friction * grady * (1 + gradx^2 + grady^2) / sqrt_y /
            sqrt(gradx^2 + grady^2 + (gradx^2 + grady^2)^2)
    end

    return dalphadt
end

"""
    simulate(dem, params, x0, y0; ux0, uy0, dt, tend, t0, curvature, atol, rtol)

Integrate the mass point model using the Tsit5 adaptive solver, saving the
state at fixed time intervals `dt`.

# Arguments
- `dem`: `DEM` with elevation and derivative interpolants.
- `params`: `MPMParams` with friction and gravity values.
- `x0`, `y0`: Initial position.
- `ux0`, `uy0`: Initial velocity components (default 0).
- `dt`: Output time step in seconds (default 1.0).
- `tend`: End time in seconds (default 300.0).
- `t0`: Start time in seconds (default 0.0).
- `curvature`: Include topographic curvature effects (default `false`).
- `atol`, `rtol`: Absolute and relative tolerances for the solver.

# Returns
Matrix with columns `[t, x, y, ux, uy, u]` where `u` is the total velocity magnitude.
"""
function simulate(
    dem::DEM,
    params::MPMParams,
    x0::Real,
    y0::Real;
    ux0::Real = 0.0,
    uy0::Real = 0.0,
    dt::Real = 1.0,
    tend::Real = 300.0,
    t0::Real = 0.0,
    curvature::Bool = false,
    atol::Float64 = 1e-6,
    rtol::Float64 = 1e-6,
)
    alpha0 = [float(x0), float(y0), float(ux0), float(uy0)]

    f(alpha, p, t) = MPM_rhs(alpha, p, dem; curvature)

    prob = ODEProblem(f, alpha0, (float(t0), float(tend)), params)
    sol = solve(prob, Tsit5(); saveat = dt, abstol = atol, reltol = rtol)

    n = length(sol.t)
    output = Matrix{Float64}(undef, n, 6)

    for i = 1:n
        alpha = sol.u[i]
        gradx, grady, = _dem_derivatives(dem, alpha[1], alpha[2])
        u = _velocity_magnitude(alpha[3], alpha[4], gradx, grady)
        output[i, :] .= (sol.t[i], alpha[1], alpha[2], alpha[3], alpha[4], u)
    end

    return output
end

end
