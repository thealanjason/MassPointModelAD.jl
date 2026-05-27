# Test case from §3 of the model PDF.
# The piecewise DEM (slope + flat) has a discontinuity at xB, so we split the
# simulation into two smooth phases and compare each against the analytical
# solution (eqs. 28–30). No curvature, no turbulent friction (xi=Inf).

const h = 100.0   # initial height (m)
const theta = pi / 4  # slope angle (45°)
const mu = 0.3     # Coulomb friction coefficient
const g = 9.8     # gravitational acceleration (m/s²)

# Analytical values (eq. 28)
const xB = h / tan(theta)
const a_AB = g * sin(theta) - mu * g * cos(theta)
const U_TxB = sqrt(2 * a_AB * h / sin(theta))
const t_B = U_TxB / a_AB
const x_C = xB + U_TxB^2 / (2 * mu * g)

const params = MPMParams(mu, Inf, g)
Z_slope(x, y) = h - x * tan(theta)
Z_flat(x, y) = 0.0

@testset "Phase 1 — slope A→B" begin
    output = simulate(Z_slope, params, 0.0, 0.0; dt = 0.01, tend = t_B, curvature = false)

    @test output[end, 4] ≈ U_TxB rtol=1e-10   # ux at t_B (eq. 29)
    @test output[end, 2] ≈ xB rtol=1e-10   # x  at t_B (eq. 30)
end

@testset "Phase 1 + 2 — full trajectory A→B→C" begin
    # Phase 1: slope
    output1 = simulate(Z_slope, params, 0.0, 0.0; dt = 0.01, tend = t_B, curvature = false)

    x_end1 = output1[end, 2]
    ux_end1 = output1[end, 4]

    # Phase 2: flat, seeded with numerical end-state of phase 1
    tend2 = ux_end1 / (mu * g)
    output2 = simulate(
        Z_flat,
        params,
        x_end1,
        0.0;
        ux0 = ux_end1,
        dt = 0.01,
        tend = tend2,
        curvature = false,
    )

    @test output2[end, 2] ≈ x_C rtol=1e-10   # x at t_C (eq. 28)
end
