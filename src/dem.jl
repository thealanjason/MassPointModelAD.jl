struct DEM{I}
    x_coords::Vector{Float64}
    y_coords::Vector{Float64}
    fz::I
    fgrad_x::I
    fgrad_y::I
    fgrad_xx::I
    fgrad_xy::I
    fgrad_yy::I
end

function _np_gradient_1d(f::AbstractVector, h::Float64)
    n = length(f)
    g = Vector{Float64}(undef, n)
    g[1] = (-3f[1] + 4f[2] - f[3]) / (2h)
    for i = 2:(n-1)
        g[i] = (f[i+1] - f[i-1]) / (2h)
    end
    g[n] = (f[n-2] - 4f[n-1] + 3f[n]) / (2h)
    return g
end

function _gradient_along_x(Z::Matrix{Float64}, h::Float64)
    nx, ny = size(Z)
    G = Matrix{Float64}(undef, nx, ny)
    for j = 1:ny
        G[:, j] = _np_gradient_1d(view(Z, :, j), h)
    end
    return G
end

function _gradient_along_y(Z::Matrix{Float64}, h::Float64)
    nx, ny = size(Z)
    G = Matrix{Float64}(undef, nx, ny)
    for i = 1:nx
        G[i, :] = _np_gradient_1d(view(Z, i, :), h)
    end
    return G
end

function _build_interpolant(x_coords, y_coords, data::Matrix{Float64})
    itp = Interpolations.interpolate(
        (x_coords, y_coords),
        data,
        Interpolations.Gridded(Interpolations.Linear()),
    )
    return Interpolations.extrapolate(itp, Interpolations.Flat())
end

function _build_dem(
    x_coords::Vector{Float64},
    y_coords::Vector{Float64},
    z_grid::Matrix{Float64},
)
    dx = x_coords[2] - x_coords[1]
    dy = y_coords[2] - y_coords[1]
    grad_x = _gradient_along_x(z_grid, dx)
    grad_y = _gradient_along_y(z_grid, dy)
    grad_xx = _gradient_along_x(grad_x, dx)
    grad_xy = _gradient_along_y(grad_x, dy)
    grad_yy = _gradient_along_y(grad_y, dy)
    return DEM(
        x_coords,
        y_coords,
        _build_interpolant(x_coords, y_coords, z_grid),
        _build_interpolant(x_coords, y_coords, grad_x),
        _build_interpolant(x_coords, y_coords, grad_y),
        _build_interpolant(x_coords, y_coords, grad_xx),
        _build_interpolant(x_coords, y_coords, grad_xy),
        _build_interpolant(x_coords, y_coords, grad_yy),
    )
end

"""
    dem_from_func(Z, x_coords, y_coords) -> DEM

Build a `DEM` from a callable `Z(x, y)` evaluated on a regular grid.
Useful for testing with analytical surfaces without needing a raster file.
Derivatives are computed with central differences (np.gradient convention).
"""
function dem_from_func(Z, x_coords::AbstractVector, y_coords::AbstractVector)
    x = collect(Float64, x_coords)
    y = collect(Float64, y_coords)
    z_grid = Float64[Z(xi, yj) for xi in x, yj in y]
    return _build_dem(x, y, z_grid)
end

"""
    load_dem(filepath) -> DEM

Load an ESRI ASCII raster and return a `DEM` with six bilinear interpolants:
elevation Z and its five partial derivatives (grad_x, grad_y, grad_xx, grad_xy, grad_yy).
Derivatives are computed with central differences (np.gradient convention).
"""
function load_dem(filepath::AbstractString)
    open(filepath) do io
        ncols = parse(Int, split(readline(io))[2])
        nrows = parse(Int, split(readline(io))[2])
        xllcorner = parse(Float64, split(readline(io))[2])
        yllcorner = parse(Float64, split(readline(io))[2])
        cellsize = parse(Float64, split(readline(io))[2])
        readline(io)  # NODATA_value line (ignored)

        raw = Matrix{Float64}(undef, nrows, ncols)
        for r = 1:nrows
            raw[r, :] = parse.(Float64, split(readline(io)))
        end

        # ESRI stores rows north→south; flip so row 1 = min y
        raw = raw[end:-1:1, :]

        x_coords = Float64[xllcorner + (c - 1) * cellsize for c = 1:ncols]
        y_coords = Float64[yllcorner + (r - 1) * cellsize for r = 1:nrows]

        # z_grid[ix, iy] = Z at (x_coords[ix], y_coords[iy])
        z_grid = collect(transpose(raw))

        return _build_dem(x_coords, y_coords, z_grid)
    end
end

"""
    save_dem(dem, filepath; nodata_value=-9999.0)

Write a `DEM` to an ESRI ASCII raster file.  Only the elevation Z is stored;
derivatives are recomputed by `load_dem` when the file is read back.
"""
function save_dem(dem::DEM, filepath::AbstractString; nodata_value::Float64 = -9999.0)
    ncols = length(dem.x_coords)
    nrows = length(dem.y_coords)
    xllcorner = dem.x_coords[1]
    yllcorner = dem.y_coords[1]
    cellsize = dem.x_coords[2] - dem.x_coords[1]

    z_grid = [dem.fz(x, y) for x in dem.x_coords, y in dem.y_coords]

    open(filepath, "w") do io
        println(io, "ncols         $ncols")
        println(io, "nrows         $nrows")
        println(io, "xllcorner     $xllcorner")
        println(io, "yllcorner     $yllcorner")
        println(io, "cellsize      $cellsize")
        println(io, "NODATA_value  $nodata_value")
        for iy = nrows:-1:1
            println(io, join((z_grid[ix, iy] for ix = 1:ncols), " "))
        end
    end
end
