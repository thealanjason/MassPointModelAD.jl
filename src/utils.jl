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
