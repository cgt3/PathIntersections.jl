using PathIntersections.PresetGeometries
using PathIntersections
using SparseArrays

nx, ny = 10, 10
square_size = minimum([1/(2*nx), 1/(2*ny)])
x_coords = LinRange(-1, 1, nx)
y_coords = LinRange(-1, 1, ny)

i_mid = trunc(Int, nx / 2)
j_mid = trunc(Int, ny / 2)


region = 1
elements_by_region = [ones(nx-1, ny-1), spzeros(nx-1, ny-1)]
elements_by_region[region][i_mid, j_mid] = 0
elements_by_region[region][i_mid+1, j_mid] = 0
elements_by_region[region][i_mid+1, j_mid+1] = 0
elements_by_region[region][i_mid, j_mid+1] = 0


cutcells = PiecewiseCurve[]
sub_bounds = [(0,1), (0,1), (0,1), (0,1), (0,1), (0,1)]

# 1
stop_pts = [0, 0.25, 0.375, 0.5, 0.625, 0.75, 1]
subcurves = [
    Line((x_coords[i_mid], y_coords[j_mid]),   (x_coords[i_mid], y_coords[j_mid+1])),
    Line((x_coords[i_mid], y_coords[j_mid+1]), (-square_size, y_coords[j_mid+1])),
    Line((-square_size, y_coords[j_mid+1]),    (-square_size, -square_size)),
    Line((-square_size, -square_size),         (x_coords[i_mid+1], -square_size)),
    Line((x_coords[i_mid+1], -square_size),    (x_coords[i_mid+1], y_coords[j_mid])),
    Line((x_coords[i_mid+1], y_coords[j_mid]), (x_coords[i_mid], y_coords[j_mid])),
]
# cutcells[region][i_mid,j_mid] = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
push!(cutcells, PiecewiseCurve(stop_pts, subcurves, sub_bounds))

# 2
stop_pts = [0, 0.125, 0.25, 0.375, 0.5, 0.75, 1]
subcurves = [
    Line((x_coords[i_mid+1], y_coords[j_mid]),   (x_coords[i_mid+1], -square_size)),
    Line((x_coords[i_mid+1], -square_size),      (square_size, -square_size)),
    Line((square_size, -square_size),            (square_size, y_coords[j_mid+1])),
    Line((square_size, y_coords[j_mid+1]),       (x_coords[i_mid+2], y_coords[j_mid+1])),
    Line((x_coords[i_mid+2], y_coords[j_mid+1]), (x_coords[i_mid+2], y_coords[j_mid])),
    Line((x_coords[i_mid+2], y_coords[j_mid]),   (x_coords[i_mid+1], y_coords[j_mid])),
]
# cutcells[region][i_mid+1,j_mid] = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
push!(cutcells, PiecewiseCurve(stop_pts, subcurves, sub_bounds))

# 3
stop_pts = [0, 0.125, 0.375, 0.625, 0.75, 0.875, 1]
subcurves = [
    Line((x_coords[i_mid+1], square_size),       (x_coords[i_mid+1], y_coords[j_mid+2])),
    Line((x_coords[i_mid+1], y_coords[j_mid+2]), (x_coords[i_mid+2], y_coords[j_mid+2])),
    Line((x_coords[i_mid+2], y_coords[j_mid+2]), (x_coords[i_mid+2], y_coords[j_mid+1])),
    Line((x_coords[i_mid+2], y_coords[j_mid+1]), (square_size, y_coords[j_mid+1])),
    Line((square_size, y_coords[j_mid+1]),       (square_size, square_size)),
    Line((square_size, square_size),             (x_coords[i_mid+1], square_size)),
]
# cutcells[region][i_mid+1,j_mid+1] = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
push!(cutcells, PiecewiseCurve(stop_pts, subcurves, sub_bounds))


# 4
stop_pts = [0, 0.25, 0.5, 0.625, 0.75, 0.875, 1]
subcurves = [
    Line((x_coords[i_mid], y_coords[j_mid+1]),    (x_coords[i_mid], y_coords[j_mid+2])),
    Line( (x_coords[i_mid], y_coords[j_mid+2]),   (x_coords[i_mid+1], y_coords[j_mid+2])),
    Line( (x_coords[i_mid+1], y_coords[j_mid+2]), (x_coords[i_mid+1], square_size) ),
    Line( (x_coords[i_mid+1], square_size),       (-square_size, square_size) ),
    Line( (-square_size, square_size),            (-square_size, y_coords[j_mid+1]) ),
    Line( (-square_size, y_coords[j_mid+1]),      (x_coords[i_mid], y_coords[j_mid+1]) ),
]
# cutcells[region][i_mid,j_mid+1] = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
push!(cutcells, PiecewiseCurve(stop_pts, subcurves, sub_bounds))