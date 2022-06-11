# Author: Christina G Taylor
# Date:   April 2022
# Description:
#   Code for identifying intersections between a multi-dimensional Cartesian 
#   grid and parameterized paths. Intersection points are solved using a 
#   bracketted secant method. 
#
# Missing Functionality:
#  1) Non-simple intersections are not yet fully supported; some may be
#     captured, others may not depending on the step size (ds) used.
#
#  2) Lines of intersection will return multiple points of intersection
#     all ds distance from one another. The endpoints may not be returned.
#
#  3) Curves with sufficiently varying arc length per step may result in
#     missing intersections
#
# Inputs/User-defined objects:
#   coords                : An array of arrays of coordiantes. Each subarray
#                           contains the Cartesian mesh coordinates in that
#                           dimension. Each subarray must be in increasing
#                           order and contain the endpoints of the domain.
#                            Indexing: coords[dim][coord index]
#   curves                : An array of functions. Each function parameterizes a
#                           single object. Functions are called using:
#                            pt = curve[curve index](s, params[curve index])
#   ds                    : Array of step sizes to use when searching each curve.
#                           Must have same dimension as 'curves'.
#   arc_tol                : Array of tolerances used to find intersection points.
#                           The Euclidean distance from the approximate intersection
#                           point and the true intersection point is guaranteed to
#                           be less than this value. Must have same dimension as 
#                           'curves'.
# corner_tol               : The tolerance used to determine if a given dimension
#                           participated in a previously found intersection.
#
# Outputs/Results:
#   intersections         : A vector of (s,pt,dim,I) tuples describing where the
#                           curves intersect the mesh where
#                               s = the approximate parameter value where the curve
#                                   intersects the mesh
#                              pt = the point returned by the curve at s
#                             dim = tuple of booleans indicating which dimension(s)
#                                   the intersection occured in where
#                                      true = that dim participated in the intersection
#                                     false = that dim did not participate
#                               I = a tuple of indices into the coordinate
#                                   arrays corresponding to the intersected value
#                                   (if that dimension participated in the 
#                                   intersection) or the closest coord value less
#                                   than the intersection point
#
# #########################################################################################



module PathIntersections

using ForwardDiff
using LinearAlgebra
using StructArrays

const DEFAULT_DS = 0.0
const DEFAULT_ARC_TOL = 1e-12
const DEFAULT_CORNER_TOL = 1e-12
const DEFAULT_CONTINUITY_TOL = 1e-12
const DEFAULT_REF_DOMAIN = (-1,1)
const DEFAULT_NORMALIZATION = false
const DEFAULT_BINARY_REGIONS = false

include("helper_functions.jl")

# For defining geometeries
export ConstFunction
export PiecewiseFunction
export PiecewiseCurve
export ds_by_num_steps
export ds_by_arc_length
include("PiecewiseCurve.jl")

# For preset geometeries
export PresetGeometries
include("PresetGeometries.jl")

# For finding mesh-curve intersections
export MeshIntersection
mutable struct MeshIntersection{T_param, T_point, T_participants, T_indices}
    s::T_param              # The s-value where the intersection occurred
    pt::T_point             # The approximate intersection point
    dim::T_participants     # Whether each dimension was involved in the intersection
    indices::T_indices            # The indices of the closest lower mesh bound or the       
                            # boundaries involved in the intersection
end

export find_mesh_intersections
include("find_mesh_intersections.jl")

# For mapping quadrature rules to curves
export map_line_quadrature
export map_line_quadrature_nonconsecutive
include("map_line_quadrature.jl")

export define_regions
export get_cutcell_nodes
include("define_regions.jl")

export is_contained
include("is_contained.jl")

end # module