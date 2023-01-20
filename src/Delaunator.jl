module Delaunator

include("interface.jl")
export triangulate, basictriangulation, update!, triangles, points, inhull

include("quicksort.jl")
include("geometry.jl")
export isinfinite, dualcell, firstpoint, lastpoint, segments

include("clipping.jl")
export clippedpoly, clippedpoly!, margin_bbox

include("algorithm.jl") # nothing is exported here... 

include("iterators.jl")
export neighbors, edges, edgelines, hullpoly

include("helpers.jl")
export PointsFromMatrix

end # module
