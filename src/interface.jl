
function getX(p)
    p[1]
end
function getY(p)
    p[2]
end 
function point(::Type{FloatType}, points, i::Integer) where FloatType 
    return FloatType(getX(points[i])), FloatType(getY(points[i]))
end 
function rawpoint(points, i::Integer)
    return getX(points[i]),getY(points[i])
end 

abstract type AbstractDelaunatorData end

struct BasicTriangulation{IntType,FloatType,PointsType} <: AbstractDelaunatorData
    triangles::Vector{Tuple{IntType,IntType,IntType}}
    halfedges::Vector{IntType}
    points::PointsType
    _triangles::Vector{IntType}
    _minxy::Tuple{FloatType,FloatType}
    _maxxy::Tuple{FloatType,FloatType}
end 

""" Create a function to simplify type management. """
function _basictriangulation(triangles, halfedges, points, _minxy, _maxxy)
    IntType = eltype(halfedges)
    FloatType = eltype(_minxy)
    PointsType = typeof(points)
    ntris = length(triangles)
    _ptridata = Base.unsafe_convert(Ptr{IntType}, triangles)
    _triangles = unsafe_wrap(Array, _ptridata, 3*ntris)
    return BasicTriangulation{IntType,FloatType,PointsType}( 
        triangles, halfedges, points, _triangles, _minxy, _maxxy)
end 


function Base.show(io::IO, t::BasicTriangulation)
    T = length(t.triangles)
    N = length(t.points)
    print(io, "BasicTriangulation: ($T triangles, $N points)")
end

struct Triangulation{IntType,FloatType,PointsType} <: AbstractDelaunatorData
    triangles::Vector{Tuple{IntType,IntType,IntType}}
    halfedges::Vector{IntType}
    points::PointsType
    _triangles::Vector{IntType}

    _minxy::Tuple{FloatType,FloatType}
    _maxxy::Tuple{FloatType,FloatType}
    hull::Vector{IntType}
    hullindex::Vector{IntType}
    index::Vector{IntType} 
    circumcenters::Vector{Tuple{FloatType,FloatType}}
    raystart::Vector{Tuple{FloatType,FloatType}}
    rayend::Vector{Tuple{FloatType,FloatType}}
end 

function _triangulation(triangles, halfedges, points, 
            minxy, maxxy, hull, hullindex, index, circumcenters, 
            raystart, rayend)
    IntType = eltype(halfedges)
    FloatType = eltype(minxy)
    PointsType = typeof(points)
    ntris = length(triangles)
    _ptridata = Base.unsafe_convert(Ptr{IntType}, triangles)
    _triangles = unsafe_wrap(Array, _ptridata, 3*ntris)
    return Triangulation{IntType,FloatType,PointsType}(triangles, halfedges,
        points, _triangles, minxy, maxxy, hull, hullindex, index, circumcenters, 
        raystart, rayend)
end             

"""    
    triangles(t::AbstractDelaunatorData)

Return the point indices for each triangle in the triangulation.     
"""
function triangles(t::AbstractDelaunatorData)
    t.triangles
end

"""    
    points(t::AbstractDelaunatorData)

Return the point coordinates that were given as input to the algorithm. 
Note that changing these does not dynamically change the triangulation. 
"""
function points(t::AbstractDelaunatorData)
    t.points
end

import Base.eachindex 
"""    
    eachindex(t::AbstractDelaunatorData)

Return the indicies of each point in the dataset.
"""
function eachindex(t::AbstractDelaunatorData)
    1:length(points(t))
end

"""
    inhull(t::Triangulation, i::Integer)

Return the index of the vertex in the list of hull vertices (t.hull) 
if it's in the hull, or zero if the vertex is not in the hull. 
"""
function inhull(t::Triangulation, i::Integer)
    return t.hullindex[i]
end 

"""    
    minpt, maxpt = bounds(t::AbstractDelaunatorData)

Return the coordinate bounds on the points. All of the points (as of the computation of the algorithm)
lie within minpt <= pt <= maxpt. 
"""    
function bounds(t::AbstractDelaunatorData)
    (t._minxy,t._maxxy)
end 

inttype(t::AbstractDelaunatorData) = eltype(eltype(triangles(t)))
floattype(t::AbstractDelaunatorData) = eltype(eltype(bounds(t)))

function Base.show(io::IO, t::AbstractDelaunatorData)
    T = length(triangles(t))
    N = length(t.points)
    print(io, "$N-point triangulation [($(t._minxy[1]),$(t._minxy[2])) - ($(t._maxxy[1]),$(t._maxxy[2]))] $T triangles ($(inttype(t)), $(floattype(t))), ")
end

function Base.show(io::IO, ::MIME"text/plain", t::AbstractDelaunatorData)
    println(io, t)
    tris = triangles(t) 
    T = length(tris)
    I, J = T > 10 ? (5, T-4) : (T, T+1)
    lines = [["  └─$(tris[i])" for i in 1:I]
            (T > 10 ? ["  ⋮"] : [])
            ["  └─$(tris[i])" for i in J:T]]
    print(io, join(lines, "\n"))
end



struct TriangulationTemporaries{IntType,FloatType}
    hullPrev::Vector{IntType} # edge to prev edge
    hullNext::Vector{IntType} # edge to next edge 
    hullTri::Vector{IntType} # edge to adjacent triangle
    hullHash::Vector{IntType} # angular edge hash
    edgeStack::Vector{IntType}
    ids::Vector{IntType}
    dists::Vector{FloatType}
end 

function Base.show(io::IO, t::TriangulationTemporaries)
    n = length(t.hullPrev) 
    print(io, "TriangulationTemporaries ($(eltype(t.hullPrev)), $(eltype(t.dists))) for $n points")
end


""" Create a wrapper to simplify type management. """ 
function _temporaries(hullStart, hullSize, hullPrev, hullNext, hullTri, hullHash, edgeStack, ids, dists)
    edgeStack[1] = hullStart # save the location, a bit hacky, <sigh>
    edgeStack[2] = hullSize 
    IntType = eltype(hullStart)
    FloatType = eltype(dists)
    return TriangulationTemporaries{IntType,FloatType}(
        hullPrev, hullNext, hullTri, hullHash, edgeStack, ids, dists)
end 

function _allocate_cdata(::Type{IntType}, ::Type{FloatType}, npoints::Integer) where {IntType,FloatType}
    n = npoints

    # temporary arrays for tracking the edges of the advancing convex hull
    hullPrev = Vector{IntType}(undef,n) # edge to prev edge
    hullNext = Vector{IntType}(undef,n)# edge to next edge
    hullTri = Vector{IntType}(undef,n) # edge to adjacent triangle
    hullHash = Vector{IntType}(undef,n) # angular edge hash
    edgeStack = Vector{IntType}(undef,500)

    # temporary arrays for sorting points
    ids = Vector{IntType}(undef,n)
    dists = Vector{FloatType}(undef,n)

    return _temporaries(0, 0, hullPrev, hullNext, hullTri, hullHash, edgeStack, ids, dists)
end




"""
    bt, cdata = basictriangulation(IntType=Int32,] [FloatType=Float64,] [points];[sizehint=length(points),] [tol])

Allocate a basic triangulation structure and associated compute data in order
to implement the Deluantor method. 

**This method does not actaully compute a triangulation, but only allocates
the data. See `triangulate` or `update!` for the computational methods.**

The triangulation data structure 
----------------
These data structures are explained at https://mapbox.github.io/delaunator/ 
(but here, all the indices have been modified a little).

- `triangles`: a length T array of 3 tuples, where each tuple is a triangle
- `halfedges`: The halfedge index for the edge in the triangle array. The halfedges for 
    triangle t are 3(t-1)+1, 3(t-1)+2, 3(t-2)+3. Each halfedge index gives the entry
    of the other halfedge. 
- `points` a copy of the input set of points
"""
basictriangulation(; kwargs...) = basictriangulation([]; kwargs...)
basictriangulation(points; kwargs...) = basictriangulation(Int32, points; kwargs...)
basictriangulation(::Type{IntType}, points; kwargs...) where IntType = basictriangulation(Int32, Float64, points; kwargs...)
function basictriangulation(::Type{IntType}, ::Type{FloatType}, points;
    sizehint::Integer=length(points)) where {IntType,FloatType}
    
    npts = max(sizehint, length(points), 128)

    maxTriangles = max(2 * npts - 5, 0)
    tridata =  Vector{Tuple{IntType,IntType,IntType}}(undef,maxTriangles)
    halfedges =  Vector{IntType}(undef,3*maxTriangles)

    cdata = _allocate_cdata(IntType, FloatType,npts)
    bt = _basictriangulation(tridata, halfedges, points, 
            (zero(FloatType),zero(FloatType)), (one(FloatType),one(FloatType)))
    return bt, cdata 
end 

"""
    bt, cdata = update!(points, bt, cdata)

Reuse the memory and arrays to update the triangulation. Note that the resulting
return values may be new as this routine can allocate new arrays if needed to handle
the updated set of points. This implements the key step of the Delaunator method. 
"""
function update!(points, bt::BasicTriangulation, cdata::TriangulationTemporaries; tol=eps(floattype(bt)))
    IntType = inttype(bt)
    FloatType = floattype(bt) 
    n = length(points)
    maxTriangles = max(2 * n - 5, 0)
    # extract the arrays and check their size...

    minsize(v,n) = length(v) >= n ? v : resize!(v, n) 

    tridata = minsize(bt.triangles, maxTriangles)
    _ptridata = Base.unsafe_convert(Ptr{IntType}, tridata)
    triangles = unsafe_wrap(Array, _ptridata, 3*maxTriangles)

    halfedges = minsize(bt.halfedges, 3*maxTriangles)
    fill!(@view(halfedges[1:3*maxTriangles]), -1)

    # temporary arrays for tracking the edges of the advancing convex hull
    hashSize = ceil(Int,sqrt(n))
    hullPrev = minsize(cdata.hullPrev,n)
    hullNext = minsize(cdata.hullNext,n)
    hullTri = minsize(cdata.hullTri,n)
    hullHash = minsize(cdata.hullHash,n)
    edgeStack = minsize(cdata.edgeStack,500) # TODO Revisit 
    fill!(@view(hullHash[1:n]), -1)

    ids = minsize(cdata.ids, n)
    dists = minsize(cdata.dists, n)

    # initialize all the ids
    copyto!(@view(ids[1:n]), 1:n)

    minxy, maxxy = _get_bounds(FloatType, points)
    cxy = (getX(minxy)+getX(maxxy))/2, (getY(minxy)+getY(maxxy))/2

    collinear, i0, i1, i2 = _find_seeds(points, cxy)

    if collinear
        # order collinear points by dx (or dy if all x are identical)
        # and return the list as a hull
        for i in eachindex(points)
            cxv = (point(FloatType,points,i)[1] - point(FloatType,points,1)[1])
            dists[i] = !iszero(cxv) ? cxv : (point(FloatType,points,i)[2] - point(FloatType,points,1)[2])
        end
        quicksort(ids, dists, 1, n)
        resize!(tridata, 0)
        resize!(halfedges, 0)
        hullStart = 1 
        hullSize = n
    else
        hullStart, hullSize, ntriangles = _delaunator!(
            points,
            (i0,i1,i2),
            triangles,
            halfedges,
            hullPrev, 
            hullNext,
            hullTri,
            hullHash,
            hashSize,
            edgeStack,
            ids, 
            dists,
            tol
        )
        # resize to the appropriate size 
        resize!(tridata, ntriangles)
        resize!(halfedges, 3*ntriangles)
    end 

    # put everything back together :) 
    return (_basictriangulation(tridata, halfedges, points, minxy, maxxy), 
            _temporaries(hullStart, hullSize, hullPrev, hullNext, hullTri, 
                hullHash, edgeStack, ids, dists))
end 

""" 
    hullvertices!(hull, bt, cdata)
After an [`update!`](@ref) call, we can use the data structures
to get a simple list of vertices on the convex hull. This routine
will run `push!(hull, v)` for each vertex on the convex hull. 
The order will be in order around the hull. This returns
the hull variable. 

This is part of the advanced interface. See 
[`triangulate`](@ref) which returns a more complete data
structure including the hull information automatically. 

## Sample calls
```julia-repl
julia> using Delaunator, StableRNGs, GeometryBasics
julia> bt, cdata = basictriangulation(pts)
julia> hull = hullvertices!(Vector{Int}, bt, cdata )
```
"""
function hullvertices!(hull, bt::BasicTriangulation, cdata::TriangulationTemporaries)
    n = length(bt.points) 
    if length(triangles(bt)) == 0 
        # this is a collinear array
        # the ids array should be sorted by distance 
        # from the first point. We want to add 
        # any point that isn't a duplicate.. 
        _dists = cdata.dists
        _ids = cdata.ids
        d0 = -Inf
        for i in 1:n
            id = _ids[i]
            if _dists[id] > d0 # if we have strictly large distance 
                push!(hull, id) 
                d0 = _dists[id]
            end
        end
        #hull = resize!(hull, j-1)
    else 
        hullNext = cdata.hullNext
        e = cdata.edgeStack[1] # hull start is edgeStack[1] 
        hullSize = cdata.edgeStack[2] # hull Size is edgeStack[2] 
        for i = 1:hullSize
            push!(hull, e)
            #hull[i] = e
            e = hullNext[e]
        end
    end 
    return hull 
end 


""" Store the first time a point occurs in halfedges into the index. """
function index_halfedges!(index, halfedges, triangles; init=true)
    if init 
        fill!(index, 0)
    end
    for i in eachindex(triangles)
        # this gives priority to exterior half-edges... 
        endpoint = triangles[(i-1)%3 == 2 ? i-2 : i+1] # this is the next halfedge
        if halfedges[i] == -1 || index[endpoint] == 0 
            index[endpoint] = i 
        end
    end 
    return index
end 

#=
bt, cdata = basictriangulation(points; [maxpoints=Integer]) # initialize data structures 
bt, cdata = update!(bt, points, cdata) # after the points have been changed, may incur allocations
h = gethull(bt, cdata)
index = index_halfedges(bt, cdata)
circumcenters!(cc, bt) # compute circumcenters 
=#

""" Get the largest distance along x or y dimension. """
function _max_dimension(bt)
    minxy,maxxy = bounds(bt) 
    return max(maxxy[1]-minxy[2], maxxy[2]-minxy[2])
end 

"""
    circumcenters!(array, bt; [collinearthresh=0])

Write information on the circumcenters into array. Array must 
have length(bt.triangles) allocated. `collinearthresh` is an 
option that 
"""
function circumcenters!(array, t::AbstractDelaunatorData; 
    #collineartol=_max_dimension(t)/sqrt(eps(eltype(array)))
    collineartol=0
    )
    points = t.points
    tris = triangles(t) 
    FloatType = eltype(eltype(array))
    rx,ry = begin 
        if length(tris) > 0 
            point(FloatType, points, tris[1][1]) 
        else
            zero(FloatType), zero(FloatType)
        end
    end
    
    for i in eachindex(tris)
        t1,t2,t3 = tris[i]
        x1,y1 = point(FloatType, points, t1)
        x2,y2 = point(FloatType, points, t2)
        x3,y3 = point(FloatType, points, t3)
        array[i] = truncatedcircumcenter(x1,y1,x2,y2,x3,y3,rx,ry,collineartol)
    end 
    return array
end 

"""
    raystart, rayend = rays(t)

This returns arrays that are indexed by the _index of_ a point on the convex hull.
So to get the infinite rays associated with the nearest point cell use:
```
function rays_for_point(t)
    rs, re = rays(t) # only need to compute this one 
    hullindex = inhull(t, i)
    if hullindex > 0 
        return rs[hullindex],re[hullindex]
    else
        return 
    end 
end 
```    
"""
function rays(::Type{FloatType}, hull, points) where FloatType
    #d3-delaunay/Voronoi.js code 
    #=
    let h = hull[hull.length - 1];
      let p0, p1 = h * 4;
      let x0, x1 = points[2 * h];
      let y0, y1 = points[2 * h + 1];
      vectors.fill(0);
      for (let i = 0; i < hull.length; ++i) {
        h = hull[i];
        p0 = p1, x0 = x1, y0 = y1;
        p1 = h * 4, x1 = points[2 * h], y1 = points[2 * h + 1];
        vectors[p0 + 2] = vectors[p1] = y0 - y1;
        vectors[p0 + 3] = vectors[p1 + 1] = x1 - x0;
      }
    =#
    nhull = length(hull) 
    raystart = Vector{Tuple{FloatType,FloatType}}(undef, nhull)
    rayend = Vector{Tuple{FloatType,FloatType}}(undef, nhull)
    p1 = hull[end]
    hi1 = nhull # index of last point 
    x1,y1 = points[p1] # last point in the hull 
    for (hi,h) in enumerate(hull)  # walk through the hull in order now..
        hi0,x0,y0 = hi1,x1,y1
        hi1 = hi  
        x1,y1 = points[h]
        # this is the outward orthogonal direction to the line on the hull. 
        rayend[hi0] = raystart[hi1] = (y0 - y1, x1 - x0) 
    end 
    return raystart, rayend 
end 

"""
    triangulate([Int32,] [FloatType=Float64,] points; [tol=eps(FloatType),])

Computes a triangulation of a set of points using the Delaunator algorithm.
This is designed for quick graphics applications and speed rather than exact computational geometry.

Inputs
------
- `points` is any type that has integer indexing and length supported. In addition, `p = points[i]` should 
    be a type where `p[1]` and `p[2]` are the x, y coordinates of p. Or you need to define the functions 
    `Delaunator.getX(p), Delaunator.getY(p)` for your own type p. 

    If you wish to use a a matrix to give the point information, [`PointsFromMatrix`](@ref)

- `tol` is used to determine when points are sufficiently close not to include

Return value
------------
A triangulation, with methods to explore edges, hull points, dual cells.

See also [`basictriangulation`](@ref)
"""
triangulate(points; kwargs...) = triangulate(Float64, points; kwargs... )
triangulate(::Type{FloatType}, points; kwargs...) where FloatType = triangulate(Int32, FloatType, points; kwargs...)
function triangulate(::Type{IntType}, ::Type{FloatType}, 
            points; tol=eps(FloatType), rcollineartol=0,
            ) where {IntType <: Signed, FloatType}
    n = length(points) 
    bt, cdata = basictriangulation(IntType, FloatType, points)
    bt, cdata = update!(points, bt, cdata; tol)

    hull = Vector{IntType}()
    hullvertices!(hull, bt, cdata)
    hullindex = fill!(Vector{IntType}(undef, n), 0)
    for (i,h) in enumerate(hull)
        hullindex[h] = i
    end 

    index = index_halfedges!(Vector{IntType}(undef, n), bt.halfedges, bt._triangles)

    ntris = length(bt.triangles)
    ccs = Vector{Tuple{FloatType,FloatType}}(undef, ntris)
    collineartol = _max_dimension(bt)*rcollineartol
    circumcenters!(ccs, bt; collineartol)

    raystart,rayend = rays(FloatType, hull, points) 
    
    _triangulation(triangles(bt), bt.halfedges, points, 
        bt._minxy, bt._maxxy, hull, hullindex, index, ccs, raystart, rayend)
end 


