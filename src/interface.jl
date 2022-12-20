
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
    print(io, "$T triangles, $N points")
end

struct Triangulation{IntType,FloatType,PointsType} <: AbstractDelaunatorData
    triangles::Vector{Tuple{IntType,IntType,IntType}}
    halfedges::Vector{IntType}
    points::PointsType
    _triangles::Vector{IntType}

    _minxy::Tuple{FloatType,FloatType}
    _maxxy::Tuple{FloatType,FloatType}
    hull::Vector{IntType}
    index::Vector{IntType} 
    circumcenters::Vector{Tuple{FloatType,FloatType}}
    rays::Vector{Vector{Tuple{FloatType,FloatType}}}
end 

function _triangulation(triangles, halfedges, points, 
            minxy, maxxy, hull, index, circumcenters, rays)
    IntType = eltype(halfedges)
    FloatType = eltype(_minxy)
    PointsType = typeof(points)
    ntris = length(triangles)
    _ptridata = Base.unsafe_convert(Ptr{IntType}, triangles)
    _triangles = unsafe_wrap(Array, _ptridata, 3*ntris)
    return Triangulation{IntType,FloatType,PointsType}(triangles, halfedges,
        points, _triangles, minxy, maxxy, hull, index, circumcenters, rays)
end             

"""    
    triangles(t::AbstractDelaunatorData)

Return the point indices for each triangle in the triangulation.     
"""
function triangles(t::AbstractDelaunatorData)
    t.triangles
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


""" Create the right returntype from a delaunator call. """
function _returntype(trituples, halfedges, hull, points, index)
    PointsType = typeof(points) 
    IntType = eltype(halfedges)

    # rewrap with a flat pointer
    ntris = length(trituples)
    _ptridata = Base.unsafe_convert(Ptr{IntType}, trituples)
    _triangles = unsafe_wrap(Array, _ptridata, 3*ntris)

    return Triangulation{IntType, Float64, PointsType}(
            trituples, halfedges,  points,  _triangles, (0.0,0.0), (1.0,1.0), hull, index, 
            Tuple{Float64,Float64}[], Tuple{Float64,Float64}[])
end 

function Base.show(io::IO, t::AbstractDelaunatorData)
    T = length(triangles(t))
    N = length(points(t))
    print(io, "$T triangles ($(inttype(t)), $(floattype(t))), $N points")
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
"""
basictriangulation(; kwargs...) = basictriangulation([]; kwags...)
basictriangulation(points; kwargs...) = basictriangulation(Int32, points; kwargs...)
basictriangulation(::Type{IntType}, points; kwargs...) where IntType = basictriangulation(Int32, Float64, points; kwargs...)
function basictriangulation(::Type{IntType}, ::Type{FloatType}, points;
    sizehint::Integer=length(points)) where {IntType,FloatType}
    
    npts = max(sizehint, length(points), 128)

    maxTriangles = max(2 * npts - 5, 0)
    tridata =  Vector{Tuple{IntType,IntType,IntType}}(undef,maxTriangles)
    halfedges =  Vector{IntType}(undef,3*maxTriangles)

    cdata = _allocate_cdata(IntType, FloatType,npts)
    bt = _basictriangulation(tridatay, halfedges, points, 
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
    maxTriangles = max(2 * npts - 5, 0)
    # extract the arrays and check their size...

    minsize(v,n) = length(v) <= n ? v : resize!(v, n) 

    tridata = minsize(bt.triangles, maxTriangles)
    _ptridata = Base.unsafe_convert(Ptr{IntType}, tridata)
    triangles = unsafe_wrap(Array, _ptridata, 3*maxTriangles)

    halfedges = minsize(bt.halfedges, 3*maxTriangles)
    fill!(@view(_halfedges[1:3*maxTriangles]), -1)

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
        for i in eachindex(coords)
            cxv = (point(FloatType,coords,i)[1] - point(FloatType,coords,1)[1])
            _dists[i] = !iszero(cxv) ? cxv : (point(FloatType,coords,i)[2] - point(FloatType,coords,1)[2])
        end
        quicksort(_ids, _dists, 1, n)
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
        d0 = -Inf
        for i in 1:n
            id = _ids[i]
            if _dists[id] > d0 # if we have strictly large distance 
                push!(hull, id) 
                hull[j] = id
                j += 1
                d0 = _dists[id]
            end
        end
        hull = resize!(hull, j-1)
    else 
        e = cdata.edgeStack[1] # hull start is edgeStack[1] 
        hullSize = cdata.edgeStack[2] # hull Size is edgeStack[2] 
        for i = 1:hullSize
            push!(hull, e)
            hull[i] = e
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
    circumcenters!(array, bt; [collinearthresh=sqrt(eps(eltype(array)))])

Write information on the circumcenters into array. Array must 
have length(bt.triangles) allocated. `collinearthresh` is an 
option that 
"""
function circumcenters!(array, t::AbstractDelaunatorData; 
    collineartol=_max_dimension(t)/sqrt(eps(eltype(array))))
    points = t.points
    tris = triangles(t) 
    FloatType = eltype(array)
    rx,ry = point(FloatType, points, tris[1][1])
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
    triangulate([Int32,] [FloatType=Float64,] points; [tol=eps(FloatType),] 
        [rcollineartol=sqrt(eps(FloatType))])
"""
triangulate(points; kwargs...) = delaunator!(Float64, points; kwargs... )
triangulate(::Type{FloatType}, points) where FloatType = triangulate(Int32, FloatType, points; kwargs...)
function triangulate(::Type{IntType}, ::Type{FloatType}, 
            points; tol=eps(FloatType), rcollineartol=sqrt(eps(FloatType))
            ) where {IntType <: Signed, FloatType}
    n = length(points) 
    bt, cdata = basictriangulation(IntType, FloatType, points)
    bt, cdata = update!(points, bt, cdata; tol)
    hull = Vector{IntType}()
    hullvertices!(hull, bt, cdata)
    index = index_halfedges!(Vector{IntType}(undef, n), bt.halfedges, bt._triangles)
    ntris = length(bt.triangles)
    ccs = Vector{Tuple{FloatType,FloatType}}(undef, ntris)
    collinearthresh=_max_dimension(t)*collineartol
    circumcenters!(ccs, bt; collinearthresh)
    # TODO fix the rays 
    rays = Vector{Tuple{FloatType,FloatType}}(undef, 0) 
    _triangulation(triangles(bt), bt.halfedges, points, 
        bt._minxy, bt._maxxy, hull, index, ccs, rays)
end 