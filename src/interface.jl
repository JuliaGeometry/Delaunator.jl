

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
    print(io, "$T triangles ($(eltype(eltype(T)))), $N points")
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



"""
    bt, cdata = basictriangulation(IntType=Int32,] [FloatType=Float64,] [points];[maxpoints=length(points),] [tol])

Allocate a basic triangulation structure and associated compute data in order
to implement the Deluantor method. 
"""
function basictriangulation(points;maxpoints::Integer=length(points),tol=eps(FloatType))
end 
#=
bt, cdata = basictriangulation(points; [maxpoints=Integer]) # initialize data structures 
bt, cdata = update!(bt, points, cdata) # after the points have been changed, may incur allocations
h = gethull(bt, cdata)
index = index_halfedges(bt, cdata)
circumcenters!(cc, bt) # compute circumcenters 
=#