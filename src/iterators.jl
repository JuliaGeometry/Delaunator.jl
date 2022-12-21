
# code to work with the voronoi / nearest point tesselation 


function _nextedge(i::Integer)
    (i-1)%3 == 2 ? i-2 : i+1
  end 

import Base.iterate, Base.eltype, Base.IteratorSize
struct PointNeighborIterator{TriType,IntType}
    t::TriType
    start::IntType
end
Base.IteratorSize(::Type{PointNeighborIterator{TriType,IntType}}) where {TriType,IntType} = Base.SizeUnknown()
Base.eltype(::Type{PointNeighborIterator{TriType,IntType}}) where {TriType,IntType} = IntType
Base.isdone(it::PointNeighborIterator, state) = state == it.start 
function Base.iterate(it::PointNeighborIterator, state=nothing) 
    # we are done if state == it.start again, after the first iteration
    if state == it.start 
        return nothing
    elseif state !== nothing && state < 0 
        return (edge_from_index(it.t, -state)[2], it.start)
    else
        if state === nothing 
            state = it.start
        end 
        outgoing = _nextedge(state)
        incoming = it.t.halfedges[outgoing]
        nextstate = incoming 
        if nextstate == -1
            nextstate = -outgoing
        end
        return (edge_from_index(it.t, state)[1], nextstate)
    end 
end

function triangle_from_index(t::Triangulation, i::Integer)
    return ((i-1) ÷ 3)+1 
end

function edge_from_index(t::Triangulation, i::Integer)
    ti, vi = divrem(i-1, 3)
    #tj, vj = divrem(t.halfedges[i]-1, 3)
    #return t.triangles[ti+1][vi+1], t.triangles[tj+1][vj+1]
    return t.triangles[ti+1][vi+1], t.triangles[ti+1][(vi+1) % 3 == 0 ? 1 : vi + 2]
end
  

function neighbors(t::Triangulation, i::Integer)
    return PointNeighborIterator{typeof(t),eltype(t.halfedges)}(t, t.index[i])
end 

function edges(t::Triangulation)
    return ((i,j) for i in 1:length(t.points) for j in neighbors(t,i) if (i<j))
end 

"""
    triangles(t, i)

Given a Triangulation and a point index `i`, return an iterator
over the indices of triangles that include point i. These are returned
in counter-clockwise order. 
TODO Make this an iterator... (right now, it's just an array...)
"""
function triangles(t::Triangulation, i::Integer)
    rval = inttype(t)[] 
    start = t.index[i]
    state = nothing
    while state != start && state != -1 
        if state === nothing
            state = start
        end 
        outgoing = _nextedge(state)
        incoming = t.halfedges[outgoing]
        nextstate = incoming 
        push!(rval, triangle_from_index(t, state))
        state = nextstate 
    end 
    return rval
end 

"""
    edgelines(t::Triangulation)

Return a generator that can be used with Makie's linesegments function to 
display the edges of the triangulation. Each edge is only drawn once.

# Example
```
using GLMakie
t = triangulate(rand(StableRNG(1), Point2f, 10))
linesegments(collect(edgelines(rval)))
```
"""
function edgelines(t::Triangulation)
    return ((rawpoint(t.points, e[1]), rawpoint(t.points,e[2])) for e in edges(t) )
end 

function hull(t::Triangulation)
    return t.hull
end 

"""
    hullpoly(t)

Return the coordinates of the convex hull suitable for plotting as a polygon. 

# Example
```julia-repl
julia> t = triangulate(rand(StableRNG(1), Point2f, 10))
julia> f = scatter(t.points)
julia> poly!(f.axis,collect(hullpoly(rval)),color=:transparent, strokewidth=1)
julia> f
```
"""
function hullpoly(t::Triangulation)
    return (rawpoint(t.points, p) for p in hull(t))
end 


