

"""
    InfinitePolygon

This is the type we use to represent dual cells / Voronoi cells.
It's a possibly infinite polygon stored as an array of points with
optional head and tail rays. 

Methods

- [`segments`](@ref) to get the line segments for the polygon.
"""    
struct InfinitePolygon{ArrayType,ElType}
    points::ArrayType
    head::Tuple{ElType,ElType}
    tail::Tuple{ElType,ElType}
end

import Base.eltype
eltype(::Type{InfinitePolygon{ArrayType,ElType}}) where {ArrayType, ElType} = Tuple{ElType,ElType}

function Base.show(io::IO, p::InfinitePolygon)
    N = length(p.points)
    extrastr = isinfinite(p) ? " with incoming ray $(p.head) and outgoing ray $(p.tail)" : ""
    print(io, "$N-point polygon $extrastr")
end

"""
    infinitepoly(pts, head, tail) 

Create a possibly unbounded polygon. The points should be listed in
counter-clockwise order. The head is the incoming halfray, the
tail is the departing halfray. These can be set to zero to 
create a finite polygon which is just a connected set of points.
Infinite polygons arise in the dual of a trianglation for those
edges on the convex hull.

p = infinitepoly(pts, head, tail)
isfinite(p) # return true if head/tail are both zero
"""

function infinitepoly(pts, head, tail)
    ElType = eltype(head) # this gets the float type... 
    ArrayType = typeof(pts) # we dont' care what this is and it could be an iterator... 
    return InfinitePolygon{ArrayType, ElType}(pts, head, tail)
end 

"""
    isfinite(poly)
    isinfinite(poly)

Test if the polygon is infinite or finite. 
"""
:isfinite, :isinfinite 

import Base.isfinite 
function isfinite(poly::InfinitePolygon)
    if poly.head == (0,0) && poly.tail == (0,0)
        return true
    else
        # TODO, check if the head/tail rays intersect... 
        return false 
    end 
end

isinfinite(poly) = !isfinite(poly)

firstpoint(poly) = first(poly.points)
function lastpoint(poly)
    p0 = first(poly.points)
    for p in poly.points
        p0 = p 
    end 
    return p0 
end 

"""
"""
function _isright(pa, pb, pt)
    x1,y1 = pa
    x2,y2 = pb
    x,y = pt
    return (x2 - x1)*(y - y1) < (y2 - y1)*(x - x1)
end 
function _isleft(pa, pb, pt)
    return orient(pa...,pb...,pt...)
end 



import Base.contains
"""
    contains(p::InfinitePolygon, pt)
    
Test if the infinite polygon contains a point.
    The point type pt must be able to be iterated to a pair 
    of numbers 
"""    
function contains(p::InfinitePolygon, pt)
    if isfinite(p)
        p0 = firstpoint(p)
        for p1 in p.points
            # note that dist really computes distance squared...
            if p0 == p1 || dist(p0...,p1...) <= eps(eltype(eltype(p)))  # skip over the first point, or repeated points 
                continue
            end 
            if _isleft(p0, p1, pt) == false
                return false
            end
            p0 = p1 
        end 
        # p0 is lastpoint now...
        p1 = firstpoint(p)
        if _isleft(p0, p1, pt) == false
            return false
        end 
    else # handle the infinite case... 
        # need to test if the polygon with infinite rays contains a point. 
        # for each edge, we need to test if the point is on the right side of the line.
        p0 = firstpoint(p)
        p0 = p0.+p.head # using the .+ makes it work for a combo of points and tuples
        for p1 in p.points
            if p0 == p1 || dist(p0...,p1...) <= eps(eltype(eltype(p)))  # skip over the first point, or repeated points 
                continue
            end 
            if _isleft(p0, p1, pt) == false
                return false
            end
            p0 = p1
        end 
        # p0 is last point now... 
        p1 = p0 .+ p.tail 
        if _isleft(p0, p1, pt) == false
            return false
        end
    end 
    return true
end 

import Base.iterate, Base.eltype, Base.IteratorSize
struct SegmentsIterator{PolyType,FloatType}
    p::PolyType
    dist::FloatType 
end
Base.IteratorSize(::Type{SegmentsIterator{PolyType,FloatType}}) where {PolyType,FloatType} = Base.SizeUnknown()
Base.eltype(::Type{SegmentsIterator{PolyType,FloatType}}) where {PolyType,FloatType} = Tuple{Tuple{FloatType,FloatType},Tuple{FloatType,FloatType}}
# algorithm
# for an infinite poly
#   set p0 = p0 + head
#   then walk through all points in p.points (checking if p1 == p0, or small distance)
#   return (p0,p1), (p1, nextiter)
#   once the nextiter is done...
#   return p0, p0+tail # using p0 here means that the sequence of points hsould be closed
# for a finite poly
#   set p0 = p0
#   then walk through all the points in p.points (checking if p1 == p0, or small distance)
#   return (p0, p1), (p1, nextiter)
#   once the nextiter is done
#   return p0, p1 = firstpoint(p) 
function Base.iterate(si::SegmentsIterator, state=nothing)
    p = si.p 
    if state === nothing
        piter = iterate(p.points)
        if piter === nothing # then there are no points
            return nothing
        end

        p0 = firstpoint(p) 
        if isinfinite(p) 
            p0 = p0 .+ p.head # TODO, adjust so that it's always at least dist away... 
        end
        piter = iterate(p.points)
        lastpoint = false
    else
        p0, piter, lastpoint = state
    end 
    # at this point, p0 and piter are both initialized

    while piter !== nothing
        p1,pstate = piter
        if p1 != p0 && dist(p0...,p1...) > si.dist 
            return ((p0,p1), (p1, iterate(p.points, pstate), lastpoint))
        end
        piter = iterate(p.points, pstate)
    end 

    # at this point, piter === nothing 
    if lastpoint # once last point is set, we have handled the last end-of-stream point 
        return nothing # we are done
    else
        lastpoint = true 
        if isfinite(p)
            p1 = firstpoint(p) 
        else
            p1 = p0 .+ p.tail # no need to adjust, since we guarantee return. 
        end 
        return ((p0,p1), (p1, piter, lastpoint))
    end 
end 
Base.isdone(it::SegmentsIterator, state) = state[3] == true 

"""
    segments(p::InfinitePolygon; [dist = eps()])

Return an iterator over linesegments involved in the polygon.
If the polygon is infinite, then this will not be closed and the
first two points will be along the incoming ray and the last two
points will be along the outgoing ray (each will be an arbitrary 
length along this ray)

If the polygon is finite, then it will be closed.
"""
function segments(p::InfinitePolygon; dist::Real = eps(eltype(eltype(p))))
    FloatType = eltype(firstpoint(p))
    return SegmentsIterator{typeof(p),FloatType}(p, dist)
end 

"""
    dualcell(t, i)
    dualcell(t, centers, i)

Return the finite or infinite polygon description of the dual cell
to a given point index in the triangulation. The dual cell is the 
Voronoi cell to a Delaunay triangulation. 

The default usage uses the circumcenters of the Delaunay triangles 
as the coordinates of the Voroni vertices. However, 
you can override this by giving another collection of centers. 
The number of centers must be equal to the number 
of triangles. 
"""    
dualcell(t::Triangulation, i::Integer) = dualcell(t, t.circumcenters, i)
function dualcell(t::Triangulation, centers, i::Integer)
    FT = floattype(t) 
    if ((hi = inhull(t,i))>0)
        raystart = t.raystart[hi]
        rayend = t.rayend[hi]
    else
        raystart = (zero(FT),zero(FT))
        rayend = (zero(FT),zero(FT))
    end 
    return infinitepoly((centers[i] for i in triangles(t, i)), 
            raystart, rayend)
end 

"""
    cellarea(p)

Compute the area of a possibly infinite polygon.
"""
function cellarea(p::InfinitePolygon)
    if isinfinite(p)
        return typemax(eltype(eltype(p)))
    # Area of a line is zero
    elseif length(p.points) < 3
        return zero(eltype(eltype(p)))
    else
        # Compute the area of the polygon using the shoelace formula using iterators
        area = zero(eltype(eltype(p)))
        previous, ite = iterate(p.points)
        while true
            el = iterate(p.points, ite)
            if isnothing(el)
                break
            end
            current, ite = el
            area += (previous[1] + current[1]) * (previous[2] - current[2])
            previous = current
        end
        current = first(p.points)
        return 0.5 * abs(area + (previous[1] + current[1]) * (previous[2] - current[2]))
    end
end

# monotonically increases with real angle, but doesn't need expensive trigonometry
function pseudoAngle(dx, dy)
    p = dx / (abs(dx) + abs(dy))
    return (dy > 0 ? 3 - p : 1 + p) / 4 # [0..1]
end

function dist(ax, ay, bx, by)
    dx = ax - bx
    dy = ay - by
    return dx * dx + dy * dy
end

# return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
function orientIfSure(px, py, rx, ry, qx, qy)
    l = (ry - py) * (qx - px)
    r = (rx - px) * (qy - py)
    return abs(l - r) >= 3.3306690738754716e-16 * abs(l + r) ? l - r : 0
end

# a more robust orientation test that's stable in a given triangle (to fix robustness issues)
function orient(rx, ry, qx, qy, px, py)
    #return (orientIfSure(px, py, rx, ry, qx, qy) ||
    #    orientIfSure(rx, ry, qx, qy, px, py) ||
    #    orientIfSure(qx, qy, px, py, rx, ry)) < 0;
    return orientIfSure(px, py, rx, ry, qx, qy) < 0 ||
           orientIfSure(rx, ry, qx, qy, px, py) < 0 ||
           orientIfSure(qx, qy, px, py, rx, ry) < 0
end

function inCircle(ax, ay, bx, by, cx, cy, px, py)
    dx = ax - px
    dy = ay - py
    ex = bx - px
    ey = by - py
    fx = cx - px
    fy = cy - py

    ap = dx * dx + dy * dy
    bp = ex * ex + ey * ey
    cp = fx * fx + fy * fy

    return dx * (ey * cp - bp * fy) -
           dy * (ex * cp - bp * fx) +
           ap * (ex * fy - ey * fx) < 0
end

function circumradius(ax, ay, bx, by, cx, cy)
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    d = 0.5 / (dx * ey - dy * ex)

    x = (ey * bl - dy * cl) * d
    y = (dx * cl - ex * bl) * d

    return x * x + y * y
end

function circumcenter(ax, ay, bx, by, cx, cy)
    dx = bx - ax
    dy = by - ay
    ex = cx - ax
    ey = cy - ay

    bl = dx * dx + dy * dy
    cl = ex * ex + ey * ey
    d = 0.5 / (dx * ey - dy * ex)

    x = ax + (ey * bl - dy * cl) * d
    y = ay + (dx * cl - ex * bl) * d

    return (x, y)
end

"""
    truncatedcircumcenter

For a nearly collinear triangle, then the circumcenter can be
off at a point near infinity. Since the goal of this library
is not computational geometry, a pragmatic choice is to truncate
these wildly divergent near infinite circumcenters. 

function from d3-delaunay / Voronoi.js
"""
function truncatedcircumcenter(x1::FloatType, y1::FloatType, 
        x2::FloatType, y2::FloatType, 
        x3::FloatType, y3::FloatType, 
        rx::FloatType, ry::FloatType,  
        collineartol=sqrt(eps(FloatType))) where FloatType
    
        #=
    tris = t.triangles
    points = t.points
    t1,t2,t3 = tris[i]
    x1,y1 = point(FloatType, points, t1)
    x2,y2 = point(FloatType, points, t2)
    x3,y3 = point(FloatType, points, t3)
    =#

    #=
    # TODO, make this relative? 
    dx = x2 - x1
    dy = y2 - y1
    ex = x3 - x1
    ey = y3 - y1
    #ab = (dx * ey - dy * ex) * 2
    ab = -Delaunator.orientIfSure(x2,y2,x3,y3,x1,y1)*2

    if abs(ab) < collineartol
        # degenerate case (collinear diagram)
        # almost equal points (degenerate triangle)
        # the circumcenter is at the infinity, in a
        # direction that is:
        # 1. orthogonal to the halfedge.
        a = 1/FloatType(collineartol)
        # 2. points away from the center; since the list of triangles starts
        # in the center, the first point of the first triangle
        # will be our reference
        # In Julia, we just pass this reerence point as an option. 
        a *= sign((rx - x1) * ey - (ry - y1) * ex)
        x = (x1 + x3) / 2 - a * ey
        y = (y1 + y3) / 2 + a * ex
    else
        d = 1 / ab
        bl = dx * dx + dy * dy
        cl = ex * ex + ey * ey
        x = x1 + (ey * bl - dy * cl) * d
        y = y1 + (dx * cl - ex * bl) * d
    end 
    return (x,y)
    =#
    return circumcenter(x1,y1,x2,y2,x3,y3)
end 
