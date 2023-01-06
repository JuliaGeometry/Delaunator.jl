

"""
    InfinitePolygon

fields
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
function isfinite(poly)
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

"""
import Base.iterate, Base.eltype, Base.IteratorSize
struct SegmentsIterator{PolyType,StartType}
    p::PolyType
    p0::StartType
end
Base.IteratorSize(::Type{SegmentsIterator{PolyType,IntType}}) where {PolyType,IntType} = Base.SizeUnknown()
Base.eltype(::Type{SegmentsIterator{PolyType}}) where {PolyType} = eltype(PolyType)
#Base.isdone(it::SegmentsIterator, state) = state
#Base.length()

The logic here is: we cache the start point in the state.

function Base.iterate(it::SegmentsIterator, state=nothing) 
    if # check done...
    else
        if state === nothing
            p0 = firstpoint(it.p)
            state = 
    # we are done if state == it.start again, after the first iteration
    if state == it.start || state == -1
        return nothing
    else
        if state === nothing 
            state = it.start
        end 
        if state == 0 
            return nothing
        end
        outgoing = _nextedge(state)
        incoming = it.t.halfedges[outgoing]
        nextstate = incoming 
        return (triangle_from_index(it.t, state), nextstate)
    end 
end

function Base.length(it::TriangleNeighborIterator) 
    len = 0 
    start = it.start 
    if start > 0
        len += 1
        halfedges = it.t.halfedges 
        state = start
        outgoing = _nextedge(state)
        state = halfedges[outgoing]
        while (state != -1 && state != start)
            outgoing = _nextedge(state)
            state = halfedges[outgoing]
            len += 1
        end 
    end 
    return len 
end 
"""

"""
    segments(p::InfinitePolygon; )

Return an iterator over linesegments involved in the polygon.
If the polygon is infinite, then this will not be closed.
If the polygon is finite, then it will be closed.
"""
function segments(p::InfinitePolygon)
    if isfinite(p)
    else
    end 
end 

"""
    dualcell(t, i)
    dualcell(t, centers, i)

Return the finite or infinite polygon description 
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

# allow us to use the bounding box api for generic tuple bounding boxes... 
import Base.maximum 
origin(bbox::Tuple{Tuple{FT,FT},Tuple{FT,FT}}) where FT = bbox[1]
maximum(bbox::Tuple{Tuple{FT,FT},Tuple{FT,FT}}) where FT = bbox[2]

function _inside_bbox(pt, lowerleft, upperright)
    if pt[1] < lowerleft[1]
        return false
    elseif pt[1] > upperright[1]
        return false
    end
    if pt[2] < lowerleft[2]
        return false
    elseif pt[2] > upperright[2]
        return false
    end
    return true
end 

"""
    clippedpoly(p::InfinitePolygon, bbox)

returns nothing if the poly is entirely outside the bounding box.
Otherwise, return a set of points that represent the infinite
polygon clipped to the bounding box. 
"""
function clippedpoly(p::InfinitePolygon, bbox)
    if isfinite(p)
        return _clip_finite(p, bbox)
    else
        return _clip_infinite(p, bbox)
    end 
end 

#=
function _clip_finite(p, bbox)
    lowerleft = origin(bbox)
    upperright = maximum(bbox)
    
    for p in poly.points
        if _inside_bbox(p, lowerleft, upperright)

        else
            lastoutside = p 
        end 
    end 
end 
=#


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

