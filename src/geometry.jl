
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
    # TODO, make this relative? 
    dx = x2 - x1
    dy = y2 - y1
    ex = x3 - x1
    ey = y3 - y1
    ab = (dx * ey - dy * ex) * 2

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
end 

