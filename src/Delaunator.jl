module Delaunator

using OffsetArrays

const EPSILON = 2^-52 # eps()
const EDGE_STACK = OffsetVector{UInt32}(undef,0:511)

mutable struct DelaunatorData
    _triangles::OffsetVector{UInt32}
    _halfedges::OffsetVector{UInt32}
end



function delaunator!(points)
    n = length(points)
    coords = OffsetVector{Float64}(undef, 0:(n*2))

    @show n
    @show coords
    for i = 0:n-1
        p = points[i+1]
        @show p
        coords[2 * i] = defaultGetX(p)
        coords[2 * i + 1] = defaultGetY(p)
    end
    @show "1"
    #return new Delaunator(coords);

    n = div(n,2) # coords.length >> 1
    #if (n > 0 && typeof coords[0] !== 'number') throw new Error('Expected coords to contain numbers.')

    # arrays that will store the triangulation graph
    maxTriangles = max(2 * n - 5, 0)
    @show typeof(maxTriangles),  maxTriangles
    _triangles =  OffsetVector{UInt32}(undef,0:(maxTriangles * 3)-1)
    _halfedges = OffsetVector{UInt32}(undef,0:(maxTriangles * 3)-1)

    # temporary arrays for tracking the edges of the advancing convex hull
    _hashSize = ceil(Int,sqrt(n))
    hullPrev = OffsetVector{UInt32}(undef,0:n-1) # edge to prev edge
    hullNext = OffsetVector{UInt32}(undef,0:n-1)# edge to next edge
    hullTri = OffsetVector{UInt32}(undef,0:n-1) # edge to adjacent triangle
    hullHash = fill!(OffsetVector{Int32}(undef,0:n-1), -1) # angular edge hash

    # temporary arrays for sorting points
    _ids = OffsetVector{UInt32}(undef,0:n-1)
    _dists = OffsetVector{Float64}(undef,0:n-1)

    #this.update()

    #n = n/2 # coords.length >> 1 # n/2

    # populate an array of point indices; calculate input data bbox
    minX = Inf
    minY = Inf
    maxX = -Inf
    maxY = -Inf

    for i = 0:n-1
        x = coords[2 * i]
        y = coords[2 * i + 1]
        x < minX && (minX = x)
        y < minY && (minY = y)
        x > maxX && (maxX = x)
        y > maxY && (maxY = y)
        _ids[i] = i
    end
    cx = (minX + maxX) / 2
    cy = (minY + maxY) / 2
    @show cx, cy

    minDist = Inf
    i0 = 0
    i1 = 0
    i2 = 0

    # pick a seed point close to the center
    for i = 0:n-1
        d = dist(cx, cy, coords[2 * i], coords[2 * i + 1])
        @show d
        if d < minDist
            i0 = i
            @show minDist
            minDist = d
        end
    end
    i0x = coords[2 * i0]
    i0y = coords[2 * i0 + 1]
    @show i0x,i0y
    minDist = Inf

    # find the point closest to the seed
    for i = 0:n-1
        i == i0 && continue
        @show minDist
        d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1])
        if d < minDist && d > 0
            i1 = i
            minDist = d
        end
    end
    i1x = coords[2 * i1]
    i1y = coords[2 * i1 + 1]

    @show i1x, i1y

    minRadius = Inf

    # find the third point which forms the smallest circumcircle with the first two
    for i = 1:n-1
        i == i0 || i == i1 && continue
        r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1])
        if r < minRadius
            i2 = i
            minRadius = r
        end
    end
    i2x = coords[2 * i2]
    i2y = coords[2 * i2 + 1]

    if minRadius == Inf
        # order collinear points by dx (or dy if all x are identical)
        # and return the list as a hull
        for i = 0:n-1
            cxv = (coords[2 * i] - coords[0])
            _dists[i] = !iszero(cxv) ? cxv : (coords[2 * i + 1] - coords[1])
        end
        quicksort(_ids, _dists, 0, n - 1)
        hull = OffsetVector{UInt32}(undef,0:n-1)
        j = 0
        d0 = -Inf
        for i = 0:n-1
            id = _ids[i]
            if _dists[id] > d0
                hull[j] = id
                j += 1
                d0 = _dists[id]
            end
        end
        hull = resize!(hull, j+1)
        triangles = OffsetVector{UInt32}(undef, 0:0)
        halfedges = OffsetVector{UInt32}(undef, 0:0)
        return hull # TODO
    end

    # swap the order of the seed points for counter-clockwise orientation
    if orient(i0x, i0y, i1x, i1y, i2x, i2y)
        i = i1
        x = i1x
        y = i1y
        i1 = i2
        i1x = i2x
        i1y = i2y
        i2 = i
        i2x = x
        i2y = y
    end

    center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
    _cx = center[1]
    _cy = center[2]

    for i = 0:n-1
        _dists[i] = dist(coords[2 * i], coords[2 * i + 1], _cx, _cy)
    end

    # sort the points by distance from the seed triangle circumcenter
    quicksort(_ids, _dists, 0, n - 1)

    # set up the seed triangle as the starting hull
    _hullStart = i0
    hullSize = 3

    hullNext[i0] = hullPrev[i2] = i1
    hullNext[i1] = hullPrev[i0] = i2
    hullNext[i2] = hullPrev[i1] = i0

    hullTri[i0] = 0
    hullTri[i1] = 1
    hullTri[i2] = 2

    fill!(hullHash,-1)
    hullHash[_hashKey(i0x, i0y,_cx,_cy,_hashSize)] = i0
    hullHash[_hashKey(i1x, i1y,_cx,_cy,_hashSize)] = i1
    hullHash[_hashKey(i2x, i2y,_cx,_cy,_hashSize)] = i2

    trianglesLen = _addTriangle(_triangles, _halfedges, 0, i0, i1, i2, -1, -1, -1)

    xp = 0.0
    yp = 0.0
    for k = 0:length(_ids)-1
        i = _ids[k]
        x = coords[2 * i]
        y = coords[2 * i + 1]

        # skip near-duplicate points
        if k > 0 && abs(x - xp) <= EPSILON && abs(y - yp) <= EPSILON
            continue
        end
        xp = x
        yp = y

        # skip seed triangle points
        if i === i0 || i === i1 || i === i2
            continue
        end

        # find a visible edge on the convex hull using edge hash
        start = 0
        key = _hashKey(x, y, _cx, _cy, _hashSize)
        @show key
        for j = 0:_hashSize-1
            start = hullHash[(key + j) % _hashSize]
            start != -1 && start != hullNext[start] && break
        end

        start = hullPrev[start]
        e = start
        q = hullNext[e]
        while !orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])
            e = q
            if e == start
                e = -1
                break
            end
            q = hullNext[e]
        end

        e == -1 && continue # likely a near-duplicate point; skip it

        # add the first triangle from the point
        t = _addTriangle(e, i, hullNext[e], -1, -1, hullTri[e])

        # recursively flip triangles from the point until they satisfy the Delaunay condition
        hullTri[i] = _legalize(t + 2)
        hullTri[e] = t # keep track of boundary triangles on the hull
        hullSize += 1

        # walk forward through the hull, adding more triangles and flipping recursively
        n = hullNext[e]
        q = hullNext[n]
        while orient(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1])
            t = _addTriangle(n, i, q, hullTri[i], -1, hullTri[n])
            hullTri[i] = _legalize(t + 2)
            hullNext[n] = n # mark as removed
            hullSize -= 1
            n = q
            q = hullNext[n]
        end

        # walk backward from the other side, adding more triangles and flipping
        if e == start
            q = hullPrev[e]
            while orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])
                t = _addTriangle(q, i, e, -1, hullTri[e], hullTri[q])
                _legalize(t + 2)
                hullTri[q] = t
                hullNext[e] = e # mark as removed
                hullSize -= 1
                e = q
                q = hullPrev[e]
            end
        end

        # update the hull indices
        _hullStart = hullPrev[i] = e
        hullNext[e] = hullPrev[n] = i
        hullNext[i] = n

        # save the two new edges in the hash table
        hullHash[_hashKey(x, y,_cx,_cy,_hashSize)] = i
        hullHash[_hashKey(coords[2 * e], coords[2 * e + 1],_cx,_cy,_hashSize)] = e
    end

    hull = OffsetVector{UInt32}(undef,0:hullSize-1)
    for i = 0:hullSize-1
        e = _hullStart
        hull[i] = e
        e = hullNext[e]
    end

    # trim typed triangle mesh arrays
    return (resize!(_triangles, trianglesLen), resize!(_halfedges, trianglesLen))
end


function _hashKey(x, y, cx, cy, _hashSize)
    return floor(Int, pseudoAngle(x - cx, y - cy) * _hashSize) % _hashSize
end

function _legalize(a, triangles, halfedges, coords)

    i = 0
    ar = 0

    #recursion eliminated with a fixed-size stack
    while true
        b = halfedges[a]

         # if the pair of triangles doesn't satisfy the Delaunay condition
         # (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
         # then do the same check/flip recursively for the new pair of triangles
         #
         #           pl                    pl
         #          /||\                  /  \
         #       al/ || \bl            al/    \a
         #        /  ||  \              /      \
         #       /  a||b  \    flip    /___ar___\
         #     p0\   ||   /p1   =>   p0\---bl---/p1
         #        \  ||  /              \      /
         #       ar\ || /br             b\    /br
         #          \||/                  \  /
         #           pr                    pr
         #
        a0 = a - a % 3
        ar = a0 + (a + 2) % 3

        if b == -1 #convex hull edge
            i == 0 && break
            i -= 1
            a = EDGE_STACK[i]
            continue
        end

        b0 = b - b % 3
        al = a0 + (a + 1) % 3
        bl = b0 + (b + 2) % 3

        p0 = triangles[ar]
        pr = triangles[a]
        pl = triangles[al]
        p1 = triangles[bl]

        illegal = inCircle(coords[2 * p0], coords[2 * p0 + 1],
                           coords[2 * pr], coords[2 * pr + 1],
                           coords[2 * pl], coords[2 * pl + 1],
                           coords[2 * p1], coords[2 * p1 + 1])

        if illegal
            triangles[a] = p1
            triangles[b] = p0

            hbl = halfedges[bl]

            # edge swapped on the other side of the hull (rare); fix the halfedge reference
            if hbl == -1
                e = _hullStart
                # CAUTION do-while
                if hullTri[e] == bl
                    hullTri[e] = a
                else
                    e = hullPrev[e]
                    while e != _hullStart
                        if this.hullTri[e] == bl
                            this.hullTri[e] = a
                            break
                        end
                        e = hullPrev[e]
                    end
                end
            end
            _link(halfedges, a, hbl)
            _link(halfedges, b, halfedges[ar])
            _link(halfedges, ar, bl)

            br = b0 + (b + 1) % 3

            # don't worry about hitting the cap: it can only happen on extremely degenerate input
            if i < EDGE_STACK.length
                EDGE_STACK[i] = br
                i += 1
            end
        else
            i == 0 && break
            i -= 1
            a = EDGE_STACK[i]
        end
    end

    return ar
end


function _link(_halfedges, a, b)
    @show a, b
    _halfedges[a] = b == -1 ? typemax(UInt32) : b
    b != -1 && (_halfedges[b] = a)
end


# add a new triangle given vertex indices and adjacent half-edge ids
function _addTriangle(_triangles, _halfedges, t, i0, i1, i2, a, b, c)

    _triangles[t] = i0
    _triangles[t + 1] = i1
    _triangles[t + 2] = i2

    _link(_halfedges, t, a)
    _link(_halfedges, t + 1, b)
    _link(_halfedges, t + 2, c)

    return t + 3
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

function quicksort(ids, dists, left, right)
    if right - left <= 20
        i = i = left + 1
        while i <= right
            temp = ids[i]
            tempDist = dists[temp]
            j = i - 1
            while j >= left && dists[ids[j]] > tempDist
                ids[j + 1] = ids[j]
                j -= 1
            end
            ids[j + 1] = temp
            i += 1
        end
    else
        median = (left + right) >> 1
        i = left + 1
        j = right
        swap(ids, median, i)
        dists[ids[left]] > dists[ids[right]] && swap(ids, left, right)
        dists[ids[i]] > dists[ids[right]] && swap(ids, i, right)
        dists[ids[left]] > dists[ids[i]] && swap(ids, left, i)

        temp = ids[i]
        tempDist = dists[temp]
        while true
            i += 1
            while dists[ids[i]] < tempDist
                i += 1
            end
            j -= 1
            while dists[ids[j]] > tempDist
                j -= 1
            end
            j < i && break
            swap(ids, i, j)
        end
        ids[left + 1] = ids[j]
        ids[j] = temp

        if (right - i + 1 >= j - left)
            quicksort(ids, dists, i, right)
            quicksort(ids, dists, left, j - 1)
        else
            quicksort(ids, dists, left, j - 1)
            quicksort(ids, dists, i, right)
        end
    end
end

function swap(arr, i, j)
    tmp = arr[i]
    arr[i] = arr[j]
    arr[j] = tmp
end

function defaultGetX(p)
    return p[1]
end

function defaultGetY(p)
    return p[2]
end


end # module
