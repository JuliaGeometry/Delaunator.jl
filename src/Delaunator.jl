module Delaunator

using OffsetArrays

const EPSILON = 2^-52 # eps()
const EDGE_STACK = OffsetVector{UInt32}(undef,0:511)



function delaunator(coords) {
    const n = length(points.length)
    const coords = Vector{Float64}(undef, n * 2)

    for (let i = 0; i < n; i++) {
        const p = points[i];
        coords[2 * i] = defaultgetX(p);
        coords[2 * i + 1] = defaultgetY(p);
    }

    return new Delaunator(coords);

    const n = coords.length >> 1
    if (n > 0 && typeof coords[0] !== 'number') throw new Error('Expected coords to contain numbers.')

    this.coords = coords

    # arrays that will store the triangulation graph
    const maxTriangles = Math.max(2 * n - 5, 0)
    this._triangles = new Uint32Array(maxTriangles * 3)
    this._halfedges = new Int32Array(maxTriangles * 3)

    # temporary arrays for tracking the edges of the advancing convex hull
    this._hashSize = Math.ceil(Math.sqrt(n));
    this._hullPrev = new Uint32Array(n) # edge to prev edge
    this._hullNext = new Uint32Array(n) # edge to next edge
    this._hullTri = new Uint32Array(n) # edge to adjacent triangle
    this._hullHash = new Int32Array(this._hashSize).fill(-1) # angular edge hash

    # temporary arrays for sorting points
    this._ids = new Uint32Array(n)
    this._dists = new Float64Array(n)

    this.update()
end


function update(coords, hullPrev, hullNext, hullTri, hullHash)
    n = coords.length >> 1

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

    minDist = Inf
    i0 = 0
    i1 = 0
    i2 = 0

    # pick a seed point close to the center
    for i = 0:n-1
        d = dist(cx, cy, coords[2 * i], coords[2 * i + 1])
        if d < minDist
            i0 = i
            minDist = d
        end
    end
    i0x = coords[2 * i0]
    i0y = coords[2 * i0 + 1]

    minDist = Inf

    # find the point closest to the seed
    for i = 0:n-1
        i == i0 && continue
        d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1])
        if d < minDist && d > 0
            i1 = i
            minDist = d
        end
    end
    i1x = coords[2 * i1]
    i1y = coords[2 * i1 + 1]

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
        for (let i = 0; i < n; i++)
            this._dists[i] = (coords[2 * i] - coords[0]) || (coords[2 * i + 1] - coords[1]);
        end
        quicksort(this._ids, this._dists, 0, n - 1)
        hull = new Uint32Array(n)
        let j = 0;
        for (let i = 0, d0 = -Infinity; i < n; i++) {
            const id = this._ids[i];
            if (this._dists[id] > d0) {
                hull[j++] = id;
                d0 = this._dists[id];
            }
        }
        this.hull = hull.subarray(0, j);
        this.triangles = new Uint32Array(0);
        this.halfedges = new Uint32Array(0);
        return;
    end

    # swap the order of the seed points for counter-clockwise orientation
    if (orient(i0x, i0y, i1x, i1y, i2x, i2y))
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
        this._dists[i] = dist(coords[2 * i], coords[2 * i + 1], _cx, _cy);
    end

    # sort the points by distance from the seed triangle circumcenter
    quicksort(_ids, _dists, 0, n - 1);

    # set up the seed triangle as the starting hull
    _hullStart = i0
    hullSize = 3

    hullNext[i0] = hullPrev[i2] = i1
    hullNext[i1] = hullPrev[i0] = i2
    hullNext[i2] = hullPrev[i1] = i0

    hullTri[i0] = 0;
    hullTri[i1] = 1;
    hullTri[i2] = 2;

    hullHash.fill(-1);
    hullHash[this._hashKey(i0x, i0y)] = i0;
    hullHash[this._hashKey(i1x, i1y)] = i1;
    hullHash[this._hashKey(i2x, i2y)] = i2;

    this.trianglesLen = 0;
    this._addTriangle(i0, i1, i2, -1, -1, -1);

    for (let k = 0, xp, yp; k < this._ids.length; k++) {
        const i = this._ids[k];
        const x = coords[2 * i];
        const y = coords[2 * i + 1];

        // skip near-duplicate points
        if (k > 0 && Math.abs(x - xp) <= EPSILON && Math.abs(y - yp) <= EPSILON) continue;
        xp = x
        yp = y

        // skip seed triangle points
        if (i === i0 || i === i1 || i === i2) continue;

        // find a visible edge on the convex hull using edge hash
        let start = 0;
        for (let j = 0, key = this._hashKey(x, y); j < this._hashSize; j++) {
            start = hullHash[(key + j) % this._hashSize];
            if (start !== -1 && start !== hullNext[start]) break;
        }

        start = hullPrev[start];
        let e = start, q;
        while (q = hullNext[e], !orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
            e = q;
            if (e === start) {
                e = -1;
                break;
            }
        }
        if (e === -1) continue; // likely a near-duplicate point; skip it

        // add the first triangle from the point
        let t = this._addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

        // recursively flip triangles from the point until they satisfy the Delaunay condition
        hullTri[i] = this._legalize(t + 2);
        hullTri[e] = t; // keep track of boundary triangles on the hull
        hullSize++;

        // walk forward through the hull, adding more triangles and flipping recursively
        let n = hullNext[e];
        while (q = hullNext[n], orient(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1])) {
            t = this._addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
            hullTri[i] = this._legalize(t + 2);
            hullNext[n] = n; // mark as removed
            hullSize--;
            n = q;
        }

        // walk backward from the other side, adding more triangles and flipping
        if (e === start) {
            while (q = hullPrev[e], orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
                t = this._addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
                this._legalize(t + 2);
                hullTri[q] = t;
                hullNext[e] = e; // mark as removed
                hullSize--;
                e = q;
            }
        }

        // update the hull indices
        this._hullStart = hullPrev[i] = e;
        hullNext[e] = hullPrev[n] = i;
        hullNext[i] = n;

        // save the two new edges in the hash table
        hullHash[this._hashKey(x, y)] = i;
        hullHash[this._hashKey(coords[2 * e], coords[2 * e + 1])] = e;
    }

    this.hull = new Uint32Array(hullSize);
    for (let i = 0, e = this._hullStart; i < hullSize; i++) {
        this.hull[i] = e;
        e = hullNext[e];
    }

    // trim typed triangle mesh arrays
    this.triangles = this._triangles.subarray(0, this.trianglesLen);
    this.halfedges = this._halfedges.subarray(0, this.trianglesLen);
}


function _hashKey(x, y, cx, cy, _hashSize)
    return floor(pseudoAngle(x - _cx, y - _cy) * _hashSize) % _hashSize
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
            if hbl === -1
                e = _hullStart
                # CAUTION do-while
                if _hullTri[e] == bl
                    _hullTri[e] = a
                else
                    e = _hullPrev[e]
                    while e != _hullStart
                        if this._hullTri[e] == bl
                            this._hullTri[e] = a
                            break
                        end
                        e = this._hullPrev[e]
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
    _halfedges[a] = b
    b !== -1 && _halfedges[b] = a
end


# add a new triangle given vertex indices and adjacent half-edge ids
function _addTriangle(_triangles, _halfedges, t, i0, i1, i2, a, b, c)

    _triangles[t] = i0
    _triangles[t + 1] = i1
    _triangles[t + 2] = i2

    _link(_halfedges, t, a)
    _link(_halfedges, t + 1, b)
    _link(_halfedges, t + 2, c)

    nothing
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
    return (orientIfSure(px, py, rx, ry, qx, qy) ||
        orientIfSure(rx, ry, qx, qy, px, py) ||
        orientIfSure(qx, qy, px, py, rx, ry)) < 0
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
                #ids[j + 1] = ids[j--]
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
            # TODO
            #do i++; while (dists[ids[i]] < tempDist);
            #do j--; while (dists[ids[j]] > tempDist);
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
