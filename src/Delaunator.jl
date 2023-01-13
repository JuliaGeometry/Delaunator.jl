module Delaunator

include("interface.jl")
export triangulate, basictriangulation, update!, triangles, points, inhull

include("quicksort.jl")
include("geometry.jl")
export isinfinite, dualcell, firstpoint, lastpoint, segments

include("clipping.jl")
export clippedpoly, clippedpoly!

include("algorithm.jl") # nothing is exported here... 

include("iterators.jl")
export neighbors, edges, edgelines, hullpoly

#include("clipping.jl")

"""
    delaunator([IntType=Int32,] [FloatType=Float64,] points; [tol=eps(FloatType])])

Computes a Delaunay triangulation of a set of points using the Delaunator algorithm.
This is designed for quick graphics applications and speed rather than exact computational geometry.

Inputs
------
- `points` is any type that has integer indexing and length supported. In addition, `p = points[i]` should
be a type where `p[1]` and `p[2]` are the x, y coordinates of p. Or you need to define the functions 
`Delaunator.getX(p), Delaunator.getY(p)` for your own type p. 

- `tol` is used to determine when points are sufficiently close not to include.

Return value
------------
a named tuple with fields 
- `triangles`: a length T array of 3 tuples, where each tuple is a triangle
- `halfedges`: The halfedge index for the edge in the triangle array. The halfedges for 
    triangle t are 3(t-1)+1, 3(t-1)+2, 3(t-2)+3. Each halfedge index gives the entry
    of the other halfedge. 
- `hull`: The sequence of point indexes that define the convex hull of the input points. 
- `points` a copy of the input set of points
- `index`: a point to halfedge/triangle index. This gives the leftmost halfedge index,
   to make it feasible to iterate over the    

These data structures are explained at https://mapbox.github.io/delaunator/ 
(but here, all the indices have been modified a little).
"""


function triangle_edges(d::Triangulation, i::Integer)
    return (3*(i-1)+1,3*(i-1)+2,3*(i-1)+3)
end 

#=
delaunator(points) = delaunator(Int32, Float32, points)

# pickup the FloatType from dists... 
delaunator!(triangles, halfedges, hull, points) = 
delaunator!(dists, triangles, halfedges, hull, points)
function delaunator(::Type{IntType}, ::Type{FloatType}, points) where {IntType, FloatType}
#    _triangles = 
end 

function delaunator!(::Type{FloatType}, triangles, halfedges, hull, points) where {FloatType}
IntType = promote_type(eltype(triangles), eltype(halfedges), eltype(hull))

end 
=#

delaunator!(coords) = delaunator!(Float64, coords)
delaunator!(::Type{FloatType}, coords) where FloatType = delaunator!(Int32, FloatType, coords)

function delaunator!(::Type{IntType}, ::Type{FloatType}, coords; tol=eps(FloatType)) where {IntType <: Signed, FloatType}
    n = length(coords)

    # arrays that will store the triangulation graph
    maxTriangles = max(2 * n - 5, 0)
    #@show typeof(maxTriangles),  maxTriangles
    _tridata =  Vector{Tuple{IntType,IntType,IntType}}(undef,maxTriangles)
    #_halfedges =  Vector{Tuple{IntType,IntType,IntType}}(undef,maxTriangles)
    #_triangles =  Vector{IntType}(undef,3*maxTriangles)
    #_triangles = reinterpret(IntType, _tridata)
    _ptridata = Base.unsafe_convert(Ptr{IntType}, _tridata)
    _triangles = unsafe_wrap(Array, _ptridata, 3*maxTriangles)
    _halfedges =  fill!(Vector{IntType}(undef,3*maxTriangles),-1)
    _edgeindex = fill!(Vector{IntType}(undef, n), 0)

    # temporary arrays for tracking the edges of the advancing convex hull
    _hashSize = ceil(Int,sqrt(n))
    hullPrev = Vector{IntType}(undef,n) # edge to prev edge
    hullNext = Vector{IntType}(undef,n)# edge to next edge
    hullTri = Vector{IntType}(undef,n) # edge to adjacent triangle
    hullHash = fill!(Vector{IntType}(undef,n), -1) # angular edge hash
    edgeStack = Vector{IntType}(undef,500)

    legalize = (t,_hullStart)->_legalize(FloatType, t, _triangles, _halfedges, coords, edgeStack, hullPrev, hullTri, _hullStart)
    # temporary arrays for sorting points
    _ids = Vector{IntType}(undef,n)
    _dists = Vector{FloatType}(undef,n)

    

    #this.update()

    #n = n/2 # coords.length >> 1 # n/2

    # populate an array of point indices; calculate input data bbox
    minX = Inf
    minY = Inf
    maxX = -Inf
    maxY = -Inf

    for i in 1:n
        x,y = point(FloatType,coords, i)
        x < minX && (minX = x)
        y < minY && (minY = y)
        x > maxX && (maxX = x)
        y > maxY && (maxY = y)
        _ids[i] = i
    end
    cx = (minX + maxX) / 2
    cy = (minY + maxY) / 2
    #@show cx, cy

    minDist = Inf
    i0 = 1
    i1 = 1
    i2 = 1

    # pick a seed point close to the center
    for i in 1:n
        d = dist(cx, cy, point(FloatType,coords,i)...)
        #@show d
        if d < minDist
            i0 = i
            #@show minDist
            minDist = d
        end
    end
    i0x,i0y = point(FloatType,coords,i0)
    #@show i0x,i0y
    minDist = Inf

    # find the point closest to the seed
    for i in 1:n
        i == i0 && continue
        #@show minDist
        d = dist(i0x, i0y, point(FloatType,coords,i)...)
        if d < minDist && d > 0
            i1 = i
            minDist = d
        end
    end
    i1x,i1y = point(FloatType,coords,i1)

    #@show i1x, i1y

    minRadius = Inf

    # find the third point which forms the smallest circumcircle with the first two
    #@show i0, i1
    for i in 2:n
        i == i0 || i == i1 && continue
        r = circumradius(i0x, i0y, i1x, i1y, point(FloatType,coords,i)...)
        if r < minRadius
            i2 = i
            minRadius = r
        end
    end
    i2x,i2y = point(FloatType,coords,i2)

    if minRadius == Inf
        # order collinear points by dx (or dy if all x are identical)
        # and return the list as a hull
        for i in eachindex(coords)
            cxv = (point(FloatType,coords,i)[1] - point(FloatType,coords,1)[1])
            _dists[i] = !iszero(cxv) ? cxv : (point(FloatType,coords,i)[2] - point(FloatType,coords,1)[2])
        end
        quicksort(_ids, _dists, 1, n)
        hull = Vector{IntType}(undef,n)
        j = 1
        d0 = -Inf
        for i in 1:n
            id = _ids[i]
            if _dists[id] > d0
                hull[j] = id
                j += 1
                d0 = _dists[id]
            end
        end
        hull = resize!(hull, j-1)
        
        return _returntype(
                resize!(_tridata, 0), 
                resize!(_halfedges, 0), 
                hull, coords, _edgeindex)
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

    _cx,_cy = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);

    for i = 1:n
        _dists[i] = dist(point(FloatType,coords,i)...,  _cx, _cy)
    end

    # sort the points by distance from the seed triangle circumcenter
    quicksort(_ids, _dists, 1, n)

    # set up the seed triangle as the starting hull
    _hullStart = i0
    hullSize = 3

    hullNext[i0] = hullPrev[i2] = i1
    hullNext[i1] = hullPrev[i0] = i2
    hullNext[i2] = hullPrev[i1] = i0

    hullTri[i0] = 1
    hullTri[i1] = 2
    hullTri[i2] = 3

    fill!(hullHash,-1)
    hullHash[_hashKey(i0x, i0y,_cx,_cy,_hashSize)] = i0
    hullHash[_hashKey(i1x, i1y,_cx,_cy,_hashSize)] = i1
    hullHash[_hashKey(i2x, i2y,_cx,_cy,_hashSize)] = i2

    trianglesLen = _addTriangle(_triangles, _halfedges, 1, i0, i1, i2, -1, -1, -1)

    xp = 0.0
    yp = 0.0
    for k = 1:length(_ids)
        i = _ids[k]
        x,y = point(FloatType,coords,i)

        # skip near-duplicate points
        if k > 1 && abs(x - xp) <= tol && abs(y - yp) <= tol
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
        #@show key
        for j = 1:_hashSize
            start = hullHash[((key + j) % _hashSize)+1]
            start != -1 && start != hullNext[start] && break
        end

        start = hullPrev[start]
        e = start
        q = hullNext[e]
        while !orient(x, y, point(FloatType,coords,e)..., point(FloatType,coords,q)...)
            e = q
            if e == start
                e = -1
                break
            end
            q = hullNext[e]
        end

        e == -1 && continue # likely a near-duplicate point; skip it

        # add the first triangle from the point
        trianglesLen = _addTriangle(_triangles, _halfedges, trianglesLen, e, i, hullNext[e], -1, -1, hullTri[e])
        t = trianglesLen-3

        # recursively flip triangles from the point until they satisfy the Delaunay condition
        hullTri[i] = legalize(t + 2, _hullStart)
        hullTri[e] = t # keep track of boundary triangles on the hull
        hullSize += 1

        # walk forward through the hull, adding more triangles and flipping recursively
        n = hullNext[e]
        q = hullNext[n]
        while orient(x, y, point(FloatType,coords,n)..., point(FloatType,coords,q)...)
            trianglesLen = _addTriangle(_triangles, _halfedges, trianglesLen, n, i, q, hullTri[i], -1, hullTri[n])
            t = trianglesLen-3
            hullTri[i] = legalize(t + 2, _hullStart)
            hullNext[n] = n # mark as removed
            hullSize -= 1
            n = q
            q = hullNext[n]
        end

        # walk backward from the other side, adding more triangles and flipping
        if e == start
            q = hullPrev[e]
            while orient(x, y, point(FloatType,coords,q)...,  point(FloatType,coords,e)...)
                trianglesLen = _addTriangle(_triangles, _halfedges, trianglesLen, q, i, e, -1, hullTri[e], hullTri[q])
                t = trianglesLen-3
                legalize(t + 2, _hullStart)
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
        hullHash[_hashKey(point(FloatType,coords,e)...,_cx,_cy,_hashSize)] = e
    end

    hull = Vector{IntType}(undef,hullSize)
    e = _hullStart
    for i = 1:hullSize
        hull[i] = e
        e = hullNext[e]
    end

    # trim typed triangle mesh arrays
    # return (triangles=resize!(_triangles, (trianglesLen-1)), # ÷ 3), 
    #     halfedges = resize!(_halfedges, trianglesLen-1), hull)
    index_halfedges!(_edgeindex, @view(_halfedges[1:trianglesLen-1]), @view(_triangles[1:trianglesLen-1]); init=false)

    return _returntype(
                resize!(_tridata, (trianglesLen-1) ÷ 3), 
                resize!(_halfedges, trianglesLen-1), 
                hull, coords, _edgeindex)
end


function _hashKey(x, y, cx, cy, _hashSize)
    return (floor(Int, pseudoAngle(x - cx, y - cy) * _hashSize) % _hashSize)+1
end

function _legalize(::Type{FloatType}, a, triangles, halfedges, coords, EDGE_STACK, hullPrev, hullTri, _hullStart) where {FloatType}

    i = 1
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
        #a0 = a - a % 3
        a0 = (a-1) - ((a-1)%3)+ 1 # translate to 1 based indexes
        ar = (a0-1) + ((a-1) + 2) % 3 + 1 # translate to 1 based indexes

        if b == -1 #convex hull edge
            i == 1 && break
            i -= 1
            a = EDGE_STACK[i]
            continue
        end

        b0 = (b-1) - ((b-1) % 3) + 1
        al = (a0-1) + ((a-1) + 1) % 3 + 1
        bl = (b0-1) + ((b-1) + 2) % 3 + 1

        p0 = triangles[ar]
        pr = triangles[a]
        pl = triangles[al]
        p1 = triangles[bl]

        illegal = inCircle(point(FloatType,coords,p0)..., 
                           point(FloatType,coords,pr)..., 
                           point(FloatType,coords,pl)..., 
                           point(FloatType,coords,p1)...)

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
                        if hullTri[e] == bl
                            hullTri[e] = a
                            break
                        end
                        e = hullPrev[e]
                    end
                end
            end

            _link(halfedges, a, hbl)
            _link(halfedges, b, halfedges[ar])
            _link(halfedges, ar, bl)
            #@assert(_validate_halfedges(halfedges) == true)

            br = (b0-1) + ((b-1) + 1) % 3 + 1

            # don't worry about hitting the cap: it can only happen on extremely degenerate input
            if i < length(EDGE_STACK)
                EDGE_STACK[i] = br
                i += 1
            end
        else
            i == 1 && break
            i -= 1
            a = EDGE_STACK[i]
        end
    end

    return ar
end


function _link(_halfedges, a, b)
    #println("Link: $a, $b")
    #_halfedges[a] = b == -1 ? typemax(UInt32) : b
    _halfedges[a] = b
    if b != -1 
        _halfedges[b] = a
    end 
end



# add a new triangle given vertex indices and adjacent half-edge ids
function _addTriangle(_triangles, _halfedges, t, i0, i1, i2, a, b, c)

    _triangles[t] = i0
    _triangles[t + 1] = i1
    _triangles[t + 2] = i2

    _link(_halfedges, t, a)
    _link(_halfedges, t + 1, b)
    _link(_halfedges, t + 2, c)
    #@assert(_validate_halfedges(_halfedges) == true)

    return t + 3
end



end # module
