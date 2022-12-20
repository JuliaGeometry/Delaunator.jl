""" Calculate upper and lowerbounds on the coordinates in points. """
function _get_bounds(::Type{FloatType}, points) where FloatType
    minX = Inf
    minY = Inf
    maxX = -Inf
    maxY = -Inf

    for i in eachindex(points)
        x,y = point(FloatType,points, i)
        x < minX && (minX = x)
        y < minY && (minY = y)
        x > maxX && (maxX = x)
        y > maxY && (maxY = y)
    end

    return (minX,minY),(maxX,maxY)
end 

function _find_seeds(points, c::Tuple)
    FloatType = eltype(c)
    n = length(points)

    cx, cy = getX(c), getY(c)

    minDist = Inf
    i0 = 1
    i1 = 1
    i2 = 1
    

    # pick a seed point close to the center
    for i in 1:n
        d = dist(cx, cy, point(FloatType,points,i)...)
        if d < minDist
            i0 = i
            minDist = d
        end
    end
    i0x,i0y = point(FloatType,points,i0)
    minDist = Inf

    # find the point closest to the seed
    for i in 1:n
        i == i0 && continue
        d = dist(i0x, i0y, point(FloatType,points,i)...)
        if d < minDist && d > 0
            i1 = i
            minDist = d
        end
    end
    i1x,i1y = point(FloatType,points,i1)

    minRadius = Inf
    for i in 2:n
        i == i0 || i == i1 && continue
        r = circumradius(i0x, i0y, i1x, i1y, point(FloatType,points,i)...)
        if r < minRadius
            i2 = i
            minRadius = r
        end
    end

    return minRadius == Inf, i0, i1, i2
end 

function _delaunator!(
    coords,
    seeds::Tuple,
    _triangles,
    _halfedges,
    hullPrev, 
    hullNext,
    hullTri,
    hullHash,
    _hashSize::Integer,
    edgeStack,
    _ids, 
    _dists::AbstractVector{FloatType},
    tol
) where {FloatType}
    #FloatType = eltype(_dists)
    n = length(coords)
    legalize = (t,_hullStart)->_legalize(FloatType, t, _triangles, _halfedges, coords, edgeStack, hullPrev, hullTri, _hullStart)

    i0, i1, i2 = seeds
    i0x,i0y = point(FloatType,coords,i0)
    i1x,i1y = point(FloatType,coords,i1)
    i2x,i2y = point(FloatType,coords,i2)

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
    for k = 1:n
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
    
    return _hullStart, hullSize, trianglesLen ÷ 3
end 

