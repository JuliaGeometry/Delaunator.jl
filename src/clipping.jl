"""
    bbox = margin_bbox(t, [margin])
    bbox = margin_bbox(t, xmargin, ymargin)

Returns a bounding box (bbox) for the triangulation after applying a margin.
The default value of margin is 0.05. 
"""
margin_bbox(t::AbstractDelaunatorData, margin::Real=0.05) = margin_bbox(t, margin, margin)
function margin_bbox(t::AbstractDelaunatorData, xmargin::Real, ymargin::Real) 
    xsize, ysize = t._maxxy[1] - t._minxy[1], t._maxxy[2] - t._minxy[1]
    if abs(xsize) == 0 
        xsize = one(typeof(xsize))
    end 
    if abs(ysize) == 0 
        ysize = one(typeof(ysize))
    end 
    return t._minxy[1]-xmargin*xsize, t._minxy[2]-ymargin*ysize, t._maxxy[1]+xmargin*xsize, t._maxxy[2]+ymargin*ysize
end 

function canonicalize_bbox(p, bbox)
    FloatType = eltype(eltype(p))
    FloatType.(bbox)
end 

"""
    clippedpoly(p::InfinitePolygon, bbox; [closed=true])
    clippedpoly!(pts, p::InfinitePolygon, bbox) 
    
returns an empty array if the poly is entirely outside the bounding box.
Otherwise, return a set of points that represent edges of the infinite
polygon clipped to the bounding box. The set of points will 
be closed (where the first point is equal to the last) if the
closed=true option is set. The set of points _may_ be closed
even if this isn't set. Using closed=true results in simpler
behavior. This is not an option on the mutating version, see below
for how to get this functionality.  

The mutating version will update the pts array by using

    - `push!(pts, <newpt>)` 
    - `last(pts)`
    - `isempty(pts)`
    
It will return the input type pts. 

Example code
============
```
# generate polygon regions for all of the dualcells 
# clipped to a 5% expansion of the point bounding box
# as a list of NaN separated paths, with all 
# polygons closed. 
using GeometryBasics
t = triangulate(rand(Point2f, 15))
ppts = Point2f[]
for i in eachindex(t)
    ind = lastindex(ppts)
    clippedpoly!(ppts, dualcell(t, i), margin_bbox(t, 0.05))
    # check if the polygon was closed... 
    if lastindex(ppts) > ind # then we added point
        if ppts[ind+1] != ppts[end] # check if we are closed
            push!(ppts, ppts[ind+1]) # close the polygon
        end 
    end 
    push!(ppts, (NaN,NaN)) # add the NaN separator 
end 
```

See also [`dualcell`](@ref). 
"""
:clippedpoly, :clippedpoly!
function clippedpoly(p::InfinitePolygon, bbox; closed::Bool=true)
    pts = clippedpoly!(eltype(p)[], p, bbox) 
    if closed && length(pts) > 0 && pts[begin] != pts[end]
        push!(pts, pts[begin]) # close the polygon
    end 
    return pts
end 

function clippedpoly!(pts, p::InfinitePolygon, bbox)
    bbox = canonicalize_bbox(pts, bbox)
    if isfinite(p)
        return _clip_finite!(pts, p, bbox)
    else
        return _clip_infinite!(pts, p, bbox)
    end 
    return pts 
end 

"""
Compute intersections between the bbox and a point
and ray combo. Or on the line between two points. 
For a line intersection, use tmin=0, tmax=1.
For a ray intersection, use tmin=0, tmax=Inf. 

This will return the origin point (pt) if there 
are no intersections with the bbox. This will return
a duplicated point if there is only a single intersection.

These return a coordinate that is guaranteed to be on the bbox,
unless the return value is pt. 
""" 
function bbox_intersection(pt, ray, bbox; tmin=0, tmax=Inf)
    # we can have at most two intersections...

    xmin,ymin,xmax,ymax = bbox 
    vx,vy = ray
    x0,y0 = pt 
    t1 = Inf
    t2 = Inf 
    c = Inf 

    make_negative_inf(x) = x < 0 ? Inf : x 

    t1 = Inf
    t2 = Inf
    x1,y1,x2,y2 = Inf,Inf,Inf,Inf

    # try the line on ymin
    #handle_dimsion()
    #t1,t2,x1,y1,x1,y2 = 

    t = (ymin - y0) / vy 
    if t >= tmin && t <= tmax  
        # check that this is valid...
        x,y = x0 + t*vx, ymin
        if xmin <= x <= xmax
            if t1 < Inf
                # assign to t2 
                t2 = t
                x2,y2 = x,y
            else
                t1 = t
                x1,y1 = x,y
            end
        end
    end

    t = (ymax - y0) / vy 
    if t >= tmin && t <= tmax  
        # check that this is valid...
        x,y = x0 + t*vx, ymax
        if xmin <= x <= xmax
            if t1 < Inf
                # assign to t2 
                t2 = t
                x2,y2 = x,y
            else
                t1 = t
                x1,y1 = x,y
            end
        end
    end

    t = (xmin - x0) / vx
    if t >= tmin && t <= tmax  
        # check that this is valid...
        x,y = xmin, y0 + t*vy
        if ymin <= y <= ymax 
            if t1 < Inf
                # assign to t2 
                t2 = t
                x2,y2 = x,y
            else
                t1 = t
                x1,y1 = x,y
            end
        end
    end

    t = (xmax - x0) / vx
    if t >= tmin && t <= tmax  
        # check that this is valid...
        x,y = xmax, y0 + t*vy
        if ymin <= y <= ymax 
            if t1 < Inf
                # assign to t2 
                t2 = t
                x2,y2 = x,y
            else
                t1 = t
                x1,y1 = x,y
            end
        end
    end

    if t1 > t2 
        # swap
        x1,y1,x2,y2 = x2,y2,x1,y1
        t1,t2 = t2,t1 
    end 
    if t1 == Inf
        return (x0,y0),(x0,y0), t1, t2
    elseif t2 == Inf
        return (x1,y1),(x1,y1), t1, t2
    else
        return (x1,y1),(x2,y2), t1, t2
    end 
end 

function bbox_intersect_line(pa, pb, bbox)
    return bbox_intersection(pa, pb .- pa, bbox; tmin=0, tmax=1)
end

function bbox_intersect_dir(pa, pb, bbox)
    return bbox_intersection(pa, pb .- pa, bbox; tmin=0, tmax=Inf)
end 

# based on 
# https://observablehq.com/@mbostock/to-infinity-and-back-again
function _regioncode(x,y,xmin,ymin,xmax,ymax)
    # we need <= / >= here because
    # we combine regioncode / edgecode 
    rval = 0 
    if x <= xmin 
        rval = 1
    elseif x >= xmax
        rval = 2
    end
    if y <= ymin
        rval |= 4
    elseif y >= ymax
        rval |= 8 
    end
    return rval 
end 

_regioncode(p, lowerleft, upperright) = _regioncode(p...,lowerleft...,upperright...)
_regioncode(p, bbox) = _regioncode(p...,bbox...)

"""
Walk the corners by their region codes to find points
on the bounding box to add as we move from outside
the bbox back inside...

```
  # This is the order of codes in counter-clockwise order. 
  # Top-Left        Top         Top-Right
  # 1001     <-     1000    <-  1010
  #   |                           ^
  #   v                           |
  # 0001                         0010
  # Left                         Right 
  #   |                           ^
  #   v                           |  
  # 0101     ->    0100     ->  0110
  # Bottom Left   Bottom       Bottom right 
  # 
  # bottom right means x >= xmax (0010) and y <= ymin (0100)
```  

"""
function _add_bbox_points(pts, c0, c1, p, bbox)
    xmin,ymin,xmax,ymax = bbox 
    while c0 != c1
        if c0 == 0b0101 # this is bottom left... -> move to bottom
            c = (xmin, ymin) # test bottom left 
            c0 = 0b0100 # move to bottom
        elseif c0 == 0b0100 # bottom -> move to bottom-right
            c0 = 0b0110
            continue
        elseif c0 == 0b0110 # bottom right -> right 
            c = (xmax, ymin) # test bottom right 
            c0 = 0b0010
        elseif c0 == 0b0010 # right -> top right
            # and 
            c0 = 0b1010
            continue 
        elseif c0 == 0b1010 # top right -> top
            c = (xmax,ymax) # test if top right is in the poly
            c0 = 0b1000
        elseif c0 == 0b1000 # top  -> top left 
            c0 = 0b1001
            continue 
        elseif c0 == 0b1001 # top left -> left 
            # test if top left is in the poly
            c = (xmin, ymax)
            c0 = 0b0001  
        elseif c0 == 0b0001 # left -> bottom left
            # and test if bottom left in poly
            c0 = 0b0101
            continue 
        end
        if contains(p, c)
            _add_pt(pts, c) 
        end 
    end
end 

function _add_pt(pts, p1)
    if !isempty(pts)
        if last(pts) == p1
            return nothing 
        end
    end
    push!(pts, p1)
    return nothing
end 

function _test_and_add_bbox(pts, p, bbox)
    xmin,ymin,xmax,ymax = bbox 
    if contains(p, (0.5*xmin+0.5*xmax,0.5*ymin + 0.5*ymax)) 
        # then we need to add everything!
        _add_pt(pts, (xmin,ymin))
        _add_pt(pts, (xmax,ymin))
        _add_pt(pts, (xmax,ymax))
        _add_pt(pts, (xmin,ymax))
    end 
end

function _handle_segment(pts,p,bbox,p1,p2,c1,c2,lastoutsidecode,closedoutside)
    addpts = false 
    entryregion = 0 
    if c1 == 0 && c2 == 0 
        # this is all inside
        _add_pt(pts, p1) # this will check if it's a duplicate... 
        _add_pt(pts, p2) 
        addpts = true 
    elseif c1 != 0 && c2 == 0
        # then we come inside...
        p1b = bbox_intersect_line(p1, p2, bbox)[1]
        c1b = _regioncode(p1b, bbox)
        entryregion = c1b 
        if lastoutsidecode == 0 
            # then nothing else has been outside... 
            # so just project the outside to inside ray 
            _add_pt(pts, p1b)
            addpts = true 
        else
            # then there was another point outside...
            _add_bbox_points(pts, lastoutsidecode, c1b, p, bbox)
            _add_pt(pts, p1b)
            addpts = true 
            closedoutside = true 
        end
        _add_pt(pts, p2)
        zeropoints = false 
    elseif c1 == 0 && c2 != 0
        # then we go outside...
        # this is when we leave the region
        _add_pt(pts, p1)
        p2b = bbox_intersect_line(p1, p2, bbox)[1]
        _add_pt(pts, p2b) 
        addpts = true 
        c2b = _regioncode(p2b, bbox) 
        lastoutsidecode = c2b 
        closedoutside = false 
    else # c1 != 0 && c2 != 0 
        if c1 != c2 # then we may go through the box
            p1b,p2b,t1,t2 = bbox_intersect_line(p1, p2, bbox)
            if isfinite(t1) # then we have at least one interesction... 
                # but the return value in this case is p1b == p2b if
                # we only have one interesction, so we can just 
                # assume we have two. 
                c1b = _regioncode(p1b, bbox) 
                c2b = _regioncode(p2b, bbox) 
                _add_pt(pts, p1b)
                _add_bbox_points(pts, c1b, c2b, p, bbox)
                _add_pt(pts, p2b)
                addpts = true 
                entryregion = c1b 
                lastoutsidecode = c2b
                closedoutside = false 
            end
        end
    end 
    return addpts, lastoutsidecode, entryregion, closedoutside
end

function _clip_finite!(pts, p, bbox)
    lastoutsidecode = 0 
    closedoutside = true 
    lastcode = 0 
    anypoints = false

    for (p1,p2) in segments(p)
        c1 = _regioncode(p1, bbox)
        c2 = _regioncode(p2, bbox)
        addpts, lastoutsidecode, entryregion, closedoutside = _handle_segment(
                pts,p,bbox,p1,p2,c1,c2,lastoutsidecode,closedoutside)
        anypoints = anypoints | addpts
        lastcode = c2 
    end 

    if closedoutside == false
        # then we still need to close the polygon
        _add_bbox_points(pts, lastoutsidecode, lastcode, p, bbox)
    end 
    if anypoints == false 
        # this _could_ mean the polygon is entirely inside...
        # so test if the entire bbox is inside the polygon
        _test_and_add_bbox(pts, p, bbox)
    end 
    return pts
end 

function _clip_infinite!(pts, p, bbox)
    lastoutsidecode = 0 
    closedoutside = true 
    firstout2incode = 0 
    anypoints = false
    first = true 
    lastp1 = firstpoint(p) 
    lastp2 = firstpoint(p)

    for (p1,p2) in segments(p) 
        lastp1, lastp2 = p1, p2 # save the point 
        c1 = _regioncode(p1, bbox)
        c2 = _regioncode(p2, bbox)
        addpts = false 
        if first
            first = false 
            # this is the incoming ray...
            c2 = _regioncode(p2, bbox) 
            if c2 == 0 # then it's inside... 
                # this is reversed because we are looking at the ray
                # from p2 (inside) to p1 (direction...)
                p1b = bbox_intersect_dir(p2, p1, bbox)[1] 
                _add_pt(pts, p1b) 
                _add_pt(pts, p2)
                firstout2incode = _regioncode(p1b, bbox)
                addpts = true 
            else
                # it's possible it goes through entirely.
                # we reverse this because we are looking at the ray
                p2b,p1b,t1,t2 = bbox_intersect_dir(p2,p1,bbox) 
                if isfinite(t1) # then it touches the bbox...
                    _add_pt(pts, p1b)
                    _add_pt(pts, p2b)
                    closedoutside=false 
                    firstout2incode = _regioncode(p1b, bbox)
                    lastoutsidecode = _regioncode(p2b, bbox)
                    addpts = true 
                end
            end 
        else
            addpts, lastoutsidecode, entryregion, closedoutside = _handle_segment(
                pts,p,bbox,p1,p2,c1,c2,lastoutsidecode,closedoutside)
            if entryregion != 0 && firstout2incode == 0 
                firstout2incode = entryregion 
            end 
        end 
        anypoints = anypoints | addpts
        lastcode = c2 
    end

    # handle last ray...
    # the issue is we have already "sorta" handled it...
    # because we just used the last ray 
    # but treated it as finite. 
    # but we can just continue it in the same
    # direction but from p2 instead of p1 
    ray = lastp2 .- lastp1 # compute the p1 to p2 direction 
    ptail = lastp2  # now same ray but from p2 instead 
    ptail2 = ptail .+ ray 
    # now project this to the boundary

    # check if ptail -> ray 
    #if intersect_ray_bbox(lastp2, lastp2 .+ ray, )
    #doadd, padd = projectray_short(ptail, ray, bbox)
    pt1,pt2,t1,t2 = bbox_intersect_dir(ptail, ptail2, bbox) 
    if isfinite(t1) 
        if firstout2incode == 0 
            firstout2incode = _regioncode(pt1,bbox)
        end
        _add_pt(pts, pt1)
        ptail = pt1
    end
    if isfinite(t2) 
        _add_pt(pts, pt2) 
        ptail = pt2
    end 
    ctail = _regioncode(ptail, bbox)

    # need to add this point att then
    if firstout2incode == 0 
        # we never came in to the bbox 
        # this _could_ mean the polygon is entirely inside...
        # so test if the entire bbox is inside the polygon
        _test_and_add_bbox(pts, p, bbox)
    else
        _add_bbox_points(pts, ctail, firstout2incode, p, bbox)
    end 
    
    return pts
end