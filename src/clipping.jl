
# allow us to use the bounding box api for generic tuple bounding boxes... 
#import Base.maximum 
#origin(bbox::Tuple{Tuple{FT,FT},Tuple{FT,FT}}) where FT = bbox[1]
#maximum(bbox::Tuple{Tuple{FT,FT},Tuple{FT,FT}}) where FT = bbox[2]


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
    return (t._minxy[1]-xmargin*xsize, t._minxy[2]-ymargin*ysize), (t._maxxy[1]+xmargin*xsize, t._maxxy[2]+ymargin*ysize)
end 

function canonicalize_bbox(p, bbox)
    FloatType = eltype(eltype(p))
    FloatType.(bbox)
end 
"""
    clippedpoly(p::InfinitePolygon, bbox)
    clippedpoly!(pts, p::InfinitePolygon, bbox) 
    
returns an empty array if the poly is entirely outside the bounding box.
Otherwise, return a set of points that represent the infinite
polygon clipped to the bounding box. 

The mutating version will update the pts array by using
    - `push!(pts, <newpt>)` 
    - `last(pts)`
    - `isempty(pts)`
It will return the input type pts

Example code
============
```
# generate polygon regions for all of the dualcells 
# clipped to a 5% expansion of the point bounding box
# as a list of NaN separated paths. 
using GeometryBasics
t = triangulate(rand(Point2f, 15))
ppts = Point2f[]
for i in eachindex(t)
    clippedpoly!(ppts, dualcell(t, i), margin_bbox(t, 0.05))
    push!(ppts, NaN)
end 
```
"""
:clippedpoly, :clippedpoly!
function clippedpoly(p::InfinitePolygon, bbox)
    clippedpoly!(eltype(p)[], p, bbox) 
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


#=
https://stackoverflow.com/questions/9043805/test-if-two-lines-intersect-javascript-function

// returns true if the line from (a,b)->(c,d) intersects with (p,q)->(r,s)
function intersects(a,b,c,d,p,q,r,s) {
  var det, gamma, lambda;
  det = (c - a) * (s - q) - (r - p) * (d - b);
  if (det === 0) {
    return false;
  } else {
    lambda = ((s - q) * (r - a) + (p - r) * (s - b)) / det;
    gamma = ((b - d) * (r - a) + (c - a) * (s - b)) / det;
    return (0 < lambda && lambda < 1) && (0 < gamma && gamma < 1);
  }

  =#
function intersect_lines(a1,a2,b1,b2)
  a,b = a1 
  c,d = a2
  p,q = b1
  r,s = b2
  det = (c-a)*(s-q) - (r-p)*(d-b)
  if det == 0 
    return false
  else
    lambda = ((s - q) * (r - a) + (p - r) * (s - b)) / det;
    gamma = ((b - d) * (r - a) + (c - a) * (s - b)) / det;
    return (0 < lambda && lambda < 1) && (0 < gamma && gamma < 1)
  end 
end 

function intersect_bbox(p1, p2, bbox)
    xmin,ymin,xmax,ymax = bbox 
    if intersect_lines(p1,p2,(xmin,ymin),(xmax,ymin)) || 
        intersect_lines(p1,p2,(xmin,ymin),(xmin,ymax)) ||
        intersect_lines(p1,p2,(xmax,ymin),(xmax,ymax)) || 
        intersect_lines(p1,p2,(xmin,ymax),(xmax,ymax))
        return true
    else
        return false
    end 
end 

function dual_intersections(pt, ray, bbox) 
    xmin,ymin,xmax,ymax = bbox 
    vx,vy = ray
    x0,y0 = pt 
    t1 = Inf
    t2 = Inf 
    c = Inf 

    make_negative_inf(x) = x < 0 ? Inf : x 

    t1 = make_negative_inf( (ymin - y0) / vy )
    x1,y1 = x0 + t1*vx, ymin
    t2 = make_negative_inf( (ymax - y0) / vy )
    x2,y2 = x0 + t2*vx, ymax
    t3 = make_negative_inf( (xmin - x0) / vx ) 
    x3,y3 = xmin, y0 + t3*vy 
    t4 = make_negative_inf( (xmax - x0) / vx ) 
    x4,y4 = xmax, y0 + t4*vy 

    # need to find the smallest two points
    # that are positive (but all negative are Inf)
    # hopefully the compiler can optimize this,
    # see below for the codegen
    if t1 <= t2 <= t3 && t1 <= t2 <= t4 
        return (x1,y1), (x2,y2)
    elseif t1 <= t3 <= t2 && t1 <= t3 <= t4 
        return (x1,y1), (x3,y3)
    elseif t1 <= t4 <= t2 && t1 <= t4 <= t3 
        return (x1,y1), (x4,y4)
    elseif t2 <= t1 <= t3 && t2 <= t1 <= t4 
        return (x2,y2), (x1,y1)
    elseif t2 <= t3 <= t1 && t2 <= t3 <= t4 
        return (x2,y2), (x3,y3)
    elseif t2 <= t4 <= t1 && t2 <= t4 <= t3 
        return (x2,y2), (x4,y4)
    elseif t3 <= t1 <= t2 && t3 <= t1 <= t4 
        return (x3,y3), (x1,y1)
    elseif t3 <= t2 <= t1 && t3 <= t2 <= t4 
        return (x3,y3), (x2,y2)
    elseif t3 <= t4 <= t1 && t3 <= t4 <= t2 
        return (x3,y3), (x4,y4)
    elseif t4 <= t1 <= t2 && t4 <= t1 <= t3 
        return (x4,y4), (x1,y1)
    elseif t4 <= t2 <= t1 && t4 <= t2 <= t3 
        return (x4,y4), (x2,y2)
    elseif t4 <= t3 <= t1 && t4 <= t3 <= t2 
        return (x4,y4), (x3,y3)
    end 
end 

#=
function write_comparison()
  for t1=1:4, t2=1:4
    if t1==t2 
      continue
    end
    t3, t4 = 0, 0 
    for t=1:4
      if t == t1 || t == t2 
        continue
      end 
      t3 = t 
      break 
    end 
    for t=1:4
      if t == t1 || t == t2 || t == t3 
        continue
      end 
      t4 = t 
      break 
    end 
    println("    elseif t$(t1) <= t$(t2) <= t$(t3) && t$(t1) <= t$(t2) <= t$(t4) ")
    println("        return (x$(t1),y$(t1)), (x$(t2),y$(t2))")
  end 
end
write_comparison()  
=#    

function projectray_short(pt, ray, bbox)
    xmin,ymin,xmax,ymax = bbox 
    vx,vy = ray
    x0,y0 = pt 
    t = Inf 
    c = Inf 

    # check if it misses entirely.
    if y0 < ymin && vy < 0 
        return false, pt
    elseif y0 > ymax && vy > 0 
        return false, pt
    elseif x0 < xmin && vx < 0 
        return false, pt
    elseif x0 > xmax && vx > 0 
        return false, pt
    end
    
    c = (ymin - y0) / vy
    if c >= 0 && c < t
        xt = x0+c*vx
        if xmin <= xt <= xmax 
            t = c 
            x,y = xt, ymin
        end
    end

    c = (ymax - y0) / vy
    if c >= 0 && c < t
        xt = x0+c*vx
        if xmin <= xt <= xmax 
            t = c 
            x,y = xt, ymax
        end
    end

    c = (xmin - x0) / vx
    if c >= 0 && c < t
        yt = y0 + c*vy 
        if ymin <= yt <= ymax 
            t = c 
            x,y = xmin, yt
        end 
    end

    c = (xmax - x0) / vx
    if c >= 0 && c < t
        yt = y0 + c*vy 
        if ymin <= yt <= ymax 
            t = c 
            x,y = xmax, yt
        end 
    end
    
    return true, (x,y)
end 

function projectray(pt, ray, bbox)
    xmin,ymin,xmax,ymax = bbox 
    vx,vy = ray
    x0,y0 = pt 
    t = Inf 
    c = Inf 
    if vy < 0 
        if (y0 <= ymin)
            return false, pt
        else
            c = (ymin - y0) / vy
            if c < t 
                t = c 
                x,y = x0+c*vx, ymin
            end
        end
    else
        if (y0 >= ymax)
            return false, pt 
        else
            c = (ymax - y0) / vy 
            if c < t
                t = c
                x,y = x0+c*vx, ymax
            end 
        end
    end

    if vx < 0 # clip on left
        if (x0 <= xmin)
            return false, pt
        else
            c = (xmin - x0) / vx
            if c < t
                t = c
                x,y = xmin, y0 + c*vy 
            end
        end 
    else
        if (x0 >= xmax)
            return false, pt 
        else
            c = (xmax - x0) / vx
            if c < t
                t = c
                x,y = xmax, y0 + c*vy
            end
        end
    end

    return true, (x,y)
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
#=
function _add_bbox_points(pts, c0, c1, p, bbox)
    xmin,ymin,xmax,ymax = bbox 
    while c0 != c1
        if c0 == 0b0101 # this is bottom left... -> move to bottom
            c0 = 0b0100 # move to bottom
            continue 
        elseif c0 == 0b0100 # bottom -> move to bottom-right
            # and test if bottom right is in the poly
            c = (xmax, ymin)
            c0 = 0b0110
        elseif c0 == 0b0110 # bottom right -> right 
            c0 = 0b0010
            continue
        elseif c0 == 0b0010 # right -> top right
            # and test if top right is in the poly
            c = (xmax,ymax)
            c0 = 0b1010
        elseif c0 == 0b1010 # top right -> top
            c0 = 0b1000
            continue
        elseif c0 == 0b1000 # top  -> top left 
            # and test if top left is in the poly
            c = (xmin, ymax)
            c0 = 0b1001
        elseif c0 == 0b1001 # top left -> left 
            c0 = 0b0001  
            continue
        elseif c0 == 0b0001 # left -> bottom left
          # and test if bottom left in poly
            c = (xmin, ymin)
            c0 = 0b0101
        end
        if contains(p, c)
            _add_pt(pts, c) 
        end 
    end
end 
=#

function _add_pt(pts, p1)
    if !isempty(pts)
        if last(pts) == p1
            return nothing 
        end
    end
    push!(pts, p1)
    return nothing
end 

function _clip_finite!(pts, p, bbox)
    lastoutsidecode = 0 
    closedoutside = true 
    lastcode = 0 
    for (p1,p2) in segments(p) 
        c1 = _regioncode(p1, bbox)
        c2 = _regioncode(p2, bbox)
        lastcode = c2 
        if c1 == 0 && c2 == 0 
            # all inside
            _add_pt(pts, p1) # this will check if it's a duplicate... 
            _add_pt(pts, p2) 
        elseif c1 != 0 && c2 == 0 
            # then we come back inside... 
            p1b = projectray_short(p1, p2 .- p1, bbox)[2] 
            c1b = _regioncode(p1b, bbox)
            if lastoutsidecode == 0 
                # then nothing else has been outside... 
                # so just project the outside to inside ray 
                _add_pt(pts, p1b)
            else
                # then there was another point outside...
                _add_bbox_points(pts, lastoutsidecode, c1b, p, bbox)
                _add_pt(pts, p1b)
            end 
            _add_pt(pts, p2)
            closedoutside = true 
        elseif c1 == 0 && c2 != 0 
            # this is when we leave the region
            _add_pt(pts, p1)
            p2b = projectray_short(p1, p2 .- p1, bbox)[2] 
            _add_pt(pts, p2b) 
            c2b = _regioncode(p2b, bbox) 
            lastoutsidecode = c2b 
            closedoutside = false 
        else
            # if both are outside, don't do anything...
            # need to check if the ray intersects the box...
            if c1 != c2 # then we go through the box
                if intersect_bbox(p1, p2, bbox) 
                    # find potentailly both intersections...
                    p1b, p2b = dual_intersections(p1, p2 .- p1, bbox)
                    if isfinite(p1b[1]) && isfinite(p1b[2])
                        # at least one of these pointos must be finite...
                        c1b = _regioncode(p1b, bbox) 
                        c2b = _regioncode(p2b, bbox) 
                        _add_pt(pts, p1b)
                        _add_bbox_points(pts, c1b, c2b, p, bbox)
                        _add_pt(pts, p2b)
                    end
                    lastoutsidecode = c2b
                    closedoutside = false 
                end 
            end 
        end
    end 
    if closedoutside == false
        # then we still need to close the polygon
        _add_bbox_points(pts, lastoutsidecode, lastcode, p, bbox)
    end 
    return pts
end 

function _clip_infinite!(pts, p, bbox)
    lastoutsidecode = 0 
    closedoutside = true 
    firstout2incode = 0 
    first = true 
    lastp1 = firstpoint(p) 
    lastp2 = firstpoint(p)
    for (p1,p2) in segments(p) 
        lastp1, lastp2 = p1, p2 # save the point 
        c1 = _regioncode(p1, bbox)
        c2 = _regioncode(p2, bbox)
        if first 
            first = false 
            # this is the incoming ray...
            c2 = _regioncode(p2, bbox) 
            if c2 == 0 # then it's inside... 
                # then we need to project it.
                p1b = projectray_short(p2, p1 .- p2, bbox)[2] 
                _add_pt(pts, p1b) 
                _add_pt(pts, p2)
                firstout2incode = _regioncode(p1b, bbox)
            end
        else
            if c1 == 0 && c2 == 0 
                # all inside
                _add_pt(pts, p1) # this will check if it's a duplicate... 
                _add_pt(pts, p2) 
            elseif c1 != 0 && c2 == 0 
                # then we come back inside... 
                p1b = projectray_short(p1, p2 .- p1, bbox)[2] 
                c1b = _regioncode(p1b, bbox)
                if firstout2incode == 0 
                    firstout2incode = c1b 
                end 

                if lastoutsidecode == 0 
                    # then nothing else has been outside... 
                    # so just project the outside to inside ray 
                    _add_pt(pts, p1b)
                else
                    # then there was another point outside...
                    _add_bbox_points(pts, lastoutsidecode, c1b, p, bbox)
                    _add_pt(pts, p1b)
                end 
                _add_pt(pts, p2)
                closedoutside = true 
            elseif c1 == 0 && c2 != 0 
                # this is when we leave the region
                _add_pt(pts, p1)
                p2b = projectray_short(p1, p2 .- p1, bbox)[2] 
                _add_pt(pts, p2b) 
                c2b = _regioncode(p2b, bbox) 
                lastoutsidecode = c2b 
                closedoutside = false 
            else c1 != 0 && c2 != 0 
                if c1 != c2 # then we may through the box
                    if intersect_bbox(p1, p2, bbox) 
                        # find potentailly both intersections...
                        p1b, p2b = dual_intersections(p1, p2 .- p1, bbox)
                        if isfinite(p1b[1]) && isfinite(p1b[2])
                            # at least one of these pointos must be finite...
                            c1b = _regioncode(p1b, bbox) 
                            c2b = _regioncode(p2b, bbox) 
                            _add_pt(pts, p1b)
                            _add_bbox_points(pts, c1b, c2b, p, bbox)
                            _add_pt(pts, p2b)
                        end
                        lastoutsidecode = c2b
                        closedoutside = false 
                    end 
                end 
            end 
        end 
    end 
    # handle last ray...
    # the issue is we have already "sorta" handled it...
    # because we just used the last ray 
    # but treated it as finite. 
    # but we can just continue it in the same
    # direction but from p2 instead of p1 
    ray = lastp2 .- lastp1 # compute the p1 to p2 direction 
    ptail = lastp2  # now same ray but from p2 instead 
    # now project this to the boundar
    doadd, padd = projectray(ptail, ray, bbox)
    if doadd
        _add_pt(pts, padd)
        ptail = padd 
    end 
    
    ctail = _regioncode(ptail, bbox)
    # need to add this point att then
    if firstout2incode == 0 
        # we never came in to the bbox 
        xmin,ymin,xmax,ymax = bbox 
        if contains(p, (0.5*xmin+0.5*xmax,0.5*ymin + 0.5*ymax)) 
            # then we need to add everything!
            _add_pt(pt, (xmin,ymin))
            _add_pt(pt, (xmax,ymin))
            _add_pt(pt, (xmax,ymax))
            _add_pt(pt, (xmin,ymax))
        end 
    else
        _add_bbox_points(pts, ctail, firstout2incode, p, bbox)
    end 
    
    return pts
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
