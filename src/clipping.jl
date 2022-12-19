
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