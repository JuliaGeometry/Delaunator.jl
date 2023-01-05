
"""
    boundedcells(t [,centers]; [kwargs...])

The bounded cell diagram, or the bounded Voronoi diagram, can be computed
from any of the various centers associated with a triangulation.
Using the circumcenters (the default) gives the nearest point, or Voronoi
diagram. 

These bounded diagrams need to be computed with respect to a bounding box,
which is determined by the keyword arguments

## Optional Arguments
- `margin=0.05` the default margin to add to the bounding box
- `xmargin=margin` the xaxis margin to add
- `ymargin=margin` the yaxis margin to add
- `boundingbox=box_with_margin(t, xmargin, ymargin)` the default bounding box
Note that providing bounding box overrides all of the other optional arguments.
"""
boundedcells(t::Triangulation;kwargs...) = boundedcells(t, circumcenters(t); kwargs...)
function boundedcells(t::Triangulation, centers; 
    margin=0.05, xmargin=margin, ymargin=margin, 
    boundingbox=box_with_margin(t, xmargin, ymargin))
    
end 

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