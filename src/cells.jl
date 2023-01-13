



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
