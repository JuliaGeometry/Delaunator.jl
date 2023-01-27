# Delaunator.jl


## API Overview

The API consists of a set of types for static triangulations, and getting 
info about the triangulation. 

```
t = Triangulation(points)
triangles(t) # get an iterator over the triangles 
triangles(t)
edges(t) # get an iterator over the edges 
neighbors(t, i) # neighbors of point i
hull(t) # get an iterator over the hull vertices
inhull(t, i)
dualcell(t, i) # get a description of the dual cell
dualcell(t, centers, i) # 


# drawing methods to translate from data to polots... 
hullpoly(t)
edgelines(t) 

# methods to get work with dualcell polys
p = dualcell(t, 1)
contains(p, pt)
segments(p)
clippedpoly(p, bbox)
```

There is also a set of types for more dynamic scenarios, where you may not
want all the computed information for speed. 
```
bt, cdata = basictriangulation(points; [maxpoints=Integer]) # initialize data structures 
bt, cdata = update!(bt, points, cdata) # after the points have been changed, may incur allocations
h = gethull(bt, cdata)
index = index_halfedges(bt, cdata)
```

Planned implementations
-----------------------
```
# cell diagram methods
cells(t, bbox) #  get an iterator over the nearest point / voronoi cells given the bounding box
bd = celldiagram(t [, centers]; [margin=0.05, boundingbox=mar()])
cellarea(t, i, [bbox]) # get cell area given bbox 
cellneighbors(t, i, [bbox]) # get cell neighbors given bbox 

# searching methods
nearestpoint(t, p)
findtriangle(t, p)
```





