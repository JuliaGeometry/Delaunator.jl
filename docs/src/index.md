# Delaunator.jl


## API Overview

The API consists of a set of types for static triangulations, and getting 
info about the triangulation. 

```
t = Triangulation(points)
triangles(t) # get an iterator over the triangles 
cells(t, bbox) #  get an iterator over the nearest point / voronoi cells given the bounding box
edges(t) # get an iterator over the edges 
neighbors(t, i) # neighbors of point i
cellneighbors(t, i, [bbox]) # get cell neighbors given bbox 
hull(t) # get an iterator 
cellarea(t, i, [bbox]) # get cell area given bbox 

# searching methods
nearestpoint(t, p)
findtriangle(t, p)

# drawing methods to translate from data to polots... 
hullpoly(t)
trianglepolys(t) 
edgelines(t) 
```

There is also a set of types for more dynamic scenarios, where you may not
want all the computed information for speed. 
```
bt, cdata = basictriangulation(points; [maxpoints=Integer]) # initialize data structures 
bt, cdata = update!(bt, points, cdata) # after the points have been changed, may incur allocations
h = gethull(bt, cdata)
index = index_halfedges(bt, cdata)
circumcenters!(cc, bt) # compute circumcenters 
```




```





```@index
```

```@autodocs
Modules = [Delaunator]
```
