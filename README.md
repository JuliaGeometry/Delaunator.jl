# Delaunator

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sjkelly.github.io/Delaunator.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sjkelly.github.io/Delaunator.jl/dev)

[![Codecov](https://codecov.io/gh/sjkelly/Delaunator.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sjkelly/Delaunator.jl)
[![Coveralls](https://coveralls.io/repos/github/sjkelly/Delaunator.jl/badge.svg?branch=master)](https://coveralls.io/github/sjkelly/Delaunator.jl?branch=master)

The Delaunator algorithm computes a simple 2d triangulation of an arbitrary set of points in the plane _quickly_. This package provides a Julia implementation of the algorithm along with a number of supporting routines that operate on the Delaunator data structures. 

On a 2020 M1 Macbook Air, the Deluantor.jl package will triangulate a million points in about 500ms including computing circumcenters for the triangles for the closest site diagram. 

```
using GeometryBasics, Delaunator, StableRNGs, CairoMakie, Statistics, LinearAlgebra

# generate random points. 
ipts = randn(StableRNG(1), ComplexF64, 100)
ipts = sqrt.(abs.(ipts)).*ipts./abs.(ipts)
pts = Point2f.(zip(real(ipts),imag(ipts)))
t = triangulate(pts) # triangulate them! 
f = Figure()
ax = Axis(f[1,1])
hidespines!(ax)
hidedecorations!(ax) 
rdata = rand(StableRNG(1), length(triangles(t)))
# work over all triangles
for (i,tri) in enumerate(triangles(t))
  tripts = pts[collect(tri)]
  poly!(ax , tripts, color=[norm(mean(tripts))*(1+0.5*rdata[i])],
    colorrange=(0,3), 
    colormap=map(c->RGBAf(c.r,c.g,c.b,0.9), Makie.to_colormap(:matter)))
end
scatter!(ax, pts, color=:black, markersize=5)
f
```

Synopsis
--------
```
using Delaunator, GeometryBasics, StableRNGs
t = triangulate(rand(StableRNG(1), Point2f, 10))
@show hull(t) 
for i in eachindex(t) # iterate over each point index 
  @show inhull(t, i)
end 
for n in neighbors(t, 1) # iterate over neighbors of i
end 
for (i,j) in edges(t) # iterate over all edges in triangulation
end 

triangles(t, 1) # get the triangles that touch point 1. 

## 
using CairoMakie
f = linesegments(collect(edgelines(t)), color=:black, linewidth=0.75) # draw the edges 
poly!(f.axis, collect(hullpoly(t)), color=:transparent, strokewidth=1)

## dual cells, aka Voronoi cells (this interface is still a bit rough)
# get an infinite poly description of the points nearest cell i 
p = dualcell(t, 1) 
contains(p, (x,y)) # test if a point (x,y) is in the polygon
@show isfinite(p) # test if it's finite or infinite 
cp = clippedpoly(p, bbox) # produce a finite polygon
```

Example uses
------------
```
# Draw the key edges in the dual cells of a triangulation, aka the Voronoi diagram
pts = rand(StableRNG(1), Point2f, 15)
rval = triangulate(pts)
f = scatter(pts)
text!(f.axis, pts, text=map(i->"$i", 1:length(pts)))
for i in eachindex(rval) 
  p = Delaunator.dualcell(rval, i)
  linesegments!(f.axis, collect(segments(p)),
    xautolimits=false,yautolimits=false, color=:black)
end 
```

In the future
-------------
```
## Searching (NOT IMPLEMENTED)
findtriangle(t, pt) 
nearestpoint(t, pt) 

## nearest point cell diagram aka Voronoi diagram. 
# this will compute and store the info in a more compact way. 
# note that you can get all of this with the clippedpoly 
# functions now, except for cellarea... 
bc = boundedcells(t; margin=0.05) # get the nearest point cells
cellarea(bc, i)
cellpoly(bc, i) 
neighbors(bc, i) # get an iterator over neighbors of the bounded cells 
```

  