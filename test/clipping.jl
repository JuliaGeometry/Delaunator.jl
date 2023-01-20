@testset "clippedpoly! demo" begin 
  # generate polygon regions for all of the dualcells 
  # clipped to a 5% expansion of the point bounding box
  # as a list of NaN separated paths, with all 
  # polygons closed. 
  using GeometryBasics, StableRNGs
  t = triangulate(rand(StableRNG(1), Point2f, 15))
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
end 

@testset "clipping with all inside or outside" begin

  p = Delaunator.InfinitePolygon([(-1.0,-1.0),(1.0,-1.0),(1.0,1.0),(-1.0,1.0)], (0.0,0.0), (0.0,0.0))
  bbox = (-0.5,-0.5,0.5,0.5)
  pts = clippedpoly(p, bbox) 
  @test pts==([(-0.5,-0.5),(0.5,-0.5),(0.5,0.5),(-0.5,0.5),(-0.5,-0.5)])

  bbox = (-0.5,1.5,0.5,2.5)
  pts = clippedpoly(p, bbox) 
  @test isempty(pts)

  # this is a vertical division (if -5, it's all outside...)
  # if +5, it's all inside 
  p = Delaunator.InfinitePolygon([(-5.0,0.0)], (0.0,-1.0), (0.0,1.0))
  pts = clippedpoly(p, (-2,-2,2,2))
  @test isempty(pts) 

  p = Delaunator.InfinitePolygon([(5.0,0.0)], (0.0,-1.0), (0.0,1.0))
  pts = clippedpoly(p, (-2,-2,2,2))
  @test pts==([(-2.0,-2.0),(2.0,-2.0),(2.0,2.0),(-2.0,2.0),(-2.0,-2.0)])

  # Reverse the direction... that'll flip the inerpretation.
  p = Delaunator.InfinitePolygon([(-5.0,0.0)], (0.0,1.0), (0.0,-1.0))
  pts = clippedpoly(p, (-2,-2,2,2))
  @test pts==([(-2.0,-2.0),(2.0,-2.0),(2.0,2.0),(-2.0,2.0),(-2.0,-2.0)])

  p = Delaunator.InfinitePolygon([(5.0,0.0)], (0.0,1.0), (0.0,-1.0))
  pts = clippedpoly(p, (-2,-2,2,2))
  @test isempty(pts) 
  
  # this is a horizontal division (+5 is outside, -5 is inside now...)
  p = Delaunator.InfinitePolygon([(0.0,5.0)], (-1.0,0.0), (1.0,0.0))
  pts = clippedpoly(p, (-2,-2,2,2))
  @test isempty(pts) 

  p = Delaunator.InfinitePolygon([(0.0,-5.0)], (-1.0,0.0), (1.0,0.0))
  pts = clippedpoly(p, (-2,-2,2,2))
  @test pts==([(-2.0,-2.0),(2.0,-2.0),(2.0,2.0),(-2.0,2.0),(-2.0,-2.0)])
end 

function in_bbox(pt, bbox)
  xmin,ymin,xmax,ymax = bbox 
  return xmin <= pt[1] <= xmax && ymin <= pt[2] <= ymax 
end 


@testset "previous bugs" begin 
    # this is a horizontal ray on the x-axis... 
    p = Delaunator.InfinitePolygon([(-5.0,0.0)], (-1.0,0.0), (1.0,0.0))
    pts = clippedpoly(p, (-2,-2,2,2))
    @test pts==([(-2.0,0.0),(2.0,0.0),(2.0,2.0),(-2.0,2.0),(-2.0,0.0)])

    # this reverses the direction... so flips the side.
    p = Delaunator.InfinitePolygon([(-5.0,0.0)], (1.0,0.0), (-1.0,0.0))
    pts = clippedpoly(p, (-2,-2,2,2))
    @test pts==([(2.0,0.0),(-2.0,0.0),(-2.0,-2.0),(2.0,-2.0),(2.0,0.0)])


    ipts = randn(StableRNG(1), ComplexF64, 100)
    ipts = sqrt.(abs.(ipts)).*ipts./abs.(ipts)
    pts = Point2f.(zip(real(ipts),imag(ipts)))
    t = triangulate(pts)
    bbox = margin_bbox(t, 0.05) 
    ppts = clippedpoly(dualcell(t, 31), bbox)
    for pt in ppts 
      @test in_bbox(pt, bbox) 
    end 

end 

@testset "bbox intersections" begin
  # make sure we don't intersect ...
  @test Delaunator.bbox_intersection((5.0,0.0), (-1.0,0.0), (-2.0,-2.0,2.0,2.0); tmin=0, tmax=1) == ((5.0, 0.0), (5.0, 0.0))
  # but now we do
  @test Delaunator.bbox_intersection((5.0,0.0), (-1.0,0.0), (-2.0,-2.0,2.0,2.0); tmin=0, tmax=Inf) == ((2.0, 0.0), (-2.0, 0.0))
end