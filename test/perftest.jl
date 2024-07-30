using GeometryBasics, StableRNGs, Delaunator, BenchmarkTools
pts = rand(StableRNG(1), Point2f,  1_000_000)

@time triangulate(pts);
@time triangulate(pts);
@time triangulate(pts);
@btime triangulate($pts); 
@btime triangulate(newpts) setup=begin newpts=rand(StableRNG(1), Point2f,  1_000_000); end; 

# Before robust predicate, took 0.9 - 0.98 seconds on Intel Mac