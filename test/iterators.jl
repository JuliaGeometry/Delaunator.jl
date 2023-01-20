@testset "triangle neighbors" begin 
  t = triangulate([[0, 1], [1, 0], [0, 0], [1, 1]])
  # triangles = [(1,2,3),(1,2,4)]
  @test Set(triangles(t, 1)) == Set([1,2])
  @test Set(triangles(t, 2)) == Set([1,2])
  @test Set(triangles(t, 3)) == Set([1])
  @test Set(triangles(t, 4)) == Set([2])
end 