using Delaunator
using Pkg
using Test
using JSON, GeometryBasics

function tsort(t::Tuple{I,I,I}) where {I<:Integer}
    a,b,c = t
    if a <= b <= c
        return (a,b,c)
    elseif a <= c <= b
        return (a,c,b)
    elseif b <= a <= c
        return (b,a,c)
    elseif b <= c <= a
        return (b,c,a)
    elseif c <= a <= b
        return (c, a, b)
    else
        assert(c <= b <= a)
        return (c,b,a)
    end
end 

@testset "simple inputs" begin 
    @test tsort.(triangulate([[0, 1], [1, 0], [1, 1]]).triangles) == [(1,2,3)]
    @test tsort.(triangulate([[0, 1], [1, 0], [0, 0], [1, 1]]).triangles) == [(1,2,3),(1,2,4)]
    triangulate([[516, 661], [369, 793], [426, 539], [273, 525], [204, 694], [747, 750], [454, 390]])
end 

points = map(x->Int.(x), JSON.parsefile("fixtures/ukraine.json"))
issue13 = map(x->Float64.(x), JSON.parsefile("fixtures/issue13.json"))
issue43 = map(x->Float64.(x), JSON.parsefile("fixtures/issue43.json"))
issue44 = map(x->Float64.(x), JSON.parsefile("fixtures/issue44.json"))
robustness1 = map(x->Float64.(x), JSON.parsefile("fixtures/robustness1.json"))
robustness2 = map(x->Float64.(x), JSON.parsefile("fixtures/robustness2.json"))

include("iterators.jl")

@testset "clipping" begin 
include("clipping.jl")
end 

@testset "helpers" begin
    @testset "PointsFromMatrix" begin 
        A = reshape(1:20, 4, 5)
        pts = PointsFromMatrix(A)
        @test pts == vec(reinterpret(Tuple{Int,Int},A[1:2,:]))

        pts = PointsFromMatrix(A, 3, 4)
        @test pts == vec(reinterpret(Tuple{Int,Int},A[3:4,:]))

        @test_throws DimensionMismatch pts = PointsFromMatrix(A, 7, 8)
    end 
end

# test('triangulates plain array', (t) => {
#     const d = new Delaunator([].concat(...points));
#     t.same(d.triangles, Delaunator.from(points).triangles);
#     t.end();
# });



# test('triangulates typed array', (t) => {
#     const d = new Delaunator(Float64Array.from([].concat(...points)));
#     t.same(d.triangles, Delaunator.from(points).triangles);
#     t.end();
# });

# test('constructor errors on invalid array', (t) => {
#     /* eslint no-new: 0 */
#     t.throws(() => {
#         new Delaunator({length: -1});
#     }, /Invalid typed array length/);
#     t.throws(() => {
#         new Delaunator(points);
#     }, /Expected coords to contain numbers/);
#     t.end();
# });

# test('produces correct triangulation', (t) => {
#     validate(t, points);
#     t.end();
# });

# test('produces correct triangulation after modifying coords in place', (t) => {
#     const d = Delaunator.from(points);

#     validate(t, points, d);
#     t.equal(d.trianglesLen, 5133);

#     const p = [80, 220];
#     d.coords[0] = p[0];
#     d.coords[1] = p[1];
#     const newPoints = [p].concat(points.slice(1));

#     d.update();
#     validate(t, newPoints, d);
#     t.equal(d.trianglesLen, 5139);

#     t.end();
# });




# test('returns empty triangulation for small number of points', (t) => {
#     let d = Delaunator.from([]);
#     t.same(d.triangles, []);
#     t.same(d.hull, []);
#     d = Delaunator.from(points.slice(0, 1));
#     t.same(d.triangles, []);
#     t.same(d.hull, [0]);
#     d = Delaunator.from(points.slice(0, 2));
#     t.same(d.triangles, []);
#     t.same(d.hull, [1, 0]); // [0, 1] is also correct
#     t.end();
# });

# test('returns empty triangulation for all-collinear input', (t) => {
#     const d = Delaunator.from([[0, 0], [1, 0], [3, 0], [2, 0]]);
#     t.same(d.triangles, []);
#     t.same(d.hull, [0, 1, 3, 2]); // [2, 3, 0, 1] is also correct
#     t.end();
# });

@testset "empty triangulation from all-collinear" begin 
    pts = [Point2f(0,0), Point2f(1,0), Point2f(3,0), Point2f(2,0)]
    d = triangulate(pts)
    @test isempty(d.triangles)
    x = Int.(d.hull)
    @test Int.(collect(d.hull)) == [1, 2, 4, 3]
end

# test('supports custom point format', (t) => {
#     const d = Delaunator.from(
#         [{x: 5, y: 5}, {x: 7, y: 5}, {x: 7, y: 6}],
#         p => p.x,
#         p => p.y);
#     t.same(d.triangles, [0, 2, 1]);
#     t.end();
# });

# function orient([px, py], [rx, ry], [qx, qy]) {
#     const l = (ry - py) * (qx - px);
#     const r = (rx - px) * (qy - py);
#     return Math.abs(l - r) >= 3.3306690738754716e-16 * Math.abs(l + r) ? l - r : 0;
# }
# function convex(r, q, p) {
#     return (orient(p, r, q) || orient(r, q, p) || orient(q, p, r)) >= 0;
# }

function orient(p, r, q)
    px,py = p
    rx,ry = r
    qx,qy = q
    l = (ry - py) * (qx - px)
    r = (rx - px) * (qy - py);
    return (abs(l - r) >= 3.3306690738754716e-16) * (abs(l + r) ? l - r : 0)
end 

function convex(r, q, p)
    return (orient(p, r, q) || orient(r, q, p) || orient(q, p, r)) >= 0 
end 

function validate_halfedges(d)
    for i in eachindex(d.halfedges)
        i2 = d.halfedges[i]
        if i2 != -1 && d.halfedges[i2] != i 
            return false 
        end
    end 
    return true
end 

function validate_dualcells(t)
    for i in eachindex(t) # for each point
        p = dualcell(t, i)
        @test contains(p, t.points[i]) == true 
    end
end 

function hull_area(points, d)
    j = lastindex(d.hull)
    hull_areas = Float64[] 
    for i in eachindex(d.hull)
        x0,y0 = points[d.hull[j]] 
        x,y = points[d.hull[i]]
        push!(hull_areas, (x - x0) * (y + y0)) # do shoelace area... 
        j = i
    end 
    return kahansum(hull_areas)
end 

function triangle_area(points, d)
    #tris = reshape(d.triangles, 3, length(d.triangles) ÷ 3)
    triareas = map(d.triangles) do tri 
        ax,ay = points[tri[1]]
        bx,by = points[tri[2]]
        cx,cy = points[tri[3]]
        return abs((by - ay) * (cx - bx) - (bx - ax) * (cy - by))
    end 
    return kahansum(triareas)
end 


# Kahan and Babuska summation, Neumaier variant; accumulates less FP error
function kahansum(x)
    sum = 0.0
    err = 0.0
    for i in eachindex(x)
        k = x[i]
        m = sum + k
        err += abs(sum) >= abs(k) ? sum - m + k : k - m + sum
        sum = m
    end
    return sum + err
end

function validate_area(points, d)
    hullarea = hull_area(points, d)
    triarea = triangle_area(points, d)
    @test abs(hullarea - triarea)/abs(hullarea) <= 2*eps(Float64)
end

function validate_clippedpolys(t)
    # make a bounding box that includes everything 
    bbox = Delaunator.margin_bbox(t, 0.05)
    for i in eachindex(t) # for each point
        p = dualcell(t, i)
        cpts = clippedpoly(p, bbox)
        if last(cpts) == first(cpts)
            cpts = cpts[1:end-1]
        end 
        FT = eltype(eltype(cpts))
        # make a finite polygon with all the points
        cp = Delaunator.InfinitePolygon(cpts, 
            (zero(FT),zero(FT)),(zero(FT),zero(FT)))
        @test contains(cp, t.points[i]) == true 
    end
end 

function validate(points, d; skip_dualcell=false, skip_clipped=false)
    @test validate_halfedges(d) == true
    validate_area(points, d)
    if !skip_dualcell
        validate_dualcells(d) 
        if !skip_clipped
            validate_clippedpolys(d)
        end 
    end 
end

validate(points; kwargs...) = validate(points, triangulate(points); kwargs...)

@testset "simple inputs" begin
    pts = [[0,0],[0,1],[1,0]]
    validate(pts)
end



# test('robustness', (t) => {
#     validate(t, robustness1);
#     validate(t, robustness1.map(p => [p[0] / 1e9, p[1] / 1e9]));
#     validate(t, robustness1.map(p => [p[0] / 100, p[1] / 100]));
#     validate(t, robustness1.map(p => [p[0] * 100, p[1] * 100]));
#     validate(t, robustness1.map(p => [p[0] * 1e9, p[1] * 1e9]));
#     validate(t, robustness2.slice(0, 100));
#     validate(t, robustness2);
#     t.end();
# });

@testset "js-delaunator robustness" begin
    # TODO, fix 
    pts = robustness1
    validate(pts; skip_dualcell=true)
    validate(map(x->(x[1]/1e9, x[2]/1e9), pts), skip_dualcell=true)
    validate(map(x->(x[1]/100, x[2]/100), pts), skip_dualcell=true)
    validate(map(x->(x[1]*100, x[2]*100), pts), skip_dualcell=true)
    validate(map(x->(x[1]*1e9, x[2]*1e9), pts), skip_dualcell=true)

    pts = robustness2
    validate(pts, skip_dualcell=true)
    validate(pts[1:100],skip_dualcell=true)
end


# test('issue #11', (t) => {
#     validate(t, [[516, 661], [369, 793], [426, 539], [273, 525], [204, 694], [747, 750], [454, 390]]);
#     t.end();
# });

# test('issue #13', (t) => {
#     validate(t, issue13);
#     t.end();
# });

# test('issue #24', (t) => {
#     validate(t, [[382, 302], [382, 328], [382, 205], [623, 175], [382, 188], [382, 284], [623, 87], [623, 341], [141, 227]]);
#     t.end();
# });

# test('issue #43', (t) => {
#     validate(t, issue43);
#     t.end();
# });

# test('issue #44', (t) => {
#     validate(t, issue44);
#     t.end();
# });




@testset "js-delaunator issues" begin 

    @testset "issue #21" begin
        pts = [[516, 661], [369, 793], [426, 539], [273, 525], [204, 694], [747, 750], [454, 390]]
        validate(pts)
    end

    @testset "issue #13" begin    
        validate(issue13)
    end


    @testset "issue #24" begin
        pts = [[382, 302], [382, 328], [382, 205], [623, 175], [382, 188], [382, 284], [623, 87], [623, 341], [141, 227]]
        validate(pts)
    end

    @testset "issue #43" begin    
        validate(issue44)
    end

    @testset "issue #44" begin    
        validate(issue44)
    end

end 

import Aqua
@testset "Aqua tests" begin 
    Aqua.test_all(Delaunator)
end

if VERSION >= v"1.8"
    @testset "JET tests" begin
        import JET
        projectinfo = Pkg.project()
        deps = Pkg.dependencies()
        @assert deps[projectinfo.dependencies["JET"]].version >= v"0.7.0"
        JET.report_package(Delaunator)
    end
end

