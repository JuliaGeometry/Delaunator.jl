

struct PointsFromMatrix{T <: AbstractMatrix, I1, I2} <: AbstractArray{Tuple{eltype(T),eltype(T)}, 1}
    A::T
end 

"""
    PointsFromMatrix(A [,i1=1,i2=2])

This implicitly extracts 2d point tuples from a matrix
using row indices i1 and i2 for the coordinates. The columns
of the matrix because individual points. 

    PointsFromMatrix(A) == vec(reinterpret(Tuple{Int,Int},A[1:2,:]))

This function can be used to transparently map a matrix into a Delaunator set of
points. (There is no copying involved).

Example
-------
```
A = reshape(1:20, 4, 5)
rval = triangulate(PointsFromMatrix(A))) # uses A[1:2,:]
````
"""
function PointsFromMatrix(A::T) where T <: AbstractMatrix 
    1 in axes(A,1) || throw(DimensionMismatch("1 must be in axes(A,1) (does A have enough rows?)"))
    2 in axes(A,1) || throw(DimensionMismatch("2 must be in axes(A,1) (does A have enough rows?)"))
    return PointsFromMatrix{T, Val{1}, Val{2}}(A) 
end 

function PointsFromMatrix(A::T, i1::Integer, i2::Integer) where T <: AbstractMatrix 
    i1 in axes(A,1) || throw(DimensionMismatch("$i1 must be in axes(A,1) (does A have enough rows?)"))
    i2 in axes(A,1) || throw(DimensionMismatch("$i2 must be in axes(A,1) (does A have enough rows?)"))
    return PointsFromMatrix{T, Val{i1}, Val{i2}}(A) 
end 

pointindices(::PointsFromMatrix{T,Val{I1},Val{I2}}) where {T,I1,I2} = (I1,I2)
import Base.getindex, Base.firstindex, Base.lastindex  
function Base.getindex(pts::PointsFromMatrix, i::Integer)
    i1,i2 = pointindices(pts)
    return (pts.A[i1,i],pts.A[i2,i])
end 
Base.firstindex(pts::PointsFromMatrix) = firstindex(pts.A, 2)
Base.lastindex(pts::PointsFromMatrix) = lastindex(pts.A, 2)

import Base.size,Base.IndexStyle 
Base.size(pts::PointsFromMatrix) = (size(pts.A,2),)
Base.IndexStyle(::Type{<:PointsFromMatrix}) = IndexLinear()

