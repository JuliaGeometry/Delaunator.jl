

struct PointsFromMatrix{T <: AbstractMatrix, I1, I2} <: AbstractArray{Tuple{eltype(T),eltype(T)}, 1}
    A::T
end 

"""
    PointsFromMatrix(A [,i1=1,i2=2])

Create an interface to a matrix type A to extract points from each column of A.
The optional values i1 and i2 give the row indices of the points in the column. 
This input can be used to transparently map a matrix into a Delaunator set of
points. (There is no copying involved)
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

