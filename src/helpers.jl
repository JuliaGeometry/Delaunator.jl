
# This was used to validate halfedges for debugging, can be removed soon.
function _validate_halfedges(halfedges)
  for i in eachindex(halfedges)
      i2 = halfedges[i]
      if i2 != -1 && halfedges[i2] != i 
          @show halfedges
          return false 
      end
  end 
  return true
end 