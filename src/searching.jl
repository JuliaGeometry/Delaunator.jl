using GeometryBasics
using Delaunator


function _prevedge(i::Integer)
    (i-1) % 3 == 0 ? i+2 : i -1
end 

function findtriangle(t, p)
  current = 1
  start = current
  n = 0

  while true
      next = Delaunator._nextedge(current)
      pc = t.points[(t._triangles[current])]
      pn = t.points[(t._triangles[next])]

      o = Delaunator.orientIfSure(pc[1], pc[2],pn[1], pn[2],p[1], p[2])

      if o >= 0
          current = next
          if start == current
              break
          end 
      else 
          current = t.halfedges[current]
          if current == -1
              break 
          end 
          n = n + 1
          if (n % 2 == 1)
              current = Delaunator._nextedge(current)
          else 
              current = Delaunator._prevedge(current)
          end 
          start = current
      end 
  end
  return current == -1 ? -1 : div(current, 3) +1 
end
