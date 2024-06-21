#=
Original code by 
Jonathan Richard Shewchuk
Computer Science Division
University of California at Berkeley
Berkeley, California 94720-1776

@article{shewchuk97a,
author = {Jonathan Richard Shewchuk},
title = {Adaptive {P}recision {F}loating-{P}oint {A}rithmetic and {F}ast
  {R}obust {G}eometric {P}redicates},
journal = {Discrete \& Computational Geometry},
volume = 18,
number = 3,
pages = {305--363},
month = oct,
year = 1997
}

Originally ported to Julia by 
vchuravy Valentin Churavy
in AdaptivePredicates.jl

The MIT License (MIT)

Copyright © 2024: Valentin Churavy, and other contributors

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

function find_epsilon()
    every_other = true
    half = 0.5
    epsilon = 1.0
    splitter = 1.0
    check = 1.0
    # /* Repeatedly divide `epsilon' by two until it is too small to add to    */
    # /*   one without causing roundoff.  (Also check if the sum is equal to   */
    # /*   the previous sum, for machines that round up instead of using exact */
    # /*   rounding.  Not that this library will work on such machines anyway. */
    cond = true
    while cond
      lastcheck = check
      epsilon *= half
      if (every_other)
        splitter *= 2.0
      end
      every_other = !every_other
      check = 1.0 + epsilon
      cond = ((check != 1.0) && (check != lastcheck))
    end
    splitter += 1.0
    return epsilon, splitter
end

const epsilon, splitter = find_epsilon()
@assert epsilon == eps(1.0) / 2 
const resulterrbound = (3.0 + 8.0 * epsilon) * epsilon
const ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon
const ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon
const ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon
const o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon
const o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon
const o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon
const iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon
const iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon
const iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon
const isperrboundA = (16.0 + 224.0 * epsilon) * epsilon
const isperrboundB = (5.0 + 72.0 * epsilon) * epsilon
const isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon

robust_orient(ax,ay,bx,by,cx,cy) = orient((ax,ay),(bx,by),(cx,cy))

function orient(a, b, c)
  dir = orient2(a, b, c)
  if dir > 0.0
    return true
  else
    return false 
  end 
end

"""
Return a positive value if the points pa, pb, and pc occur
in counterclockwise order a negative value if they occur
in clockwise order and zero if they are collinear.  The
result is also a rough approximation of twice the signed
area of the triangle defined by the three points.                                                             
"""
function orient2(pa, pb, pc)::Float64
  detleft = (pa[1] - pc[1]) * (pb[2] - pc[2])
  detright = (pa[2] - pc[2]) * (pb[1] - pc[1])
  det = detleft - detright

  if (detleft > 0.0)
    if (detright <= 0.0)
      return det
    else
      detsum = detleft + detright
	end
  elseif (detleft < 0.0)
    if (detright >= 0.0) 
      return det
    else 
      detsum = -detleft - detright
	end
  else
    return det
  end

  errbound = ccwerrboundA * detsum
  if (det >= errbound) || (-det >= errbound)
    return det
  end

  return orient2adapt(pa, pb, pc, detsum)
end

function orient2adapt(pa, pb, pc, detsum::Float64)::Float64
  acx = (pa[1] - pc[1])
  bcx = (pb[1] - pc[1])
  acy = (pa[2] - pc[2])
  bcy = (pb[2] - pc[2])

  detleft, detlefttail = Two_Product(acx, bcy)
  detright, detrighttail = Two_Product(acy, bcx)

  B3, B2, B1, B0 = Two_Two_Diff(detleft, detlefttail, detright, detrighttail)
  B = (B0, B1, B2, B3)

  det = estimate(4, B)
  errbound = ccwerrboundB * detsum
  if ((det >= errbound) || (-det >= errbound))
    return det
  end

  acxtail = Two_Diff_Tail(pa[1], pc[1], acx)
  bcxtail = Two_Diff_Tail(pb[1], pc[1], bcx)
  acytail = Two_Diff_Tail(pa[2], pc[2], acy)
  bcytail = Two_Diff_Tail(pb[2], pc[2], bcy)

  if ((acxtail == 0.0) && (acytail == 0.0)
      && (bcxtail == 0.0) && (bcytail == 0.0))
    return det
  end

  errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det)
  det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail)
  if ((det >= errbound) || (-det >= errbound))
    return det
  end

  s1, s0 = Two_Product(acxtail, bcy)
  t1, t0 = Two_Product(acytail, bcx)
  u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
  u = (u0, u1, u2, u3)
  C1, C1length = fast_expansion_sum_zeroelim(4, B, 4, u, Val(8))

  s1, s0 = Two_Product(acx, bcytail)
  t1, t0 = Two_Product(acy, bcxtail)
  u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
  u = (u0, u1, u2, u3)
  C2, C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, Val(12))

  s1, s0 = Two_Product(acxtail, bcytail)
  t1, t0 = Two_Product(acytail, bcxtail)
  u3, u2, u1, u0 = Two_Two_Diff(s1, s0, t1, t0)
  u = (u0, u1, u2, u3)
  D, Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, Val(16))  
  return @inbounds D[Dlength] # originally Dlength - 1 
end


function Split(a)
  c = (splitter * a)
  abig = (c - a)
  ahi = c - abig
  alo = a - ahi
  return ahi, alo
end

function Two_Product_Tail(a, b, x)
  ahi, alo = Split(a)
  bhi, blo = Split(b)
  err1 = x - (ahi * bhi)
  err2 = err1 - (alo * bhi)
  err3 = err2 - (ahi * blo)
  y = (alo * blo) - err3
  return y
end

function Two_Product(a, b)
	x = a * b
	return x, Two_Product_Tail(a, b, x)
end

function Two_Diff_Tail(a, b, x)
  bvirt = (a - x)
  avirt = x + bvirt
  bround = bvirt - b
  around = a - avirt
  y = around + bround
  return y
end

function Two_Diff(a, b)
  x = (a - b)
  return x, Two_Diff_Tail(a, b, x)
end

function Two_Sum_Tail(a, b, x)
  bvirt = (x - a)
  avirt = x - bvirt
  bround = b - bvirt
  around = a - avirt
  y = around + bround
  return y
end

function Two_Sum(a, b) 
  x = (a + b)
  return x, Two_Sum_Tail(a, b, x)
end

function Fast_Two_Sum_Tail(a, b, x)
    bvirt = x - a
    y = b - bvirt
    return y
end

function Fast_Two_Sum(a, b)
    x = (a + b)
    return x, Fast_Two_Sum_Tail(a, b, x)
end

function Two_One_Diff(a1, a0, b)
	i, x0 = Two_Diff(a0, b)
	x2, x1 = Two_Sum(a1, i)
	return x2, x1, x0
end

function Two_Two_Diff(a1, a0, b1, b0)
  _j, _0, x0 = Two_One_Diff(a1, a0, b0)
  x3, x2, x1 = Two_One_Diff(_j, _0, b1)
  return x3, x2, x1, x0
end

Absolute(a) = ((a) >= 0.0 ? (a) : -(a))

function estimate(elen, e)
  Q = e[1]
  for eindex in 2:elen
    Q += e[eindex]
  end
  return Q
end

"""
    Sets h = e + f.  See the long version of my paper for details. 

	If round-to-even is used (as with IEEE 754), maintains the strongly
	nonoverlapping property.  (That is, if e is strongly nonoverlapping, h
	will be also.)  Does NOT maintain the nonoverlapping or nonadjacent
	properties.                  
"""
@inline function fast_expansion_sum_zeroelim(elen::Int, e, flen::Int, f, ::Val{N})::Tuple{NTuple,Int} where {N}
  h = ntuple(i -> zero(typeof(e[1])), Val(N))
  @inbounds begin
    enow = e[1]
    fnow = f[1]
    eindex = findex = 1
    if ((fnow > enow) == (fnow > -enow))
      Q = enow
      enow = e[eindex += 1]
    else
      Q = fnow
      fnow = f[findex += 1]
    end
    hindex = 0
    if ((eindex < elen) && (findex < flen)) # still < since pre-increment
      if ((fnow > enow) == (fnow > -enow))
        Qnew, hh = Fast_Two_Sum(enow, Q)
        enow = e[eindex += 1]
      else
        Qnew, hh = Fast_Two_Sum(fnow, Q)
        fnow = f[findex += 1]
      end
      Q = Qnew
      if hh != 0.0
        #h[hindex+=1] = hh
        h = Base.setindex(h, hh, hindex+=1)
      end

      while (eindex < elen) && (findex < flen)  # still < since pre-increment
        if (fnow > enow) == (fnow > -enow) 
          Qnew, hh = Two_Sum(Q, enow)
          enow = e[eindex += 1]
        else 
          Qnew, hh = Two_Sum(Q, fnow)
          fnow = f[findex += 1]
        end
        Q = Qnew
        if hh != 0.0
          #h[hindex += 1] = hh
          h = Base.setindex(h, hh, hindex+=1)
        end
      end
    end

    while eindex <= elen
      Qnew, hh = Two_Sum(Q, enow)
      eindex += 1
      # We need an extra iteration to calculate Q
      # but we don't want to access e
      if eindex <= elen
          enow = e[eindex]
      end
      Q = Qnew
      if hh != 0.0
        #h[hindex += 1] = hh
        h = Base.setindex(h, hh, hindex+=1)
      end
    end

    while findex <= flen
      Qnew, hh = Two_Sum(Q, fnow)
      findex += 1
      if findex <= flen
          fnow = f[findex]
      end
      Q = Qnew
      if hh != 0.0
        #h[hindex += 1] = hh
        h = Base.setindex(h, hh, hindex+=1)
      end
    end
    if (Q != 0.0) || (hindex == 0)
      #h[hindex += 1] = Q
      h = Base.setindex(h, Q, hindex+=1)
    end
    return h, hindex
  end
end