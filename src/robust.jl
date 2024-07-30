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
and also
Daniel VandenHeuvel

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
robust_incircle(ax,ay,bx,by,cx,cy,dx,dy) = incircle((ax,ay),(bx,by),(cx,cy),(dx,dy)) < 0

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
  @inbounds begin 
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
    return D[Dlength] # originally Dlength - 1 
  end 
end


function incircle(pa, pb, pc, pd, cache=nothing)::Float64
  @inbounds begin
      adx = pa[1] - pd[1]
      bdx = pb[1] - pd[1]
      cdx = pc[1] - pd[1]
      ady = pa[2] - pd[2]
      bdy = pb[2] - pd[2]
      cdy = pc[2] - pd[2]

      bdxcdy = bdx * cdy
      cdxbdy = cdx * bdy
      alift = adx * adx + ady * ady

      cdxady = cdx * ady
      adxcdy = adx * cdy
      blift = bdx * bdx + bdy * bdy

      adxbdy = adx * bdy
      bdxady = bdx * ady
      clift = cdx * cdx + cdy * cdy

      det = alift * (bdxcdy - cdxbdy) +
            blift * (cdxady - adxcdy) +
            clift * (adxbdy - bdxady)

      permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift +
                  (Absolute(cdxady) + Absolute(adxcdy)) * blift +
                  (Absolute(adxbdy) + Absolute(bdxady)) * clift
      errbound = iccerrboundA * permanent
      if (det > errbound) || (-det > errbound)
          return det
      end

      return _incircleadapt_round1(pa, pb, pc, pd, permanent, cache)
  end
end




@inline function Split(a)
  c = (splitter * a)
  abig = (c - a)
  ahi = c - abig
  alo = a - ahi
  return ahi, alo
end

@inline function Two_Product_Tail(a, b, x)
  ahi, alo = Split(a)
  bhi, blo = Split(b)
  err1 = x - (ahi * bhi)
  err2 = err1 - (alo * bhi)
  err3 = err2 - (ahi * blo)
  y = (alo * blo) - err3
  return y
end

@inline function Two_Product(a, b)
	x = a * b
	return x, Two_Product_Tail(a, b, x)
end

@inline function Two_Product_Presplit(a, b, bhi, blo)
  x = (a*b)
  ahi, alo = Split(a)
  err1 = x - (ahi * bhi)
  err2 = err1 - (alo * bhi)
  err3 = err2 - (ahi * blo)
  y = (alo * blo) - err3
  return x, y
end

@inline function Two_Diff_Tail(a, b, x)
  bvirt = (a - x)
  avirt = x + bvirt
  bround = bvirt - b
  around = a - avirt
  y = around + bround
  return y
end

@inline function Two_Diff(a, b)
  x = (a - b)
  return x, Two_Diff_Tail(a, b, x)
end

@inline function Two_Sum_Tail(a, b, x)
  bvirt = (x - a)
  avirt = x - bvirt
  bround = b - bvirt
  around = a - avirt
  y = around + bround
  return y
end

@inline function Two_Sum(a, b) 
  x = (a + b)
  return x, Two_Sum_Tail(a, b, x)
end

@inline function Fast_Two_Sum_Tail(a, b, x)
    bvirt = x - a
    y = b - bvirt
    return y
end

@inline function Fast_Two_Sum(a, b)
    x = (a + b)
    return x, Fast_Two_Sum_Tail(a, b, x)
end

@inline function Two_One_Diff(a1, a0, b)
	i, x0 = Two_Diff(a0, b)
	x2, x1 = Two_Sum(a1, i)
	return x2, x1, x0
end

@inline function Two_Two_Diff(a1, a0, b1, b0)
  _j, _0, x0 = Two_One_Diff(a1, a0, b0)
  x3, x2, x1 = Two_One_Diff(_j, _0, b1)
  return x3, x2, x1, x0
end

@inline function Two_One_Sum(a1, a0, b)
  _i, x0 = Two_Sum(a0, b)
  x2, x1 = Two_Sum(a1, _i)
  return x2, x1, x0
end

@inline function Two_Two_Sum(a1, a0, b1, b0)
  _j, _0, x0 = Two_One_Sum(a1, a0, b0)
  x3, x2, x1 = Two_One_Sum(_j, _0, b1)
  return x3, x2, x1, x0
end


#Absolute(a) = ((a) >= 0.0 ? (a) : -(a))
@inline function Absolute(a) 
   return abs(a)
end    

@inline function Square_Tail(a, x)
  ahi, alo = Split(a)
  err1 = x - (ahi * ahi)
  err3 = err1 - ((ahi + ahi) * alo)
  y = (alo * alo) - err3
  return y
end

@inline function Square(a)
  x = a * a
  y = Square_Tail(a, x)
  return x, y
end


@inline function estimate(elen, e)
  Q = e[1]
  for eindex in 2:elen
    Q += e[eindex]
  end
  return Q
end

@inline setindex!!(tup::Tuple, value, index) = @inbounds Base.setindex(tup, value, index)
@inline setindex!!(tup::AbstractVector, value, index) = @inbounds Base.setindex!(tup, value, index)

@inline function scale_expansion_zeroelim(elen::Int, e, b, ::Val{N})::Tuple{NTuple,Int} where {N} 
  h = ntuple(i -> zero(typeof(e[1])), Val(N))
  @inbounds begin
    bhi, blo = Split(b)
    Q, hh = Two_Product_Presplit(e[1], b, bhi, blo)
    hindex = 1
    if !iszero(hh)
        h = setindex!!(h, hh, hindex)
        hindex += 1
    end
    for eindex in 2:elen
        enow = e[eindex]
        product1, product0 = Two_Product_Presplit(enow, b, bhi, blo)
        sum, hh = Two_Sum(Q, product0)
        if !iszero(hh)
            h = setindex!!(h, hh, hindex)
            hindex += 1
        end
        Q, hh = Fast_Two_Sum(product1, sum)
        if !iszero(hh)
            h = setindex!!(h, hh, hindex)
            hindex += 1
        end
    end
    if !iszero(Q) || isone(hindex)
        h = setindex!!(h, Q, hindex)
        hindex += 1
    end
    return h, hindex - 1
end
  
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
  return fast_expansion_sum_zeroelim(elen, e, flen, f, h)
end 

@inline function fast_expansion_sum_zeroelim(elen::Int, e, flen::Int, f, h)
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
        h = setindex!!(h, hh, hindex+=1)
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
          h = setindex!!(h, hh, hindex+=1)
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
        h = setindex!!(h, hh, hindex+=1)
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
        h = setindex!!(h, hh, hindex+=1)
      end
    end
    if (Q != 0.0) || (hindex == 0)
      h = setindex!!(h, Q, hindex+=1)
    end
    return h, hindex
  end
end

function _incircleadapt_round1(pa, pb, pc, pd, permanent, cache=nothing)
  T = Float64
  @inbounds begin
    adx = pa[1] - pd[1]
    bdx = pb[1] - pd[1]
    cdx = pc[1] - pd[1]
    ady = pa[2] - pd[2]
    bdy = pb[2] - pd[2]
    cdy = pc[2] - pd[2]

    bdxcdy1, bdxcdy0 = Two_Product(bdx, cdy)
    cdxbdy1, cdxbdy0 = Two_Product(cdx, bdy)
    bc3, bc2, bc1, bc0 = Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0)
    bc = (bc0, bc1, bc2, bc3)
    axbc, axbclen = scale_expansion_zeroelim(4, bc, adx, Val(8))
    axxbc, axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, Val(16))
    aybc, aybclen = scale_expansion_zeroelim(4, bc, ady, Val(8))
    ayybc, ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, Val(16))
    adet, alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, Val(32))

    cdxady1, cdxady0 = Two_Product(cdx, ady)
    adxcdy1, adxcdy0 = Two_Product(adx, cdy)
    ca3, ca2, ca1, ca0 = Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0)
    ca = (ca0, ca1, ca2, ca3)
    bxca, bxcalen = scale_expansion_zeroelim(4, ca, bdx, Val(8))
    bxxca, bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, Val(16))
    byca, bycalen = scale_expansion_zeroelim(4, ca, bdy, Val(8))
    byyca, byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, Val(16))
    bdet, blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, Val(32))

    adxbdy1, adxbdy0 = Two_Product(adx, bdy)
    bdxady1, bdxady0 = Two_Product(bdx, ady)
    ab3, ab2, ab1, ab0 = Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0)
    ab = (ab0, ab1, ab2, ab3)
    cxab, cxablen = scale_expansion_zeroelim(4, ab, cdx, Val(8))
    cxxab, cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, Val(16))
    cyab, cyablen = scale_expansion_zeroelim(4, ab, cdy, Val(8))
    cyyab, cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, Val(16))
    cdet, clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, Val(32))

    # Julia can handle up to Val(32) without any allocations or weird behavior
    # so let's try and check if we can solve this within Val(32)... 
    # if not, then we will fall back to the full expansion... 
    abdet, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, Val(32))
    fin1, finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, Val(32))

    # check to make sure we haven't truncated too much
    if ablen < 32 || finlength < 32

      det = estimate(finlength, fin1)
      errbound = iccerrboundB * permanent
      if (det ≥ errbound) || (-det ≥ errbound)
          return det
      end

      adxtail = Two_Diff_Tail(pa[1], pd[1], adx)
      adytail = Two_Diff_Tail(pa[2], pd[2], ady)
      bdxtail = Two_Diff_Tail(pb[1], pd[1], bdx)
      bdytail = Two_Diff_Tail(pb[2], pd[2], bdy)
      cdxtail = Two_Diff_Tail(pc[1], pd[1], cdx)
      cdytail = Two_Diff_Tail(pc[2], pd[2], cdy)

      if iszero(adxtail) && iszero(bdxtail) && iszero(cdxtail) &&
        iszero(adytail) && iszero(bdytail) && iszero(cdytail)
          return det
      end

      errbound = iccerrboundC * permanent + resulterrbound * Absolute(det)
      detadd = ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) -
                                          (bdy * cdxtail + cdx * bdytail)) +
                2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
              ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) -
                                          (cdy * adxtail + adx * cdytail)) +
                2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
              ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) -
                                          (ady * bdxtail + bdx * adytail)) +
                2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx))
      det = T(det + detadd) # Had to change this to match how C handles the 2.0 multiplication with Float32
      if (det ≥ errbound) || (-det ≥ errbound)
          return det
      end
    end

    # now run the allocations ... 
    h48, h64, h1152_1, h1152_2 = _allocs_or_cache(cache)

    # if they were too short, rerun them...
    if ablen == 32 || finlength == 32
      abdetv, ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, h64)
      fin1v, finlength = fast_expansion_sum_zeroelim(ablen, abdetv, clen, cdet, h1152_1)

      # recheck termination condition...
      det = estimate(finlength, fin1v)
      errbound = iccerrboundB * permanent
      if (det ≥ errbound) || (-det ≥ errbound)
          return det
      end

      adxtail = Two_Diff_Tail(pa[1], pd[1], adx)
      adytail = Two_Diff_Tail(pa[2], pd[2], ady)
      bdxtail = Two_Diff_Tail(pb[1], pd[1], bdx)
      bdytail = Two_Diff_Tail(pb[2], pd[2], bdy)
      cdxtail = Two_Diff_Tail(pc[1], pd[1], cdx)
      cdytail = Two_Diff_Tail(pc[2], pd[2], cdy)

      if iszero(adxtail) && iszero(bdxtail) && iszero(cdxtail) &&
        iszero(adytail) && iszero(bdytail) && iszero(cdytail)
          return det
      end

      errbound = iccerrboundC * permanent + resulterrbound * Absolute(det)
      detadd = ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) -
                                          (bdy * cdxtail + cdx * bdytail)) +
                2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
              ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) -
                                          (cdy * adxtail + adx * cdytail)) +
                2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
              ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) -
                                          (ady * bdxtail + bdx * adytail)) +
                2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx))
      det = T(det + detadd) # Had to change this to match how C handles the 2.0 multiplication with Float32
      if (det ≥ errbound) || (-det ≥ errbound)
          return det
      end

    else # this means we had accuracy, but didn't have tolerance... so just copy over the final state 
      # copy fin1 to h1152_1, since 
      copyto!(h1152_1, fin1)
    end

    finnow = h1152_1
    finother = h1152_2

    if !iszero(bdxtail) || !iszero(bdytail) || !iszero(cdxtail) || !iszero(cdytail)
        adxadx1, adxadx0 = Square(adx)
        adyady1, adyady0 = Square(ady)
        aa3, aa2, aa1, aa0 = Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0)
        aa = (aa0, aa1, aa2, aa3)
    end
    if !iszero(cdxtail) || !iszero(cdytail) || !iszero(adxtail) || !iszero(adytail)
        bdxbdx1, bdxbdx0 = Square(bdx)
        bdybdy1, bdybdy0 = Square(bdy)
        bb3, bb2, bb1, bb0 = Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0)
        bb = (bb0, bb1, bb2, bb3)
    end
    if !iszero(adxtail) || !iszero(adytail) || !iszero(bdxtail) || !iszero(bdytail)
        cdxcdx1, cdxcdx0 = Square(cdx)
        cdycdy1, cdycdy0 = Square(cdy)
        cc3, cc2, cc1, cc0 = Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0)
        cc = (cc0, cc1, cc2, cc3)
    end

    if !iszero(adxtail)
        axtbc, axtbclen = scale_expansion_zeroelim(4, bc, adxtail, Val(8))
        temp16a, temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2adx, Val(16))

        axtcc, axtcclen = scale_expansion_zeroelim(4, cc, adxtail, Val(8))
        temp16b, temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, Val(16))

        axtbb, axtbblen = scale_expansion_zeroelim(4, bb, adxtail, Val(8))
        temp16c, temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, Val(16))

        temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
        temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
        finnow, finother = finother, finnow
    end
    if !iszero(adytail)
        aytbc, aytbclen = scale_expansion_zeroelim(4, bc, adytail, Val(8))
        temp16a, temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2 * ady, Val(16))

        aytbb, aytbblen = scale_expansion_zeroelim(4, bb, adytail, Val(8))
        temp16b, temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, Val(16))

        aytcc, aytcclen = scale_expansion_zeroelim(4, cc, adytail, Val(8))
        temp16c, temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, Val(16))

        temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
        temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
        finnow, finother = finother, finnow
    end
    if !iszero(bdxtail)
        bxtca, bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, Val(8))
        temp16a, temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2 * bdx, Val(16))

        bxtaa, bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, Val(8))
        temp16b, temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, Val(16))

        bxtcc, bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, Val(8))
        temp16c, temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, Val(16))

        temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
        temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
        finnow, finother = finother, finnow
    end
    if !iszero(bdytail)
        bytca, bytcalen = scale_expansion_zeroelim(4, ca, bdytail, Val(8))
        temp16a, temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2 * bdy, Val(16))

        bytcc, bytcclen = scale_expansion_zeroelim(4, cc, bdytail, Val(8))
        temp16b, temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, Val(16))

        bytaa, bytaalen = scale_expansion_zeroelim(4, aa, bdytail, Val(8))
        temp16c, temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, Val(16))

        temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
        temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
        finnow, finother = finother, finnow
    end
    if !iszero(cdxtail)
        cxtab, cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, Val(8))
        temp16a, temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2 * cdx, Val(16))

        cxtbb, cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, Val(8))
        temp16b, temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, Val(16))

        cxtaa, cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, Val(8))
        temp16c, temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, Val(16))

        temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
        temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
        finnow, finother = finother, finnow
    end
    if !iszero(cdytail)
        cytab, cytablen = scale_expansion_zeroelim(4, ab, cdytail, Val(8))
        temp16a, temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2 * cdy, Val(16))

        cytaa, cytaalen = scale_expansion_zeroelim(4, aa, cdytail, Val(8))
        temp16b, temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, Val(16))

        cytbb, cytbblen = scale_expansion_zeroelim(4, bb, cdytail, Val(8))
        temp16c, temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, Val(16))

        temp32a, temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
        temp48, temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c, temp32alen, temp32a, h48)
        finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
        finnow, finother = finother, finnow
    end

    if !iszero(adxtail) || !iszero(adytail)
        if !iszero(bdxtail) || !iszero(bdytail) || !iszero(cdxtail) || !iszero(cdytail)
            ti1, ti0 = Two_Product(bdxtail, cdy)
            tj1, tj0 = Two_Product(bdx, cdytail)
            u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
            u = (u0, u1, u2, u3)
            negate = -bdy
            ti1, ti0 = Two_Product(cdxtail, negate)
            negate = -bdytail
            tj1, tj0 = Two_Product(cdx, negate)
            v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
            v = (v0, v1, v2, v3)
            bct, bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, Val(8))

            ti1, ti0 = Two_Product(bdxtail, cdytail)
            tj1, tj0 = Two_Product(cdxtail, bdytail)
            bctt3, bctt2, bctt1, bctt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
            bctt = (bctt0, bctt1, bctt2, bctt3)
            bcttlen = 4
        else
            bct = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
            bctlen = 1
            bctt = (zero(T), zero(T), zero(T), zero(T))
            bcttlen = 1
        end

        if !iszero(adxtail)
            temp16a, temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, Val(16))
            axtbct, axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, Val(16))
            temp32a, temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2adx, Val(32))
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
            if !iszero(bdytail)
                temp8, temp8len = scale_expansion_zeroelim(4, cc, adxtail, Val(8))
                temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, Val(16))
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(cdytail)
                temp8, temp8len = scale_expansion_zeroelim(4, bb, -adxtail, Val(8))
                temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, Val(16))
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                finnow, finother = finother, finnow
            end

            temp32a, temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail, Val(32))
            axtbctt, axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, Val(8))
            temp16a, temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2adx, Val(16))
            temp16b, temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail, Val(16))
            temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
            temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(adytail)
            temp16a, temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, Val(16))
            aytbct, aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, Val(16))
            temp32a, temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2ady, Val(32))
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow

            temp32a, temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail, Val(32))
            aytbctt, aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, Val(8))
            temp16a, temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2ady, Val(16))
            temp16b, temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail, Val(16))
            temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
            temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
            finnow, finother = finother, finnow
        end
    end

    if !iszero(bdxtail) || !iszero(bdytail)
        if !iszero(cdxtail) || !iszero(cdytail) || !iszero(adxtail) || !iszero(adytail)
            ti1, ti0 = Two_Product(cdxtail, ady)
            tj1, tj0 = Two_Product(cdx, adytail)
            u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
            u = (u0, u1, u2, u3)
            negate = -cdy
            ti1, ti0 = Two_Product(adxtail, negate)
            negate = -cdytail
            tj1, tj0 = Two_Product(adx, negate)
            v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
            v = (v0, v1, v2, v3)
            cat, catlen = fast_expansion_sum_zeroelim(4, u, 4, v, Val(8))

            ti1, ti0 = Two_Product(cdxtail, adytail)
            tj1, tj0 = Two_Product(adxtail, cdytail)
            catt3, catt2, catt1, catt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
            catt = (catt0, catt1, catt2, catt3)
            cattlen = 4
        else
            cat = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
            catlen = 1
            catt = (zero(T), zero(T), zero(T), zero(T))
            cattlen = 1
        end

        if !iszero(bdxtail)
            temp16a, temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, Val(16))
            bxtcat, bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, Val(16))
            temp32a, temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2bdx, Val(32))
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
            if !iszero(cdytail)
                temp8, temp8len = scale_expansion_zeroelim(4, aa, bdxtail, Val(8))
                temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail, Val(16))
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(adytail)
                temp8, temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, Val(8))
                temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, Val(16))
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                finnow, finother = finother, finnow
            end

            temp32a, temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail, Val(32))
            bxtcatt, bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, Val(8))
            temp16a, temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2bdx, Val(16))
            temp16b, temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail, Val(16))
            temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
            temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(bdytail)
            temp16a, temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, Val(16))
            bytcat, bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, Val(16))
            temp32a, temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2bdy, Val(32))
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow

            temp32a, temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail, Val(32))
            bytcatt, bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, Val(8))
            temp16a, temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2bdy, Val(16))
            temp16b, temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail, Val(16))
            temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
            temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
            finnow, finother = finother, finnow
        end
    end

    if !iszero(cdxtail) || !iszero(cdytail)
        if !iszero(adxtail) || !iszero(adytail) || !iszero(bdxtail) || !iszero(bdytail)
            ti1, ti0 = Two_Product(adxtail, bdy)
            tj1, tj0 = Two_Product(adx, bdytail)
            u3, u2, u1, u0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
            u = (u0, u1, u2, u3)
            negate = -ady
            ti1, ti0 = Two_Product(bdxtail, negate)
            negate = -adytail
            tj1, tj0 = Two_Product(bdx, negate)
            v3, v2, v1, v0 = Two_Two_Sum(ti1, ti0, tj1, tj0)
            v = (v0, v1, v2, v3)
            abt, abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, Val(8))

            ti1, ti0 = Two_Product(adxtail, bdytail)
            tj1, tj0 = Two_Product(bdxtail, adytail)
            abtt3, abtt2, abtt1, abtt0 = Two_Two_Diff(ti1, ti0, tj1, tj0)
            abtt = (abtt0, abtt1, abtt2, abtt3)
            abttlen = 4
        else
            abt = (zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
            abtlen = 1
            abtt = (zero(T), zero(T), zero(T), zero(T))
            abttlen = 1
        end

        if !iszero(cdxtail)
            temp16a, temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, Val(16))
            cxtabt, cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, Val(16))
            temp32a, temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2cdx, Val(32))
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow
            if !iszero(adytail)
                temp8, temp8len = scale_expansion_zeroelim(4, bb, cdxtail, Val(8))
                temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail, Val(16))
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                finnow, finother = finother, finnow
            end
            if !iszero(bdytail)
                temp8, temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, Val(8))
                temp16a, temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail, Val(16))
                finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen, temp16a, finother)
                finnow, finother = finother, finnow
            end

            temp32a, temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail, Val(32))
            cxtabtt, cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, Val(8))
            temp16a, temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2cdx, Val(16))
            temp16b, temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail, Val(16))
            temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
            temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
            finnow, finother = finother, finnow
        end
        if !iszero(cdytail)
            temp16a, temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, Val(16))
            cytabt, cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, Val(16))
            temp32a, temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2cdy, Val(32))
            temp48, temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp32alen, temp32a, h48)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len, temp48, finother)
            finnow, finother = finother, finnow

            temp32a, temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail, Val(32))
            cytabtt, cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, Val(8))
            temp16a, temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2cdy, Val(16))
            temp16b, temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail, Val(16))
            temp32b, temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a, temp16blen, temp16b, Val(32))
            temp64, temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a, temp32blen, temp32b, h64)
            finother, finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len, temp64, finother)
            finnow, finother = finother, finnow
        end
    end
    return finnow[finlength]
  end
end 


function incircle_cache()
  @inbounds begin 
    # allocate h48 a little bigger to hit a 64-byte boundary... 
    total_cache = 64+64+1152+1152
    mem = Vector{Float64}(undef, total_cache)
    off = 0 
    h48 = view(mem, off .+ (1:48))
    off += 64
    h64 = view(mem, off .+ (1:64))
    off += 64 
    h1152_1 = view(mem, off .+ (1:1152))
    off += 1152
    h1152_2 = view(mem, off .+ (1:1152))
    return (h48, h64, h1152_1, h1152_2)  
  end 
end 

_allocs_or_cache(::Nothing) = incircle_cache() 
_allocs_or_cache(cache) = cache

