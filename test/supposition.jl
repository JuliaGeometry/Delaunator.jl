@testset "Supposition" begin 
  import AdaptivePredicates
  import ExactPredicates
  # The following tests compare the output to AdaptivePredicates.jl and ExactPredicates.jl 
  # these are correct for the range -142:201 as explained in the README for 
  # adaptive predicates. 
  check_range(x::Float64) = iszero(x) || exponent(abs(x)) ∈ -142:201 
  using Supposition

  fgen = Data.Floats{Float64}(infs=false, nans=false)
  r2gen = @composed _tuple(a=fgen, b=fgen) = (a, b)
  r3gen = @composed _tuple(a=fgen, b=fgen, c=fgen) = (a, b, c)
  # Compare against sign from ExactPredicates
  @check function _orient2(p=r2gen, q=r2gen, r=r2gen)
    assume!(all(check_range, (p...,q...,r...)))
    ap = sign(Delaunator.orient2(p, q, r))
    c = ExactPredicates.orient(p, q, r)
    event!("AdaptiveOrient2", ap)
    event!("ExactPredicatesOrient2d", c)
    ap == c
  end
  @check function _orient2ap(p=r2gen, q=r2gen, r=r2gen)
    assume!(all(check_range, (p...,q...,r...)))
    ap = Delaunator.orient2(p, q, r)
    c = AdaptivePredicates.orient2(p, q, r)
    event!("AdaptiveOrient2", ap)
    event!("AdaptivePredicatesOrient2d", c)
    ap == c
  end
  
  @check function _incircle(p=r2gen, q=r2gen, r=r2gen, s=r2gen)
    assume!(all(check_range, (p...,q...,r...,s...)))
    ap = sign(Delaunator.incircle(p, q, r, s))
    c = ExactPredicates.incircle(p, q, r, s)
    event!("AdaptiveIncircle", ap)
    event!("ExactPredicatesIncircle", c)
    ap == c
  end
  @check function _incircle_ap(p=r2gen, q=r2gen, r=r2gen, s=r2gen)
    assume!(all(check_range, (p...,q...,r...,s...)))
    ap = Delaunator.incircle(p, q, r, s)
    c = AdaptivePredicates.incircle(p, q, r, s)
    event!("AdaptiveOrient2", ap)
    event!("AdaptivePredicatesIncircle", c)
    ap == c
  end
  

end 