#=

Functions to check TaylorMap/VectorField sanity during operations 
involving possible two or more TaylorMaps, specifically, ensuring that types 
match. checkstates is provided for only a states check (spin and stochastic), 
and checkinplace is provided for an in-place operation (!) where promotion 
may be occurring (first argument is destination). All out-of-place operators 
must internally call in-place operators, so there is no equivalent 
out-of-place version.

The use of these functions is preferred to using "where {S,T,U,V}" in the 
function definition itself because these functions give descriptive errors, 
rather than just "function not found". Furthermore, the function headers 
become way too verbose when using "where". All sanity checks are type stable, 
and the checks should be optimized away in the JIT compilation if valid.

=#
@inline checkspin(stuff...) = all(x->isnothing(x.q), stuff) || all(x->!isnothing(x.q), stuff) || error("Atleast one map/vector field includes spin while others do not")
@inline checkstochastic(maps::TaylorMap...) = all(x->isnothing(x.s), maps) || all(x->!isnothing(x.s), maps) || error("Atleast one map includes stochasticity while others do not")
@inline function checkvarsparams(stuff...)
  NV = nvars(first(stuff))
  NN = ndiffs(first(stuff))
  return all(stuff) do m
    nvars(m) == NV || error("Atleast one argument has TPSA with different number of variables")
    ndiffs(m) == NN || error("Atleast one argument has TPSA with different number of parameters")
  end
end

@inline checkstates(stuff...) = checkvarsparams(filter(x->(x isa Union{TaylorMap,VectorField}), stuff)...) && checkspin(filter(x->(x isa Union{TaylorMap,VectorField}), stuff)...) && checkstochastic(filter(x->(x isa TaylorMap), stuff)...)

# checkinplace is preferred to using the "where {S,..}.." syntax as this 
# gives descriptive errors rather than just "function not found"
@inline function checkinplace(m::Union{TaylorMap,VectorField}, stuff...; internal_promotion=false)
  checkstates(m, stuff...)
  # Checks that the output map has all types properly promoted
  # or NOT promoted if unneccessary
  maps = filter(x->(x isa TaylorMap), stuff)
  mapsvfs = filter(x->(x isa Union{TaylorMap,VectorField}), stuff)
  nums = filter(x->(x isa Number), stuff)
  eltypes = map(x->typeof(x), nums) # scalars only affect x and Q, not x0 or E in FPP

  xtypes = map(x->eltype(x.x), mapsvfs)
  xnumtypes = map(x->TI.numtype(x),xtypes)

  if m isa TaylorMap
    x0types = map(x->eltype(x.x0), maps)
    if internal_promotion
      outx0type = promote_type(x0types..., xnumtypes...) # reference orbit in composition is affected by orbital part
      eltype(m.x0) == outx0type || error("Output $(typeof(m)) reference orbit type $(eltype(m.x0)) must be $outx0type")
    else
      all(t->t==eltype(m.x0), x0types) || error("All reference orbit eltypes must be the same (output eltype is $(eltype(m.x0)))")
    end
  end

  if internal_promotion
    outxtype = promote_type(xtypes..., eltypes...)
    eltype(m.x) == outxtype || error("Output $(typeof(m)) orbital ray type $(eltype(m.x)) must be $xtype")
  else
    all(t->t==eltype(m.x), xtypes) || error("All orbital ray eltypes must be the same (output eltype is $(eltype(m.x)))")
  end

  if !isnothing(m.q)
    qtypes = map(x->eltype(x.q), mapsvfs)
    if internal_promotion
      outqtype = promote_type(qtypes..., eltypes...)
      eltype(m.q) == outqtype || error("Output $(typeof(m)) quaternion type $(eltype(m.q)) must be $outqtype")
    else
      all(t->t==eltype(m.q), qtypes) || error("All quaternion eltypes must be the same (output eltype is $(eltype(m.q)))")
    end
  end

  # Part of the promotion is stochasticity:
  # the output map must include stochasticity if any input includes stochasticity:
  if m isa TaylorMap && !isnothing(m.s)
    stypes = map(x->eltype(x.s), maps)
    if internal_promotion
      outtype = promote_type(xnumtypes..., stypes...)
      eltype(m.s) == outtype || error("Output $(typeof(m)) stochastic matrix type $(eltype(m.s)) must be $outtype")
    else
      all(t->t==eltype(m.s), stypes) || error("All stochastic matrix eltypes must be the same (output eltype is $(eltype(m.s)))")
    end
  end

  return true
end

