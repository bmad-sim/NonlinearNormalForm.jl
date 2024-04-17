# --- log ---
logpb!(na::Cint, ma::Vector{Ptr{RTPSA}}, mb::Vector{Ptr{RTPSA}}, mc::Vector{Ptr{RTPSA}}) = (@inline; mad_tpsa_logpb!(na, ma, mb, mc))
logpb!(na::Cint, ma::Vector{Ptr{CTPSA}}, mb::Vector{Ptr{CTPSA}}, mc::Vector{Ptr{CTPSA}}) = (@inline; mad_ctpsa_logpb!(na, ma, mb, mc))


"""


All arguments must have correct types. No promotion here yet. I am tired. Maybe soon

### Keyword arguments
- `'work_low` -- Tuple of 3 vectors of length nv with type equal to the lowtype of `F`.

"""
function log!(F::VectorField{T,U}, m2::DAMap{S,T,U,V}, m1::DAMap{S,T,U,V}; work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_log_work_low(F)) where {S,T,U,V}
  nv = numvars(F)

  Fx_low = work_low[1]
  m2x_low = work_low[2]
  m1x_low = work_low[3]

  @assert length(Fx_low) >= nv "Incorrect length for work_low[1] = Fx_low; received $(length(Fx_low)), should be >=$nv"
  @assert length(m2x_low) >= nv "Incorrect length for work_low[1] = m2x_low; received $(length(m2x_low)), should be >=$nv"
  @assert length(m1x_low) >= nv "Incorrect length for work_low[1] = m1x_low; received $(length(m1x_low)), should be >=$nv"
  @assert eltype(outx_low) == lowtype(outT) "Incorrect eltype of work_low[1] = outx_low. Received $(eltype(Fx_low)), should be $(lowtype(T))"
  @assert eltype(m2x_low) == lowtype(outT) "Incorrect eltype of work_low[2] = m2x_low. Received $(eltype(m2x_low)), should be $(lowtype(T))"
  @assert eltype(m1x_low) == lowtype(outT) "Incorrect eltype of work_low[3] = m1x_low. Received $(eltype(m1x_low)), should be $(lowtype(T))"
 
end


"""
    logpb(mf::Vector{<:Union{TPS,ComplexTPS}}, mi::Vector{<:Union{TPS,ComplexTPS}}=vars(first(F)))

Given a final map `mf` and initial map `mi`, this function calculates the vector field `F`
such that `mf=exp(F⋅∇)mi`. If `mi` is not provided, it is assumed to be the identity.

```julia-repl
julia> d = Descriptor(2,10); x = vars()[1]; p = vars()[2];

julia> time = 0.01; k = 2; m = 0.01;

julia> h = p^2/(2m) + 1/2*k*x^2;

julia> hf = getvectorfield(h);

julia> map = exppb(-time*hf);

julia> logpb(map) == -time*hf
true
```
"""
function log(mf::Vector{<:Union{TPS,ComplexTPS}}, mi::Vector{<:Union{TPS,ComplexTPS}}=vars(first(mf)))
  desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(mf[1].tpsa).d))
  if length(mf) != desc.nv || length(mi) != desc.nv
    error("Vector length != number of variables in the GTPSA")
  end
  ma1, mb1 = promote(mf, mi)
  m1 = map(t->t.tpsa, ma1) 
  m2 = map(t->t.tpsa, mb1) 
  mc = zero.(ma1)
  m3 = map(t->t.tpsa, mc)
  GC.@preserve ma1 mb1 logpb!(Cint(length(mf)), m1, m2, m3)  
  return mc
end
