function show(io::IO, m::Union{Probe{TPS},TaylorMap})
  println(io, typeof(m),":")
  extralines=0
  if GTPSA.show_header
    println(io, "-----------------------")
    desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(first(m.v).tpsa).d))
    show_GTPSA_info(io, desc)
    println(io, "-----------------------")
    extralines = 2 + (desc.nv > 0 ? 2 : 0) + (desc.np > 0 ? 2 : 0)
  end
  GTPSA.show_map(io, m.x0 + m.v, extralines, true)
end

function show(io::IO, m::Probe{Float64})
  println(io, typeof(m), ":")
  for i in LinearIndices(m.v)
    @printf(io, "%-3s %23.16le\n", "$(i):", m.v[i])
  end
end
