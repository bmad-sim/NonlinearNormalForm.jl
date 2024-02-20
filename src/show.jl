function show(io::IO, m::Union{Probe,TaylorMap})
  println(io, typeof(m),":")
  extralines=0
  if GTPSA.show_header
    println(io, "-----------------------")
    desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(first(m.x).tpsa).d))
    show_GTPSA_info(io, desc)
    println(io, "-----------------------")
    extralines = 2 + (desc.nv > 0 ? 2 : 0) + (desc.np > 0 ? 2 : 0)
  end
  println(io, "Coordinate System Origin:")
  for i in LinearIndices(m.x)
    @printf(io, "%-3s %23.16le\n", "$(i):", m.x0[i])
  end
  println(io, "\nOrbital Ray:")
  if eltype(m.x) == TPS
    GTPSA.show_map(io, m.x, 1+length(m.x)+extralines, true)
  else
    for i in LinearIndices(m.x)
      @printf(io, "%-3s %23.16le\n", "$(i):", m.x[i])
    end
  end
end